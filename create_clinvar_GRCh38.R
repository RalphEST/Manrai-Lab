# Arjun Manrai
# Ralph Estanboulieh
# NHLBI BDC

# libraries
library(RSQLite)
library(tidyverse)
library(dplyr)
library(DBI)
library(odbc)
library(RMariaDB)
library(rlist)
library(rjson)
library(httr)
library(jsonlite)
library(rvest)
setwd(dirname(rstudioapi::getSourceEditorContext()$path))


# helper functions
filter.geneinfo<- function(panel, df) {
  # filters INFO column of df passed for strings in panel
  panel <- paste0(panel,':')
  s  <- paste0('GENEINFO=', panel, collapse="|")
  df <- df %>%
    filter(grepl(s,INFO))
  return(df)
}

annotate.clnsig <- function(info) {
  # isolates and returns clinical significance into an array
  pattern <- "CLNSIG=[^;]*;"
  clnsig <- str_extract(info,pattern)
  clnsig <- substr(clnsig,8,nchar(clnsig)-1)
  return(clnsig)
}

onehot.clnsig <- function(df) {
  # adds 4 columns (on for VUS, Pathogenic, Likely pathogenic, and Conflicting) in a one-hot encoding fashion
  df <- df %>%
    mutate(VUS = ifelse(CLNSIG=='Uncertain_significance',1,0)) %>%
    mutate(Pathogenic = ifelse(CLNSIG=='Pathogenic' | CLNSIG=='Pathogenic/Likely_pathogenic',1,0)) %>%
    mutate(LikelyPathogenic = ifelse(CLNSIG=='Likely_pathogenic',1,0)) %>%
    mutate(Conflicting = ifelse(CLNSIG=='Conflicting_interpretations_of_pathogenicity',1,0))
  return(df)
}

# main
# download GRCh38 from https://www.ncbi.nlm.nih.gov/genome/guide/human/
src.directory <- ''
clinvar0 <- read_tsv(paste0(src.directory,'GRCh38_latest_clinvar.vcf'),
                 comment='#',
               col_names = c('CHROM','POS','ID','REF',
                             'ALT','QUAL','FILTER','INFO'))
panel <- c('MYH7','TNNT2','MYBPC3','TPM1','MYL3','TNNI3','MYL2','ACTC1')

clinvar <- filter.geneinfo(panel,clinvar0)


clinvar$CLNSIG <- annotate.clnsig(clinvar$INFO)
clinvar <- onehot.clnsig(clinvar)
clinvar <- clinvar %>%
  mutate_all(type.convert) %>%
  mutate_if(is.double,as.integer)

# summarize
# clinvar %>% group_by(CLNSIG) %>% summarise(n=n()) %>% arrange(desc(n))

# store as SQLite db
remove_existing_db <- 0 # warning this will delete previous db!
if (remove_existing_db) {
  file.remove("../../data/prebuilt-dbs/clinvar_grch38.db")
}
conn <- dbConnect(SQLite(), dbname = "clinvarSQL.db")
dbWriteTable(conn, "clinvar", clinvar)
dbListTables(conn) # double check table was created


dbDisconnect(conn)
con <- dbConnect(RMariaDB::MariaDB(), group = "")

# connecting to and querying the UCSC Genome Browser via RMariaDB

# creating connection with the database
con <- dbConnect(
  drv = RMariaDB::MariaDB(), 
  username = "genome",
  password = NULL, 
  host = "genome-mysql.soe.ucsc.edu", 
  port = 3306,
  db = "hg38"
)


# manually find a gene
findGene <- function(chr, pos) {
  x = dbGetQuery(con, paste0("SELECT curated.chrom, curated.txStart, curated.txEnd, link.name as geneSymbol \
                 FROM ncbiRefSeqCurated curated, ncbiRefSeqLink link\
                 WHERE curated.name=link.id \
                       AND curated.chrom ='chr",chr,"'
                       AND curated.txStart < ",pos,"
                       AND curated.txEnd >",pos))
  return(x)
}

# our checking function: adds a "check" column (TRUE if position in panel gene, FALSE otherwise)
refSeq.check = function(df, panel) {
  df$check_refSeq = FALSE # initialize "check" column
  check <- paste(panel, collapse = " OR link.name LIKE ") # prepare query
  ranges <- dbGetQuery(con, paste0("SELECT curated.chrom, curated.txStart, curated.txEnd, link.name as geneSymbol \
                 FROM ncbiRefSeqCurated curated, ncbiRefSeqLink link\
                 WHERE curated.name=link.id \
                       AND (link.name LIKE", check, ")"))
  ranges <- unique(ranges) # make unique
  for (row in 1:nrow(ranges)) {
    chr = str_sub(ranges[row, "chrom"], 4, -1)
    start = ranges[row, "txStart"]
    end = ranges[row, "txEnd"]
    df$check_refSeq[df$CHROM == chr & df$POS >= start & df$POS <= end] = TRUE
  }
  return(df)
}

# evaluate the positions
panel_ <- c("'MYH7'","'TNNT2'","'MYBPC3'","'TPM1'","'MYL3'","'TNNI3'","'MYL2'","'ACTC1'")
clinvar_ <- refSeq.check(clinvar0, panel_)
sum(clinvar_$check_refSeq) # 6526

# this is different from the GENEINFO parsing method, let's compare both methods

# another checking function: adds a "check" column (TRUE if position in panel gene, FALSE otherwise) 
GENEINFO.check<- function(panel, df) {
  # filters INFO column of df passed for strings in panel
  panel <- paste0(panel,':')
  s  <- paste0('GENEINFO=', panel, collapse="|")
  df$check_GENEINFO = FALSE
  df$check_GENEINFO[grepl(s,df$INFO)] = TRUE
  return(df)
}

clinvar_ <- GENEINFO.check(panel, clinvar_)
sum(clinvar_$check_GENEINFO) # 6406 

# let's compare more directly and select the rows where at least one of the checks is true
atLeastOne = clinvar_ %>% filter(check_GENEINFO | check_refSeq)

# we need to name the genes of each of the 6000+ positions

# in order to do this efficiently (SQL queries are very slow), we use the fact that 
# the table is ordered to bound the number of such queries by the number of distinct genes 
# in the data set. Essentially, we keep track of the gene we are in with a positional variable, 
# and every time we 'leave' a gene locus, we run a new SQL query

# helper function: check whether a position is within a range
gene.range <- function(pos, range) {
    # pos = c(chrom, pos)
    # range = c(chrom, start, end)
    return(paste0('chr', pos[1])==range[1] && pos[2]>=range[2] && pos[2]<=range[3])
}

# the gene naming function
refSeqInfo.gene <- function(df) {
    curr_gene = c('chr1',1,0) # initialize current gene
    curr_gene_name = "GENE" # initialize current gene name
    for (row in 1:nrow(df)) {
      pos = c(df[row, "CHROM"], df[row, "POS"])
      if (!gene.range(pos, curr_gene)) {
        print(paste0("row: ",row))
        query = paste0("SELECT curated.chrom, curated.txStart, curated.txEnd, link.name as geneSymbol
                   FROM ncbiRefSeqCurated curated, ncbiRefSeqLink link
                   WHERE curated.name=link.id
                         AND curated.chrom='chr",pos[1],"'
                         AND curated.txStart < ",pos[2],"
                         AND curated.txEnd > ",pos[2])
        query = dbGetQuery(con, query)
        if (nrow(query)>0){
          curr_gene = c(query$chrom[1], query$txStart[1], query$txEnd[1])
          curr_gene_name = query$geneSymbol[1]
        } else {
          curr_gene_name = "NO_GENE"
        }
      }
      df[row, "refSeqGene"] = curr_gene_name
      # name using GENEINFO
      df[row, "GENEINFOGene"] = str_extract(df[row, "INFO"], "(?<=GENEINFO=)(.*?)(?=\\:)")
    }
    return(df)
  }

atLeastOne <- refSeqInfo.gene(atLeastOne)

sum(is.na(atLeastOne$GENEINFOGene)) # 251 rows without GENEINFO

onlyInRefSeq <- atLeastOne %>% filter(check_refSeq & !check_GENEINFO) # 251 rows
onlyInGENEINFO <- atLeastOne %>% filter(!check_refSeq & check_GENEINFO) # 131 rows (where refSeq disagrees with GENEINFO)


# >> MedGen Concept ID 
# HCM                     | C0007194
# Primary Familial HCM    | C0949658
# Familial HCM 1          | C3495498
# Familial HCM 2          | C1861864
# Familial HCM 8          | C1837471
# >> HPO
# HCM                     | HP:0001638
# Asym. Septal Hypertrophy| HP:0001670
# Conc. Hypert. Cardiom.  | HP:0005157
# Apic. Hypert. Cardiom.  | HP:0031992
# >> ORHPANET
# HCM                     | ORPHA:217569
# Rare Fam. Dis. with HCM | ORPHA:99739
# Non-fam. HCM            | ORPHA:217598
# >> SNOMED_CT
# HCM                     | 233873004
# Fetal HCM               | 472324001
# HCM with other dis.     | 445336009
# HCM with obstr.         | 195020003
# Hyper. Obstr. Cardiom.  | 45227007
# Primary HCM             | 700065003
# >>OMIM
# Infantile HCM           | 500006

patterns = list("HPO"="(?<=Human_Phenotype_Ontology\\:HP\\:)(.*?)(?=[\\|\\,\\;])", 
                "SNOMED_CT"="(?<=SNOMED_CT\\:)(.*?)(?=[\\|\\,\\;])",
                "MedGen"="(?<=MedGen\\:)(.*?)(?=[\\|\\,\\;])",
                "Orphanet"="(?<=ORPHA)(.*?)(?=[\\|\\,\\;])",
                "OMIM"="(?<=OMIM\\:)(.*?)(?=[\\|\\,\\;])",
                "MONDO"="(?<=MONDO\\:MONDO\\:)(.*?)(?=[\\|\\,\\;])")

extract.codes <- function (df, pattern) {
  codes = c()
  for(row in 1:nrow(df)) {
    s = df[row, "INFO"]
    extr_codes = str_extract_all(s, pattern)
    codes = c(codes, 
              unlist(extr_codes))
  }
  return(unique(codes))
}

hpo = extract.codes(atLeastOne, 
                    patterns$HPO)
snomed_ct = extract.codes(atLeastOne, 
                          patterns$SNOMED_CT)
medgen = extract.codes(atLeastOne, 
                       patterns$MedGen)
mondo = extract.codes(atLeastOne, 
                      patterns$MONDO)
orphanet = extract.codes(atLeastOne, 
                         patterns$Orphanet)
omim = extract.codes(atLeastOne, 
                     patterns$OMIM)


# interpreting OMIM codes
omim_meanings = tibble(entry = omim)
for(row in 1:nrow(omim_meanings)){
  res = GET("https://api.omim.org/api/entry/search", 
            query = list(search = omim_meanings[row, "entry"], 
                         apiKey = "J-YQ2DurSO2YVHYi0GY8HQ", 
                         format = "json",
                         include = "clinicalSynopsis",
                         retrieve = "clinicalSynopsis"))
  res = rawToChar(res$content)
  res = fromJSON(res)
  titles = res$omim$searchResponse$clinicalSynopsisList$clinicalSynopsis$preferredTitle
  omim_meanings[row, "Clinical Synopsis"] = paste(titles, collapse = " | ")
}

# interpreting MedGen codes
medgen_meanings = tibble(entry = medgen)

for(row in 1:nrow(medgen_meanings)){
  print(row)
  scrape <- read_html(paste0("https://www.ncbi.nlm.nih.gov/medgen/", medgen_meanings[row, "entry"]))
  medgen_meanings[row, "MedGen Title"] = scrape %>% html_nodes(".MedGenTitleText") %>% html_text()
}

# interpreting HPO codes
hpo_meanings = tibble(entry = hpo)

for(row in 1:nrow(hpo_meanings)){
  res = GET(paste0("http://hpo.jax.org/api/hpo/term/HP:", hpo_meanings[row, "entry"]))
  res = rawToChar(res$content)
  res = fromJSON(res)
  name = res$details$name
  hpo_meanings[row, "Clinical Synopsis"] = ifelse(is.null(name), "NO RESULT", name)
}

# interpreting MONDO codes
mondo_meanings = tibble(entry = mondo)

for(row in 1:nrow(mondo_meanings)){
  res = GET(paste0("https://api.monarchinitiative.org/api/bioentity/MONDO:", mondo_meanings[row, "entry"]))
  res = rawToChar(res$content)
  res = fromJSON(res)
  name = res$label
  mondo_meanings[row, "Label"] = ifelse(is.null(name), "NO RESULT", name)
}

# interpreting Orphanet codes
orphanet_meanings = tibble(entry = orphanet)

for(row in 1:nrow(orphanet_meanings)){
  res = GET(paste0("https://api.monarchinitiative.org/api/bioentity/ORPHA:", orphanet_meanings[row, "entry"]))
  res = rawToChar(res$content)
  res = fromJSON(res)
  name = res$label
  orphanet_meanings[row, "Label"] = ifelse(is.null(name), "NO RESULT", name)
}


# interpreting SNOMED_CT codes
snomed_ct_meanings = tibble(entry = snomed_ct)

for(row in 1:nrow(snomed_ct_meanings)){
  scrape <- read_html(paste0("https://snomedbrowser.com/Codes/Details/", snomed_ct_meanings[row, "entry"]))
  scrape = scrape %>% html_node("#selected-concept") %>% html_text()
  name = str_extract(scrape, "(?<=Name:)(.*?)(?=\\\r)")
  snomed_ct_meanings[row, "Name"] = ifelse(is.null(name), "NO RESULT", name)
}

# turn everything into .CSV files for manual review
write_csv(hpo_meanings, "hpo_codes.csv" )
write_csv(omim_meanings, "omim_meanings.csv" )
write_csv(mondo_meanings, "mondo_meanings.csv" )
write_csv(snomed_ct_meanings, "snomed_ct_meanings.csv" )
write_csv(medgen_meanings, "medgen_meanings.csv" )
write_csv(orphanet_meanings, "orphanet_meanings.csv" )

# select wanted codes HERE!
hpo_wanted = c("1639", "1638", "1712", "5171")
omim_wanted = 
medgen_wanted = medgen[1:5]
orphanet_wanted = orphanet[1:5]
snomed_ct_wanted = snomed_ct[1:5]
mondo_wanted = mondo[1:5]

hpo_pattern = paste0("Human_Phenotype_Ontology:HP:(", paste(hpo_wanted, collapse = "|"), ")")
omim_pattern = paste0("OMIM:(", paste(omim_wanted, collapse = "|"), ")")
medgen_pattern = paste0("MedGen:(", paste(medgen_wanted, collapse = "|"), ")")
orphanet_pattern = paste0("ORPHA(", paste(orphanet_wanted, collapse = "|"), ")")
snomed_ct_patterd = paste0("SNOMED_CT:(", paste(snomed_ct_wanted, collapse = "|"), ")")
mondo_pattern = paste0("MONDO:(", paste(mondo_wanted, collapse = "|"), ")")

global_pattern = paste(hpo_pattern, 
                        omim_pattern, 
                        medgen_pattern, 
                        orphanet_pattern, 
                        snomed_ct_patterd, 
                        mondo_pattern, sep = "|")

hcm <- function (df, pattern) {
  for(row in 1:nrow(df)) {
    print(row)
    df[row, "HCM"] = str_detect(df[row, "INFO"], pattern)
  }
  return(df)
}

hcm_clinvar = hcm(atLeastOne, hpo_pattern)
