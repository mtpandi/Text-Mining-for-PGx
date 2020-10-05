#prepare the gold standard
library(dplyr)
library(pbapply)

#combine ctd with pharmgkb
ctd <- read.delim('CTD.tsv', header = T, sep = '\t', stringsAsFactors = F)
ctd[ctd =='']<-NA
ctd <- ctd %>% subset(!is.na(DrugBankIDs))
ctd$ChemicalID<- pblapply(ctd$ChemicalID, function(x) gsub("MESH:", "", x)) 
ctd <- ctd %>%  mutate_if(is.list, as.character)

#pharmgkb 
pgkb <-  read.delim('chemicals.tsv', header = T, sep = '\t', stringsAsFactors = F)
pgkb[pgkb=='']<-NA

pgkb <- pgkb %>% subset(Type!='Drug Class') %>% subset(!is.na(External.Vocabulary)) %>%
  select(PharmGKB.Accession.Id, Type, Cross.references)


getdrugbank <- function(df){
  
  term <- 'DrugBank:'
  id <- df[['PharmGKB.Accession.Id']]
  text <- df[['Cross.references']]
  
  if(grepl(term, text)){
    
    
    x<- strsplit(text, ',') %>% unlist 
    x <- x[grep(term, x)]
    
    x<- gsub(term, "", x)
    
    
    res <- data.frame(DrugBankIDs = x, PharmGKBID = id, stringsAsFactors = F) 
    
    
  }else{
    
    res <- data.frame(DrugBankIDs = 'None',  PharmGKBID = id, stringsAsFactors = F)
  } 
  return(res)
}

ids <- pbapply(pgkb, 1, getdrugbank) %>% bind_rows() %>% subset(DrugBankIDs!='None')

ctd <- left_join(ctd, ids, by='DrugBankIDs')
ctd <- left_join(ctd, pgkb[,c('PharmGKB.Accession.Id', 'Type')], by=c('PharmGKBID'='PharmGKB.Accession.Id'))

ctd <- ctd %>% subset(Type %in% c('Drug', 'Drug,Biological Intermediate', 'Prodrug',
                                  'Drug,Ion', 'Drug,Metabolite'))

#in addition, same chemicals were removed, after it was noted 
#that they are mostly noise (highlighted in the 'ChemicalstoFilter.xlsx)
noise <- readxl::read_xlsx('noisychemicals.xlsx')
noise <- noise %>% select(ChemicalID) %>% unlist() %>% unname()

ctd <- ctd %>% subset(!(ChemicalID %in% noise)) %>% select(ChemicalID, PharmGKBID)
rm(noise, ids, pgkb)

#load gold standards
pos <- readxl::read_excel('relationships_true_positives.xlsx')
neg <- readxl::read_excel('relationships_true_negatives_v2.xlsx')


#map entity id 2 and pharmgkb id
pos <- left_join(pos, ctd, by=c('Entity2_id'= 'PharmGKBID'))
neg <- left_join(neg, ctd, by=c('Entity2_id'= 'PharmGKBID'))

#keep only columns of interest 
pos <- pos %>% 
  subset(!is.na(ChemicalID)) %>%
  subset(!is.na(PMIDs)) %>%
  select(Entity1_name, ChemicalID, PMIDs) %>%
  mutate(Class = 'corr')


neg <- neg %>% 
  subset(!is.na(ChemicalID)) %>%
  subset(!is.na(PMIDs)) %>%
  select(Entity1_name, ChemicalID, PMIDs) %>%
  mutate(Class = 'nocorr')

#join the pos and neg dfs
gs.set <- rbind(pos, neg)


#keep a list of the corresponding pmids
pmids.set <- gs.set %>% select(PMIDs) %>% 
  mutate(PMIDs = strsplit(PMIDs, ';')) %>%
  tidyr::unnest(PMIDs) %>% unique() 

write.csv(pmids.set, 'pmidsset.csv', row.names = FALSE)

write.csv(gs.set, 'gs.set.csv', row.names = FALSE)