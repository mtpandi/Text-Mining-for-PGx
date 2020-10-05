#how the f to keep only drugs
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
                                  'Drug,Ion', 'Drug,Metabolite')) %>% select(ChemicalID)



xlsx::write.xlsx(ctd, 'ChemicalstoFilter.xlsx', row.names = F)

# this file needs to be manually curated and re-uploaded to furter 
# remove chemical that mostly add noise 

noise <- readxl::read_xlsx('noisychemicals.xlsx')
noise <- noise %>% select(ChemicalID) %>% unlist() %>% unname()

ctd <- ctd %>% subset(!(ChemicalID %in% noise))
write.csv(ctd, 'FilterChemicals.csv', row.names = F) #this file will be used for filtering

