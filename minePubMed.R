library(easyPubMed)
library(dplyr)
library(RCurl)
library(pbapply)
library(data.table)
library(tidyr)
library(stringr)

source('pubtator_v3.R')

source('TermNormalization.R')

paperids <- function(inpquery){
  
  #function needed
  get_pmcids <- function(x) {
    res<-data.frame(PMID=character(), PMCID=character())
    json_file <- RJSONIO::fromJSON(getURL(paste("https://www.ncbi.nlm.nih.gov/pmc/utils/idconv/v1.0/?ids=",
                                                x, "&idtype=pmid&format=json&tool=my_tool&email=my_email@example.com" ,
                                                sep = "")))
    for(i in 1:length(json_file$records)){
      if("pmcid" %in% names(json_file$records[[i]])){
        pmcid<- json_file$records[[i]]['pmcid'] %>% as.character()
      }else{
        pmcid <-NA
      }
      pmid<-json_file$records[[i]]['pmid'] %>% as.character()
      ires<-data.frame(pmcid,pmid, stringsAsFactors = FALSE)
      res<-rbind(res, ires)
    }
    return(res)
  }
  
  #access pubmed
  my_query <- get_pubmed_ids(inpquery)
  
  # Download by 5000-item batches
  batches <- seq(from = 0, to = my_query$Count, by = 5000)
  abstracts <- pblapply(batches,  function(i) {
    fetch_pubmed_data(my_query, format = "uilist", retmax = 5000, retstart = i)  
  })
  
  pmids <- pblapply(abstracts, as.integer) %>% unlist
  
  #convert to pmcid 
  pmid_start <- seq(from = 1, to = length(pmids), by = 200)
  ids <- pblapply(pmid_start, function(i){
    pmid_end <- ifelse((i+199<=length(pmids)),i+199, length(pmids))
    search<-paste0(pmids[i:pmid_end], collapse = ",")
    get_pmcids(search)
  })
  data <- do.call("rbind", ids) 
  
  return(data)
  
}


# Query pubmed and fetch many results
query<-'((pharmacogen*[Text Word]) AND ("humans"[MeSH Terms])) NOT (Review[ptyp])'
data<-paperids(query)



#Get PubTaror annotations
annots<-pubtator_anns(data[is.na(data$pmcid),'pmid'], data[,'pmcid']) 


#keep only the ones I need
annots <- annots %>% subset(type %in% c('Chemical', 'Mutation', 'Gene'))

#'mask' urls
 
maskurls <- function(df){
  
  #split the df
  ok <- df %>% subset(!(grepl('(http[^ ]*)|(www\\.[^ ]*)', TEXT)))
  tofix <- df %>% subset((grepl('(http[^ ]*)|(www\\.[^ ]*)', TEXT)))
  
  #create an invetory of urls
  aa <- lapply(tofix$TEXT, function(x) str_extract_all(x, '(http[^ ]*)|(www\\.[^ ]*)'))
  aa <- aa %>% unlist() %>% unique() 
  aa <- str_replace_all(aa, '[[:punct:]]+$', '') 
  aa <- unique(aa)
  
  resurls <- data.frame('original'=aa, stringsAsFactors = F)
  resurls$new <- lapply(resurls$original, function(x) str_replace_all(x, '\\.', '\\_'))
  resurls$new <- as.character(resurls$new)
  
  #need a function to search through the TEXT column and wherever she finds a term from 
  for(i in 1:nrow(resurls)){
    old <- resurls[i, 'original']
    new <- resurls[i, 'new']
    
    tofix <- tofix %>% 
      mutate(TEXT=ifelse(grepl(old, TEXT,  fixed = T),
                         str_replace_all(TEXT, stringr::fixed(old), new), TEXT))
  }
  
  df <- rbind(ok, tofix)
  return(df)
}

annots <- maskurls(annots)

#fix the issue with periods (affects the way sentences are created)
#before splitting the sentences I must replace any . found in the the TEXT_ann column to avoid splitting 
#the replacement must be made inside the TEXT as well as inside the TEXT_ann columns

maskper <- function(data){
  
  ok <- data %>%
    subset(!grepl('p\\.|c\\.|n\\.|m\\.|r\\.|g\\.|v\\.|V\\.', TEXT_ann)) %>%
    mutate(ann=TEXT_ann)
  
  dots <- data %>%
    subset(grepl('p\\.|c\\.|n\\.|m\\.|r\\.|g\\.|v\\.|V\\.', TEXT_ann)) 
  
  unique.TEXT <- dots %>% select(TEXT) %>% unique %>% unlist 
  df.list <- list()
  
  for(i in 1:length(unique.TEXT)){
    
    df <- dots %>% subset(TEXT == unique.TEXT[[i]])  %>%
      mutate(ann=TEXT_ann,
             TEXT_ann=ifelse(grepl('p\\.|c\\.|n\\.|m\\.|r\\.|g\\.|v\\.|V\\.', TEXT_ann), 
                             str_replace_all(TEXT_ann, c('p\\.' = 'p_',
                                                         'c\\.' = 'c_',
                                                         'n\\.'='n_',
                                                         'm\\.'='m_',
                                                         'r\\.'='r_',
                                                         'g\\.'='g_',
                                                         'v\\.'='v_',
                                                         'V\\.'='V_'))
                             , TEXT_ann))  
    
    text <- df[1,'TEXT']
    tochange <- df %>%
      distinct(TEXT_ann, .keep_all = T)
    
    for(j in 1:nrow(tochange)){
      old <- tochange[j, 'ann']
      new <- tochange[j, 'TEXT_ann']
      text <- str_replace_all(text, stringr::fixed(old), new)
    }
    
    df <- df %>% mutate(TEXT = text)
    df.list[[i]] <- df
  }
  
  newdots <- data.table::rbindlist(df.list) %>% as.data.frame()
  
  res <- rbind(ok, newdots)
  return(res)
  
  
}

annots <- maskper(annots)

#also, Fig., Figs., fig. and figs. etc must also be replaced 
annots <- annots %>%
  mutate(TEXT = ifelse(grepl('Fig\\.|fig\\.|FIG\\.|Figs\\.|figs\\.|FIGS\\.|\\bi\\.e\\.|\\be\\.g\\.|\\bcf\\.|\\bS\\.E\\.|\\bal\\.|\\bvs\\.|\\betc\\.|\\bS\\.D\\.|\\=\\.|\\=\\s\\.', TEXT), 
                       str_replace_all(TEXT, c('Fig\\.'='Fig',
                                               'fig\\.'='fig',
                                               'FIG\\.'='FIG',
                                               'Figs\\.'='Figs',
                                               'figs\\.'='figs',
                                               'FIGS\\.'='FIGS',
                                               '\\bi\\.e\\.'='ie',
                                               '\\be\\.g\\.'='eg',
                                               '\\bcf\\.'='cf',
                                               '\\bS\\.E\\.'='SE',
                                               '\\bal\\.'='al',
                                               '\\bvs\\.'='vs' ,
                                               '\\betc\\.'='etc',
                                               '\\bS\\.D\\.'='SD',
                                               '\\=\\.'='=0.',
                                               '\\=\\s\\.'='=0.')), TEXT))


#problematic gene annotations (eg p=.049 is mistakenly annotated as gene)
annots <- annots %>% subset(!(grepl('\\=\\.|\\=\\s\\.', TEXT_ann)))

#some sections are merely noise -> remove them early on the analysis
#Revome Section Types that are mostly adding noise
annots <- annots %>% 
  subset(text_Type %in% c('abstract', 'paragraph'))


#normalize the terms
norm_ann <- termnorm(annots)
norm_ann <- norm_ann %>% subset(!is.na(TERM))


#next parts include working with the text, so a good idea is to 
#remove any unwanted white-spaces with stringr::str_squish()

norm_ann$TEXT <- pblapply(norm_ann$TEXT,  function(x) stringr::str_squish(x))
norm_ann <- norm_ann %>%  mutate_if(is.list, as.character)

#Star alleles
findstars <- function(df){
  
  
  df_int <- df %>% subset((type=='Gene')&(TEXT_ann %in% c("CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A13", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", 
                                                          "CYP2D6", "CYP2E1", "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2R1", "CYP2W1", "CYP3A4",
                                                          "CYP3A5", "CYP3A7", "CYP3A43", "CYP4F2", "NUDT15", "DPYD", "HLA-A", "HLA-B", "GSTM1", "GSTP1", 
                                                          "GSTT1", "NAT1", "NAT2", "SLC15A2", "SLC22A2", "SLCO1B1", "SLCO2B1", "TPMT", "UGT1A1",
                                                          "UGT2B7", "UGT2B15", "UGT2B17", "VKORC1")))
  
  df_rest <- df %>% subset(!((type=='Gene')&(TEXT_ann %in% c("CYP1A1", "CYP1A2", "CYP1B1", "CYP2A6", "CYP2A13", "CYP2B6", "CYP2C8", "CYP2C9", "CYP2C19", 
                                                             "CYP2D6", "CYP2E1", "CYP2F1", "CYP2J2", "CYP2R1", "CYP2S1", "CYP2R1", "CYP2W1", "CYP3A4",
                                                             "CYP3A5", "CYP3A7", "CYP3A43", "CYP4F2", "NUDT15", "DPYD", "HLA-A", "HLA-B", "GSTM1", "GSTP1", 
                                                             "GSTT1", "NAT1", "NAT2", "SLC15A2", "SLC22A2", "SLCO1B1", "SLCO2B1", "TPMT", "UGT1A1",
                                                             "UGT2B7", "UGT2B15", "UGT2B17", "VKORC1"))))
  
  fstar <- function(df){
    
    found <- str_locate_all(as.character(df['TEXT']), 
                            fixed(paste(as.character(df['TEXT_ann']), "*", sep=""))) %>% 
      as.data.frame()
    
    if(nrow(found) != 0){
      
      if(!grepl('\\(|\\)', as.character(df['TEXT_ann']))){
        
        foundt <- str_extract_all(as.character(df['TEXT']),
                                  paste(as.character(df['TEXT_ann']), "\\*[^. ]*", sep="")) %>% 
          unlist %>% unname
        
        if(nrow(found)==length(foundt)){
          
          found.res <- data.frame(found, text_found = foundt, stringsAsFactors = F)
          
          found.res <- found.res %>%
            mutate(difference = abs(start - as.numeric(df['offset']))) %>%
            filter(difference <3)
          
          if(nrow(found.res) != 0){
            
            res <- as.character(found.res['text_found'])
            
          }else{
            res <- 'No'
          }
          
        }else{
          res <- 'No'
        }
        
      }else{
        res<-'No'
      }
      
    }else{
      res <- 'No'
    }
    
    return(res)
    
    
  }
  
  df_int$isit <- apply(df_int, 1, fstar) 
  
  df_int <- df_int %>%
    mutate(type = ifelse(isit=='No', type, 'StarAllele'),
           TEXT_ann = ifelse(isit=='No', TEXT_ann, isit)) %>%
    select(-isit)
  
  df_genes <- df_int %>% subset(type=='Gene')
  
  df_sa <- df_int %>% subset(type=='StarAllele')
  
  rm(df_int)
  
  df_rest <- bind_rows(df_rest, df_genes)
  
  rm(df_genes)
  
  splitkeepfirst <- function(text){
    x <- strsplit(text, '[,);/(+&]') %>% unlist
    x <- x[1] 
    return(x)
  }
  
  df_sa$TEXT_ann <- lapply(df_sa$TEXT_ann, splitkeepfirst) %>% as.character()
  
  #remove drugs
  df_sa$TEXT_ann <- lapply(df_sa$TEXT_ann, function(x) gsub('-[a-z]+', '', x)) %>%
    as.character()
  
  
  #If last character is colon, remove it
  df_sa$TEXT_ann <- lapply(df_sa$TEXT_ann, function(x) gsub('\\:$', '', x)) %>% as.character()
  df_sa$TEXT_ann <- lapply(df_sa$TEXT_ann, function(x) gsub('\\-$', '', x)) %>% as.character()
  
  #df_sa <- df_sa %>% mutate(TERM= paste(TERM, '_SA', sep=""), TERM)
  
  res.df <- bind_rows(df_rest, df_sa)
  
  return(res.df)
  
  
}
norm_ann <- findstars(norm_ann)

#For Star Alleles -> replace id with TEXT_ann
norm_ann <- norm_ann %>% mutate(TERM = ifelse(type=='StarAllele', TEXT_ann, TERM))


#Get the sentences (from TEXT) for the annotates that are type==Gene|Mutation|Chemical|StarAllele,
#using TEXT_ann in the regex

find.sents <- function(ann, text, offset){
  
  ann <- as.character(ann)
  text <- as.character(text)
  offset <- as.integer(offset)
  
  #HGVS exceptions:c m n r p g 
  sents.text <- grep(ann , 
                     unlist(strsplit(text, 
                                     '(?<=\\.)\\s(?=[A-Z])|(?<=\\.)(?=[A-Z])|(?<=\\.)\\s(?=[\\(])|(?<=\\.)(?=[\\(])|(?<=\\.)(?=rs|hsa|mRNA|cDNA|rRNA|tRNA|shRNA|lncRNA|siRNA|mtDNA|miRNA|p|c|n|m|r|g)|(?<=\\.)\\s(?=rs|hsa|mRNA|cDNA|rRNA|tRNA|shRNA|lncRNA|siRNA|mtDNA|miRNA|p|c|n|m|r|g)|(?<=[[:alpha:]\\)]\\.)\\s(?=[0-9])|(?<=[[:alpha:]\\)]\\.)(?=[0-9])|(?<=[[:alpha:]\\)\\s]\\.)\\s(?=[0-9])|(?<=[[:alpha:]\\)\\s]\\.)(?=[0-9])', 
                                     perl=TRUE)), 
                     value=TRUE, fixed = TRUE) %>%
    unique()
  
  if(length(sents.text)>1){
    
    sent.loc <- str_locate(text, fixed(sents.text)) %>% as.data.frame()
    
    sentences <- sent.loc %>%
      mutate(sent = sents.text) 
    
    sentence <- sentences %>%
      filter( (as.integer(start)<offset) & ((as.integer(end)+1)>offset) ) %>%
      select(sent) %>%  unlist() %>% as.character()
    
  }else{
    sentence <- sents.text %>% as.character()
  }
  
  return(sentence)
}

norm_ann$sentence <- pbapply(norm_ann, 1, 
                             function(x){find.sents(x['TEXT_ann'], x['TEXT'], x['offset'])})

norm_ann$sentence <- norm_ann$sentence %>% as.character()
norm_ann$sentence[norm_ann$sentence== "character(0)"] <- NA
norm_ann <- norm_ann %>% subset(!(is.na(sentence)))


#create group based on sentence and filter those that are in the same sentence 
sentences <- norm_ann %>%
  select(pmid, pmcid, journal, authors, type, TEXT_ann, identifier, sentence, Section_Type, TERM) %>%
  subset(type != 'Gene') %>%
  group_by(pmid, sentence) %>%
  mutate(typecomb = paste(type, collapse=',')) %>% 
  subset(grepl('Chemical', typecomb)&grepl('Mutation|StarAllele', typecomb)) %>%
  select(-typecomb) %>%
  ungroup() 

#fix duplicates resulting from Variants
#unique combinations of pmid-sentence-identifier (for mutations only)
s1 <- sentences %>% subset(type=='Mutation') %>% 
  distinct(pmid, identifier, sentence, .keep_all = T)
s2 <- sentences %>% subset(type!='Mutation')
sentences <- bind_rows(s1,s2)
rm(s1,s2)

#keep some metadata for the sentences
sents_metadata <- sentences %>%
  select(pmid, pmcid, journal, authors, sentence, Section_Type) %>%
  distinct(pmid, sentence, .keep_all = T) 

#take all the possible combinations
unique.sents <- sentences %>% select(sentence) %>% unique %>% unlist 
df.list <- list()

for(i in 1:length(unique.sents)){
  
  df <- sentences %>% subset(sentence == unique.sents[[i]])
  
  chemids <- df %>% subset(type=='Chemical') %>% select(TEXT_ann, identifier, TERM) %>% unique()
  colnames(chemids) <- c('chemical', 'ChemicalID', 'NormChemical')
  varids <- df %>% subset(type == "Mutation" | type == "StarAllele") %>% select(TEXT_ann, identifier, TERM) %>% unique()
  colnames(varids) <- c('variant', 'VariantID', 'NormVariant')
  
  df.combs<- expand.grid(chemical = filter(df, type == "Chemical")$TEXT_ann,
                         variant = filter(df, type == "Mutation" | type == "StarAllele")$TEXT_ann) %>%
    mutate(sentence = as.character(df[1,'sentence']),
           pmid = as.character(df[1,'pmid'])) %>%
    mutate_if(is.factor, as.character)
  
  
  df.combs <- left_join(df.combs, chemids, by='chemical')  %>% distinct()
  df.combs <- left_join(df.combs, varids, by='variant')  %>% distinct()
  
  df.list[[i]] <- df.combs
}

finalsentences <- data.table::rbindlist(df.list) %>% 
  select(pmid, chemical, ChemicalID, NormChemical, variant, VariantID, NormVariant, sentence) 

rm(df.list, unique.sents, df, df.combs, i, chemids, varids)

#want to break the df in 2
#one with sentences with only 1 pair (1 combination of pmid-sentence)
#one with sentences with multiple pairs (N combinationS of pmid-sentence)

one.comb <- finalsentences%>%
  group_by(pmid, sentence) %>%
  filter(n() == 1) %>%
  ungroup

mult.comb <- finalsentences%>%
  group_by(pmid, sentence) %>%
  filter(n() > 1) %>%
  ungroup


#save as xlsx files
xlsx::write.xlsx2(as.data.frame(one.comb),'1pairSentences.xlsx', row.names=F)

xlsx::write.xlsx2(as.data.frame(mult.comb),'NpairSentences.xlsx', row.names=F)

#Those sentenses need to get curated as following:
# 1 -> the context of the sentence implies a definite PGx association
# 0 -> the context of the sentence does not imply a definite PGx association ->
# eg, either there is no association, or it is unclear whether a PGx association exists
