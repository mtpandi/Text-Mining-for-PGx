pubtator_anns <- function(pmids, pmcids){
  
  library(httr)
  library(jsonlite)
  library(pbapply)
  library(stringr)
  library(reutils)
  library(gsubfn)
  library(dplyr)
  library(data.table)
  
  prefix <- "https://www.ncbi.nlm.nih.gov/research/pubtator-api/publications/export/biocjson?"
  
  
  jsonf_pmid <- function(df){
    
    #split input pmids in batches of 100 and save in list x
    x <- split(df, ceiling(seq_along(pmids)/100))
    x <- pblapply(x, function(y) paste0(y, collapse = ","))
    
    #get output as json and separate objects
    json_file <- pblapply(x, function(y) httr::GET(paste0(prefix, "pmids=", y, sep = "")))
    aa <- pblapply(json_file, function(y) httr::content(y, as = "text", encoding = 'UTF-8'))
    test <- pblapply(aa, function(y) strsplit(y, "\n")[[1]]) 
    test <- test %>% unlist() %>% as.vector()
    
    #make sure the jsons are good quality
    test <- test[which((lapply(test, jsonlite::validate) %>% unlist) == TRUE)]
    test <- test[!grepl('accessions\\\"\\: \\[\\]', test)]
    test <- test[!grepl('^$', test)]
    return(test)
    
  }
  
  jsonf_pmcid <- function(df){
    
    #split in batches of 100
    x<- split(df, ceiling(seq_along(pmcids)/100))
    x <- pblapply(x, function(y) paste0(y, collapse = ","))
    
    #get the json output and separate its' objects
    json_file <- pblapply(x, function(y) httr::GET(paste0(prefix, "pmcids=", y, sep = "")))
    aa <- pblapply(json_file, function(y) httr::content(y, as = "text", encoding = 'UTF-8'))
    test <- pblapply(aa, function(y) strsplit(y, "\n")[[1]]) %>% unlist() %>% as.vector()
    
    #keep the good quality and useful jsons
    test <- test[which((lapply(test, jsonlite::validate) %>% unlist) == TRUE)]
    test <- test[!grepl('accessions\\\"\\: \\[\\]', test)]
    test <- test[!grepl('^$', test)]
    
    return(test)
  }

  
  element.exists <- function(var, element1, element2) {
    tryCatch({
      if(length(var[[element1]][[element2]][[2]]) > -1)
        return(T)
    }, error = function(e) {
      return(F)
    })
  }

  getanns <- function(x){
    if(length(x) != 0 ) {
      locations <- rbindlist(x[["locations"]], fill=T) %>% as.data.frame()
      
      annotations <- data.frame(ID = x[["infons"]]["identifier"],
                                TYPE = x[["infons"]]["type"],
                                TEXT_ann = x[["text"]], 
                                OFFSET = locations["offset"],
                                LENGTH = locations["length"],
                                stringsAsFactors = F)
      
    }else {
      annotations <- data.frame(ID = 'None',
                                TYPE = 'None',
                                TEXT_ann = 'None', 
                                OFFSET = 'None',
                                LENGTH = 'None',
                                stringsAsFactors = F)
    }
    
    return(annotations)
    
  }
  
  pmid_anns <- function(x){
    
    data <- jsonlite::fromJSON(x)
    
    #extract info
    #infons 
    text <- data.frame(TEXT = data[["passages"]][["text"]], 
                       OFFSET= as.numeric(data[["passages"]][["offset"]]),
                       TYPE = data[["passages"]][["infons"]]["type"],
                       SECTION = data[["passages"]][["infons"]]["section"],
                       stringsAsFactors = F)
    
    #annotations
    if(element.exists(data, "passages", "annotations")) {
      
      accessions<-data[['accessions']]
      accessions <- data.frame(do.call('rbind', 
                                       strsplit(as.character(accessions),'@',fixed=TRUE)),
                               stringsAsFactors = F)
      colnames(accessions) <- c('Type', 'ids')
      
      if( any(grepl('chemical', accessions))&
          any(grepl("gene|mutation", accessions))&
          any(accessions[accessions[['Type']]=='species', 'ids']=='9606')){
        
        anns<- lapply(data[["passages"]][["annotations"]], getanns)  
        
      }else{
        anns<-NULL
      }
      
      
    } else {
      
      anns <- NULL
    }
    
    if(length(anns) != 0){
      #anns is a list of dfs with the annotations, 1 for each passage but, we can't see anywhere 
      #the passage it refers
      
      for(i in 1:nrow(text)){
        anns[[i]]['TEXT'] <- text[['TEXT']][i] %>% as.character()
        anns[[i]]['text_offset'] <- text[['OFFSET']][i] %>% as.numeric()
        anns[[i]]['text_Type'] <- text[['type']][i] %>% as.character()
        anns[[i]]['Section'] <- text[['section']][i] %>% as.character()
      }
      
      pmid_res<-rbindlist(anns, fill = T) %>%
        mutate(offset = offset - text_offset, 
               pmid = ifelse( length(as.character(data[["pmid"]]))!=0, as.character(data[["pmid"]]), 'None'),
               pmcid = ifelse( length(as.character(data[["pmcid"]]))!=0, as.character(data[["pmcid"]]), 'None'),
               journal = ifelse(length(as.character(data[["passages"]][["infons"]][1,"journal"]))!=0, as.character(data[["passages"]][["infons"]][1,"journal"]), 'No data'),
               #data["journal"] %>% as.character() just the journal
               year =  ifelse(length(as.character(data[["year"]]))!=0, as.character(data[["year"]]), 'None') ,
               authors =  ifelse(length(data[["authors"]])!=0, paste(data[["authors"]], sep=",", collapse = ","), 'None'), 
               Section_Type = ifelse(Section=='Title', 'TITLE', 'ABSTRACT')) %>%
        select(c(pmid, pmcid, journal, year, authors, TEXT_ann, identifier, type, Section, offset, length, TEXT , text_Type, Section_Type)) %>% 
        as.data.frame()
      
    }else{
      pmid_res <- data.frame(pmid = 'None', 
                             pmcid = 'None', 
                             journal = 'None', 
                             year = 'None', 
                             authors = 'None', 
                             TEXT_ann = 'None', 
                             identifier = 'None', 
                             type = 'None', 
                             Section = 'None', 
                             offset = 'None', 
                             length = 'None', 
                             TEXT = 'None', stringsAsFactors = F)
    }
  
   
    return(pmid_res)
    
  }
  
  pmcid_anns <- function(x){
    
    data <- jsonlite::fromJSON(x)
    
    text <- data.frame(TEXT = data[["passages"]][["text"]], 
                       OFFSET= as.numeric(data[["passages"]][["offset"]]),
                       TYPE = data[["passages"]][["infons"]]["type"],
                       SECTION = data[["passages"]][["infons"]]["section"],
                       SECTION_TYPE = data[["passages"]][["infons"]]["section_type"],
                       stringsAsFactors = F)
    
   
    if(element.exists(data, "passages", "annotations")) {
      
      accessions<-data[['accessions']]
      accessions <- data.frame(do.call('rbind', 
                                       strsplit(as.character(accessions),'@',fixed=TRUE)),
                               stringsAsFactors = F)
      colnames(accessions) <- c('Type', 'ids')
      
      if( any(grepl('chemical', accessions))&
          any(grepl("gene|mutation", accessions))&
          any(accessions[accessions[['Type']]=='species', 'ids']=='9606')){
        
        anns<- lapply(data[["passages"]][["annotations"]], getanns)  
        
      }else{
        anns<-NULL
      }
      
    } else {
      
      anns <- NULL
    }

    
  if(length(anns)!=0){
    for(i in 1:nrow(text)){
      anns[[i]]['TEXT'] <- text[['TEXT']][i] %>% as.character()
      anns[[i]]['text_offset'] <- text[['OFFSET']][i] %>% as.numeric()
      anns[[i]]['text_Type'] <- text[['type']][i] %>% as.character()
      anns[[i]]['Section'] <- text[['section']][i] %>% as.character()
      anns[[i]]['Section_Type'] <- text[['section_type']][i] %>% as.character()
    }
    
    pmcid_res<-rbindlist(anns, fill = T) %>%
      mutate(offset = offset - text_offset, 
             pmid = ifelse( length(as.character(data[["pmid"]]))!=0, as.character(data[["pmid"]]), 'None'),
             pmcid = ifelse( length(as.character(data[["pmcid"]]))!=0, as.character(data[["pmcid"]]), 'None'),
             journal = ifelse(length(as.character(data[["passages"]][["infons"]][1,"journal"]))!=0, as.character(data[["passages"]][["infons"]][1,"journal"]), 'None'),
             #data["journal"] %>% as.character() just the journal
             year =  ifelse(length(as.character(data[["year"]]))!=0, as.character(data[["year"]]), 'None') ,
             authors =  ifelse(length(data[["authors"]])!=0, paste(data[["authors"]], sep=",", collapse = ","), 'None')) %>%
      select(c(pmid, pmcid, journal, year, authors, TEXT_ann, identifier, type, Section, offset, length, TEXT , text_Type, Section_Type)) %>% 
      as.data.frame()
    
    
  }else{
    pmcid_res <- data.frame(pmid = 'None', 
                           pmcid = 'None', 
                           journal = 'None', 
                           year = 'None', 
                           authors = 'None', 
                           TEXT_ann = 'None', 
                           identifier = 'None', 
                           type = 'None', 
                           Section = 'None', 
                           offset = 'None', 
                           length = 'None', 
                           TEXT = 'None', 
                           text_Type = 'None',
                           Section_Type = 'None', stringsAsFactors = F)
  }
      
   
    return(pmcid_res)
  }
  
  
  
  if(all(is.na(pmcids))) {
    
    message("No PMCIDs found, procceding only with PMIDs\n")
    
    message("Fetching json files")
    test <- jsonf_pmid(pmids)
    
    message("Extracting annotations")
    pmid_res <- pblapply(test,pmid_anns)
    
    pmid_res[sapply(pmid_res, is.null)] <- NULL
    
    message("Finishing up")
    
    results <- rbindlist(pmid_res, fill=T)  %>% as.data.frame()
    results[results=="None"] <- NA
    #drop unecessary annotations 
    results <- results[!((results[['Section_Type']] %in% c("COMP_INT", "REF", "ABBR", "KEYWORD", "AUTH_CONT", "ACK_FUND"))|(is.na(results[['TEXT_ann']]))|(results[['type']] %in% c("Disease","CellLine"))|(is.na(results[['type']]))),]
    

  } else {
    
    ##############pmids
    
    message("Processing both PMCIDs and PMIDs\n")
    
    message("1. PMIDs\n")
    
    message("Fetching json files")
    test_pmid <- jsonf_pmid(pmids)
    
    message("Extracting annotations")
    pmid_res <- pblapply(test_pmid, pmid_anns)
    
    pmid_res[sapply(pmid_res, is.null)] <- NULL
    
    message("Finishing up...\n")
    pmid_res <- rbindlist(pmid_res, fill = T)  %>% as.data.frame()
    
    
    
    ##############pmcids
    
    message("2. PMCIDs\n")
    
    message("Fetching json files")
    test_pmcid <- jsonf_pmcid(pmcids)
    
    message("Extracting annotations")
    pmcid_res <- pblapply(test_pmcid, pmcid_anns)
    
    pmcid_res[sapply(pmcid_res, is.null)] <- NULL
    
    message("Finishing up...\n")
    pmcid_res <- rbindlist(pmcid_res, fill = T)  %>% as.data.frame()
    
    
    ###################merge pmids-pmcids
    message("Combining PMIDs-PMCIDs...\n")
    results <- rbindlist(list(pmid_res, pmcid_res), fill=T) %>% as.data.frame()
    results[results=="None"] <- NA
    #drop unecessary annotations 
    results <- results[!((results[['Section_Type']] %in% c("COMP_INT", "REF", "ABBR", "KEYWORD", "AUTH_CONT", "ACK_FUND"))|(is.na(results[['TEXT_ann']]))|(results[['type']] %in% c("Disease","CellLine"))|(is.na(results[['type']]))),]
    
    

    #filter for species=9606 & Gene|Mutation & chemical 
    
    
  }
  
  return(results)
  
  
}