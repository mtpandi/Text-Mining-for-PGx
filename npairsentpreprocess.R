text.preproc <- function(df){
  
  library(dplyr)
  library(stringr)
  library(pbapply)


  
  # #1. replace the entities of interest
  # entity_replace <- function(sentence, variant, chemical){
  #   sentence <- stringr::str_replace_all(sentence, stringr::fixed(variant), 'ClassVariant')
  #   sentence <- stringr::str_replace_all(sentence, stringr::fixed(chemical), 'ClassChemical')
  #   sentence <- gsub('\\b\\d+\\s', '', sentence)
  #   return(sentence)
  # }
  # 
  # df$sentence <- pbmapply(entity_replace, df$sentence, df$variant , df$chemical)
  
  #2. convert to lower
  df$sentence <- pbsapply(df$sentence, tolower)
  
  #3. remove urls
  df <- df %>%
    mutate(sentence = ifelse(grepl('(http[^ ]*)|(www\\_[^ ]*)', sentence), 
                             stringr::str_replace_all(sentence, '(http[^ ]*)|(www\\_[^ ]*)', ''), 
                             sentence))
  
  #4. remove numbers and punctuation
  df$sentence <- lapply(df$sentence, function(x) gsub('[[:punct:]0-9]', ' ', x))
  df$sentence <- as.character(df$sentence)
  
  #5. remove words with 3 characters or less
  removesmall <- function(x){
    x <- as.character(x)
    x <- strsplit(x, ' ') %>% unlist()
    x <- ifelse(!grepl('\\bno\\b|\\bnot\\b|\\bnor\\b|\\badr\\b|\\bae\\b|\\baes\\b', x), gsub('\\b\\w{1,3}\\b','', x),x )
    x <- paste0(x, collapse = ' ')
    return(x)
  }
  
  df$sentence <- lapply(df$sentence, removesmall)
  df$sentence <- as.character(df$sentence)
  #6. remove stopwords
  
  stop <- read.csv('customstopwords.csv', stringsAsFactors = FALSE) %>% unlist %>% unname %>% as.character()
  
  df$sentence <- qdap::rm_stopwords(df$sentence, stopwords = stop, separate=FALSE)
  
  
  #make sure no unwanted spaces exist
  df$sentence <- pblapply(df$sentence,  function(x) stringr::str_squish(x))
  df <- df %>%  mutate_if(is.list, as.character)
  
  
  return(df)
}
