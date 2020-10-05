library(dplyr)
library(caret)
library(text2vec)
source('sentpreprocess.R')

#load the data
one.pair <- readxl::read_xlsx('1pairSentences_curated.xlsx')  %>%
  subset(is.na(Class)) %>%
  select(-Class)

vocab <- readRDS('onepair_vocabulary.RDS')


#prepare the data using word2vec
test <- text.preproc(one.pair) %>% 
  mutate(pairid = row_number())

# #should I filter with the same words???
# filterms <- '\\bassociat|\\bmodulat|\\binfluenc|\\balter|\\baffect|\\bcorrelat|\\bsignifican|\\bpredictive|\\bpredictor|\\bincreas|\\bdecreas|\\breduc|\\binteract|\\bimpact|\\bresponse|\\btolerance|\\bresistan|marker\\b|markers\\b|\\btoxicity|\\badverse|\\bAE\\b|\\bADR\\b|\\bADRs\\b|\\bAEs\\b'
# 
# 
# test <- test %>%
#   subset(grepl(filterms, sentence, ignore.case = TRUE, perl = TRUE))



#tokenization
it_test <- itoken(test$sentence, 
                  tokenizer = word_tokenizer,
                  id = test$pairid,
                  progressbar = TRUE)

# creating vocabulary and document-term matrix
dtm_test <- create_dtm(it_test, vocab_vectorizer(vocab))
dtm <- Matrix::as.matrix(dtm_test)



# create a function to get metrics for each classifier
classifmetr <- function(classifiermodel){
  #classify the dtm, add to test df and compare with the gold standard
  test$PredictedClass <- predict(classifiermodel, 
                                 dtm, 
                                 type = 'raw')
  
  #gold standards
  gs.set <- read.csv('gs.set.csv')
  
  #filter the gold standard sets, to keep only the same pmids with the test set
  pmids <- test %>% select(pmid) %>% unique() %>% unlist() %>% unname() %>% paste0(collapse = '|')
  gs.set <- gs.set %>% subset(grepl(pmids, PMIDs))
  
  #filter again to keep only the pairs appearing both in the test set
  testpairs <- test %>% select(NormVariant, ChemicalID) %>% distinct(.keep_all = TRUE)
  gs.set <- inner_join(gs.set, testpairs, by=c("Entity1_name" = "NormVariant", "ChemicalID"))
  
  #true positives by gold standard
  gsp <- gs.set %>% subset(Class == 'corr') %>% distinct(.keep_all = TRUE) %>% select(Entity1_name, ChemicalID)
  
  #true negatives by gold standard
  gsn <- gs.set %>% subset(Class == 'nocorr') %>% distinct(.keep_all = TRUE) %>% select(Entity1_name, ChemicalID)
  
  
  #pairs predicted as positive/negative
  pred.pos <- test %>% subset(PredictedClass == 'corr') %>%
    select(NormVariant, ChemicalID) %>%
    distinct(.keep_all = TRUE)
  
  pred.neg <- test %>% subset(PredictedClass == 'nocorr') %>%
    select(NormVariant, ChemicalID) %>%
    distinct(.keep_all = TRUE)
  
  #keep the commons only as nocorr
  commons <- inner_join(pred.pos, pred.neg, by=c("NormVariant", "ChemicalID")) 
  pred.pos <- anti_join(pred.pos, commons, by=c("NormVariant", "ChemicalID"))
  
  #true positives (pred.pos found in gsp)
  tp <- inner_join(pred.pos, gsp, by=c("NormVariant"="Entity1_name", "ChemicalID"))
  
  #true negatives (pred.neg found in gsn)
  tn <- inner_join(pred.neg, gsn, by=c("NormVariant"="Entity1_name", "ChemicalID"))
  
  #false positives (pred.pos found in gsn)
  fp <- inner_join(pred.pos, gsn, by=c("NormVariant"="Entity1_name", "ChemicalID"))
  
  #false negatives (pred.neg  found in gsp)
  fn <- inner_join(pred.neg, gsp, by=c("NormVariant"="Entity1_name", "ChemicalID"))
  
  
  #metrics
  accuracy <- (nrow(tp) + nrow(tn))/(nrow(tp)+nrow(fp)+nrow(tn)+nrow(fn))
  precision <- nrow(tp)/(nrow(tp)+nrow(fp))
  specificity <- nrow(tn)/(nrow(tn)+nrow(fp))
  recall <- nrow(tp)/(nrow(tp)+nrow(fn)) #sensitivity
  npv <- nrow(tn)/(nrow(tn)+nrow(fn)) #negative predictive value = precision for the negative class
  
  res <- data.frame(Accuracy = accuracy,
                    Precision = precision, 
                    Specificity = specificity,
                    Recall = recall,
                    NPV =npv)
  
  return(res)
}

# xgboost
xgb.res <- classifmetr(readRDS('1pair_xgboosttune.RDS'))

# svm
svm.res <- classifmetr(readRDS('1pair_svmtune.RDS'))

# glmnet
glmnet.res <- classifmetr(readRDS('1pair_glmnettune.RDS'))