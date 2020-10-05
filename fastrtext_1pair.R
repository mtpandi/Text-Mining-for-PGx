library(fastrtext)
library(dplyr)
source('sentpreprocess.R')

#load the data
one.pair <- readxl::read_xlsx('1pairSentences_curated.xlsx')  %>% 
  na.omit(Class) 

#filter the sentences - will that improve performance?
filterms <- '\\bassociat|\\bmodulat|\\binfluenc|\\balter|\\baffect|\\bcorrelat|\\bsignifican|\\bpredictive|\\bpredictor|\\bincreas|\\bdecreas|\\breduc|\\binteract|\\bimpact|\\bresponse|\\btolerance|\\bresistan|marker\\b|markers\\b|\\btoxicity|\\badverse|\\bAE\\b|\\bADR\\b|\\bADRs\\b|\\bAEs\\b'


data <- one.pair %>%
  subset(grepl(filterms, sentence, ignore.case = TRUE, perl = TRUE))




#text classification
data <- text.preproc(data) %>% 
  mutate(Class = ifelse(Class==1, 'corr', 'nocorr'),
         Class = paste0("__label__", Class)) %>%
  as.data.frame()



kfoldcv <- function(k, df){
  
  metrics <- list()
  
  ids <- caret::createFolds(df$Class, k, list=FALSE)
  
  for(i in 1:k){
    
    train <- df[ids!=i,] 
    test <- df[ids==i,]
    
    #fasttext
    #fasttext
    model <- build_supervised(documents = train$sentence,
                              targets = train$Class,
                              model_path = 'my_model',
                              loss = "ns",
                              label = "__label__",
                              dim = 200, #size of vector
                              minCount = 5,
                              ws = 2, #window size
                              lr = 0.1, 
                              epoch = 50, 
                              thread = 10,
                              wordNgrams = 2)
    
    model_call <- load_model(model)
    
    
    predictions <- predict(model_call, test$sentence)
    
    pred_df <- data.frame(pred = as.factor(names(unlist(predictions))), true = as.factor(test$Class))
    
    metrics[[i]] <- caret::confusionMatrix(pred_df$pred, pred_df$true, positive='__label__corr')
    
    
    
    
  }
  
  cvdf <- list()
  for(i in 1:length(metrics)){
    cvdf[[i]] <- data.frame("Accuracy" = metrics[[i]][['overall']][['Accuracy']] %>% as.numeric(),
                            "Kappa" = metrics[[i]][['overall']][['Kappa']] %>% as.numeric(),
                            "AccuracyLower" = metrics[[i]][['overall']][['AccuracyLower']] %>% as.numeric(),
                            "AccuracyUpper" = metrics[[i]][['overall']][['AccuracyUpper']] %>% as.numeric(),
                            "Sensitivity" = metrics[[i]][['byClass']][['Sensitivity']] %>% as.numeric(),
                            "Specificity" = metrics[[i]][['byClass']][['Specificity']] %>% as.numeric(),
                            "Pos Pred Value" = metrics[[i]][['byClass']][['Pos Pred Value']] %>% as.numeric(),
                            "Neg Pred Value" = metrics[[i]][['byClass']][['Neg Pred Value']] %>% as.numeric(),
                            "Precision" = metrics[[i]][['byClass']][['Precision']] %>% as.numeric(),
                            "Recall" = metrics[[i]][['byClass']][['Recall']] %>% as.numeric(),
                            "F1" = metrics[[i]][['byClass']][['F1']] %>% as.numeric(),
                            "Prevalence" = metrics[[i]][['byClass']][['Prevalence']] %>% as.numeric(),
                            "Detection Rate" = metrics[[i]][['byClass']][['Detection Rate']] %>% as.numeric(),
                            "Detection Prevalence" = metrics[[i]][['byClass']][['Detection Prevalence']] %>% as.numeric(),
                            "Balanced Accuracy" = metrics[[i]][['byClass']][['Balanced Accuracy']] %>% as.numeric()
    )
  }
  
  cvdf <- as.data.frame(t(bind_rows(cvdf)))
  cvdf <- cvdf %>%
    tibble::rownames_to_column('metric') %>%
    mutate(cvmetrics = rowMeans(select(cvdf, starts_with("V")), na.rm = T)) %>%
    tibble::column_to_rownames('metric') %>%
    as.data.frame()
  
  return(cvdf)
  
}


classifres <- kfoldcv(10, data)


write.csv(classifres, 'onepair_fastrtext_10cv.csv', row.names = TRUE)


#fasttext
set.seed(17071994)
model <- build_supervised(documents = data$sentence,
                          targets = data$Class,
                          model_path = 'my_model',
                          loss = "ns",
                          label = "__label__",
                          dim = 200, #size of vector
                          minCount = 5,
                          ws = 2, #window size
                          lr = 0.1, 
                          epoch = 50, 
                          thread = 10,
                          wordNgrams = 2)

model_call <- load_model(model)



#load test data 
test <- readxl::read_xlsx('1pairSentences_curated.xlsx')  %>%
  subset(is.na(Class)) %>%
  select(-Class)

#prepare the data 
test <- text.preproc(test) 

test <- test %>%
  subset(grepl(filterms, sentence, ignore.case = TRUE, perl = TRUE))




test$PredictedClass <- predict(model_call, sentences = test$sentence) %>% unlist() %>% names


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
pred.pos <- test %>% subset(PredictedClass == '__label__corr') %>%
  select(NormVariant, ChemicalID) %>%
  distinct(.keep_all = TRUE)

pred.neg <- test %>% subset(PredictedClass == '__label__nocorr') %>%
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


