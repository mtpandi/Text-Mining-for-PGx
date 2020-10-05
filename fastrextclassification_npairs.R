library(fastrtext)
library(dplyr)
source('npairsentpreprocess.R')

#load the data
n.pairs <- readxl::read_xlsx('NpairSentences_curated.xlsx') %>%
  na.omit(Class) 

#preprocessing
#1. replace the entities of interest
sentence.list <- n.pairs %>% distinct(sentence) %>% unlist %>% unname

entity_replace <- function(df.sent){
  
  variant.list <- df.sent %>% select(variant) %>% unique() %>% unlist %>% unname
  chemical.list <- df.sent %>% select(chemical) %>% unique() %>% unlist %>% unname
  
  sent <- df.sent$sentence[1]
  for(i in 1:length(variant.list)){
    sent <- stringr::str_replace_all(sent, stringr::fixed(variant.list[i]), 'ClassVariant')
  }
  
  for(i in 1:length(chemical.list)){
    sent <- stringr::str_replace_all(sent, stringr::fixed(chemical.list[i]), 'ClassChemical')
  }
  
  sent <- gsub('\\b\\d+\\s', '', sent)
  
  df.sent <- df.sent %>% 
    mutate(sentence = sent)
  return(df.sent)
}

sent.dfs.list <- list()
for(i in 1:length(sentence.list)){
  inp.df <- n.pairs %>% subset(sentence == sentence.list[i])
  sent.dfs.list[[i]] <- entity_replace(inp.df)
}

data <- bind_rows(sent.dfs.list)
rm(i, inp.df, sent.dfs.list, sentence.list)

data <- text.preproc(data) %>% as.data.frame()

data <- data %>%
  group_by(sentence) %>%
  summarise(Class = mean(Class)) %>%
  mutate(Class = case_when(Class ==0 ~ 'nocorr',
                           Class == 1 ~'corr', 
                           TRUE ~ 'both'),
         Class = paste0("__label__", Class))


#filter the sentences - will that improve performance?
filterms <- '\\bassociat|\\bmodulat|\\binfluenc|\\balter|\\baffect|\\bcorrelat|\\bsignifican|\\bpredictive|\\bpredictor|\\bincreas|\\bdecreas|\\breduc|\\binteract|\\bimpact|\\bresponse|\\btolerance|\\bresistan|marker\\b|markers\\b|\\btoxicity|\\badverse|\\bAE\\b|\\bADR\\b|\\bADRs\\b|\\bAEs\\b'


data <- data %>%
  subset(grepl(filterms, sentence, ignore.case = TRUE, perl = TRUE)) 



#text classification


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
                            "corr_Sensitivity" = metrics[[i]][['byClass']][['Class: __label__corr','Sensitivity']] %>% as.numeric(),
                            "corr_Specificity" = metrics[[i]][['byClass']][['Class: __label__corr','Specificity']] %>% as.numeric(),
                            "corr_Pos Pred Value" = metrics[[i]][['byClass']][['Class: __label__corr','Pos Pred Value']] %>% as.numeric(),
                            "corr_Neg Pred Value" = metrics[[i]][['byClass']][['Class: __label__corr','Neg Pred Value']] %>% as.numeric(),
                            "corr_Precision" = metrics[[i]][['byClass']][['Class: __label__corr','Precision']] %>% as.numeric(),
                            "corr_Recall" = metrics[[i]][['byClass']][['Class: __label__corr','Recall']] %>% as.numeric(),
                            "corr_F1" = metrics[[i]][['byClass']][['Class: __label__corr','F1']] %>% as.numeric(),
                            "corr_Prevalence" = metrics[[i]][['byClass']][['Class: __label__corr','Prevalence']] %>% as.numeric(),
                            "corr_Detection Rate" = metrics[[i]][['byClass']][['Class: __label__corr','Detection Rate']] %>% as.numeric(),
                            "corr_Detection Prevalence" = metrics[[i]][['byClass']][['Class: __label__corr','Detection Prevalence']] %>% as.numeric(),
                            "corr_Balanced Accuracy" = metrics[[i]][['byClass']][['Class: __label__corr', 'Balanced Accuracy']] %>% as.numeric(),
                            "nocorr_Sensitivity" = metrics[[i]][['byClass']][['Class: __label__nocorr','Sensitivity']] %>% as.numeric(),
                            "nocorr_Specificity" = metrics[[i]][['byClass']][['Class: __label__nocorr','Specificity']] %>% as.numeric(),
                            "nocorr_Pos Pred Value" = metrics[[i]][['byClass']][['Class: __label__nocorr','Pos Pred Value']] %>% as.numeric(),
                            "nocorr_Neg Pred Value" = metrics[[i]][['byClass']][['Class: __label__nocorr','Neg Pred Value']] %>% as.numeric(),
                            "nocorr_Precision" = metrics[[i]][['byClass']][['Class: __label__nocorr','Precision']] %>% as.numeric(),
                            "nocorr_Recall" = metrics[[i]][['byClass']][['Class: __label__nocorr','Recall']] %>% as.numeric(),
                            "nocorr_F1" = metrics[[i]][['byClass']][['Class: __label__nocorr','F1']] %>% as.numeric(),
                            "nocorr_Prevalence" = metrics[[i]][['byClass']][['Class: __label__nocorr','Prevalence']] %>% as.numeric(),
                            "nocorr_Detection Rate" = metrics[[i]][['byClass']][['Class: __label__nocorr','Detection Rate']] %>% as.numeric(),
                            "nocorr_Detection Prevalence" = metrics[[i]][['byClass']][['Class: __label__nocorr','Detection Prevalence']] %>% as.numeric(),
                            "nocorr_Balanced Accuracy" = metrics[[i]][['byClass']][['Class: __label__nocorr','Balanced Accuracy']] %>% as.numeric(),
                            "both_Sensitivity" = metrics[[i]][['byClass']][['Class: __label__both','Sensitivity']] %>% as.numeric(),
                            "both_Specificity" = metrics[[i]][['byClass']][['Class: __label__both','Specificity']] %>% as.numeric(),
                            "both_Pos Pred Value" = metrics[[i]][['byClass']][['Class: __label__both','Pos Pred Value']] %>% as.numeric(),
                            "both_Neg Pred Value" = metrics[[i]][['byClass']][['Class: __label__both','Neg Pred Value']] %>% as.numeric(),
                            "both_Precision" = metrics[[i]][['byClass']][['Class: __label__both','Precision']] %>% as.numeric(),
                            "both_Recall" = metrics[[i]][['byClass']][['Class: __label__both','Recall']] %>% as.numeric(),
                            "both_F1" = metrics[[i]][['byClass']][['Class: __label__both','F1']] %>% as.numeric(),
                            "both_Prevalence" = metrics[[i]][['byClass']][['Class: __label__both','Prevalence']] %>% as.numeric(),
                            "both_Detection Rate" = metrics[[i]][['byClass']][['Class: __label__both','Detection Rate']] %>% as.numeric(),
                            "both_Detection Prevalence" = metrics[[i]][['byClass']][['Class: __label__both','Detection Prevalence']] %>% as.numeric(),
                            "both_Balanced Accuracy" = metrics[[i]][['byClass']][['Class: __label__both','Balanced Accuracy']] %>% as.numeric()
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


classifres <- kfoldcv(5, data)

write.csv(classifres, 'fastrtext_npairs5cv.csv', row.names = TRUE)



