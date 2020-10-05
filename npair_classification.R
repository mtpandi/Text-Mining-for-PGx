library(caret)
library(text2vec)
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
                           TRUE ~ 'both'))


#filter the sentences - will that improve performance?
filterms <- '\\bassociat|\\bmodulat|\\binfluenc|\\balter|\\baffect|\\bcorrelat|\\bsignifican|\\bpredictive|\\bpredictor|\\bincreas|\\bdecreas|\\breduc|\\binteract|\\bimpact|\\bresponse|\\btolerance|\\bresistan|marker\\b|markers\\b|\\btoxicity|\\badverse|\\bAE\\b|\\bADR\\b|\\bADRs\\b|\\bAEs\\b'


data <- data %>%
  subset(grepl(filterms, sentence, ignore.case = TRUE, perl = TRUE)) 



#prepare the data using word2vec
data <- data %>%
  mutate(pairid = row_number())  %>%
  as.data.frame()


#tokenization
it_data <- itoken(data$sentence, 
                  tokenizer = word_tokenizer,
                  id = data$pairid,
                  progressbar = TRUE)

# creating vocabulary and document-term matrix
vocab <- create_vocabulary(it_data) %>%
  prune_vocabulary(term_count_min = 3,
                   term_count_max = 400) 

dtm_data <- create_dtm(it_data, vocab_vectorizer(vocab))
dtm <- Matrix::as.matrix(dtm_data)

saveRDS(vocab, 'npair_vocabulary.RDS')
###############

#classification
control <- trainControl(method="cv", number=5, verboseIter = TRUE, savePredictions="final")
seed <- 170794

metric <- "Accuracy"

getmetrics <- function(df){
  
  metrics <- df %>%
    select(pred, obs, Resample) %>%
    group_by(Resample) %>%
    group_map(~ caret::confusionMatrix(.x$pred, .x$obs, positive='corr'))
  
  cvdf <- list()
  for(i in 1:length(metrics)){
    cvdf[[i]] <- data.frame("Accuracy" = metrics[[i]][['overall']][['Accuracy']] %>% as.numeric(),
                            "Kappa" = metrics[[i]][['overall']][['Kappa']] %>% as.numeric(),
                            "AccuracyLower" = metrics[[i]][['overall']][['AccuracyLower']] %>% as.numeric(),
                            "AccuracyUpper" = metrics[[i]][['overall']][['AccuracyUpper']] %>% as.numeric(),
                            "corr_Sensitivity" = metrics[[i]][['byClass']][['Class: corr','Sensitivity']] %>% as.numeric(),
                            "corr_Specificity" = metrics[[i]][['byClass']][['Class: corr','Specificity']] %>% as.numeric(),
                            "corr_Pos Pred Value" = metrics[[i]][['byClass']][['Class: corr','Pos Pred Value']] %>% as.numeric(),
                            "corr_Neg Pred Value" = metrics[[i]][['byClass']][['Class: corr','Neg Pred Value']] %>% as.numeric(),
                            "corr_Precision" = metrics[[i]][['byClass']][['Class: corr','Precision']] %>% as.numeric(),
                            "corr_Recall" = metrics[[i]][['byClass']][['Class: corr','Recall']] %>% as.numeric(),
                            "corr_F1" = metrics[[i]][['byClass']][['Class: corr','F1']] %>% as.numeric(),
                            "corr_Prevalence" = metrics[[i]][['byClass']][['Class: corr','Prevalence']] %>% as.numeric(),
                            "corr_Detection Rate" = metrics[[i]][['byClass']][['Class: corr','Detection Rate']] %>% as.numeric(),
                            "corr_Detection Prevalence" = metrics[[i]][['byClass']][['Class: corr','Detection Prevalence']] %>% as.numeric(),
                            "corr_Balanced Accuracy" = metrics[[i]][['byClass']][['Class: corr', 'Balanced Accuracy']] %>% as.numeric(),
                            "nocorr_Sensitivity" = metrics[[i]][['byClass']][['Class: nocorr','Sensitivity']] %>% as.numeric(),
                            "nocorr_Specificity" = metrics[[i]][['byClass']][['Class: nocorr','Specificity']] %>% as.numeric(),
                            "nocorr_Pos Pred Value" = metrics[[i]][['byClass']][['Class: nocorr','Pos Pred Value']] %>% as.numeric(),
                            "nocorr_Neg Pred Value" = metrics[[i]][['byClass']][['Class: nocorr','Neg Pred Value']] %>% as.numeric(),
                            "nocorr_Precision" = metrics[[i]][['byClass']][['Class: nocorr','Precision']] %>% as.numeric(),
                            "nocorr_Recall" = metrics[[i]][['byClass']][['Class: nocorr','Recall']] %>% as.numeric(),
                            "nocorr_F1" = metrics[[i]][['byClass']][['Class: nocorr','F1']] %>% as.numeric(),
                            "nocorr_Prevalence" = metrics[[i]][['byClass']][['Class: nocorr','Prevalence']] %>% as.numeric(),
                            "nocorr_Detection Rate" = metrics[[i]][['byClass']][['Class: nocorr','Detection Rate']] %>% as.numeric(),
                            "nocorr_Detection Prevalence" = metrics[[i]][['byClass']][['Class: nocorr','Detection Prevalence']] %>% as.numeric(),
                            "nocorr_Balanced Accuracy" = metrics[[i]][['byClass']][['Class: nocorr','Balanced Accuracy']] %>% as.numeric(),
                            "both_Sensitivity" = metrics[[i]][['byClass']][['Class: both','Sensitivity']] %>% as.numeric(),
                            "both_Specificity" = metrics[[i]][['byClass']][['Class: both','Specificity']] %>% as.numeric(),
                            "both_Pos Pred Value" = metrics[[i]][['byClass']][['Class: both','Pos Pred Value']] %>% as.numeric(),
                            "both_Neg Pred Value" = metrics[[i]][['byClass']][['Class: both','Neg Pred Value']] %>% as.numeric(),
                            "both_Precision" = metrics[[i]][['byClass']][['Class: both','Precision']] %>% as.numeric(),
                            "both_Recall" = metrics[[i]][['byClass']][['Class: both','Recall']] %>% as.numeric(),
                            "both_F1" = metrics[[i]][['byClass']][['Class: both','F1']] %>% as.numeric(),
                            "both_Prevalence" = metrics[[i]][['byClass']][['Class: both','Prevalence']] %>% as.numeric(),
                            "both_Detection Rate" = metrics[[i]][['byClass']][['Class: both','Detection Rate']] %>% as.numeric(),
                            "both_Detection Prevalence" = metrics[[i]][['byClass']][['Class: both','Detection Prevalence']] %>% as.numeric(),
                            "both_Balanced Accuracy" = metrics[[i]][['byClass']][['Class: both','Balanced Accuracy']] %>% as.numeric()
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

# SVM Linear
set.seed(seed)
fit.svm <- train(dtm, y=as.factor(data$Class), method="svmLinear", metric=metric,  trControl=control)

svm.cv <- getmetrics(fit.svm$pred)


# xgboost
set.seed(seed)
xgboost.tune <- expand.grid(nrounds = c(50, 100, 150, 200),
                            max_depth = seq(4,12,2),
                            eta = seq(0.1, 0.5, 0.1),
                            gamma = 0,
                            colsample_bytree = 1,
                            min_child_weight = 1,
                            subsample = c(0.7, 0.75, 1))

fit.xgboost <- train(dtm, y=as.factor(data$Class), method="xgbTree", 
                     metric=metric, 
                     tuneGrid = xgboost.tune ,
                     trControl=control)

xgboost.cv <- getmetrics(fit.xgboost$pred)

#glmnet
set.seed(seed)
fit.glmnet <- train(dtm, y=as.factor(data$Class), method="glmnet", 
                    metric=metric , trControl=control)

glmnet.cv <- getmetrics(fit.glmnet$pred)

#save the models

saveRDS(fit.glmnet, 'npairs_glmnettune.RDS')
saveRDS(fit.xgboost, 'npairs_xgboosttune.RDS')
saveRDS(fit.svm, 'npairs_svmtune.RDS')


#graphs

#load fastertext cv resultls
fastrtext.cv <- read.csv('fastrtext_npairs5cv.csv')
library(tidyverse)
fastrtext.cv <- fastrtext.cv %>% remove_rownames %>% column_to_rownames(var="X")


library(ggplot2)
library(reshape2)


as.data.frame(t(svm.cv)) %>%
  tibble::rownames_to_column( "CVFold") %>%
  mutate(CVFold = str_replace_all(CVFold, '^V', 'Fold0'), 
         Method= 'LinearSVM') -> aa

as.data.frame(t(glmnet.cv)) %>%
  tibble::rownames_to_column( "CVFold") %>%
  mutate(CVFold = str_replace_all(CVFold, '^V', 'Fold0'),
         Method = 'glmnet') -> bb

as.data.frame(t(xgboost.cv)) %>%
  tibble::rownames_to_column( "CVFold") %>%
  mutate(CVFold = str_replace_all(CVFold, '^V', 'Fold0'),
         Method = 'xgboost') -> cc

as.data.frame(t(fastrtext.cv)) %>%
  tibble::rownames_to_column( "CVFold") %>%
  mutate(CVFold = str_replace_all(CVFold, '^V', 'Fold0'),
         Method = 'fastText') -> dd

boxplot_data <- rbind(aa,bb,cc, dd) %>%
  mutate_if(is.numeric, ~round(.,2))



# 4 graphs
# 1 with the accuracy and kappa 
# 3 with the metrics by class
accuracydata <-  boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(Accuracy, AccuracyLower, AccuracyUpper, Kappa, Method) %>%
  reshape2::melt(id.vars = c("Method"))

corrdata <- boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(starts_with('corr'), Method) %>%
  reshape2::melt(id.vars = c("Method")) %>%
  mutate(variable = gsub('^corr_', '', variable),
         variable = gsub('Sensitivity', 'Sensitivity/Recall', variable),
         variable = gsub('Precision', 'Precision/Pos.Pred.Value', variable)) %>%
  subset(variable != 'Recall') %>%
  subset(variable != 'Pos.Pred.Value')

nocorrdata <- boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(starts_with('nocorr'), Method) %>%
  reshape2::melt(id.vars = c("Method")) %>%
  mutate(variable = gsub('^nocorr_', '', variable),
         variable = gsub('Sensitivity', 'Sensitivity/Recall', variable),
         variable = gsub('Precision', 'Precision/Pos.Pred.Value', variable)) %>%
  subset(variable != 'Recall') %>%
  subset(variable != 'Pos.Pred.Value')

bothdata <- boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(starts_with('both'), Method) %>%
  reshape2::melt(id.vars = c("Method")) %>%
  mutate(variable = gsub('^both_', '', variable),
         variable = gsub('Sensitivity', 'Sensitivity/Recall', variable),
         variable = gsub('Precision', 'Precision/Pos.Pred.Value', variable)) %>%
  subset(variable != 'Recall') %>%
  subset(variable != 'Pos.Pred.Value')


# total accuracies graph
ggplot(accuracydata, aes(x=variable, y=value, fill = Method, label = value))+
  facet_grid(Method ~ .)+
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values=c( "#FED880", '#789DB7', "#688F4E", "#85678C"))+
  coord_flip() -> plotcv

ggsave('model_accuracies_npairs.png',
       plotcv,
       device = "png",
       dpi = 600,
       scale = 2)



#also a graph with the total classifier metrics
boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(starts_with('corr'), starts_with('nocorr'), starts_with('both'), Method)  %>%
  reshape2::melt(id.vars = c("Method")) %>%
  mutate(Class = case_when(grepl('^corr\\_', variable) ~ 'Correlated',
                           grepl('^nocorr\\_', variable) ~ 'noCorrelated',
                           TRUE ~ 'Both'),
         variable = gsub('^corr_|^nocorr_|^both_', '', variable),
         variable = gsub('Sensitivity', 'Sensitivity/Recall', variable),
         variable = gsub('Precision', 'Precision/Pos.Pred.Value', variable)) %>%
  subset(variable != 'Recall') %>%
  subset(variable != 'Pos.Pred.Value') %>%
  spread(Class, value) %>%
  mutate(Total = Both*0.15 + Correlated*0.50+ noCorrelated*0.35) %>%
  subset(variable %in% c('Balanced.Accuracy', 'F1', 'Precision/Pos.Pred.Value', 
                         'Sensitivity/Recall', 'Specificity')) %>%
  rename(Metrics=variable) %>%
  reshape2::melt(id.vars = c('Method', 'Metrics')) -> total

ggplot(total, aes(x=Metrics, y=value, fill = Method, label = value))+
  facet_grid(variable ~ .)+
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values=c( "#FED880", '#789DB7', "#688F4E", "#85678C"))+
  coord_flip() -> npairplot

ggsave('metrics_npairs.tiff',
       npairplot,
       device = 'tiff',
       dpi = 600,
       width=15, height=5, 
       compression = "lzw")