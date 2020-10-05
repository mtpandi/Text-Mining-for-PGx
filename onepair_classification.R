library(caret)
library(text2vec)
library(dplyr)
source('sentpreprocess.R')

#load the data
one.pair <- readxl::read_xlsx('1pairSentences.xlsx')  %>% 
  na.omit(Class) 


#filter the sentences - will that improve performance?
filterms <- '\\bassociat|\\bmodulat|\\binfluenc|\\balter|\\baffect|\\bcorrelat|\\bsignifican|\\bpredictive|\\bpredictor|\\bincreas|\\bdecreas|\\breduc|\\binteract|\\bimpact|\\bresponse|\\btolerance|\\bresistan|marker\\b|markers\\b|\\btoxicity|\\badverse|\\bAE\\b|\\bADR\\b|\\bADRs\\b|\\bAEs\\b'


data <- one.pair %>%
  subset(grepl(filterms, sentence, ignore.case = TRUE, perl = TRUE))



#prepare the data using word2vec
data <- text.preproc(data) %>% as.data.frame()
data <- data %>%
  mutate(pairid = row_number(),
         Class = ifelse(Class==1, 'corr', 'nocorr'))


#tokenization
it_data <- itoken(data$sentence, 
                   tokenizer = word_tokenizer,
                   id = data$pairid,
                   progressbar = TRUE)

# creating vocabulary and document-term matrix
vocab <- create_vocabulary(it_data) %>%
  prune_vocabulary(term_count_min = 3,
                   term_count_max = 300) 

dtm_data <- create_dtm(it_data, vocab_vectorizer(vocab))
dtm <- Matrix::as.matrix(dtm_data)

saveRDS(vocab, 'onepair_vocabulary.RDS')



###############

#classification
control <- trainControl(method="cv", number=10, verboseIter = TRUE, savePredictions="final")
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


#save the models tuned
saveRDS(fit.glmnet, '1pair_glmnettune.RDS')
saveRDS(fit.xgboost, '1pair_xgboosttune.RDS')
saveRDS(fit.svm, '1pair_svmtune.RDS')


#load the fastrtext results to plot all of them
fastrtext.cv <- read.csv('onepair_fastrtext_10cv.csv')
library(tidyverse)
fastrtext.cv <- fastrtext.cv %>% remove_rownames %>% column_to_rownames(var="X")

#plots
library(ggplot2)
library(reshape2)
library(stringr)

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


#accuracies and kappa
accuracies_data <- boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(Accuracy, AccuracyLower, AccuracyUpper, Kappa, Method) %>%
  reshape2::melt(id.vars = c("Method")) 

#rest metrics
metrics_data <- boxplot_data %>%
  subset(CVFold == 'cvmetrics') %>%
  select(-c(CVFold, AccuracyLower, AccuracyUpper, Kappa)) %>%
  reshape2::melt(id.vars = c("Method")) %>% 
  mutate(variable = gsub('Sensitivity', 'Sensitivity/Recall', variable),
         variable = gsub('Precision', 'Precision/Pos.Pred.Value', variable)) %>%
  subset(variable != 'Recall') %>%
  subset(variable != 'Pos.Pred.Value')

ggplot(accuracies_data, aes(x=variable, y=value, fill = Method, label = value))+
  facet_grid(Method ~ .)+
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values=c( "#FED880", '#789DB7', "#688F4E", "#85678C"))+
  coord_flip() -> plotacc

ggplot(metrics_data, aes(x=variable, y=value, fill = Method, label = value))+
  facet_grid(Method ~ .)+
  geom_bar(stat = "identity") +
  geom_text(size = 3, position = position_stack(vjust = 0.5))+
  theme(axis.text.x=element_blank(), axis.title.x = element_blank(), axis.ticks.x = element_blank(), axis.title.y = element_blank())+
  scale_fill_manual(values=c( "#FED880", '#789DB7', "#688F4E", "#85678C"))+
  coord_flip() -> plotmetr


#save the graphs
ggsave('model_accuracies_1pair.png',
       plotacc,
       device = "png",
       dpi = 600,
       scale = 2)

ggsave('model_metrics_1pair.tiff',
       plotmetr,
       device = 'tiff',
       dpi = 600,
       scale = 2, 
       compression = "lzw")

