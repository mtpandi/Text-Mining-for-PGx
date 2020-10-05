#create a custom list of stopwords
library(dplyr)

pmed <- c("a", "about", "again", "all", "almost", "also", "always", "among", "an", "and", 
          "another", "any", "are", "as", "at", "be", "because", "been", "before", "being", "between", 
          "both", "but", "by", "can", "could", "did", "do", "does", "done", "due", "during", 
          "each", "either", "enough", "especially", "etc", "for", "from", "further", "had", 
          "has", "have", "having", "here", "how", "however", "i", "if", "in", "into", "is", "it", "its", 
          "itself", "just", "kg", "km", "made", "mainly", "make", "may", "mg", "might", "ml", "mm", "most", 
          "mostly", "must", "nearly",  "obtained", "of", "often", "on", "our", 
          "overall", "perhaps", "pmid", "pmcid", "quite", "rather", "really", "regarding", "seem", "seen", 
          "several", "should", "since", "so", "some", "such", "than", 
          "that", "the", "their", "theirs", "them", "then", "there", "therefore", "these", "they", "this", 
          "those", "through", "thus", "to", "upon", "use", "used", "using", "various", "very", "was", "we", 
          "were", "what", "when", "which", "while", "with", "within", "without", "would")

extras <-    c("patient","patients", "polymorphism","polymorphisms",  "genotype", "genotypes", "study", "studies",
               "allele", "alleles", "dose", "doses","gene", "genes", "variant", "variants" , "carrier", "carriers",
               "table", "tables", "analysis", "analyses",  "single", "homozygous", "cohort",
               "cohorts", "figures", "nucleotides", "populations", "types", "individuals", "snps", "africans",
               "figure", "common", "figs", "nucleotide", "baseline", "population", "americans",
               "type", "individual", "plasma", "snp",  "african", "american", "therapy", "therapies",
               "wide", "glucose", "chinese", "subject", "subjects", "wild", "daily", "transporter",
               "transporters", "europeans", "blacks", "metabolites", "participants", "chromosomes",
               "european", "heterozygous", "blood", "cotinine", "group","groups", "intronic", "platelet",
               "systolic", "black", "metabolite", "participant", "chromosome", "diastolic", "missense",
               "panels", "proteins", "receptors", "whites", "promoters", "smokers", "enzymes", "indians",
               "oral", "panel",  "serum", "south", "ancestry", "anesthesia", "opioid", "protein", 
               "receptor", "white", "coding", "homozygote","homozygotes", "administration", "myopathy", "ototoxicity", 
               "promoter", "smoker", "supplementary", "upstream", "vitro", "vivo", "silico", "enzyme",
               "months", "recipients", "cells", "structures", "haplotypes", "conformations", "blockers",
               "indian", "intergenic", "japanese", "fold", "metabolizer","metabolizers", "meta", "muscle", "month", "players",
               "recipient", "receiving", "renal", "sample", "smoking", "three", "trial", "albumin", 
               "ergogenic", "encoding","children", "cell", "structure", "metastatic", "gwas", "drinking",
               "estrogen", "ethnic", "caucasian","catalytic", "cardiovascular", "background", "alpha", "beta",
               "placebo", "liver", "haplotype", "diabetes", "cycling", "consumption", "conformation", "mutants",
               "blocker", "transplant", "hypertension", "binding", "player", "gender", "channel", "bound",
               "weight",  "versus", "basketball", "cancer", "chemotherapy", "mutant", "mutation","mutations", "clinical")

pmed <- append(pmed, extras)
# pmed <- append(pmed, 
#                stopwords::stopwords(source= 'snowball') %>% as.data.frame() %>%
#                  subset(!(`.` %in% pmed)) %>%
#                  subset(!grepl('[[:punct:]]', `.`)) %>%
#                  filter(`.` != 'cannot') %>%
#                  unlist %>%
#                  as.character())


#gene names
genes <- read.csv('genes.csv', header = T, sep=',') 
names(genes) <- 'symbol'


genes$symbol <- lapply(genes$symbol, function(x) gsub('[[:punct:]0-9]', ' ', x))
genes$symbol <- as.character(genes$symbol)

genes$symbol  <- lapply(genes$symbol, function(x) strsplit(x, '\\s')[[1]][1])
genes$symbol <- as.character(genes$symbol)

genes <- genes %>%
  mutate(symbol = tolower(symbol)) %>%
  distinct() %>%
  subset(nchar(symbol) > 3) %>%
  unlist() %>%
  as.character()

swords <- append(pmed, genes)
write.csv(swords, 'customstopwords.csv', row.names = F)