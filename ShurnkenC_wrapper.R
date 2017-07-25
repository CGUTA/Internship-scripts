library(pamr)
sample.data <- read.table("Consensus_Mix_Matrix_20170209_GPL_570_CUP_data_.txt", row.names = 1, sep="\t")
sample.info <- fread("cup_usable_samples.csv")
no_nas <- is.na(sample.info$tissue_simple) | 
  sample.info$tissue_simple %in% c("adipose",
                                   "esophagus",
                                   "eye",
                                   "gF-fallopian",
                                   "gF-vagina",
                                   "gM-testicles",
                                   "peritoneum")
no_nas <- is.na(sample.info$tissue_simple) |
  sample.info$tissue_type %in% c("cirrotic_tissue",
                                 "metaplastic_tissue",
                                 "normal_histology",
                                 "cancer")
#khan.data <- pamr.from.excel("khan.txt", 65, sample.labels=TRUE)
no_purifiable <- read.csv("non_purifiables.csv")
no_nas <- no_purifiable$x | is.na(sample.info$purity_label)
khan.data <- list(x = as.matrix(sample.data[!no_nas,]),
                  y = sample.info$tissue_simple[!no_nas],
                  #y = sample.info$gender[!no_nas],
                  genenames = NULL,
                  geneid = rownames(sample.data),
                  samplelabels = sample.info$X.1[!no_nas],
                  batchlabels = NULL
)
khan.train <- pamr.train(khan.data)
khan.results<- pamr.cv(khan.train, khan.data)
khan.results