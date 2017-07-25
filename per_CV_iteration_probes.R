###Scripts for doing the variance subsets taking labs into account
library(parallel)
library(class)
library(plyr)
raw_data <- read.table("Desktop/DatasetCarlos_GPL570_mRNA_RMA_CleanedIdentifiers_QCed.txt", sep="\t", row.names = 1, header = T)
label_data <- read.csv2("cup_usable_samples.csv", row.names = 1)
cross_validation_groups <- label_data$lab_batch
probes <- sapply( levels(cross_validation_groups), function(holdup_group) {
  print(holdup_group)
  training_selector <- cross_validation_groups != holdup_group
  return( names(sort(apply(raw_data[training_selector,],2,mad), decreasing = T)) ) # for PCAs scores
  #return( names(sort(apply(raw_data[,training_selector],1,mad), decreasing = T)[1:5000]) ) # for probes
})

## Output file
probes <- read.csv2("probes.csv", row.names = 1)