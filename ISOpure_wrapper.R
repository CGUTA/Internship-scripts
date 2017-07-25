library(data.table)
library(ggplot2)
sample.info <- fread("cup_usable_samples.csv")
sample.data <- fread("DatasetCarlos_GPL570_mRNA_RMA_CleanedIdentifiers_QCed.txt")
tissue.to.avoid <- sample.info$tissue_simple %in% c("adipose", #these tissues dont have more than one lab
                                                    "esophagus",
                                                    "eye",
                                                    "gF-fallopian",
                                                    "gF-vagina",
                                                    "gM-testicles",
                                                    "peritoneum"
) |
  sample.info$incidental_tissue %in% c("abdominal", #these tissues dont have healthys
                                       "bone/lung",
                                       "peritoneum",
                                       "small_intestine",
                                       "skin/parotid",
                                       "urethra",
                                       "soft_tissue",
                                       "soft_tissue/lymphoid_node/brain/adrenal",
                                       "adipose_tissue_omentum"
  ) |
  is.na(sample.info$tissue_simple) # CUPs
usable.sample.info <- sample.info[!tissue.to.avoid, with = T]
set.seed(42)
n.labs <- length( unique( usable.sample.info$lab_batch))
numbers <- 1:n.labs
#for.deconvolution <- sample(numbers, floor(n.labs/2))
#for.validation <- setdiff(numbers, for.deconvolution)
for.validation <- sample(numbers, floor(n.labs/2))
for.deconvolution <- setdiff(numbers, for.deconvolution)
usable.sample.info[, splits := .GRP, by = lab_batch]
deconvolution.sample.info <- usable.sample.info[splits %in% for.deconvolution, with=T]
sample.data <- fread("DatasetCarlos_GPL570_mRNA_RMA_CleanedIdentifiers_QCed.txt")
relevant.probes <- sample.data[, c("V1", deconvolution.sample.info$X.1), with = F ][, .(mads = apply(.SD, 1, mad), V1 = V1), .SDcols = 2:1209][order(mads, decreasing = T)][1:1000, V1]
deconvolution.sample.data <- sample.data[, c("V1", deconvolution.sample.info$X.1), with = F ][V1 %in% relevant.probes, with = T]
rm(sample.data)
####################################################
###start ISOPURE protocol
library(ISOpureR)
tumor.selector <- deconvolution.sample.info$incidental_tissue == "liver" &
  deconvolution.sample.info$tissue_simple == "colon"
normal.selector <- deconvolution.sample.info$tissue_type == "normal_histology" &
  deconvolution.sample.info$tissue_simple == "liver"

tumor.selector <- deconvolution.sample.info[which(tumor.selector)]
normal.selector <- deconvolution.sample.info[which(normal.selector)]

tumordata <- as.matrix(deconvolution.sample.data[, deconvolution.sample.data[, tumor.selector$X.1], with =F])
normaldata <- as.matrix(deconvolution.sample.data[, deconvolution.sample.data[, normal.selector$X.1], with =F])

if (ncol(normaldata) >= ncol(tumordata)){
  cutoff <- ncol(tumordata) -1
} else{
  cutoff <- ncol(normaldata)
}
counter <- 1
sapply(c(1,3,5,7), function(cutoff){
  ptm <- proc.time()
  ISOpureS1model <- ISOpure.step1.CPE(2^tumordata, 2^normaldata[,1:3])
  proc.time() - ptm
})
ptm <- proc.time()
ISOpureS1model <- ISOpure.step1.CPE(2^tumordata, 2^normaldata[,1:cutoff])
proc.time() - ptm
ptm <- proc.time()
ISOpureS2model <- ISOpure.step2.PPE(2^tumordata, 2^normaldata[,1:cutoff], ISOpureS1model)
proc.time() - ptm
purified.tumors <- t(log2(ISOpureS2model$cc_cancerprofiles))
impure.tumors <- t(tumordata)