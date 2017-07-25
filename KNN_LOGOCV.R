library(data.table)
library(ggplot2)

###FUN#########################################################################
bad.tissues <- function() {
  sample.table <- fread("cup_usable_samples.csv")
  tissue.to.avoid <- sample.table$tissue_simple %in% c("adipose", #these tissues dont have more than one lab
                                                       "esophagus",
                                                       "eye",
                                                       "gF-fallopian",
                                                       "gF-vagina",
                                                       "gM-testicles",
                                                       "peritoneum"
  ) |
    sample.table$incidental_tissue %in% c("abdominal", #these tissues dont have healthys
                                          "bone/lung",
                                          "peritoneum",
                                          "small_intestine",
                                          "skin/parotid",
                                          "urethra",
                                          "soft_tissue",
                                          "soft_tissue/lymphoid_node/brain/adrenal",
                                          "adipose_tissue_omentum"
    ) |
    is.na(sample.table$tissue_simple) # CUPs
  return(tissue.to.avoid)
}
non.purifiable.tissues <- function() {
  sample.table <- fread("cup_usable_samples.csv")
  tissue.to.avoid <- bad.tissues()
  non.purifiable <- sample.table$tissue_simple %in% c("bone", #these tisues dont have enough reliable helthys
                                                      "galbladder",
                                                      "gF-vulva",
                                                      "gM-penis",
                                                      "gM-prostate",
                                                      "gM-testicles",
                                                      "heart",
                                                      "nasopharingleal",
                                                      "oral", #issue with ISOpure
                                                      "pancreas"
  ) |
    sample.table$lab_batch %in% c("spiller", # these labs have batch effects that place them far away from the majority
                                  "milind",
                                  "okayama",
                                  "pedersen",
                                  "gudjoson",
                                  "yaoY",
                                  "stephan"
    ) |
    sample.table$tissue_type %in% c("metastatic", # these have to be purified in a different manner
                                    "cirrotic_tissue"
    ) | bad.tissues()
  sample.table[!non.purifiable & tissue_type=="cancer", n_strength:= .N, by= .(adjacent_tumor,tissue_simple), with = T]
  few.tissues <- (sample.table$n_strength <=9 & sample.table$tissue_type=="cancer")
  final.selector <- non.purifiable | few.tissues
  return(final.selector)
}
load.data <- function(dir,tissue.to.avoid) {
  sample.table <- fread(dir)
  return(sample.table[!tissue.to.avoid, with = T])
}

run.knn.curve <- function(labels,cross.validation.groups,row.cutoff) {
  curve_error <- sapply(1:10, function(i) {
    no_cores <- 3
    cluster <- makeCluster(no_cores, type="FORK")
    ptm <- proc.time()
    tally <- parSapply(cluster, levels(cross.validation.groups), function(holdup_group) {
      point_data <- sample.data[,filters[1:row.cutoff,get(holdup_group)], with=F] ### here ends editing
      training_selector <- cross.validation.groups != holdup_group
      
      prediction.vector <- run.knn(
        labelst <- labels[ training_selector ],
        labelsv <- labels[ !training_selector ],
        train.points <- point_data[ training_selector, ],
        test.points <- point_data[ !training_selector, ],
        step = i*2 + 1
      )
      
      names(prediction.vector) <- row.names(point_data[ !training_selector, ])
      return(prediction.vector)
    }, USE.NAMES=F)
    stopCluster(cluster)
    ### do something with predicition vector
  })
}

run.knn <- function(labelst,labelsv,train.points,test.points,step) {
  predictions <- knn(train = train.points, 
                     test = test.points,
                     cl = labelst, 
                     k = step
  )
  hits <- predictions == labelsv
  return(hits)
}

################################################################################
#INIT
tissue.to.avoid <- bad.tissues()
sample.data <- load.data("relevant_probe_matrix.csv",tissue.to.avoid)
sample.info <- load.data("cup_usable_samples.csv",tissue.to.avoid)
set.seed(42)
###############################################################################
##KNN VALIDATION CURVES
filters <- fread("probes.csv", row.names = 1) ## file with top 1000probes personalized for cv-iteration


validation.error.values <- sapply(c(50,100,250,500,1000), function(i) {
  predictions <- run.knn.curve(
    labels <- sample.info$tissue_simple,
    cross.validation.groups <- sample.info$lab_batch,
    filters <- filters,
    row.cutoff <- i
  )
})sample.info$tissue_simple