ICA_predict <- function(ICAS, sample, filter){
  current <- as.data.frame(matrix(0, ncol = 0, nrow = 636))
  current["ICAs"] <- sample*ICAS$directions
  current["likely"] <- ICAS$affected
  current["AUC"] <- ICAS$auc_all_samples
  dt <- as.data.table(current[filter,])
  #print(dt[,sum(ICAs), by=likely][order(V1, decreasing = T),])
  n <- dt[,sum(ICAs), by=likely][V1 > 0][order(V1, decreasing = T),][1]
  #print(dt)
  #return(dt[,mean(ICAs), by=likely][,likely])
  return(n$likely)
  #dt
  #predicted <- dt[order(ICAs, decreasing = T),likely[2]]
  #cum_weigths <- dt[, sum(ICAs), by=likely][,V1]
  #return(cum_weigths)
  #return(dt$likely)
  #return(dt[, sum(ICAs), by=likely][,V1])
  #predicted <- dt[,mean(ICAs), by=likely][,V1]
  #all <- dt[, abs(ICAs)]
  #return( shannon_diversity(abs(cum_weigths)))
}
shannon_diversity <- function(x){
  return( log( 1/prod(x ^ x) ))
}

Get_index_of_bests <- function(ICAS, top, mask){
  #return(ICAS[order(1-auc_all_samples)][, .SD[1:top,], by=affected][,affected])
  if (is.na(top)){
    return(ICAS[!(affected %in% mask),V1])
  }
  return(
    ICAS[order(1-auc_all_samples)][, .SD[1:top,], by=affected][!(affected %in% mask),V1]
  )
}
n_top <- NA
predictions <- sapply(1:7651, function(x){
  gender <- sample.info$gender[x]
  if(is.na(gender)){
    if(transposed[x,ICAV201] < 0.01){
      gender <- "male"
    } else{
      gender <- "female"
    }
  }
  if(gender == "female"){
    gender_mask <- c("gM-prostate",
                     "gM-testicle",
                     "gM-penis")
  } else{
    gender_mask <- c("gF-fallopian",
                     "gF-uterus_cervix",
                     "gF-breast",
                     "gF-vulva",
                     "gF-vagina",
                     "gF-ovary",
                     "gF-uterus")
  }
  if(sample.info[x,tissue_type] == "normal_histology"){
    prediction <- ICA_predict(healthy_ICAS,
                              sample,
                              Get_index_of_bests(healthy_ICAS,#[discr_cancer_version==F], # only use ICAS that predict healthy version
                                                 n_top,
                                                 c("lymphoid_node","heart", gender_mask)
                              )
    )
  } else if(sample.info[x,tissue_type] == "cancer"){
    prediction <- ICA_predict(cancer_ICAS,
                              sample,
                              Get_index_of_bests(cancer_ICAS,
                                                 n_top,
                                                 c("lymphoid_node","heart", gender_mask)
                              )
    )
  } else if(sample.info[x,tissue_type] == "metastatic" | sample.info[x,tissue_type] == "cancer_unknown_primary"){
    prediction <- ICA_predict(cancer_ICAS,
                              sample,
                              Get_index_of_bests(cancer_ICAS,
                                                 n_top,
                                                 c(gender_mask, sample.info[x,incidental_tissue]) # mask incidental tissue
                              )
    )
  } else {
    prediction <- NA
  }
  print(x)
  return(prediction)
})

n_top <- 3
predictions <- sapply(1:7651, function(x){
  sample <- t(transposed[x,-(1)])
  prediction <- ICA_predict(cnICAs,
                            sample,
                            cnICAs[order(p_value_all_samples)][2:n_top,V1]
  )
  print(x)
  return(prediction)
})