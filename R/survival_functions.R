#'@title Fit CoxPH model with cluster labels and plot the KM curves
#' Fit a CoxPH model with the cluster labels obtained by hierarchical clustering of the computed JSD matrix
#'
#' @param surv is a data.frame with a column named "SampleID" denoting the patient IDs, a column named 
#' "Survival.in.days" denoting the patient survival, a column named "Censored" denoting the censoring status, and
#' a column named "Cluster" denoting the cluster label of thge patients found using JSD based hierarchical clusetring
#' 
#' @return The function coxPH_plot fits a CoxPH model with the cluster labels and plots the KM curves specific to the 
#' clusters accompanied by the estimated Hazard Ratio (HR) and the p-value for testing its significance 
#' @export

coxPH_plot = function(surv){
  ggsobj <- ggsurvplot(fit = survminer::surv_fit(Surv(`Survival.in.days`,Censored) ~ Cluster, data = surv), 
                       xlab = "Days", ylab = "Proportion alive",risk.table = T)
  
  res.cox<-coxph(Surv(`Survival.in.days`,Censored) ~ Cluster, data = surv)
  
  p1 <- KMplot(ggsobj, summary(res.cox)$conf.int[2], summary(res.cox)$coefficients[5])
  #print(p1)
  
  return(list(plot = p1, HR = summary(res.cox)$conf.int[2], pvalue = summary(res.cox)$coefficients[5]))
}


#' @export
coxme_model = function(surv, JSD_mat){
  exp_dist = exp(-JSD_mat)
  res.coxme<-coxme(Surv(`Survival.in.days`,Censored) ~ (1|SampleID), data = surv,
                   varlist = coxmeMlist(list(exp_dist),rescale = FALSE,
                   pdcheck = TRUE, positive = TRUE))
  LRT <- 2*(res.coxme$loglik[2] - res.cox$loglik[1])
  LRT_pvalue <- 1-pchisq(LRT,1)
  return(list(Variance_estimate = res.coxme$vcoef, LRT = LRT, LRT_pvalue = LRT_pvalue))
}