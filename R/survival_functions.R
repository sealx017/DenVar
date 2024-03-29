#'@title Fit CoxPH (and Coxme) model with cluster labels and plot the KM curves
#' Fit a CoxPH model and a Coxme model with the cluster labels obtained by hierarchical clustering of the computed JSD matrix
#'
#' @param surv is a data.frame with a column named "SampleID" denoting the patient IDs, a column named
#' "Survival.in.days" denoting the patient survival, a column named "Censored" denoting the censoring status, and
#' a column named "Cluster" denoting the cluster label of thge patients found using JSD based hierarchical clusetring
#'
#' @return The function, coxPH_plot fits a CoxPH model with the cluster labels and plots the KM curves specific to the
#' clusters accompanied by the estimated Hazard Ratio (HR) and the p-value for testing its significance. The function,
#' coxme_model fits a CoxPH model with a random effects capturing the characteristics of the JSD matrix (Coxme model).
#' @export
KMplot<-function(ggsurvobject, HR, pval)
{

  print(ggsurvobject$plot + ggplot2::annotate(
    "text",
    x = Inf, y = Inf,
    vjust = 1, hjust = 1,
    label = paste0("HR = ", round(HR, 4), "\n p < ", round(pval,4)),
    size = 5))
}

#' @export
coxPH_plot = function(surv, name = "Non-recurrence probability"){

  ggsobj <- ggsurvplot(fit = survminer::surv_fit(Surv(Time,Censored) ~ Cluster, data = surv),
                       xlab = "Time", ylab = name, risk.table = T)

  res.cox<-coxph(Surv(Time,Censored) ~ Cluster, data = surv)

  p1 <- KMplot(ggsobj, summary(res.cox)$conf.int[2], summary(res.cox)$coefficients[5])
  #print(p1)

  return(list(plot = p1, HR = summary(res.cox)$conf.int[2], pvalue = summary(res.cox)$coefficients[5]))
}


#' @export
coxme_model = function(surv, JSD_mat){
  exp_dist = exp(-JSD_mat)
  res.coxme<-coxme(Surv(Time,Censored) ~ (1|ID), data = surv,
                   varlist = coxmeMlist(list(exp_dist),rescale = FALSE,
                   pdcheck = TRUE, positive = TRUE))
  LRT <- 2*(res.coxme$loglik[2] - res.coxme$loglik[1])
  LRT_pvalue <- 1-pchisq(LRT,1)
  return(list(Variance_estimate = res.coxme$vcoef, LRT = LRT, LRT_pvalue = LRT_pvalue))
}
