
######## Varying observations size 


run_experiment_n<-function(brm_obj,generator,b=0.2,nu=100,levels=c(0.9),alpha=0.05,num_obs=c(16,32),K=4,stat_tests=uniformity_tests,seed=444,seed_test=555,summary_stats=list(mean=mean,var=var)) {
  
#args:
  # brm_obj: stan model object
  #level:([0,1]) central interval level for coverage
  #alpha:( [0,1])  significance level for estimate conf interval of the estimated coverage proportion
  #num_obs : (vector) varying number of observations
  # K : (int) number of simulations 
  # stat_tests: (list) list of uniformity test statistic
  #seed: key for reproducible data
  #seed_test: key for reproducible test data.

  #Output:( list) : realized statistics (list of df), emp_coverage+ conf interval ( list) accross simulations. 

  
coverage_results<-list()
pit_statistic_results<-list()
pvalues_results <- list()
pit_summary_results <- list()
# fitted_model <- list()

for (n in num_obs) {

   pit_mat<-matrix(numeric(0), ncol=n, nrow=0)
   test_pit_mat<-matrix(numeric(0), ncol=n,nrow = 0)
   loo_pit_mat<-matrix(numeric(0), ncol = n, nrow =0)
   covered_pp<-matrix(NA, nrow = K, ncol = length(levels))  #numeric(K)
  covered_loo<-matrix(NA, nrow = K, ncol = length(levels))
  
  colnames(covered_pp) <- paste0("level_", levels)
  colnames(covered_loo) <- paste0("level_", levels)
  
  test_data<-generator(N=n, seed=seed_test+n,beta=b,df=nu)         # generate test data
  
  # repeat accross datasets 
  
  for (k in 1:K) {
    
     sim_data<-generator(N=n,seed=seed+k+n,beta=b,df=nu)             # generate observed data
        
     # fit the model to the data
     fit<-update(brm_obj, newdata = sim_data)
     # Coverage
    # test-coverage state

    y_pred_test<-posterior_predict(fit,newdata = test_data[1])
    
    f_vec1 <- Vectorize(function(level) {
       test_interval <- predictive_interval(y_pred_test,prob=level)
       mean(test_data$y>=test_interval[, 1] & test_data$y <= test_interval[, 2])})    #coverage_state(y_pred_test,test_data$y,level=level)

         covered_pp[k,]<-f_vec1(levels)

   # covered_pp[k,]<-sapply(levels, function(level) {
    #   test_interval <- predictive_interval(y_pred_test,prob=level)
     #  mean(test_data$y>=test_interval[, 1] & test_data$y <= test_interval[, 2])})     #coverage_state(y_pred_test,test_data$y,level=level)

    # loo-coverage state                                                                #loo::psis(-log_lik(fit), cores = 4)
    psis <-loo(fit, save_psis = TRUE, cores = 1)$psis_object
    
     f_vec2 <- Vectorize( function(level) {
      loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     covered_loo[k,]<-f_vec2(levels)


   # covered_loo[k,]<-sapply(levels,function(level) {
    #   psis <-loo(fit, save_psis = TRUE, cores = 4)$psis_object                                  #loo::psis(-log_lik(fit), cores = 4)
     # loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      #mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     
         # Compute pit with different approaches (pit, test_pit, loo_pit)
         pit_vals<-compute_all_pit(fit,psis,sim_data,test_data)
         # Store 
         pit_mat<-rbind(pit_mat,pit_vals$pit)
         test_pit_mat<-rbind(test_pit_mat,pit_vals$test_pit)
         loo_pit_mat<-rbind(loo_pit_mat,pit_vals$loo_pit)
    rm(fit)
  gc() }
  
  # Empirical coverage :
  emp_test_coverage<-as.data.frame(t( apply(covered_pp,2,function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})))
  
 # summary_cov_pp<-lapply(as.data.frame(covered_pp),function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})
   # emp_test_coverage <-as.data.frame(do.call(rbind, summary_cov_pp))
    emp_test_coverage$level <-levels # names(summary_cov_pp)


  emp_loo_coverage<-as.data.frame(t(apply(covered_loo,2,function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})))
 
  #summary_cov_loo<-lapply(as.data.frame(covered_loo), function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})
   #emp_loo_coverage <- as.data.frame(do.call(rbind, summary_cov_loo))
   emp_loo_coverage$level <- levels #names(summary_cov_loo)

  # pit vals summary
  summary_pit<-as.data.frame( sapply(summary_stats, function(f) apply(pit_mat, 2, f) ))
  summary_test_pit<- as.data.frame( sapply(summary_stats, function(f) apply(test_pit_mat,2, f) ))
  summary_loo_pit<- as.data.frame( sapply(summary_stats, function(f) apply(loo_pit_mat,2, f) ))

  # collect information about pit's vals uniformity : uniformity test statistic
  
  #decision-based critical value:  compute sample statistics to compare agains a critical value 
  stat_pit<- as.data.frame( sapply(stat_tests, function(f) apply(pit_mat, 1, f) ))
  stat_test_pit<- as.data.frame( sapply(stat_tests, function(f) apply(test_pit_mat, 1, f) ))
  stat_loo_pit<- as.data.frame( sapply(stat_tests, function(f) apply(loo_pit_mat, 1, f) ))

  #pvaue-based decision: evaluate p-values at specified significance level (alpha)
  pvals_pit <- compute_pvalues(stat_pit,n,stat_tests)         
  pvals_test_pit<-compute_pvalues(stat_test_pit,n,stat_tests)
  pvals_loo_pit <- compute_pvalues(stat_loo_pit,n,stat_tests)

  # map to  the corresponding number of observations used
  pit_summary_results[[paste0("num_obs_",n)]] <- list(pit=summary_pit,test_pit=summary_test_pit,loo_pit=summary_loo_pit)
  pvalues_results[[paste0("num_obs_",n)]] <-list(pit=pvals_pit,test_pit=pvals_test_pit,loo_pit=pvals_loo_pit) 
  pit_statistic_results[[paste0("num_obs_",n)]]<-list(pit=stat_pit,test_pit=stat_test_pit,loo_pit=stat_loo_pit)
  coverage_results[[paste0("num_obs_",n)]]<-list(test_coverage=emp_test_coverage ,loo_coverage=emp_loo_coverage) 
  #data.frame( rbind( test_coverage = emp_test_coverage, loo_coverage = emp_loo_coverage ))
 # fitted_model[[paste0("num_obs_",n)]] <- fit
}
  return( list(pvalues=pvalues_results,statistics=pit_statistic_results,coverage=coverage_results,summary=pit_summary_results) )
}






###################### control non linearity  ( varying beta in data generating mechanism  based Hilbert space GP )



run_experiment_b<-function(brm_obj,generator=simul_hsgp,b=c(0.2,0.4),nu=100,levels=c(0.9),alpha=0.05,n=32,K=4,stat_tests=uniformity_tests,seed=444,seed_test=555,summary_stats=list(mean=mean,var=var)) {
  
#args:
  # brm_obj: stan model object
  #level:([0,1]) central interval level for coverage
  #alpha:( [0,1])  significance level for estimate conf interval of the estimated coverage proportion
  #num_obs : (vector) varying number of observations
  # K : (int) number of simulations 
  # stat_tests: (list) list of uniformity test statistic
  #seed: key for reproducible data
  #seed_test: key for reproducible test data.

  #Output:( list) : realized statistics (list of df), emp_coverage+ conf interval ( list) accross simulations. 

  
coverage_results<-list()
pit_statistic_results<-list()
pvalues_results <- list()
pit_summary_results <- list()
# fitted_model <- list()

  for ( i in seq_along(b) ) {
    
    beta <- b[i]

   pit_mat<-matrix(numeric(0), ncol=n, nrow=0)
   test_pit_mat<-matrix(numeric(0), ncol=n,nrow = 0)
   loo_pit_mat<-matrix(numeric(0), ncol = n, nrow =0)
   covered_pp<-matrix(NA, nrow = K, ncol = length(levels))  #numeric(K)
   covered_loo<-matrix(NA, nrow = K, ncol = length(levels))
  
  colnames(covered_pp) <- paste0("level_", levels)
  colnames(covered_loo) <- paste0("level_", levels)
  
 
   test_data<-generator(N=n, seed=seed_test+i,beta=beta,df=nu)         # generate test data
  # repeat accross datasets 
  
  for (k in 1:K) {
    
     sim_data<-generator(N=n,seed=seed+k+i,beta=beta,df=nu)             # generate observed data
        
     # fit the model to the data
     fit<-update(brm_obj, newdata = sim_data)
     # Coverage
    # test-coverage state

    y_pred_test<-posterior_predict(fit,newdata = test_data[1])
    
    f_vec1 <- Vectorize(function(level) {
       test_interval <- predictive_interval(y_pred_test,prob=level)
       mean(test_data$y>=test_interval[, 1] & test_data$y <= test_interval[, 2])})    #coverage_state(y_pred_test,test_data$y,level=level)

         covered_pp[k,]<-f_vec1(levels)

   # covered_pp[k,]<-sapply(levels, function(level) {
    #   test_interval <- predictive_interval(y_pred_test,prob=level)
     #  mean(test_data$y>=test_interval[, 1] & test_data$y <= test_interval[, 2])})     #coverage_state(y_pred_test,test_data$y,level=level)

    # loo-coverage state                                                                #loo::psis(-log_lik(fit), cores = 4)
    psis <-loo(fit, save_psis = TRUE, cores = 1)$psis_object
    
     f_vec2 <- Vectorize( function(level) {
      loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     covered_loo[k,]<-f_vec2(levels)


   # covered_loo[k,]<-sapply(levels,function(level) {
    #   psis <-loo(fit, save_psis = TRUE, cores = 4)$psis_object                                  #loo::psis(-log_lik(fit), cores = 4)
     # loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      #mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     
         # Compute pit with different approaches (pit, test_pit, loo_pit)
         pit_vals<-compute_all_pit(fit,psis,sim_data,test_data)
         # Store 
         pit_mat<-rbind(pit_mat,pit_vals$pit)
         test_pit_mat<-rbind(test_pit_mat,pit_vals$test_pit)
         loo_pit_mat<-rbind(loo_pit_mat,pit_vals$loo_pit)

    rm(fit)
    gc()    }
  
  # Empirical coverage :
  emp_test_coverage<-as.data.frame(t(apply(covered_pp,2,function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})))
  
 # summary_cov_pp<-lapply(as.data.frame(covered_pp),function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})
   # emp_test_coverage <-as.data.frame(do.call(rbind, summary_cov_pp))
    emp_test_coverage$level <-levels # names(summary_cov_pp)


  emp_loo_coverage<-as.data.frame(t(apply(covered_loo,2,function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})))
 
  #summary_cov_loo<-lapply(as.data.frame(covered_loo), function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})
   #emp_loo_coverage <- as.data.frame(do.call(rbind, summary_cov_loo))
   emp_loo_coverage$level <- levels #names(summary_cov_loo)

  # pit vals summary
  summary_pit<-as.data.frame( sapply(summary_stats, function(f) apply(pit_mat,2, f) ))
  summary_test_pit<- as.data.frame( sapply(summary_stats, function(f) apply(test_pit_mat,2, f) ))
  summary_loo_pit<- as.data.frame( sapply(summary_stats, function(f) apply(loo_pit_mat,2, f) ))

  # collect information about pit's vals uniformity : uniformity test statistic
  
  #decision-based critical value:  compute sample statistics to compare agains a critical value 
  stat_pit<- as.data.frame( sapply(stat_tests, function(f) apply(pit_mat, 1, f) ))
  stat_test_pit<- as.data.frame( sapply(stat_tests, function(f) apply(test_pit_mat, 1, f) ))
  stat_loo_pit<- as.data.frame( sapply(stat_tests, function(f) apply(loo_pit_mat, 1, f) ))

  #pvaue-based decision: evaluate p-values at specified significance level (alpha)
  pvals_pit <- compute_pvalues(stat_pit,n,stat_tests)         
  pvals_test_pit<-compute_pvalues(stat_test_pit,n,stat_tests)
  pvals_loo_pit <- compute_pvalues(stat_loo_pit,n,stat_tests)

  # map to  the corresponding number of observations used
  pit_summary_results[[paste0("beta_",beta)]] <- list(pit=summary_pit,test_pit=summary_test_pit,loo_pit=summary_loo_pit)
  pvalues_results[[paste0("beta_",beta)]] <-list(pit=pvals_pit,test_pit=pvals_test_pit,loo_pit=pvals_loo_pit) 
  pit_statistic_results[[paste0("beta_",beta)]]<-list(pit=stat_pit,test_pit=stat_test_pit,loo_pit=stat_loo_pit)
  coverage_results[[paste0("beta_",beta)]]<-list(test_coverage=emp_test_coverage ,loo_coverage=emp_loo_coverage) 
  #data.frame( rbind( test_coverage = emp_test_coverage, loo_coverage = emp_loo_coverage ))
 # fitted_model[[paste0("num_obs_",n)]] <- fit
}
  return( list(pvalues=pvalues_results,statistics=pit_statistic_results,coverage=coverage_results,summary=pit_summary_results) )
}





################# controlling the degree of freedom ( student_t model)

run_experiment_nu<-function(brm_obj,generator=simul_stt,b=0.2,nu=c(1,5,10),levels=c(0.9),alpha=0.05,n=32,K=4,stat_tests=uniformity_tests,seed=444,seed_test=555,summary_stats=list(mean=mean,var=var)) {
  
#args:
  # brm_obj: stan model object
  #level:([0,1]) central interval level for coverage
  #alpha:( [0,1])  significance level for estimate conf interval of the estimated coverage proportion
  #num_obs : (vector) varying number of observations
  # K : (int) number of simulations 
  # stat_tests: (list) list of uniformity test statistic
  #seed: key for reproducible data
  #seed_test: key for reproducible test data.

  #Output:( list) : realized statistics (list of df), emp_coverage+ conf interval ( list) accross simulations. 

  
coverage_results<-list()
pit_statistic_results<-list()
pvalues_results <- list()
pit_summary_results <- list()
# fitted_model <- list()

  for ( i in seq_along(nu) ) {
    
    df <- nu[i]

   pit_mat<-matrix(numeric(0), ncol=n, nrow=0)
   test_pit_mat<-matrix(numeric(0), ncol=n,nrow = 0)
   loo_pit_mat<-matrix(numeric(0), ncol = n, nrow =0)
   covered_pp<-matrix(NA, nrow = K, ncol = length(levels))  #numeric(K)
   covered_loo<-matrix(NA, nrow = K, ncol = length(levels))
  
  colnames(covered_pp) <- paste0("level_", levels)
  colnames(covered_loo) <- paste0("level_", levels)
  
 
   test_data<-generator(N=n, seed=seed_test+i,beta=b,df=df)         # generate test data
  # repeat accross datasets 
  
  for (k in 1:K) {
    
     sim_data<-generator(N=n,seed=seed+k+i,beta=b, df=df)             # generate observed data
        
     # fit the model to the data
     fit<-update(brm_obj, newdata = sim_data)
     # Coverage
    # test-coverage state

    y_pred_test<-posterior_predict(fit,newdata = test_data[1])
    
    f_vec1 <- Vectorize(function(level) {
       test_interval <- predictive_interval(y_pred_test,prob=level)
       mean(test_data$y>=test_interval[, 1] & test_data$y <= test_interval[, 2])})    #coverage_state(y_pred_test,test_data$y,level=level)

         covered_pp[k,]<-f_vec1(levels)

   # covered_pp[k,]<-sapply(levels, function(level) {
    #   test_interval <- predictive_interval(y_pred_test,prob=level)
     #  mean(test_data$y>=test_interval[, 1] & test_data$y <= test_interval[, 2])})     #coverage_state(y_pred_test,test_data$y,level=level)

    # loo-coverage state                                                                #loo::psis(-log_lik(fit), cores = 4)
    psis <-loo(fit, save_psis = TRUE, cores = 1)$psis_object
    
     f_vec2 <- Vectorize( function(level) {
      loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     covered_loo[k,]<-f_vec2(levels)


   # covered_loo[k,]<-sapply(levels,function(level) {
    #   psis <-loo(fit, save_psis = TRUE, cores = 4)$psis_object                                  #loo::psis(-log_lik(fit), cores = 4)
     # loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      #mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     
         # Compute pit with different approaches (pit, test_pit, loo_pit)
         pit_vals<-compute_all_pit(fit,psis,sim_data,test_data)
         # Store 
         pit_mat<-rbind(pit_mat,pit_vals$pit)
         test_pit_mat<-rbind(test_pit_mat,pit_vals$test_pit)
         loo_pit_mat<-rbind(loo_pit_mat,pit_vals$loo_pit)

    rm(fit)
    gc()    }
  
  # Empirical coverage :
  emp_test_coverage<-as.data.frame(t(apply(covered_pp,2,function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})))
  
 # summary_cov_pp<-lapply(as.data.frame(covered_pp),function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})
   # emp_test_coverage <-as.data.frame(do.call(rbind, summary_cov_pp))
    emp_test_coverage$level <-levels # names(summary_cov_pp)


  emp_loo_coverage<-as.data.frame(t(apply(covered_loo,2,function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})))
 
  #summary_cov_loo<-lapply(as.data.frame(covered_loo), function(cov_vec) {margin_coverage_ci(cov_vec, alpha = alpha)})
   #emp_loo_coverage <- as.data.frame(do.call(rbind, summary_cov_loo))
   emp_loo_coverage$level <- levels #names(summary_cov_loo)

  # pit vals summary
  summary_pit<-as.data.frame( sapply(summary_stats, function(f) apply(pit_mat,2, f) ))
  summary_test_pit<- as.data.frame( sapply(summary_stats, function(f) apply(test_pit_mat,2, f) ))
  summary_loo_pit<- as.data.frame( sapply(summary_stats, function(f) apply(loo_pit_mat,2, f) ))

  # collect information about pit's vals uniformity : uniformity test statistic
  
  #decision-based critical value:  compute sample statistics to compare agains a critical value 
  stat_pit<- as.data.frame( sapply(stat_tests, function(f) apply(pit_mat, 1, f) ))
  stat_test_pit<- as.data.frame( sapply(stat_tests, function(f) apply(test_pit_mat, 1, f) ))
  stat_loo_pit<- as.data.frame( sapply(stat_tests, function(f) apply(loo_pit_mat, 1, f) ))

  #pvaue-based decision: evaluate p-values at specified significance level (alpha)
  pvals_pit <- compute_pvalues(stat_pit,n,stat_tests)         
  pvals_test_pit<-compute_pvalues(stat_test_pit,n,stat_tests)
  pvals_loo_pit <- compute_pvalues(stat_loo_pit,n,stat_tests)

  # map to  the corresponding number of observations used
  pit_summary_results[[paste0("df_",df)]] <- list(pit=summary_pit,test_pit=summary_test_pit,loo_pit=summary_loo_pit)
  pvalues_results[[paste0("df_",df)]] <-list(pit=pvals_pit,test_pit=pvals_test_pit,loo_pit=pvals_loo_pit) 
  pit_statistic_results[[paste0("df_",df)]]<-list(pit=stat_pit,test_pit=stat_test_pit,loo_pit=stat_loo_pit)
  coverage_results[[paste0("df_",df)]]<-list(test_coverage=emp_test_coverage ,loo_coverage=emp_loo_coverage) 
  #data.frame( rbind( test_coverage = emp_test_coverage, loo_coverage = emp_loo_coverage ))
 # fitted_model[[paste0("num_obs_",n)]] <- fit
}
  return( list(pvalues=pvalues_results,statistics=pit_statistic_results,coverage=coverage_results,summary=pit_summary_results) )
}

