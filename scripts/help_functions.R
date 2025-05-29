

###  helping functions ####


##############  pit values computation

compute_all_pit<-function(brm_model,psis,sim_data,test_data){
      # get hypothetical replicated data
      y_rep <-posterior_predict(brm_model)
      # compute pit
      y<-sim_data$y
      obs_pit<-pit(y_rep,y) 
      # compute loo-pit
     # psis_2 <- loo(brm_model, save_psis = TRUE, cores = 4)$psis_object
      lw<-weights(psis)      # normalized  log weighs from psis scheme
      y_rep <-posterior_predict(brm_model,newdata=sim_data[1])
      loo_pit<-pit(y_rep,y,weights=lw,log = TRUE)
      # compute pit with test data
      y_pred_test <- posterior_predict(brm_model,newdata = test_data[1]) #test_data[1]
      y_test<-test_data$y
      test_pit<-pit(y_pred_test,y_test)
      results<-list(pit=obs_pit,test_pit=test_pit, loo_pit=loo_pit)
     }





# estimate the critical value under the null hypothesis for a given test statistic via Monte Carlo simulation.

compute_critical_val<-function(statistic=kolmogorov_test,sample_size=50, n_iterations = 10000, alpha = 0.05){
  #sample_size:   sample size
  # n_iterations:  number of Monte Carlo replicates
  # alpha:         significance level

  #  approximates the sampling distribution of the test statistic under the null.
  T_vals <- replicate(n_iterations, { #repeat 
    sample<- runif(sample_size)    # Simulate data under the null hypothesis (Uniform[0,1])
    statistic(sample)              # Calculate the test statistic for the simulated data
  })
  
#  Estimate the critical value for the upper bound of the accepted region
  crit_val <- quantile(T_vals, 1 - alpha)  # 95th percentile for a 5% significance level
  return(crit_val)
  }





# estimate p-value with MC simulation  under the null hypothesis.
mc_simulation <- function(statistic=kolmogorov_test,sample_size=50, n_iterations =100000) {
      T_vals <- replicate(n_iterations, { #repeat 
      sample<- runif(sample_size)    # Sample  under the null hypothesis (Uniform[0,1])
      statistic(sample)              # evaluate the test statistic
     })
  }



compute_pvalues <- function(obs_stat_df,size,stat_tests=uniformity_tests){
  #obs_stat (df) : observed statistics over experiment for different Tests.

     # Vectorized computation of  p(T>=t) via MC simulation

     # Samples dist of different test the under the null ( df ) 
     T_null_df <- as.data.frame(sapply(stat_tests, function(f)  mc_simulation(f,size))) 

        pvals <- as.data.frame( vapply(seq_along(T_null_df), function(j) {
           sapply(obs_stat_df[[j]], function(t) mean(T_null_df[[j]] >= t))
          }, FUN.VALUE = numeric(nrow(obs_stat_df)))
        )
        names(pvals) <- names(obs_stat_df) 
            return(pvals)
}





rejection_rates<-function(Results,obs_sizes=NULL,nu_s=NULL,bs=NULL,N=NULL,stat_tests=uniformity_tests,alpha=0.05,n_iter=10000) {

  # rstat_list (list) : list of list of df containing realized statistics for diff test
  # obs_sizes (vector): vector containing diff number of observations size.
  rstat_list <- Results$statistics
  rej_rates <- list()

  for (key in c("pit","test_pit","loo_pit")) {
    rate_mat<-matrix(numeric(0), ncol=length(stat_tests), nrow=0)

    if( !is.null(obs_sizes) ) {
      
    for( n in obs_sizes) {
      df <- rstat_list[[paste0("num_obs_",n)]][[key]]
       # estimate critical vals via MC simulation for each test
       critical_values<-sapply(stat_tests, function(stat) compute_critical_val(stat,n,n_iter,alpha))
       # identify sample whose realized statistic exceed the critical value ( logical matrix)
       reject_state_mat<-sweep(df,2,critical_values,FUN=">")
       # rejection rate as a proportion of rejected samples for each test statistics
       rate<-colMeans(reject_state_mat)

      rate_mat <- rbind(rate_mat,rate)
      rownames(rate_mat)[nrow(rate_mat)] <- paste0("N_",n)

    } }

    if (!is.null(bs)) {

    for( beta in bs) {
      df <- rstat_list[[paste0("beta_",beta)]][[key]]
       # estimate critical vals via MC simulation for each test
       critical_values<-sapply(stat_tests, function(stat) compute_critical_val(stat,N,n_iter,alpha))
       # identify sample whose realized statistic exceed the critical value ( logical matrix)
       reject_state_mat<-sweep(df,2,critical_values,FUN=">")
       # rejection rate as a proportion of rejected samples for each test statistics
       rate<-colMeans(reject_state_mat)

      rate_mat <- rbind(rate_mat,rate)
      rownames(rate_mat)[nrow(rate_mat)] <- paste0("beta_",beta)

    }
 }

     if ( !is.null(nu_s) ) {

    for( nu in nu_s) {
      df <- rstat_list[[paste0("df_",nu)]][[key]]
       # estimate critical vals via MC simulation for each test
       critical_values<-sapply(stat_tests, function(stat) compute_critical_val(stat,N,n_iter,alpha))
       # identify sample whose realized statistic exceed the critical value ( logical matrix)
       reject_state_mat<-sweep(df,2,critical_values,FUN=">")
       # rejection rate as a proportion of rejected samples for each test statistics
       rate<-colMeans(reject_state_mat)

      rate_mat <- rbind(rate_mat,rate)
      rownames(rate_mat)[nrow(rate_mat)] <- paste0("df_",nu)

    }
 }
  
  rej_rates[[key]] <- as.data.frame(rate_mat)

}
  
 
  return(rej_rates)  #,crit_vals=critical_values
  
}



##############  Posterior predictive Coverage:

# Posterior predictive coverage ( coverage probability for test data )
coverage_state <- function(samples_df, observed, level = 0.95) {
  
          stopifnot(ncol(samples_df) == length(observed))
          alpha <- 1 - level
          # Compute lower and upper bounds for 1-alpha central interval each variable: 95% pointwise interval 
          lower <- apply(samples_df, 2, quantile, probs = alpha / 2)
          upper <- apply(samples_df, 2, quantile, probs = 1 - alpha / 2)
  
          #if  simultanuous coverage:  whether all observed values fall within their respective intervals (simultaneous)
          # sim_covered <- all(observed >= lower & observed <= upper)
          # marginal coverage: Logical vector: TRUE if observed[i] is within interval [lower[i], upper[i]]
  
          coverage <-mean( observed >= lower & observed <= upper)
          return(coverage)  #
        }


# estimate the confidence interval of the empirical coverage 
margin_coverage_ci <- function(coverage_vec, alpha = 0.05) {
          # coverage_vec (vector): empirical coverage accross simulations
          # alpha: significance level for the conf interval of the estimated coverate proportion  
             
           n <- length(coverage_vec)
          # confidence interval of mean  empirical coverage  by normal approximation ( for large n: CLT)
           mean_coverage <- mean(coverage_vec)
          # var_coverage <- var(coverage_vec)
           std_err <- sd(coverage_vec) / sqrt(n)

           z <- qnorm(1 - alpha/2)  # Critical value from standard normal
           margin <- z * std_err

           lower <-max(0, mean_coverage - margin)
           upper <-min(1, mean_coverage + margin)
             return(c(coverage= mean_coverage, lower= lower, upper= upper))
             # return(list(mean_coverage= mean_coverage, lower_ci= lower, upper_ci = upper))
            }


