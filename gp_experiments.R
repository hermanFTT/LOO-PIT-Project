library(brms)
library(ggplot2)
library(tidybayes)
library(cmdstanr)
library(tictoc)
library(posterior)
library(bayesplot)
library(modeest)
library(ddst) # adaptive Neyman's smooth goodness of fit test 
library(ggplot2)
library(tidyr)

#install.packages("tictoc")
#options(brms.backend = "cmdstanr", mc.cores = 1)
#install.packages("modeest")

## Track configurations

seed<-978205      # seed for reproductibility 
seed_test<-677222
K<-4               # repeat experiments K times
num_obs<-c(10,20)  # for different  number of observations



###  helping functions ####


##############  pit values computation

# stan model


compute_all_pit<-function(gp_model,fit_obj,sim_data,test_data){
      # get hypothetical replicated data
      y_rep <- fit_obj$draws(format = "matrix", variable = "f")
      # compute pit
      y<-sim_data$y
      obs_pit<-pit(y_rep,y) 
       # compute loo-pit
      lw <- weights(fit_obj$loo(save_psis=TRUE,cores=4)$psis_object)              # normalized  log weighs from psis scheme
      loo_pit<-pit(y_rep,y,weights=lw,log = TRUE)
      # compute pit with test data
      preds <- gp_model$generate_quantities(fit_obj, data = test_data)            # set "seed" later
      y_pred_test=preds$draws(format="matrix",variable="f")
      y_test<-test_data$y
      test_pit<-pit(y_pred_test,y_test)                                          # new input x for test_data ?? 
      results<-list(pit=obs_pit,test_pit=test_pit, loo_pit=loo_pit)
     }




##############  Statistical uniformity tests to assess pit vals uniformity

# Watson's test
watson_u2_test <- function(sample) {
  n <- length(sample)
  U <- sort(sample)
  i <- 1:n
  U_bar <- mean(U)
  W2 <- sum((U - (2 * i - 1) / (2 * n))^2) + 1 / (12 * n)
  U2_stat <- W2 - n * (U_bar - 0.5)^2
  return(U2_stat)
}

# Kolmogorov Sminorv test
kolmogorov_test <-function(sample){
  ks.test(sample,"punif")$statistic
}

 # T1_test ( mean absolute deviation test for uniformity)
 T1_test <- function(sample) {
  n <- length(sample)
  U<-sort(sample)
 U_expected <- (1:n) / (n + 1)
  T_stat <- mean(abs(U-U_expected))
  return(T_stat)
  }

 # Adaptive Neyman smooth goodness of fit test
 neyman_test<-function(sample) {
    ddst.uniform.test(sample,compute.p = TRUE)$statistic
 }




uniformity_tests<-list( ks_stat=kolmogorov_test,u2_stat=watson_u2_test,t1_stat=T1_test)



# estimate the critical value under the null hypothesis for a given test statistic via Monte Carlo simulation.

compute_critical_val<-function(statistic=kolmogorov_test,sample_size=50, n_iterations = 100000, alpha = 0.05){
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

    ## compute rejections rates for different test stattistics

rejection_rates<-function(df,n=20,stat_tests=uniformity_tests,alpha=0.05,n_iter=100000 ){

  # estimate critical vals via MC simulation for each test
  critical_values<-sapply(stat_tests, function(stat) compute_critical_val(stat,n,n_iter,alpha))
  # identify sample whose realized statistic exceed the critical value ( logical matrix)
  reject_state_mat<-sweep(df,2,critical_values,FUN=">")
  # rejection rate as a proportion of rejected samples for each test statistics
  rate<-colMeans(reject_state_mat)
  return(list(rej_rate=rate,crit_vals=critical_values))
  
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

              return(c(emp_coverage= mean_coverage, lower_ic = lower, upper_ic = upper))
            }



##############  Data generating mechanisms & model


# simulate data to fit GP models ( Hilbert space based GP approximation)

# load the stan file
hsgp_model<- cmdstan_model("HSGP_simulate.stan")
# writeLines(readLines(hsgp_model))


simul_hsgp <- function(model=hsgp_model,bound=c(0,1),N=100,ls=0.2,sig=1,noise=0.3,M_f=160,c_f=3,inter_f=0.4,c_ls=0.2,c_sig=1,seed=111222){

       x1=seq(bound[1],bound[2],length.out=N) # univariate input
       set.seed(seed)
       beta_f=rnorm(M_f)
       # input to past to the generator
      data_list<- list(x=x1, # input ( univariate )
                          N=N, # number of observations
                          c_f=c_f, # factor c of basis functions for GP for f1
                          M_f=M_f,  # number of basis functions for GP for f1
                          sigma_f=sig, # magnitude
                          lengthscale_f=ls, #lengthscale
                       beta_f=beta_f,    # coefficient of the basis functions
                       sigman=noise,           # noise level
                       intercept_f=inter_f)

  # simulate 
      sim_hsgp_model<- hsgp_model$sample(data=data_list, fixed_param=TRUE,
                                chains=1, iter=1, iter_sampling=1)
  # draws
      y_sim=as.vector(unlist( as_draws_df(sim_hsgp_model)[,1:N]))
      sim_data <- list(N=length(x1),x=x1,y=y_sim,c_f=c_f,M_f=M_f,lengthscale_f=c_ls,sigma_f=c_sig)

 return( sim_data)
}


  # visualize
#sim_data <- simul_hsgp(seed=42)
#plot(sim_data$x1,sim_data$y)




##############  experiments  pipeline 

# GP model
#sim_data <- simul_hsgp(N=40,seed=42)

gp_model <-  cmdstan_model("gp_model.stan")   # load and compile the gp_model 

#fit <- gp_model$sample(data=sim_data,parallel_chains =4,show_messages = FALSE)


########################## experiment pipeline

run_experiment<-function(stan_model=gp_model, level=0.95,alpha=0.05,num_obs=c(10,20), K=4,stat_tests=uniformity_tests,seed=978205,seed_test=677222,seed_mc=42,summary_stats=list(mean=mean,var=var)) {
  
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
fitted_model <- list()
for (n in num_obs) {

   pit_mat<-matrix(numeric(0), ncol=n, nrow=0)
   test_pit_mat<-matrix(numeric(0), ncol=n,nrow = 0)
   loo_pit_mat<-matrix(numeric(0), ncol = n, nrow =0)
   covered_pp<-numeric(K)
   covered_loo<-numeric(K)

  
   test_data <- simul_hsgp(N=n,seed=seed_test+n)         # generate test data
  
  # repeat accross datasets 
  
  for (k in 1:K) {
    
     sim_data <- simul_hsgp(N=n,seed=seed+n+k)           # generate observed data
        
     # fit the model to the data
     fit <- stan_model$sample(data=sim_data,parallel_chains =4,show_messages = FALSE)
     # Coverage
    # test-coverage state
     preds <- gp_model$generate_quantities(fit, data = test_data)            # set "seed" later
     y_pred_test=preds$draws(format="matrix",variable="f")
     covered_pp[k]<-coverage_state(y_pred_test,test_data$y,level=level)
     # loo-coverage state
     psis <- loo::psis(-log_lik(fit), cores = 4)
     loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
     covered_loo[k]<-all(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])
     
         # Compute pit with different approaches (pit, test_pit, loo_pit)
         pit_vals<-compute_all_pit(fit,sim_data,test_data)
         # Store 
         pit_mat<-rbind(pit_mat,pit_vals$pit)
         test_pit_mat<-rbind(test_pit_mat,pit_vals$test_pit)
         loo_pit_mat<-rbind(loo_pit_mat,pit_vals$loo_pit)
         }
  #  Empirical coverage :
  emp_test_coverage <-margin_coverage_ci(covered_pp, alpha = alpha)
  emp_loo_coverage <-margin_coverage_ci(covered_loo, alpha = alpha)

  # pit vals summary
  summary_pit<-as.data.frame( sapply(summary_stats, function(f) apply(pit_mat, 1, f) ))
  summary_test_pit<- as.data.frame( sapply(summary_stats, function(f) apply(test_pit_mat, 1, f) ))
  summary_loo_pit<- as.data.frame( sapply(summary_stats, function(f) apply(loo_pit_mat, 1, f) ))

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
  coverage_results[[paste0("num_obs_",n)]]<-data.frame( rbind( test_coverage = emp_test_coverage, loo_coverage = emp_loo_coverage ))
  fitted_model[[paste0("num_obs_",n)]] <- fit
}
  return( list(pvalues=pvalues_results,statistics=pit_statistic_results,coverage=coverage_results,summary=pit_summary_results,model=fitted_model) )
}





##### trial

Results<-run_experiment(K=50,num_obs=c(20,40))

print(Results)

Rates<-rejection_rates(Results$statistics$num_obs_20$pit,n=20)

print(Rates)





########################################################### draft #######################################################################

## stan model
yrep <- fit$draws(format = "matrix", variable = "f")
log_lik <- fit$draws(format = "matrix", variable = "log_lik")

#compute weight
pp_check(sim_data$y, yrep, fun="pit_ecdf")
lw <- weights(fit$loo(save_psis=TRUE)$psis_object)

# predict for new data : use new input x1 ? 
test_data<- simul_hsgp(N=40,seed=111)
preds <- gp_model$generate_quantities(fit, data = test_data,seed=222)

y_pred_test=preds$draws(format="matrix",variable="f")


dim(y_pred_test[,-((ncol(y_pred_test)-2):ncol(y_pred_test))])
y_pred_test[1,]

plot(test_data$x,y_pred_test[3,-((ncol(y_pred_test)-2):ncol(y_pred_test))])
plot(sim_data$x,sim_data$y)

y_pred_test
