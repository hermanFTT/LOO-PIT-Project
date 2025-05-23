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

compute_all_pit<-function(brm_model,sim_data,test_data){
      # get hypothetical replicated data
      y_rep <-posterior_predict(brm_model)
      # compute pit
      y<-sim_data$y
      obs_pit<-pit(y_rep,y) 
       # compute loo-pit
      lw<-weights(psis(brm_model))      # normalized  log weighs from psis scheme
      loo_pit<-pit(y_rep,y,weights=lw,log = TRUE)
      # compute pit with test data
      y_pred_test <- posterior_predict(brm_model,newdata = test_data[1]) #test_data[1]
      y_test<-test_data$y
      test_pit<-pit(y_pred_test,y_test)
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

# Generate synthetic  data
generate_data  <- function(N = 10, seed = 978205, df=100, std = 0.2, beta1 = 1.5,resp_type="continuous") {
 
      if(resp_type=="continuous") {
      # continuous data with student_t noise
      set.seed(seed)
      x1 <- rnorm(N) # Covariate
      noise <-std*rt(N,df=df) #student_t  #std*rnorm(N)  # Gaussian noise
      y <- beta1 * x1 + noise      # Response variable
      }
      else if(resp_type=="count") {
        # count data from poisson dist
        set.seed(seed)
         x1 <- rnorm(N)  
        lambda<-exp(beta1*x1)
        y<-rpois(N,lambda)      
      }
      data.frame(x1 = x1, y = y)
    }


# simulate data to fit GP models ( Hilbert space based GP approximation)
 # load the stan file
hsgp_model<- cmdstan_model("HSGP_simulate.stan")
# writeLines(readLines(hsgp_model))

simul_hsgp <- function(model=hsgp_model,bound=c(0,1),N=100,ls=0.2,sig=1,noise=0.3,M_f=160,c_f=3,inter_f=0.4,seed=111222){

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
      sim_data <- data.frame(x1=x1,y=y_sim) #,N=length(x1),c_f=c_f,M_f=M_f,lengthscale_f=c_ls,sigma_f=c_sig)

 return( sim_data)
}


  # visualize
#sim_data <- simul_hsgp(seed=42)
#plot(sim_data$x1,sim_data$y)




##############  experiments  pipeline 

# dummy data & model compilation
sim_data<-generate_data(N=20)
# choose model
  # Gaussian noise
lin_model_fit <- brm(y ~ x1, data =sim_data,

                              prior=c(prior(normal(0,2), class="b"),
                                      prior(normal(0,0.3),class="Intercept"),
                                      prior(normal(0,0.5),class="sigma") ),
                              silent = 2, refresh = 0, cores=4)

# student t
 st_fit <- brm(
  y ~ x1, 
    data = sim_data,
    family = student(), # control degrees of freedom 
    prior = c(
      prior(normal(0, 2), class = "b"),
      prior(normal(0, 0.3), class = "Intercept"),
      prior(normal(0, 0.5), class = "sigma"),
      prior(constant(20),class="nu")),           # control the degree of freedom 
    silent = 2, refresh = 0, cores = 4
  ) 

# GP model
sim_data <- simul_hsgp(N=50,seed=42)
#gp_model <-  cmdstan_model("gp_model.stan")
#fit <- gp_model$sample(data=sim_data)

  gp_fit <- brm(
  y ~  gp(x1), # control linearity wiht gp hyperparams 
  data = sim_data,
  prior = c(
   # prior(normal(0, 2), class = "b", coef="x1"),                   # linear effect 
    prior(normal(0, 0.3), class = "Intercept"),                    # Prior for intercept
    prior(normal(0, 0.4), class = "sigma"),                        # Prior for residual noise
       prior(normal(1, 1e-6), class = "sdgp",coef="gpx1"),        # control  GP magnitude (std dev
    prior(normal(0.2,1e-6), class = "lscale",coef="gpx1")         #control GP lengthscale
  ),
  silent = 2, refresh = 0, cores = 4
  )


########################## experiment pipeline

run_experiment<-function(brm_obj=st_fit, level=0.95,alpha=0.05,num_obs=c(10,20), K=4,stat_tests=uniformity_tests,seed=978205,seed_test=677222,summary_stats=list(mean=mean,var=var)) {
  
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
 
  test_data<-generate_data(N=n, seed=seed_test+n)         # generate test data
  
  # repeat accross datasets 
  
  for (k in 1:K) {
    
     sim_data<-generate_data(N=n,seed=seed+k+n)             # generate observed data
        
     # fit the model to the data
     fit<-update(brm_obj, newdata = sim_data)
     # Coverage
     # test-coverage state
     y_pred_k<-posterior_predict(fit)
     covered_pp[k]<-coverage_state(y_pred_k,test_data$y,level=level)
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

Results<-run_experiment(K=500,num_obs=c(20,50,100,500))
save(Results, file = "results_studentt_df_20.RData")
Resumts <- load("results_studentt_df_80.RData")
print(Results)

Rates<-rejection_rates(Results$statistics$num_obs_50$test_pit,n=50)

print(Rates)


df=Results$pvalues$num_obs_500$pit


 df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()

########################################################### draft #######################################################################

## stan model
yrep <- fit$draws(format = "matrix", variable = "y_rep")
log_lik <- fit$draws(format = "matrix", variable = "log_lik")

#compute weight
pp_check(y, yrep, fun="pit_ecdf")
w_to_use <- weights(fit$loo(save_psis=TRUE)$psis_object)

# predict for new data
preds <- model$generate_quantities(fit$draws(), data = new_data)





########### Questions & discussions

 # Average accross datasets and store the results
     # pit_vals_df<-data.frame(pit=colMeans(pit_mat),test_pit=colMeans(test_pit_mat),loo_pit=colMeans(loo_pit_mat))
     # map to  the corresponding number of observations used
     # PIT_results[[paste0("num_obs_",n)]]<-pit_vals_df


 ### Friday 25 summary:
  # Information on how uniform the pit-vals are:
     # * Use uniformity tests: KS, T1 and U1 : check Graphical uniformity test paper and code
  # ** discuss with Aki about estimation the critical value ( threshold) for rejection ( samples exceding that value are rejected )
  # ** 
  
# LOO-coverage (a bit  ambiguous) repeat simulations with new dataset and new left out observation (true ) ??? (
    # check "loo_predict" from brms to help with computation using the pointwise loo predictive samples. 
  
  
        
# do we need to transform pit vals as it is done with for  power analysis in GUT paper ?
# Find the threshold for rejection rate:  compute critical value for chosen test (on uniform) and use as a rejection for samples exceeding that value ??
# compare base on the p-value ??
# long story short: compare test statistic to critical value or evaluate p-value computed from the test statistic ?





###### another ways to proceed 

# critical vals computation
compute_critical_val <- function(statistic=kolmogorov_test,sample_size=50, n_iterations = 100000, alpha = 0.05) {

 #generate random numbers from the null dist ( uniform distribution U(0, 1) ). For each sample, the test s tatistic T is calculated and stored in ascending order. The critical value t is estimated by using the T -value placed in the position "[n_iterations(1-alpha)] + 1".
  
  # approximate the sampling distribution of the test statistic under the null hypothesis. 
  T_values <- replicate(n_iterations, {
    sample<- runif(sample_size)    # Simulate data under the null hypothesis (Uniform[0,1])
    statistic(sample)              # Calculate the test statistic for the simulated data
  })

  # Sort T values  ( order statistic)
  T_values <- sort(T_values)
  # Calculate the index for the (1 - alpha) quantile
  crit_val_index <- ceiling((1 - alpha) * n_iterations)
  # Get the critical value from sorted test statistics
  critical_value <- T_values[crit_val_index]
  return(critical_value)
}




# p-values computation
compute_pvals <- function(obs_stats,size,T,n_iter=1000){
  #obs_stat (df) : observed statistics over experiment for different Tests.

     # Vectorized computation of  p(T>=t) via MC simulation

  # Samples dist of different test the under the null ( df )
  
     pvals<- vapply(obs_stats, function(t) mean( mc_simulation(T,size,n_iter) >= t ), FUN.VALUE = numeric(1) ) 

       # names(pvals) <- names(obs_stat_df) 
            return(pvals)
      }
    ## compute rejections rates for different test stattistics


  #pvals_pit  as.data.frame( sapply(stat_tests, function(T) compute_pvals( apply(pit_mat, 1, T) ),n,T ))
  #pvals_test_pit  as.data.frame( sapply(stat_tests, function(T) compute_pvals( apply(test_pit_mat, 1, T) ),n,T ))
  #pvals_loo_pit  as.data.frame( sapply(stat_tests, function(T) compute_pvals( apply(loo_pit_mat, 1, T) ),n,T ))
