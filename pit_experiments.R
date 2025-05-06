library(brms)
library(ggplot2)
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
      y_test<-test_data$y
     # y_rep_t <-posterior_predict(brm_model,ndraws=length(y_test))
      test_pit<-pit(y_rep,y_test)
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
          # simultanuous coverage:  whether all observed values fall within their respective intervals (simultaneous)
           sim_covered <- all(observed >= lower & observed <= upper)
           return(sim_covered)
          #(if  marginal coverage: Logical vector: TRUE if observed[i] is within interval [lower[i], upper[i]]
          #( inside <- observed >= lower & observed <= upper
          #return(inside)  # Logical vector))
        }


# extimate the confidence interval of the empirical coverage 
simultaneous_coverage_ci <- function(coverage_vec, alpha = 0.05) {
          # coverage_vec: logical vector of coverage statments accross simulations
          # alpha: significance level for the conf interval of the estimated proportion  
              #stopifnot(is.logical(coverage_vec))
              n <- length(coverage_vec)
              # confidence interval of the proportion ( empirical coverage) by normal approximation ( for large n: CLT)
              p_hat <- mean(coverage_vec)
              # 1-alpha/2-quantile of normal
              z <- qnorm(1 - alpha / 2)
              # Standard error
              se <- sqrt(p_hat * (1 - p_hat) / n)
              # Confidence interval
              lower <- max(0, p_hat - z * se)
              upper <- min(1, p_hat + z * se)
              return(c(emp_coverage= p_hat, lower_ic = lower, upper_ic = upper))
            }




##############  Data generating mechanism & model

# Generate synthetic Continuous data
generate_data  <- function(N = 10, seed = 978205, std = 0.3, beta1 = 1.5,resp_type="continuous") {
 
  if(resp_type=="continuous") {
  # continuous data with Gaussian noise
  set.seed(seed)
  x1 <- rnorm(N)               # Covariate
  noise <- std * rnorm(N)      # Gaussian noise
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




##############  experiments  pipeline 

# dummy data & model compilation
sim_data<-generate_data(N=2)
lin_model_fit <- brm(y ~ x1, data =sim_data,

                              prior=c(prior(normal(0,2), class="b"),
                                      prior(normal(0,0.3),class="Intercept"),
                                      prior(normal(0,0.5),class="sigma") ),
                              silent = 2, refresh = 0, cores=4)


run_experiment<-function(brm_obj=lin_model_fit, level=0.95,alpha=0.05,num_obs=c(10,20), K=4,stat_tests=uniformity_tests,key=978205,seed_test=677222) {
  
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

for (n in num_obs) {

   pit_mat<-matrix(numeric(0), ncol=n, nrow=0)
   test_pit_mat<-matrix(numeric(0), ncol=n,nrow = 0)
   loo_pit_mat<-matrix(numeric(0), ncol = n, nrow =0)
   sim_covered_pp<-numeric(K)
   sim_covered_loo<-numeric(K)
 
  test_data<-generate_data(N=n, seed=seed_test+n)         # generate test data
  
  # repeat accross datasets 
  
  for (k in 1:K) {
    
     sim_data<-generate_data(N=n,seed=seed+k+n)             # generate observed data
        
     # fit the model to the data
     fit<-update(brm_obj, newdata = sim_data)
     # Coverage
     # test-coverage state
     y_pred_k<-posterior_predict(fit)
     sim_covered_pp[k]<-coverage_state(y_pred_k,test_data$y,level=level)
     # loo-coverage state
     psis <- loo::psis(-log_lik(fit), cores = 4)
     loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
     sim_covered_loo[k]<-all(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])
     
         # Compute pit with different approaches (pit, test_pit, loo_pit)
         pit_vals<-compute_all_pit(fit,sim_data,test_data)
         # Store 
         pit_mat<-rbind(pit_mat,pit_vals$pit)
         test_pit_mat<-rbind(test_pit_mat,pit_vals$test_pit)
         loo_pit_mat<-rbind(loo_pit_mat,pit_vals$loo_pit)
         }
  #  Empirical coverage :
  emp_test_coverage <-simultaneous_coverage_ci(sim_covered_pp, alpha = alpha)
  emp_loo_coverage <-simultaneous_coverage_ci(sim_covered_loo, alpha = alpha)
   
  # collect information about pit's vals uniformity : uniformity test statistic
  # compute sample statistics to compare agains a critical value 
  stat_pit<- as.data.frame( sapply(stat_tests, function(f) apply(pit_mat, 1, f) ))
  stat_test_pit<- as.data.frame( sapply(stat_tests, function(f) apply(test_pit_mat, 1, f) ))
  stat_loo_pit<- as.data.frame( sapply(stat_tests, function(f) apply(loo_pit_mat, 1, f) ))

  # map to  the corresponding number of observations used
  pit_statistic_results[[paste0("num_obs_",n)]]<-list(pit=stat_pit,test_pit=stat_test_pit,loo_pit=stat_loo_pit)
  coverage_results[[paste0("num_obs_",n)]]<-data.frame( rbind( test_coverage = emp_test_coverage, loo_coverage = emp_loo_coverage ))
  #list(test_coverage=emp_test_coverage,loo_coverage=emp_loo_coverage)

}
  return( list(statistics=pit_statistic_results,coverage=coverage_results) )
}





##### trial

Results<-run_experiment(K=100,num_obs=c(30,60,100))

print(Results)

Rates<-rejection_rates(Results$statistics$num_obs_50$test_pit,n=50)

print(Rates)










############################################################## draft #######################################################################


# visual inspection
df=Results$statistics$num_obs_50$loo_pit

df_long <- pivot_longer(df, cols = everything(), names_to = "variable", values_to = "value")
crit_df <- data.frame( variable = names(df), critical =Rates$crit_vals)

ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 20, fill = "steelblue", color = "white") +
  facet_wrap(~ variable, scales = "free") +
    geom_vline( data = crit_df, aes(xintercept = critical, color = "Critical value"),
    linetype = "dashed" ) +
  scale_color_manual(name = "", values = c("Critical value" = "red")) +
  theme_minimal() +
  theme(legend.position = "top")








############## 


model <- function(sim_data,family="gaussian",nu=NULL,sd_gp=1,ls_gp=0.5){
  
  if(family=="gaussian") {
    
    fit<- brm(y ~ x1, data=sim_data,
        prior=c(prior(normal(0,2), class="b"),
                prior(normal(0,0.3),class="Intercept"),
                prior(normal(0,0.5),class="sigma") ),
        silent = 2, refresh = 0, cores=4)
  }
  
  if(family=="student"){
  fit <- brm(
  y ~ x1, 
    data = sim_data,
    family = student(), # control degrees of freedom 
    prior = c(
      prior(normal(0, 2), class = "b"),
      prior(normal(0, 0.3), class = "Intercept"),
      prior(normal(0, 0.5), class = "sigma"),
      prior(constant(3),class="nu")),           # control the degree of freedom 
    silent = 2, refresh = 0, cores = 4
  ) }

  if(family=='GP'){
    
 # Linear effect of x1 + non-linear GP over x1
   fit <- brm(
  y ~ x1 + gp(x1), # control linearity wiht gp hyperparams 
  data = sim_data,
  prior = c(
    prior(normal(0, 2), class = "b", coef="x1"),                   # linear effect 
    prior(normal(0, 0.3), class = "Intercept"),                    # Prior for intercept
    prior(normal(0, 0.5), class = "sigma"),                        # Prior for residual noise
    prior(normal(1, 1e-6), class = "sdgp",coef="gpx1"),           # control  GP magnitude (std dev)
    prior(normal(0.5,1e-6), class = "lscale",coef="gpx1")         #control GP lengthscale
  ),
  silent = 2, refresh = 0, cores = 4
)
}}


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



## another way to proceed ( estimate critical value)
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

