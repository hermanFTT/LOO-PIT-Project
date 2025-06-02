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
options(brms.backend = "cmdstanr", mc.cores = 4)
#install.packages("modeest")


source("scripts/Uniformity_tests.R")
source("scripts/help_functions.R")
source("scripts/generators.R")



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
      psis_2 <- loo(brm_model, save_psis = TRUE, cores = 1)$psis_object
      lw<-weights(psis_2)      # normalized  log weighs from psis scheme
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




rejection_rates<-function(Results,obs_sizes,stat_tests=uniformity_tests,alpha=0.05,n_iter=10000) {

  # rstat_list (list) : list of list of df containing realized statistics for diff test
  # obs_sizes (vector): vector containing diff number of observations size.
  rstat_list <- Results$statistics
  rej_rates <- list()

  for (key in c("pit","test_pit","loo_pit")) {
    rate_mat<-matrix(numeric(0), ncol=length(stat_tests), nrow=0)
    
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



##############  Data generating mechanisms & model

# Generate synthetic  data
simul_stt<- function(N = 10,beta=0.4,seed=111,df=100) {
      # continuous data with student_t noise
      set.seed(seed)
      x1 <- rnorm(N) # Covariate
      # student_t noise
      y <- beta * x1 + (1-beta)*rt(N,df=df)     # Response variable
      sim_data <- data.frame(x1 = x1, y = y)
      return(sim_data)
    }


simul_count <- function(N=10, beta=0.8,seed=222,df=100) {
  # count data from poisson dist
        set.seed(seed)
         x1 <- rnorm(N)  
        lambda<-exp(beta*x1)
        y<-rpois(N,lambda)
        sim_data <- data.frame(x1 = x1, y=y)
        return(sim_data)
      }

# simulate data with GP models ( Hilbert space based GP approximation)
# load the stan file
hsgp_model<- cmdstan_model("stan_file/HSGP_simulate.stan")
# writeLines(readLines(hsgp_model))

simul_hsgp <- function(N=50,beta=0.95,seed=333,model=hsgp_model,ls=0.3,sig=1,noise=0.3,M_f=160,c_f=3,inter_f=0.4,df=100){

      # x1=seq(bound[1],bound[2],length.out=N) # univariate input bound=c(0,1),
      
        set.seed(seed)
        x1 <- rnorm(N)
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

  # simulate hsgp
      sim_hsgp_model<-model$sample(data=data_list, fixed_param=TRUE,
                                chains=1, iter=1, iter_sampling=1,show_messages=FALSE)
  ## data generating mechanism
  gp<-as.vector( sim_hsgp_model$draws(format="matrix", variable="y"))
    # linear term + non linearity correction 
  y_sim <- beta*x1+(1-beta)*gp
  sim_data <- data.frame(x1=x1,y=y_sim)

 return(sim_data)
}


  # visualize
#sim_data <- simul_hsgp(seed=4,N=32,beta=0.3) #0.2-0.5
#plot(sim_data$x1,sim_data$y)





##############  experiments  pipeline 

## dummy data for model compilation

sim_data<-simul_stt(seed=4,N=128,beta=0.2,df=100) # linear + student_t
#sim_data<- simul_hsgp(seed=40,N=10,beta=0.2)  # linear + hsgp
#hist(sim_data$y,breaks=40) # visualize

##  choose model: linear model with Gaussian noie

lin_model_fit <- brm(y ~ x1, data =sim_data,

                              prior=c(prior(uniform(0,1), class="b"),
                                      prior(normal(0,0.2),class="Intercept"),
                                      prior(normal(0,0.5),class="sigma") ),
                     silent = 2, refresh = 0, cores=4)




   # inspection 
#plot(lin_model_fit)         
pp_check(lin_model_fit)
yrep <- posterior_predict(lin_model_fit)

sim_data_test <- simul_hsgp(seed=64,N=50,beta=0.3)
yrep_test <- posterior_predict(lin_model_fit,newdata =sim_data_test[1])
psis1 <- loo(lin_model_fit, save_psis = TRUE, cores = 4)$psis_object
lw <- weights(psis1)
ppc_loo_intervals(sim_data$y, yrep, psis_object = psis1,prob=0.5)


test_interval <- predictive_interval(yrep_test,prob=0.05) #
mean(sim_data_test$y >=test_interval[, 1] & sim_data_test$y <= test_interval[, 2])

loo_interval <- loo_predictive_interval(lin_model_fit, prob = 0.05, psis_object = psis1)
mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])


yrep <- posterior_predict(lin_model_fit)
#pit
s_pit <- pit(yrep,sim_data$y)
hist(s_pit,breaks=20)
# loo-pit
psis <- loo(lin_model_fit, save_psis = TRUE, cores = 1)$psis_object
lw<-weights(psis)      # normalized  log weighs from psis scheme
loo_pit<-pit(yrep,sim_data$y,weights=lw,log = TRUE)
x11()
hist(loo_pit,breaks=20)

s_pit==loo_pit


pit_loo<- function(y, yrep, lw) {
  pit <- vapply(seq_len(ncol(yrep)), function(j) {
    sel_min <- yrep[, j] < y[j]
    pit_min <- exp_log_sum_exp(lw[sel_min,j])
    sel_sup <- yrep[, j] == y[j]
    pit_sup <- pit_min + exp_log_sum_exp(lw[sel_sup,j])
    runif(1, pit_min, pit_sup)
  }, FUN.VALUE = 1)
  pmax(pmin(pit, 1), 0)
}
exp_log_sum_exp <- function(x) {
  m <- suppressWarnings(max(x))
  exp(m + log(sum(exp(x - m))))
}


pit_vals <- compute_all_pit(lin_model_fit,sim_data,test_data)






########################## experiment pipeline


model <- lin_model_fit # brm_obj
dgp <- simul_hsgp #generator

beta <- 0.2 # beta

df <- 100 #nu

levels <- seq(0.05,0.95,by=0.025) # levels

Ns <- c(16,32,64,128,256,512) # num_obs

S <- 20 # K  

key <- 444 # seed
t_key <- 555 # test_seed

Results<-run_experiment(brm_obj=model,generator=dgp,b=beta,nu=df,levels=levels,num_obs=Ns,K=S,seed=key,seed_test=t_key)









run_experiment<-function(brm_obj,generator,b=0.2,nu=100,levels=c(0.9),alpha=0.05,num_obs=c(16,32),K=4,stat_tests=uniformity_tests,seed=444,seed_test=555,summary_stats=list(mean=mean,var=var)) {
  
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
    psis <-loo(fit, save_psis = TRUE, cores = 4)$psis_object
    
     f_vec2 <- Vectorize( function(level) {
      loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     covered_loo[k,]<-f_vec2(levels)


   # covered_loo[k,]<-sapply(levels,function(level) {
    #   psis <-loo(fit, save_psis = TRUE, cores = 4)$psis_object                                  #loo::psis(-log_lik(fit), cores = 4)
     # loo_interval=loo_predictive_interval(fit, prob = level, psis_object = psis)
      #mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])})
     
         # Compute pit with different approaches (pit, test_pit, loo_pit)
         pit_vals<-compute_all_pit(fit,sim_data,test_data)
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
  coverage_results[[paste0("num_obs_",n)]]<-list(test_coverage=emp_test_coverage ,loo_coverage=emp_loo_coverage) 
  #data.frame( rbind( test_coverage = emp_test_coverage, loo_coverage = emp_loo_coverage ))
 # fitted_model[[paste0("num_obs_",n)]] <- fit
}
  return( list(pvalues=pvalues_results,statistics=pit_statistic_results,coverage=coverage_results,summary=pit_summary_results) )
}





##### trial
levels <- seq(0.05,0.95,by=0.025)

Results<-run_experiment(generator=simul_hsgp,b=0.2,nu=100,levels=c(0.9),alpha=0.05,num_obs=c(10,20), K=4 ) # seq(0.05,0.95,by=0.05) 
#save(Results, file = "results_studentt_df_20.RData")
#Results <- load("results_studentt_df_100.RData")

print(Results)

Rates<-rejection_rates(Results,Ns)

print(Rates)
#########################


 
library(latex2exp)
library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(grid)
#install.packages("patchwork")


rejection_plot<-function(Rates_list=Rates,obs_sizes=c(16,32),test="ks_stat",title="KS-test",x_lab=TRUE,y_lab=TRUE,legend=FALSE){
pit1_df <- Rates_list$pit
pit2_df <- Rates_list$test_pit
pit3_df <- Rates_list$loo_pit
# Add PIT method names
pit1_df$Method <-"pit"
pit2_df$Method <-"test_pit"
pit3_df$Method <-"loo_pit"

# Add num_obs as a column (assuming rownames are num_obs)
pit1_df$num_obs <- obs_sizes
pit2_df$num_obs <- obs_sizes
pit3_df$num_obs <- obs_sizes

# Combine all data
full_df <- bind_rows(pit1_df, pit2_df, pit3_df)

# Reshape to long format
long_df <- pivot_longer(full_df, cols = ends_with("stat"), names_to = "Test", values_to="RejectionRate") 

test_names <- unique(long_df$Test)

# Generate a list of plots, one per test

if (test  %in% as.vector(test_names)) {

   p<- ggplot(filter(long_df, Test == test), aes(x = num_obs, y = RejectionRate, color = Method)) +
    geom_line(size = 1,linetype = "dashed") +
    geom_point(size = 2) +
    # scale_y_continuous( limits = c(l,u), breaks = seq(l, u, by = 16))+
    scale_y_continuous( limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
   labs(x=NULL,y =NULL) +
     theme( text = element_text(size = 16),axis.title.x = element_text(size = 16 ), axis.title.y = element_text(size=16) )+
     theme(legend.position = "none" )+
     theme_minimal()  

  if(!is.null(title)) { p <- p+labs(title = title)+theme( plot.title = element_text( size = 18,face = "bold", hjust = 0.5 ) )}
  
  if(x_lab==TRUE) { p <- p +labs(x="Number of observations")}
  if(y_lab==TRUE) { p <- p +labs(y="Rej. rate")}
 if(legend==TRUE) { p <- p+theme(legend.position = "bottom")} 
}
return(p)

}

comb_rejection_plot <- function(Rates_list=Rates,obs_sizes=c(16,32),title=TRUE,x_lab=TRUE,row_descript=NULL,legend=TRUE){
  if(title==TRUE) {
  p1 <- rejection_plot(Rates_list,obs_sizes,test="ks_stat",title="KS-test",x_lab=FALSE,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,obs_sizes,test="t1_stat",title="T1-test",x_lab=FALSE,y_lab=FALSE)
  if(x_lab==TRUE){
  p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title="U2-test",x_lab=TRUE,y_lab=FALSE)} else{
  p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title="U2-test",x_lab=FALSE,y_lab=FALSE)} }

  else {
  p1 <- rejection_plot(Rates_list,obs_sizes,test="ks_stat",title=NULL,x_lab=FALSE,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,obs_sizes,test="t1_stat",title=NULL,x_lab=FALSE,y_lab=FALSE)
  if(x_lab==TRUE){ p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title="U2-test",x_lab=TRUE,y_lab=FALSE)}else{
   p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title=NULL,x_lab=FALSE,y_lab=FALSE)}}

  if(legend==TRUE){
    combined_plot <- p1+p2+p3+plot_layout(ncol=3,widths = c(4,4,4),heights = c(4,4,4),guides = "collect")&theme(legend.position = "bottom")} else{
    combined_plot <-  p1+p2+p3+plot_layout(ncol = 3,widths = c(4,4,4),heights = c(4,4,4),guides="collect")&theme(legend.position = "none",plot.margin = margin(0, 5, 0, 5) )}

      if(!is.null(row_descript)){
        label1 <- wrap_elements(grid::textGrob(row_descript, rot =360,just = "center",y=0.85))
        combined_plot <- label1+combined_plot+plot_layout(widths = c(0.15, 0.85)) }

      return(combined_plot)}

  
 comb_rejection_plot(obs_sizes = Ns,legend=TRUE,x_lab=TRUE,row_descript = TeX("$\\beta=0.2$") )
  



#######################
#ggsave("trial_plot.pdf", plot = final_plot, width = 7, height = 6, dpi = 300,device = cairo_pdf)
#########################################





coverage_plot<-function(Results,n=512,target="loo_coverage",type="errorbar"){
  # n (int): observation size
  # exp_results (list) : list of list (results from experiment)
  # type (string): possible vals c("errorbar","shaded","diff") 

coverage_df<- Results$coverage[[paste0("num_obs_",n)]][[target]]



if(type=="shaded") {

p <- ggplot(coverage_df, aes(x = level * 100, y =coverage * 100)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100), fill = "gray80", alpha = 0.5) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue",size=3,alpha=0.4) +
  scale_x_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  scale_y_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  labs(x = "central interval width", y = "Observed coverage") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    panel.grid.minor = element_blank() )}
  
  if(type=="errorbar"){

p <- ggplot(coverage_df, aes(x = level * 100, y =coverage * 100)) +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),width = 1,  color = "gray80",size=1) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue",size=3,alpha=0.4) +
  scale_x_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  scale_y_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  labs(x = "Central interval width", y = "Observed coverage") +
  theme_minimal() +
  theme(
    text = element_text(size = 16),
    panel.grid.minor = element_blank()
  )
 }

  if (type=="diff") {


    l <- min((coverage_df$lower-coverage_df$level)*100)
    u <- max((coverage_df$upper-coverage_df$level)*100)
    
p <- ggplot(coverage_df, aes(x = level * 100, y =(coverage-level)* 100)) +
  geom_ribbon(aes(ymin =(lower-level)*100, ymax =(upper-level)*100), fill = "gray80", alpha = 0.4) +
  geom_line(color = "black", size = 1) +
  geom_point(color = "black", size = 1) +
  geom_abline(slope = 0, intercept = 0, linetype = "solid", color = "blue",size=3,alpha=0.4) +
  scale_x_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  scale_y_continuous(limits=c(l,u), labels = percent_format(scale = 1)) +
  labs(x = "central interval width", y = "Observed coverage") +
  theme_minimal() +
  plot_layout(widths = 4, heights = 4)+
  theme(
    text = element_text(size = 16),
    panel.grid.minor = element_blank() ) }

  return(p)
 }



dual_coverages_plot <- function(Results,n=16,type="errorbar",title=TRUE){
  # Results (list): experiments results from "run_experiment()"
  # n (int): observation size
  # type ( string): plot style of credible interval of the coverage

    p1 <- coverage_plot(Results,n,"loo_coverage",type=type)
    p2 <-coverage_plot(Results,n,"test_coverage",type=type)


    label1 <- wrap_elements(grid::textGrob(TeX(paste0("$n=$",n)), rot =360,just = "center",y=0.75,  gp = gpar(fontsize = 16, fontface = "bold")))
    combined_plot <- p1 +plot_spacer()+p2+plot_layout(ncol=3,widths = c(4,0.2,4),heights = c(4,0.1,4))

    if(title==TRUE){
      col1 <- wrap_elements(grid::textGrob("LOO-coverage",just = "center",x=0.18, gp = gpar(fontsize = 16, fontface = "bold")))
      col2 <- wrap_elements(grid::textGrob("test-coverage",just = "center", x=0.37, gp = gpar(fontsize = 16, fontface = "bold")))


    final_plot<-(plot_spacer()+col1+col2)/(label1+combined_plot+plot_layout(widths = c(0.15, 0.85)))+ plot_layout(widths =c(0.1,1,1), heights = c(0.1,1), guides = "collect")

    } else { final_plot <-  label1+combined_plot+plot_layout(widths = c(0.15, 0.85)) }

    return( final_plot)
  }




#
dual_coverages_plot(Results,n=32,type="errorbar")

length(levels)










dual_coverages_plot(Results,type="shaded")


#######################
ggsave("trial_plot_2.pdf", plot = final_plot, width =13, height = 7, dpi = 300,device = cairo_pdf)
print(final_plot)
##################

########################################################### draft #######################################################################





##############



   # inspection 
#plot(lin_model_fit)         
pp_check(lin_model_fit)
yrep <- posterior_predict(lin_model_fit)

sim_data_test <- simul_hsgp(seed=64,N=50,beta=0.3)
yrep_test <- posterior_predict(lin_model_fit,newdata =sim_data_test[1])
psis1 <- loo(lin_model_fit, save_psis = TRUE, cores = 4)$psis_object
lw <- weights(psis1)
ppc_loo_intervals(sim_data$y, yrep, psis_object = psis1,prob=0.5)


test_interval <- predictive_interval(yrep_test,prob=0.05) #
mean(sim_data_test$y >=test_interval[, 1] & sim_data_test$y <= test_interval[, 2])

loo_interval <- loo_predictive_interval(lin_model_fit, prob = 0.05, psis_object = psis1)
mean(sim_data$y >=loo_interval[, 1] & sim_data$y <= loo_interval[, 2])
   



################

final_plot <- wrap_plots(plot_list, nrow = 1)
print(final_plot)
#########
df=Results$pvalues$num_obs_32$pit


 df_long <- df %>%
  pivot_longer(cols = everything(), names_to = "variable", values_to = "value")

ggplot(df_long, aes(x = value)) +
  geom_histogram(bins = 30, fill = "steelblue", color = "white") +
  facet_wrap(~ variable, scales = "free") +
  theme_bw()


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

#linear + gp
# 
