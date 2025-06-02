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
#install.packages("modeest")
#install.packages("tictoc")
options(brms.backend = "cmdstanr", mc.cores =4)




source("scripts/Uniformity_tests.R")
source("scripts/help_functions.R")
source("scripts/generators.R")
source("scripts/experiments_pipline.R")
source("scripts/plot_functions.R")



## dummy data for model compilation

#sim_data <- simul_stt(seed=41,N=10, df=100)   # student t familly :  
sim_data<- simul_hsgp(seed=48,N=32,beta=0.6)  # linear_trend + GP : beta*x+(1-beta)*gp 
  # visualize
#plot(sim_data$x1,sim_data$y)
#hist(sim_data$y,breaks=40) # visualize


##  choose model: linear model with Gaussian noie
lin_model_fit <- brm(y ~ x1, data =sim_data,
                              prior=c(prior(uniform(0,1), class="b"),
                                      prior(normal(0,0.2),class="Intercept"),
                                      prior(normal(0,0.5),class="sigma") ),
                     silent = 2, refresh = 0, cores=4)


########################## run experiment pipeline


style <- as.character(snakemake@wildcards[["style"]])
dgp <- as.character(snakemake@wildcards[["dgp"]])
N <- as.numeric(snakemake@wildcards[["n"]])
df <- as.numeric(snakemake@wildcards[["df"]])
beta <- as.numeric(snakemake@wildcards[["beta"]])

S <-as.numeric( snakemake@params[["nsims"]])
bins <- as.numeric(snakemake@params[["nbins"]])


file_1= snakemake@output[["output_file"]]

file_r1= snakemake@output[["rates_file"]]



#NB : If you not using Snakemake, please  comment out the Snakemake block above and assign values to the variables in the lines below (which are commented).



model <- lin_model_fit   # brm_obj



#S <-200                      # K  ( number of simulations ) 

key <- 444                        # seed  ( key for reproducible data )
t_key <- 555                      # test_seed
bins <- 15                        # number of bins for the pit histogram


Ns <- c(16,32,64,128,256)      # num_obs ( different observations sizes )

#beta <- 0.95             # [0,1]  (fix amount of  linearity in GP based dgp, varying N ) 

#df <- 100                # ( fix  the degree of freedom in student t based dgp, varying N )





#N <- 128                      # fix observation size

beta_s<-c(0.2,0.4,0.6,0.8,0.95)    # ( fix N,  control non linearity in GP based generating mechanism )

df_s <-c(1,5,10,20,40,80,100)       # (fixed N ,  control the degree of freedom in student t based dgp )

levels <- seq(0.05,0.95,by=0.025)  # levels ( widths of central intervals for predictive coverage)




#style <-"varying_beta"    # c("varying_N","varying_beta","varying_nu")

#dgp <-"hsgp"              # c("hsgp","student_t")    # data generating mechanism 





if(style=="varying_N") {

    if(dgp=="student_t") {
     # fixed degree of freedom ( nu ) varying num_obs.
      Results<-run_experiment_n(brm_obj=model,generator=simul_stt, b=0.2,nu=df,levels=levels,num_obs=Ns,K=S,seed=key,seed_test=t_key)
      param <- paste0(paste0(TeX("$\\nu$")," = "), df)
      
      } else if(dgp=="hsgp"){
        Results<-run_experiment_n(brm_obj=model,generator=simul_hsgp,b=beta,nu=100,levels=levels,num_obs=Ns,K=S,seed=key,seed_test=t_key)
        param <-paste0(paste0(TeX("$\\beta$")," = "), beta)
      }

  rej_rates <- rejection_rates(Results,obs_sizes=Ns)

  # plot obs_sizes vs rejection rate
  p1 <- comb_rejection_plot(rej_rates,obs_sizes=Ns,bs=NULL,nu_s=NULL,title=TRUE,x_lab="#Observations",row_descript=param,legend=TRUE)

  # plot predictive coverage for specified number of observations n in "Ns" . 
  p2 <- comb_coverages_plot(Results,n=Ns[1],beta=NULL,nu=NULL,type="errorbar",title=TRUE)
  
  # histogram plot of the averaged  pit values
  p3 <- comb_hist(Results,n=Ns[2],beta=NULL,nu=NULL,bins=bins,title=TRUE)
 }




if(style=="varying_beta") {

   ## fixed N (num_obs) control non linearity  with beta ( in hsgp based dgp ) : 


   Results<-run_experiment_b(brm_obj=model,b=beta_s,nu=100,levels=levels,n=N,K=S,seed=key,seed_test=t_key)
   rej_rates <- rejection_rates(Results,bs=beta_s,N=N)

   # plot beta vs rejection rate
   param <- paste0(paste0(TeX("$n$")," = "), N)
   p1 <- comb_rejection_plot(rej_rates,obs_sizes=NULL,bs=beta_s,nu_s=NULL,title=TRUE,x_lab=TeX("$\\beta$"),row_descript=param,legend=TRUE)

   # plot predictive coverage for specified "beta" in "beta_s" . 
  p2 <- comb_coverages_plot(Results,n=NULL,beta=beta_s[5],nu=NULL,type="diff",title=TRUE)

  # histogram plot of the averaged  pit values (specify "beta" val )
  p3 <- comb_hist(Results,n=NULL,beta=beta_s[1],nu=NULL,bins=bins,title=TRUE)

}



if(style=="varying_nu") {

    ### fixed N varying degree of freedom ( nu)  in student t based dgp.

   Results<-run_experiment_nu(brm_obj=model,b=0.2,nu=df_s,levels=levels,n=N,K=S,seed=key,seed_test=t_key)
   rej_rates <- rejection_rates(Results,nu_s=df_s,N=N)

   # Plot df vs rej_rate
   param<-paste0(paste0(TeX("$n$")," = "), N)
   p1<- comb_rejection_plot(rej_rates,obs_sizes=NULL,bs=NULL,nu_s=df_s,title=TRUE,x_lab=TeX("$\\nu$"),row_descript=param,legend=TRUE)

   # plot predictive coverage for specify degree of freedom "nu" in "df_s"
   p2<- comb_coverages_plot(Results,n=NULL,beta=NULL,nu=df_s[1],type="errorbar",title=TRUE)

    
  # histogram plot of the averaged  pit values ( specify "nu" )
  p3 <- comb_hist(Results,n=NULL,beta=NULL,nu=df_s[1],bins=bins,title=TRUE)

 }


# output files of snakemake
saveRDS(Results, file = file_1)
saveRDS(rej_rates,file=file_r1)



####

#p1
#x11()
#p2
#x11()
#p3

