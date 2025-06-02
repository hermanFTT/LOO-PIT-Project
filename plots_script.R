

library(ggplot2)
library(tidyr)
library(glue)
source("scripts/plot_functions.R")

Ns <- c(16,32,64,128,256)
beta_s<-c(0.2,0.4,0.6,0.8,0.95)    # ( fix N,  control non linearity in GP based generating mechanism )

df_s <-c(1,5,10,20,40,80,100)       # (fixed N ,  control the degree of freedom in student t based dgp )



##################################################################################################################

# varying nu
Results <- readRDS("results/student_t_varying_nu/results-varying_nu-student_t-df_9999-beta_9999-n_256.rds")
rej_rates <-readRDS("results/student_t_varying_nu/rates-varying_nu-student_t-df_9999-beta_9999-n_256.rds")
N <- 256

param<- paste0(paste0(TeX("$n$")," = "),N )  # fixed argument ( in this case the num_obs)
p1 <- comb_rejection_plot(rej_rates,obs_sizes=NULL,bs=NULL,nu_s=df_s,title=TRUE,x_lab=TeX("$\\nu$"),row_descript=param,legend=TRUE)
rej_plot <- p1+plot_annotation(caption = glue("degree of freedom vs rejection rate. DGP: student-t, fit: normal model, N= {N} observations."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 12),plot.caption.position = "plot")  )
rej_plot
ggsave(glue("results/figures/rej-stt-df-n_{N}.pdf"), plot =rej_plot, width = 8, height = 8, dpi = 300,device = cairo_pdf)



id <- 7
# plot predictive coverage for specified the degree of freedom "nu" .
p21 <- comb_coverages_plot(Results,n=NULL,beta=NULL,nu=df_s[id],type="shaded",title=TRUE)
p22 <- comb_coverages_plot(Results,n=NULL,beta=NULL,nu=df_s[id],type="diff",title=FALSE)
p2 <- wrap_plots(p21, p22, ncol = 1)+plot_annotation(caption = glue("Predictive coverage. DGP: student-t (df = {df_s[id]}), fit: normal model, N= {N} observations."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 12),plot.caption.position = "plot")  )
p2

ggsave(glue("results/figures/cov-stt-df_{df_s[id]}-n_{N}.pdf"), plot =p2, width = 8, height = 8, dpi = 300,device = cairo_pdf)



id <- 5
# histogram plot of the averaged  pit values (specify "beta" val )
a <- comb_hist(Results,n=NULL,beta=NULL,nu=df_s[id],bins=20,title=TRUE)
b <- comb_hist(Results,n=NULL,beta=NULL,nu=df_s[id+1],bins=20,title=TRUE)
c <- wrap_plots(a, b, ncol = 1)+
plot_annotation(caption = glue("Averaged pit vals over S=200 simulations. DGP: student-t (df = {df_s[id]},{df_s[id+1]}), fit: normal model, N= {N} observations."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 11),plot.caption.position = "plot")  )

ggsave(glue("results/figures/hist-stt-df_{df_s[id]},{df_s[id+1]}-n_{N}.pdf"), plot =c, width = 8, height = 8, dpi = 300,device = cairo_pdf)
















#########################################################################################################################

#varying beta
Results <- readRDS("results/hsgp_varying_beta/results-varying_beta-hsgp-df_9999-beta_9999-n_32.rds")
rej_rates <-readRDS("results/hsgp_varying_beta/rates-varying_beta-hsgp-df_9999-beta_9999-n_32.rds")

N <- 32

param<- paste0(paste0(TeX("$n$")," = "),N )  # fixed argument ( in this case the num_obs)
p1 <- comb_rejection_plot(rej_rates,obs_sizes=NULL,bs=beta_s,nu_s=NULL,title=TRUE,x_lab=TeX("$\\beta$"),row_descript=param,legend=TRUE)
rej_plot <- p1+plot_annotation(caption = glue("beta vs rejection rate. DGP: beta*x+(1-beta)*gp (hilbert space gp), fit: normal model, N= {N} observations."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 12),plot.caption.position = "plot")  )
rej_plot
ggsave(glue("results/figures/rej-hsgp-beta-n_{N}.pdf"), plot =rej_plot, width = 8, height = 8, dpi = 300,device = cairo_pdf)



id <- 5
# plot predictive coverage for specified the degree of freedom "nu" .
p21 <- comb_coverages_plot(Results,n=NULL,beta=beta_s[id],nu=NULL,type="shaded",title=TRUE)
p22 <- comb_coverages_plot(Results,n=NULL,beta=beta_s[id],nu=NULL,type="diff",title=FALSE)
p2 <- wrap_plots(p21, p22, ncol = 1)+plot_annotation(caption = glue("Predictive coverage. DGP:  beta*x+(1-beta)*gp (hilbert space gp), (beta = {beta_s[id]}), fit: normal model, N= {N} observations."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 12),plot.caption.position = "plot")  )
p2

ggsave(glue("results/figures/cov-hsgp-beta_{beta_s[id]}-n_{N}.pdf"), plot =p2, width = 8, height = 8, dpi = 300,device = cairo_pdf)



id <- 1
# histogram plot of the averaged  pit values (specify "beta" val )
a <- comb_hist(Results,n=NULL,beta=beta_s[id],nu=NULL,bins=20,title=TRUE)
b <- comb_hist(Results,n=NULL,beta=beta_s[id+1],nu=NULL,bins=20,title=TRUE)
c <- wrap_plots(a, b, ncol = 1)+
plot_annotation(caption = glue("Averaged pit vals over S=200 simulations. DGP:  beta*x+(1-beta)*gp (hilbert space gp), (beta = {beta_s[id]},{beta_s[id+1]}), fit: normal model, N= {N} observations."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 11),plot.caption.position = "plot")  )

ggsave(glue("results/figures/hist-hsgp-beta_{beta_s[id]},{beta_s[id+1]}-n_{N}.pdf"), plot =c, width = 8, height = 8, dpi = 300,device = cairo_pdf)







###########################################################################################################

## varying N
Results <- readRDS("results/student_t_varying_N/results-varying_N-student_t-df_2-beta_9999-n_9999.rds")
rej_rates <-readRDS("results/student_t_varying_N/rates-varying_N-student_t-df_2-beta_9999-n_9999.rds")
df <- 2
# plot obs_sizes vs rejection rate
param <- paste0(paste0(TeX("$\\nu$")," = "),df)
p1 <- comb_rejection_plot(rej_rates,obs_sizes=Ns,bs=NULL,nu_s=NULL,title=TRUE,x_lab="#Observations",row_descript=param,legend=TRUE)

rej_plot <- p1+plot_annotation(caption = glue("#Observations  vs rejection rate. DGP: student-t, fit: normal model, df= {df}."),theme = theme(  plot.caption = element_text( hjust = 0.5,size = 12),plot.caption.position = "plot")  )
rej_plot
ggsave(glue("results/figures/rej-stt-n-df_{df}.pdf"), plot =rej_plot, width = 8, height = 8, dpi = 300,device = cairo_pdf)



