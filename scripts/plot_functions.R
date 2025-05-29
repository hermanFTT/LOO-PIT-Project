library(latex2exp)
library(scales)
library(ggplot2)
library(tidyr)
library(dplyr)
library(patchwork)
library(grid)
#install.packages("patchwork")


rejection_plot<-function(Rates_list,obs_sizes=NULL,bs=NULL,nu_s=NULL,test="ks_stat",title="KS-test",x_lab=NULL,y_lab=TRUE,legend=FALSE){
pit1_df <- Rates_list$pit
pit2_df <- Rates_list$test_pit
pit3_df <- Rates_list$loo_pit
# Add PIT method names
pit1_df$Method <-"pit"
pit2_df$Method <-"test-pit"
pit3_df$Method <-"loo-pit"

if(!is.null(obs_sizes)) {
# Add num_obs as a column 
pit1_df$num_obs <- obs_sizes
pit2_df$num_obs <- obs_sizes
pit3_df$num_obs <- obs_sizes}

if(!is.null(bs)) {
# Add diff beta's vals  as a column 
pit1_df$num_obs <- bs
pit2_df$num_obs <- bs
pit3_df$num_obs <- bs }

if(!is.null(nu_s)) {
# Add diff degree of freedom  as a column 
pit1_df$num_obs <- nu_s
pit2_df$num_obs <- nu_s
pit3_df$num_obs <- nu_s }


# Combine all data
full_df <- bind_rows(pit1_df, pit2_df, pit3_df)

# Reshape to long format
long_df <- pivot_longer(full_df, cols = ends_with("stat"), names_to = "Test", values_to="RejectionRate") 

test_names <- unique(long_df$Test)

# Generate a list of plots, one per test

if (test  %in% as.vector(test_names)) {

   p<- ggplot(filter(long_df, Test == test), aes(x = num_obs, y = RejectionRate, color = Method,linetype=Method)) +
     geom_line(size = 1) +
     geom_point(size = 0.7) +
     scale_color_manual(values = c(
    "pit" = "#0072B2", 
    "test-pit" =  "#4CAF50",  
    "loo-pit" =  "#F94144"  
  ))+
  scale_linetype_manual(values = c(
    "test-pit" = "solid",
    "loo-pit" = "dashed",
    "pit" = "dotdash"
  ))+
    # scale_y_continuous( limits = c(l,u), breaks = seq(l, u, by = 16))+
   scale_y_continuous( limits = c(0, 1), breaks = seq(0, 1, by = 0.2))+
   labs(x=NULL,y =NULL) +
     theme( text = element_text(size = 11),axis.title.x = element_text(size = 12 ), axis.title.y = element_text(size=11) )+
     theme(legend.position = "none" )+
     theme_minimal()  

  if(!is.null(title)) { p <- p+labs(title = title)+theme( plot.title = element_text( size = 12,face = "bold", hjust = 0.5 ) )}
  
  if(!is.null(x_lab)) { p <- p +labs(x=x_lab)}
  if(y_lab==TRUE) { p <- p +labs(y="Rej. rate")}
 if(legend==TRUE) { p <- p+theme(legend.position = "bottom")} 
}
return(p)

}

comb_rejection_plot <- function(Rates_list,obs_sizes=NULL,bs=NULL,nu_s=NULL,title=TRUE,x_lab=NULL,row_descript=NULL,legend=TRUE){

  if(!is.null(obs_sizes)) {

  if(title==TRUE) {
  p1 <- rejection_plot(Rates_list,obs_sizes,test="ks_stat",title="KS-test",x_lab=NULL,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,obs_sizes,test="t1_stat",title="T1-test",x_lab=NULL,y_lab=FALSE)
  if(!is.null(x_lab)){
  p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title="U2-test",x_lab=x_lab,y_lab=FALSE)} else{
  p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title="U2-test",x_lab=NULL,y_lab=FALSE)} }

  else {
  p1 <- rejection_plot(Rates_list,obs_sizes,test="ks_stat",title=NULL,x_lab=NULL,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,obs_sizes,test="t1_stat",title=NULL,x_lab=NULL,y_lab=FALSE)
  if(!is.null(x_lab)){ p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title="U2-test",x_lab=x_lab,y_lab=FALSE)}else{
   p2 <- rejection_plot(Rates_list,obs_sizes,test="u2_stat",title=NULL,x_lab=NULL,y_lab=FALSE)}}

  if(legend==TRUE){
    combined_plot <- p1+p2+p3+plot_layout(ncol=3,widths = c(4,4,4),heights = c(4,4,4),guides = "collect")&theme(legend.position = "bottom")} else{
    combined_plot <-  p1+p2+p3+plot_layout(ncol = 3,widths = c(4,4,4),heights = c(4,4,4),guides="collect")&theme(legend.position = "none",plot.margin = margin(0, 5, 0, 5) )}

      if(!is.null(row_descript)){
        label1 <- wrap_elements(grid::textGrob(row_descript, rot =360,just = "center",y=0.85))
        combined_plot <- label1+combined_plot+plot_layout(widths = c(0.15, 0.85)) } }



    if(!is.null(bs)) {

  if(title==TRUE) {
  p1 <- rejection_plot(Rates_list,bs=bs,test="ks_stat",title="KS-test",x_lab=NULL,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,bs=bs,test="t1_stat",title="T1-test",x_lab=NULL,y_lab=FALSE)
  if(!is.null(x_lab)){
  p2 <- rejection_plot(Rates_list,bs=bs,test="u2_stat",title="U2-test",x_lab=x_lab,y_lab=FALSE)} else{
  p2 <- rejection_plot(Rates_list,bs=bs,test="u2_stat",title="U2-test",x_lab=NULL,y_lab=FALSE)} }

  else {
  p1 <- rejection_plot(Rates_list,bs=bs,test="ks_stat",title=NULL,x_lab=NULL,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,bs=bs,test="t1_stat",title=NULL,x_lab=NULL,y_lab=FALSE)
  if(!is.null(x_lab)){ p2 <- rejection_plot(Rates_list,bs=bs,test="u2_stat",title="U2-test",x_lab=x_lab,y_lab=FALSE)}else{
   p2 <- rejection_plot(Rates_list,bs=bs,test="u2_stat",title=NULL,x_lab=NULL,y_lab=FALSE)}}

  if(legend==TRUE){
    combined_plot <- p1+p2+p3+plot_layout(ncol=3,widths = c(4,4,4),heights = c(4,4,4),guides = "collect")&theme(legend.position = "bottom")} else{
    combined_plot <-  p1+p2+p3+plot_layout(ncol = 3,widths = c(4,4,4),heights = c(4,4,4),guides="collect")&theme(legend.position = "none",plot.margin = margin(0, 5, 0, 5) )}

      if(!is.null(row_descript)){
        label1 <- wrap_elements(grid::textGrob(row_descript, rot =360,just = "center",y=0.85))
        combined_plot <- label1+combined_plot+plot_layout(widths = c(0.15, 0.85)) } }


    if(!is.null(nu_s)) {

  if(title==TRUE) {
  p1 <- rejection_plot(Rates_list,nu_s=nu_s,test="ks_stat",title="KS-test",x_lab=NULL,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,nu_s=nu_s,test="t1_stat",title="T1-test",x_lab=NULL,y_lab=FALSE)
  if(!is.null(x_lab)){
  p2 <- rejection_plot(Rates_list,nu_s=nu_s,test="u2_stat",title="U2-test",x_lab=x_lab,y_lab=FALSE)} else{
  p2 <- rejection_plot(Rates_list,nu_s=nu_s,test="u2_stat",title="U2-test",x_lab=NULL,y_lab=FALSE)} }

  else {
  p1 <- rejection_plot(Rates_list,nu_s=nu_s,test="ks_stat",title=NULL,x_lab=NULL,y_lab=TRUE)
  p3 <- rejection_plot(Rates_list,nu_s=nu_s,test="t1_stat",title=NULL,x_lab=NULL,y_lab=FALSE)
  if(!is.null(x_lab)){ p2 <- rejection_plot(Rates_list,nu_s=nu_s,test="u2_stat",title="U2-test",x_lab=x_lab,y_lab=FALSE)}else{
   p2 <- rejection_plot(Rates_list,nu_s=nu_s,test="u2_stat",title=NULL,x_lab=NULL,y_lab=FALSE)}}

  if(legend==TRUE){
    combined_plot <- p1+p2+p3+plot_layout(ncol=3,widths = c(4,4,4),heights = c(4,4,4),guides = "collect")&theme(legend.position = "bottom")} else{
    combined_plot <-  p1+p2+p3+plot_layout(ncol = 3,widths = c(4,4,4),heights = c(4,4,4),guides="collect")&theme(legend.position = "none",plot.margin = margin(0, 5, 0, 5) )}

      if(!is.null(row_descript)){
        label1 <- wrap_elements(grid::textGrob(row_descript, rot =360,just = "center",y=0.85))
        combined_plot <- label1+combined_plot+plot_layout(widths = c(0.15, 0.85)) } }
  
      return(combined_plot)}







# comb_rejection_plot(legend=TRUE,x_lab=TRUE,row_descript = TeX("$\\beta=0.2$") )
  



#######################
#ggsave("trial_plot.pdf", plot = final_plot, width = 7, height = 6, dpi = 300,device = cairo_pdf)
#########################################





coverage_plot<-function(Results,n=NULL,beta=NULL,nu=NULL,target="loo_coverage",type="errorbar"){
  # n (int): observation size
  # exp_results (list) : list of list (results from experiment)
  # type (string): possible vals c("errorbar","shaded","diff") 
if(!is.null(n)) {
  coverage_df<- Results$coverage[[paste0("num_obs_",n)]][[target]]}

if(!is.null(beta)) {
  coverage_df<- Results$coverage[[paste0("beta_",beta)]][[target]] }

if(!is.null(nu)) {
  coverage_df<- Results$coverage[[paste0("df_",nu)]][[target]] }

  



if(type=="shaded") {

p <- ggplot(coverage_df, aes(x = level * 100, y =coverage * 100)) +
  geom_ribbon(aes(ymin = lower*100, ymax = upper*100), fill = "gray80", alpha = 0.6) +
  geom_line(color = "black", size = 0.8) +
  geom_point(color = "black", size = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue",size=2,alpha=0.4) +
  scale_x_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  scale_y_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  labs(x = "central interval width", y = "Observed coverage") +
  theme_minimal() +
  theme(
    text = element_text(size = 11),
    panel.grid.minor = element_blank() )}
  
  if(type=="errorbar"){

p <- ggplot(coverage_df, aes(x = level * 100, y =coverage * 100)) +
  geom_errorbar(aes(ymin = lower * 100, ymax = upper * 100),width = 1,  color = "gray80",size=1,alpha=0.6) +
  geom_line(color = "black", size = 0.8) +
  geom_point(color = "black", size = 0.4) +
  geom_abline(slope = 1, intercept = 0, linetype = "solid", color = "blue",size=2,alpha=0.4) +
  scale_x_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  scale_y_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  labs(x = "Central interval width", y = "Observed coverage") +
  theme_minimal() +
  theme(
    text = element_text(size = 11),
    panel.grid.minor = element_blank()
  )
 }

  if (type=="diff") {


    l <- min((coverage_df$lower-coverage_df$level)*100)
    u <- max((coverage_df$upper-coverage_df$level)*100)
    
p <- ggplot(coverage_df, aes(x = level * 100, y =(coverage-level)* 100)) +
  geom_ribbon(aes(ymin =(lower-level)*100, ymax =(upper-level)*100), fill = "gray80", alpha = 0.6) +
  geom_line(color = "black", size = 0.8) +
  geom_point(color = "black", size = 0.4) +
  geom_abline(slope = 0, intercept = 0, linetype = "solid", color = "blue",size=2,alpha=0.4) +
  scale_x_continuous(limits=c(0,100),breaks = c(0, 25, 50, 75, 100), labels = percent_format(scale = 1)) +
  scale_y_continuous(limits=c(l,u), labels = percent_format(scale = 1)) +
  labs(x = "central interval width", y = "Coverage diff.") +
  theme_minimal() +
 # plot_layout(widths = 4, heights = 4)+
  theme(
    text = element_text(size = 11),
    panel.grid.minor = element_blank() ) }

  return(p)
 }



comb_coverages_plot <- function(Results,n=NULL,beta=NULL,nu=NULL,type="errorbar",title=TRUE){
  # Results (list): experiments results from "run_experiment()"
  # n (int): observation size
  # type ( string): plot style of credible interval of the coverage
if(!is.null(n)){
    p1 <- coverage_plot(Results,n=n,target="loo_coverage",type=type)
    p2 <-coverage_plot(Results,n=n,target="test_coverage",type=type)
 label1 <- wrap_elements(grid::textGrob(TeX(paste0("$n=$",n)), rot =360,just = "center",y=0.75,  gp = gpar(fontsize = 12, fontface = "bold")))
}
 
if(!is.null(beta)){
    p1 <- coverage_plot(Results,beta=beta,target="loo_coverage",type=type)
    p2 <-coverage_plot(Results,beta=beta,target="test_coverage",type=type)
 label1 <- wrap_elements(grid::textGrob(TeX(paste0("$\\beta=$",beta)), rot =360,just = "center",y=0.75,  gp = gpar(fontsize = 12, fontface = "bold")))
}

 if(!is.null(nu)){
    p1 <- coverage_plot(Results,nu=nu,target="loo_coverage",type=type)
    p2 <-coverage_plot(Results,nu=nu,target="test_coverage",type=type)

 label1 <- wrap_elements(grid::textGrob(TeX(paste0("$\\nu=$",nu)), rot =360,just = "center",y=0.75,  gp = gpar(fontsize = 12, fontface = "bold")))
 }
 
 

   # label1 <- wrap_elements(grid::textGrob(TeX(paste0("$n=$",n)), rot =360,just = "center",y=0.75,  gp = gpar(fontsize = 16, fontface = "bold")))
    combined_plot <- p1 +plot_spacer()+p2+plot_layout(ncol=3,widths = c(4,0.2,4),heights = c(4,0.1,4))

    if(title==TRUE){
      col1 <- wrap_elements(grid::textGrob("loo-coverage",just = "center",x=0.22, gp = gpar(fontsize = 12, fontface = "bold")))
      col2 <- wrap_elements(grid::textGrob("test-coverage",just = "center", x=0.38, gp = gpar(fontsize = 12, fontface = "bold")))


    final_plot<-(plot_spacer()+col1+col2)/(label1+combined_plot+plot_layout(widths = c(0.15, 0.85)))+ plot_layout(widths =c(0.1,1,1), heights = c(0.1,1), guides = "collect")

    } else { final_plot <-  label1+combined_plot+plot_layout(widths = c(0.15, 0.85)) }

    return( final_plot)
  }




draw_hist<-function(Results,n=NULL,beta=NULL,nu=NULL,target="loo_pit",bins=10,title=TRUE){
  # n (int): observation size
  # exp_results (list) : list of list (results from experiment)
  # type (string): possible vals c("errorbar","shaded","diff") 
if(!is.null(n)) {
  summary_df<- Results$summary[[paste0("num_obs_",n)]][[target]]}

if(!is.null(beta)) {
  summary_df<- Results$summary[[paste0("beta_",beta)]][[target]] }

if(!is.null(nu)) {
  summary_df<- Results$summary[[paste0("df_",nu)]][[target]] }

  

if(target=="pit") {

p <-ggplot(summary_df, aes(x = mean)) +
  geom_histogram( bins = bins,color = "black",fill= "#0072B2",alpha=0.5) +    # red :"#F94144"    green : "#4CAF50"
  labs(title=NULL, x = "pit values",y = NULL ) +
   scale_x_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  theme_minimal()+
  theme( text = element_text(size = 11))
if(title==TRUE){ p <- p+labs(title = "pit")+theme(plot.title = element_text(hjust = 0.5,face = "bold",size=12)) }}

if(target=="test_pit") {

p <-ggplot(summary_df, aes(x = mean)) +
  geom_histogram( bins = bins,color = "black",fill= "#4CAF50",alpha=0.5) +    # red :"#F94144"   
  labs(title=NULL, x = "pit values",y = NULL ) +
   scale_x_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  theme_minimal()+
  theme( text = element_text(size = 11))

if(title==TRUE){ p <- p+labs(title = "test-pit")+theme(plot.title = element_text(hjust = 0.5,face="bold",size=12)) }}


if(target=="loo_pit") {

p <-ggplot(summary_df, aes(x = mean)) +
  geom_histogram( bins = bins,color = "black",fill= "#F94144",alpha=0.5) +    # red :"#F94144"   
  labs(title=NULL, x = "pit values",y = NULL ) +
  scale_x_continuous(limits=c(0,1),breaks = c(0,0.5,1))+
  theme_minimal()+
  theme( text = element_text(size = 11))

if(title==TRUE){ p <- p+labs(title = "loo-pit")+theme(plot.title = element_text(hjust = 0.5,face="bold",size=12)) }}


  return(p)
 }



comb_hist <- function(Results,n=NULL,beta=NULL,nu=NULL,bins=10,title=TRUE){
  # Results (list): experiments results from "run_experiment()"
  # n (int): observation size
  # type ( string): plot style of credible interval of the coverage
  if(!is.null(n)){
    p1 <- draw_hist(Results,n=n,beta=NULL,nu=NULL,target="pit",bins=bins,title=title)
    p2 <- draw_hist(Results,n=n,beta=NULL,nu=NULL,target="test_pit",bins=bins,title=title)
    p3 <- draw_hist(Results,n=n,beta=NULL,nu=NULL,target="loo_pit",bins=bins,title=title)
 label1 <- wrap_elements(grid::textGrob(TeX(paste0("$n=$",n)), rot =360,just = "center",y=0.85,  gp = gpar(fontsize = 16, fontface = "bold")))
}

  
if(!is.null(beta)){
    p1 <- draw_hist(Results,n=NULL,beta=beta,nu=NULL,target="pit",bins=bins,title=title)
    p2 <- draw_hist(Results,n=NULL,beta=beta,nu=NULL,target="test_pit",bins=bins,title=title)
    p3 <- draw_hist(Results,n=NULL,beta=beta,nu=NULL,target="loo_pit",bins=bins,title=title)
 label1 <- wrap_elements(grid::textGrob(TeX(paste0("$\\beta=$",beta)), rot =360,just = "center",y=0.85,  gp = gpar(fontsize = 16, fontface = "bold")))
}

 if(!is.null(nu)){
   p1 <- draw_hist(Results,n=NULL,beta=NULL,nu=nu,target="pit",bins=bins,title=title)
    p2 <- draw_hist(Results,n=NULL,beta=NULL,nu=nu,target="test_pit",bins=bins,title=title)
    p3 <- draw_hist(Results,n=NULL,beta=NULL,nu=nu,target="loo_pit",bins=bins,title=title)
 label1 <- wrap_elements(grid::textGrob(TeX(paste0("$\\nu=$",nu)), rot =360,just = "center",y=0.85,  gp = gpar(fontsize = 16, fontface = "bold")))
 }
 
 
   # label1 <- wrap_elements(grid::textGrob(TeX(paste0("$n=$",n)), rot =360,just = "center",y=0.75,  gp = gpar(fontsize = 16, fontface = "bold")))
    combined_plot <- p1 +plot_spacer()+p2+plot_spacer()+p3+plot_layout(ncol=5,widths = c(4,0.2,4,0.2,4),heights = c(4,0.1,4,0.1,4))


    final_plot<-label1+combined_plot+plot_layout(widths = c(0.15, 0.85))


    return( final_plot)
  }
