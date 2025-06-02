
##############  Data generating mechanisms & model

# Generate synthetic  data
simul_stt<- function(N = 10,beta=0.2,seed=111,df=100) {
      # continuous data with student_t noise
      set.seed(seed)
      x1 <- rnorm(N) # Covariate
      # student_t noise
      y <- beta * x1 + (1-beta)*rt(N,df=df)     # Response variable
      sim_data <- data.frame(x1 = x1, y = y)
      return(sim_data)
    }


#

simul_count <- function(N=10, beta=0.8,seed=222,df=100) {
  # count data from poisson dist
        set.seed(seed)
         x1 <- rnorm(N)  
        lambda<-exp(beta*x1)
        y<-rpois(N,lambda)
        sim_data <- data.frame(x1 = x1, y=y)
        return(sim_data)
      }

# rnbinom(n = 400, size = 1e3, mu = 10) # control overall dispersion "size" 

# simulate data with GP models ( Hilbert space based GP approximation)
 # load the stan file
hsgp_model<- cmdstan_model("stan_file/HSGP_simulate.stan")
# writeLines(readLines(hsgp_model))

simul_hsgp <- function(N=50,beta=0.4,seed=333,model=hsgp_model,ls=0.3,sig=1,noise=0.3,M_f=160,c_f=3,inter_f=0.4,df=100){

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
