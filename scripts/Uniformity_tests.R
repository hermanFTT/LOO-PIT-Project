
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
