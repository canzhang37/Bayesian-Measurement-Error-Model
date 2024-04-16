library(coda)
library(runjags)
library('Rlab')
library(truncnorm)

beta0_true <- 1
beta1_true <- -1.5
beta2_true <- -0.5
sigma_true <- 0.5
tau_true <- 5
alpha_x_true <- 2.0
beta_x_true <- 5.0

size <- c(4,8,16)
Times <- 1000

bem <- "
model{
  ## Likelihood
  for (i in 1:N){
    logY[i] ~ dnorm(beta0+x[brand_y[i]]*beta1+z[i]*beta2, 1/sigma^2)
  }
  for (k in 1:K){
    w[k] ~ dbeta((x[brand_w[k]]*tau), ((1-x[brand_w[k]])*tau))
  }
  for(j in 1:J){
    x[j] ~ dbeta(alpha_x, beta_x)
  }
  
  ## Priors
  beta0 ~ dnorm(0, 1/3^2)
  beta1 ~ dnorm(0, 1/3^2)
  beta2 ~ dnorm(0, 1/3^2)
  sigma ~ dgamma(0.01, 0.01)
  tau ~ dgamma(0.01, 0.01)
  alpha_x ~ dgamma(0.01, 0.01)
  beta_x ~ dgamma(0.01, 0.01)
}
"


require(rjags); require(runjags)

beta0_naive_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(beta0_naive_df) <- c("Lower95","Upper95", "Mean","SD")
beta1_naive_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(beta1_naive_df) <- c("Lower95","Upper95", "Mean","SD")
beta2_naive_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(beta2_naive_df) <- c("Lower95","Upper95", "Mean","SD")
sigma_naive_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(sigma_naive_df) <- c("Lower95","Upper95", "Mean","SD")
nu_naive_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(nu_naive_df) <- c("Lower95","Upper95", "Mean","SD")

beta0_bayes_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(beta0_bayes_df) <- c("Lower95","Upper95", "Mean","SD")
beta1_bayes_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(beta1_bayes_df) <- c("Lower95","Upper95", "Mean","SD")
beta2_bayes_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(beta2_bayes_df) <- c("Lower95","Upper95", "Mean","SD")
sigma_bayes_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(sigma_bayes_df) <- c("Lower95","Upper95", "Mean","SD")
nu_bayes_df <- data.frame(matrix(ncol = 4, nrow = Times))
colnames(nu_bayes_df) <- c("Lower95","Upper95", "Mean","SD")


                                                     
for (s in size){
  count <- 1
  requirement_met <- FALSE
  for (t in 1:100000){
    n_b <- 50
    set.seed((n_b+1)*t+1)
    print(n_b*s)
    print(count)
    N_brand <- n_b
    I_each <- s
    x <- rbeta(N_brand, alpha_x_true, beta_x_true)
    J <- N_brand
    
    I_list <- rep(I_each, N_brand)
    N <- sum(I_list)
    z <- rbern(N, prob = 0.4)
    
    logY_mean <- c()
    logY <- c()
    w <- c()
    brand_y<- c()
    
    m<-3
    
    for (j in 1:J){
      logY_mean <- append(logY_mean, beta0_true + x[j]*beta1_true)
      w <- append(w,rbeta(m, (x[j]*tau_true), ((1-x[j])*tau_true)))
    }
    w
    logY_mean <- rep(logY_mean,rep(s,N_brand)) + beta2_true*z

    for (i in 1:N){
      logY <- append(logY,rnorm(1,logY_mean[i], sigma_true))
    }
    brand_y <- rep(1:J,each=I_each)
    brand_w <- rep(1:J,each=m)
    K <- m * J
    
    datalist<-list(logY=logY,w=w,z=z,brand_y=brand_y,brand_w=brand_w,N=N,K=K,J=J)
    
    x_naive <- c()
    for (j in 1:J){
      x_naive <- append(x_naive,rep(mean(w[(m*(j-1)+1):(m*j)]),I_each))
    }
    naive <- lm(logY~x_naive+z)
    summary_naive <- summary(naive)
    mean_naive <- summary_naive$coefficients[,1]
    sd_naive <- summary_naive$coefficients[,2]
    L95_naive <- confint(naive)[,1]
    U95_naive <- confint(naive)[,2]

  
    beta0_naive_df[count,] <- c(L95_naive[1],U95_naive[1],mean_naive[1],sd_naive[1])
    beta1_naive_df[count,] <- c(L95_naive[2],U95_naive[2],mean_naive[2],sd_naive[2])
    beta2_naive_df[count,] <- c(L95_naive[3],U95_naive[3],mean_naive[3],sd_naive[3])

    
    write.csv(beta0_naive_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/beta0_naive_", N, ".csv", sep = ""),row.names=F)
    write.csv(beta1_naive_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/beta1_naive_", N, ".csv", sep = ""),row.names=F)
    write.csv(beta2_naive_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/beta2_naive_", N, ".csv", sep = ""),row.names=F)
    write.csv(sigma_naive_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/sigma_naive_", N, ".csv", sep = ""),row.names=F)
    write.csv(nu_naive_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/nu_naive_", N, ".csv", sep = ""),row.names=F)
    
    inits1 <- list(beta0 = 0.1, 
                   beta1 = -1, 
                   beta2 = -0.1, 
                   sigma = 0.5, 
                   tau = 100, 
                   alpha_x = 1, 
                   beta_x = 4, 
                   .RNG.name = 'base::Wichmann-Hill', 
                   .RNG.seed = n_b*Times + t)
    inits2 <- list(beta0 = 0.5, 
                   beta1 = -0.5, 
                   beta2 = -0.5, 
                   sigma = 1, 
                   tau = 200, 
                   alpha_x = 2, 
                   beta_x = 2, 
                   .RNG.name = 'base::Wichmann-Hill', 
                   .RNG.seed = (n_b+1)*Times + t)
    inits3 <- list(beta0 = 1, 
                   beta1 = -0.1, 
                   beta2 = -1, 
                   sigma = 2, 
                   tau = 400, 
                   alpha_x = 4, 
                   beta_x = 1, 
                   .RNG.name = 'base::Wichmann-Hill', 
                   .RNG.seed = (n_b+2)*Times + t)
    initlist <- list(inits1, inits2, inits3)
    
    model_attempt <- try({
      results <- run.jags(model = bem,
                          params_to_monitor <- c("beta0","beta1","beta2","sigma","tau"),
                          inits = initlist,
                          data = datalist,
                          n.chains = 3,
                          burnin= 2000,
                          sample = 4000,
                          adapt = 2000)
    }, silent = TRUE)
    
    if (inherits(model_attempt, "try-error")) {
      print(paste("Caught an error at iteration", t, ". Skipping to next iteration..."))
      next
    }
    
    summary_bayes <- summary(results)
    gelman_result <- gelman.diag(results)
    max_gelman_result <- max(gelman_result$psrf)
    if (max_gelman_result>1.1){
      next
    }
    beta0_bayes <- summary_bayes[1,c(1,3:5)]
    beta1_bayes <- summary_bayes[2,c(1,3:5)]
    beta2_bayes <- summary_bayes[3,c(1,3:5)]
    beta0_bayes_df[count,] <- beta0_bayes
    beta1_bayes_df[count,] <- beta1_bayes
    beta2_bayes_df[count,] <- beta2_bayes
    write.csv(beta0_bayes_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/beta0_bayes_", N, ".csv", sep = ""),row.names=F)
    write.csv(beta1_bayes_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/beta1_bayes_", N, ".csv", sep = ""),row.names=F)
    write.csv(beta2_bayes_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/beta2_bayes_", N, ".csv", sep = ""),row.names=F)
    write.csv(sigma_bayes_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/sigma_bayes_", N, ".csv", sep = ""),row.names=F)
    write.csv(nu_bayes_df,file=paste("C:/Users/Mayes/Desktop/Plan B/Simulation_beta/beta_1/nu_bayes_", N, ".csv", sep = ""),row.names=F)

    if(count>=Times){
      requirement_met <- TRUE
      break
    }
    count <- count + 1
  }
  if (requirement_met) {
    next
  }
}
