library(gt)
library(dplyr)
library(ggplot2)
library(Rlab)
library(ggplot2)
library(tidyr)
set.seed(42)
data <- read.csv('C:/Users/Mayes/Desktop/Plan B/Data Analysis_beta/FV_CENIC_datamerge_final_01Sep2023.csv')
data <- data[!is.na(data$Filter.Ventilation),]
data$menthol <- ifelse(data$menthol_usual_new==2,0,1)
nrow(data)
sum((data$tot_nnal_pmolperml_BSL2!=99999))
sum((data$tne_nmolperml_BSL2!=99999))
sum((data$FTND_score_SCR!=99999))
sum((data$baseline_cpd!=99999))

data_nnal <- data[which(data$tot_nnal_pmolperml_BSL2!=99999),]
data_tne <- data[which(data$tne_nmolperml_BSL2!=99999),]
data_ftnd <- data[which(data$FTND_score_SCR!=99999),]
data_cpd <- data[which(data$baseline_cpd!=99999),]


# NNAL
nnal <- data_nnal$tot_nnal_pmolperml_BSL2
log_nnal <- log(nnal)
age <- data_nnal$Demo_2_1_TEXT_SCR
sex <- ifelse(data_nnal$Demo_3_SCR==1,1,0)
race <- ifelse(data_nnal$race_3cat==0,1,0)
edu <- ifelse(data_nnal$Demo_9_SCR<=3,0,1)
duration <- age - data_nnal$Tob_A3_1_TEXT_SCR
menthol <- data_nnal$menthol
z <- cbind(age, sex,race,edu,duration,menthol)
L <- ncol(z)

mean(age)
sd(age)
mean(sex)
mean(race)
mean(menthol)
mean(edu)
mean(duration)
sd(duration)


w <- as.vector(apply(data_nnal[, (ncol(data_nnal)-2):ncol(data_nnal)], 1, c))
w <- w/100
K <- length(w)
N <- length(log_nnal)

brand_y <- data_nnal$ID.Picture.Book
J <- length(unique(brand_y))
numeric_y <- as.numeric(factor(brand_y))
numeric_w <- rep(numeric_y, each=3)

# Bayes
beta <- "
model{
  pi <- 3.141592653589793
  
  ## Likelihood
  for (i in 1:N){
    logY[i] ~ dnorm(beta0+x[brand_y[i]]*beta1+sum(beta2 * z[i,]), 1/sigma^2)
    log_lik[i] <- log(dnorm(logY[i], beta0 + x[brand_y[i]] * beta1 + sum(beta2 * z[i,]), sigma))
  }
  
  # Compute the deviance
  deviance = -2 * sum(log_lik)
  
    # Predictions for new data
  for (i in 1:N){
    pred_logY[i] <- beta0 + x[brand_y[i]] * beta1 + sum(beta2 * z[i,])
  }
  
  for (k in 1:K){
    w[k] ~ dbeta((x[brand_w[k]]*tau), ((1-x[brand_w[k]])*tau))
  }
  for(j in 1:J){
    x[j] ~ dbeta(alpha_x, beta_x)
  }
  for(l in 1:L){
    beta2[l] ~ dnorm(0, 1/3^2)
  }
  ## Priors
  beta0 ~ dnorm(0, 1/3^2)
  beta1 ~ dnorm(0, 1/3^2)
  sigma ~ dgamma(0.01, 0.01)
  tau ~ dgamma(0.01, 0.01)
  alpha_x ~ dgamma(0.01, 0.01)
  beta_x ~ dgamma(0.01, 0.01)
}
"

require(rjags); require(runjags)
datalist<-list(logY=log_nnal, w=w,z=z,brand_y=numeric_y,brand_w=numeric_w,N=N,K=K,J=J,L=L)


library(rjags)

# Assuming you've already loaded your data and the model string `beta`:
adapt <- 4000
burnin <- 4000
model <- jags.model(textConnection(beta), data=datalist, n.chains=3, n.adapt=adapt)
update(model, n.iter=burnin)

n_samples <- 8000
# samples <- coda.samples(model, variable.names=c("deviance"), n.iter=n_samples)
# 
# deviance_samples <- as.matrix(samples)[, "deviance"]
# DIC_beta <- mean(deviance_samples) + var(deviance_samples) / 2
# print(DIC_beta)





# Set the number of iterations for the post-burn-in sampling
n.iterations <- 10000  # for example, adjust based on your needs

# Sample from the posterior distribution
samples <- coda.samples(model, variable.names=c("pred_logY"), n.iter=n.iterations)

library(coda)


pred_beta <- summary(samples)


pred_beta$statistics[,1]



length(numeric_y)
length(pred_beta$statistics[,1])





w_mean <- FV_mean
x1 <- w_mean$mean[match(brand_y, w_mean$ID)]
length(x1)
y <- log_nnal
length(y)
results_naive <- lm(y ~ x1 + age + sex + race + edu + duration + menthol)
summary(results_naive)

data_simple <- data.frame(
  x1 = x1,
  age = age,
  sex = sex,
  race = race,
  edu = edu,
  duration = duration,
  menthol = menthol,
  y = y
)

new_predictions <- predict(results_naive, newdata=data_simple)


test_data <- data.frame(
  true = exp(log_nnal),
  beta = exp(pred_beta$statistics[,1]),
  simple = exp(new_predictions),
  brand = brand_y
)






library(ggplot2)  # Load the ggplot2 package


# Plot for the Beta Method
p1 <- ggplot(test_data, aes(x=true, y=beta, color=brand)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE) +
  ggtitle("Beta Method vs True Values") +
  xlab("True Values") +
  ylab("Predicted by Beta Method")

# Plot for the Simple Method
p2 <- ggplot(test_data, aes(x=true, y=simple, color=brand)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE) +
  ggtitle("Simple Method vs True Values") +
  xlab("True Values") +
  ylab("Predicted by Simple Method")

# Display the plots
print(p1)
print(p2)


# Reshape data for faceting
long_data <- reshape2::melt(test_data, id.vars=c("true", "brand"), measure.vars=c("beta", "simple"))

# Create a combined plot
p <- ggplot(long_data, aes(x=true, y=value, color=brand)) +
  geom_point() +
  geom_smooth(method='lm', se=FALSE) +
  facet_wrap(~variable, scales="free_y", ncol=1) +
  ggtitle("Comparison of Prediction Methods") +
  xlab("True Values") +
  ylab("Predicted Values")

# Display the plot
print(p)



library(reshape2)

# Reshape data to long format
long_data <- melt(test_data, id.vars=c("true", "brand"), measure.vars=c("beta", "simple"))
# Plot the data
ggplot(long_data, aes(x=true, y=value, group=brand, color=variable)) + 
  geom_line() +
  facet_wrap(~brand, scales="free") +
  labs(title="Performance Comparison Over Brands",
       x="True Values",
       y="Predicted Values",
       color="Method") +
  theme_minimal() +
  scale_color_manual(values=c("beta"="blue", "simple"="red"))




# Function to calculate RMSE
calc_rmse <- function(predictions, true_values) {
  sqrt(mean((predictions - true_values)^2))
}

# Calculate RMSE for each brand and method
rmse_data <- aggregate(cbind(beta, simple) ~ brand, data=test_data, FUN=function(x) calc_rmse(x, test_data$true))

# Reshape the data to long format
library(reshape2)
rmse_long <- melt(rmse_data, id.vars = "brand", variable.name = "method", value.name = "RMSE")

library(ggplot2)

# Plotting the data
ggplot(rmse_long, aes(x=brand, y=RMSE, fill=method)) + 
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("beta"="blue", "simple"="red")) +
  labs(title="RMSE Comparison Across Brands",
       x="Brand",
       y="RMSE",
       fill="Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))  # Rotate X labels for better readability




# Assuming rmse_long is already created with 'brand', 'method', and 'RMSE'

# Calculate the mean RMSE for each brand
rmse_long$brand <- factor(rmse_long$brand, levels=rmse_long$brand[order(-rmse_long$RMSE)])

# Create the plot
ggplot(rmse_long, aes(x=reorder(brand, RMSE, FUN=mean), y=exp(RMSE), fill=method)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.7), alpha=0.8) +
  scale_fill_manual(values=c("beta"="blue", "simple"="red")) +
  labs(title="RMSE Comparison Across Brands",
       x="Brand",
       y="RMSE",
       fill="Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  guides(fill=guide_legend(reverse=TRUE))  # Reverse the legend order if needed


# Calculate the mean RMSE for each brand
rmse_long$brand <- factor(rmse_long$brand, levels=rmse_long$brand[order(-rmse_long$RMSE)])

# Apply an exponential transformation to the RMSE values to exaggerate differences
rmse_long$exp_RMSE <- exp(rmse_long$RMSE - min(rmse_long$RMSE))  # Subtracting the min to keep values in a reasonable range

# Create the plot with exaggerated differences
ggplot(rmse_long, aes(x=reorder(brand, exp_RMSE, FUN=mean), y=exp_RMSE, fill=method)) + 
  geom_bar(stat="identity", position=position_dodge(width=0.7), alpha=0.8) +
  scale_fill_manual(values=c("beta"="blue", "simple"="red")) +
  labs(title="Exaggerated RMSE Comparison Across Brands",
       x="Brand",
       y="Exaggerated RMSE (exp)",
       fill="Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1)) +
  guides(fill=guide_legend(reverse=TRUE))



# Assuming you have the 'test_data' dataframe as previously defined

# Install necessary packages if not already installed
if (!require(reshape2)) install.packages("reshape2")
if (!require(ggplot2)) install.packages("ggplot2")

# Load the libraries
library(reshape2)
library(ggplot2)

# Function to calculate MAE
calc_mae <- function(predictions, true_values) {
  mean(abs(predictions - true_values))
}

# Calculate MAE for each brand and method
mae_data <- aggregate(cbind(beta, simple) ~ brand, data=test_data, FUN=function(x) calc_mae(x, test_data$true))

# Reshape the data to long format for ggplot2
mae_long <- melt(mae_data, id.vars = "brand", variable.name = "method", value.name = "MAE")

# Plotting the MAE data
ggplot(mae_long, aes(x=reorder(brand, MAE, FUN=mean), y=MAE, fill=method)) +
  geom_bar(stat="identity", position=position_dodge(width=0.7)) +
  scale_fill_manual(values=c("beta"="blue", "simple"="red")) +
  labs(title="MAE Comparison Across Brands",
       x="Brand",
       y="Mean Absolute Error (MAE)",
       fill="Method") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle=90, hjust=1)) +
  guides(fill=guide_legend(title="Prediction Method"))

# Display the plot
print(ggplot(mae_long, aes(x=reorder(brand, MAE, FUN=mean), y=MAE, fill=method)) +
        geom_bar(stat="identity", position=position_dodge(width=0.7)) +
        scale_fill_manual(values=c("beta"="blue", "simple"="red")) +
        labs(title="MAE Comparison Across Brands",
             x="Brand",
             y="Mean Absolute Error (MAE)",
             fill="Method") +
        theme_minimal() +
        theme(axis.text.x = element_text(angle=90, hjust=1)) +
        guides(fill=guide_legend(title="Prediction Method")))






w_mean <- data_nnal$Filter.Ventilation/100
results_naive <- lm(log_nnal ~ w_mean + age + sex + race + edu + duration + menthol)
summary(results_naive)
rss <- sum(residuals(results_naive)^2)
aic_val <- AIC(results_naive)
bic_val <- BIC(results_naive)
# 
univ <- "
model{
  ## Likelihood
  for (i in 1:N){
    logY[i] ~ dnorm(beta0+x[brand_y[i]]*beta1+sum(beta2 * z[i,]), 1/sigma^2)
    log_lik[i] <- log(dnorm(logY[i], beta0 + x[brand_y[i]] * beta1 + sum(beta2 * z[i,]), sigma))
  }
  
  # Compute the deviance
  deviance <- -2 * sum(log_lik)
  
  for (k in 1:K){
    w_star[k] ~ dnorm(x_star[brand_w[k]], 1/nu^2)
  }
  for(j in 1:J){
    x_star[j] ~ dnorm(mu_x, 1/sigma_x^2)
    x[j] <- ilogit(x_star[j])
  }
  for(l in 1:L){
    beta2[l] ~ dnorm(0, 1/4^2)
  }
  ## Priors
  mu_x ~ dnorm(0, 1/4^2)
  sigma_x ~ dgamma(0.01, 0.01)
  beta0 ~ dnorm(0, 1/4^2)
  beta1 ~ dnorm(0, 1/4^2)
  sigma ~ dgamma(0.01, 0.01)
  nu ~ dgamma(0.01, 0.01)
}
"

w_star <- log(w / (1 - w))
datalist<-list(logY=log_nnal,w_star=w_star,z=z,brand_y=numeric_y,brand_w=numeric_w,N=N,K=K,J=J,L=L)
adapt <- 4000
burnin <- 4000
model <- jags.model(textConnection(univ), data=datalist, n.chains=3, n.adapt=adapt)
update(model, n.iter=burnin)

n_samples <- 8000
samples <- coda.samples(model, variable.names=c("deviance"), n.iter=n_samples)

deviance_samples <- as.matrix(samples)[, "deviance"]
DIC_univ <- mean(deviance_samples) + var(deviance_samples) / 2
print(DIC_univ)

# post_mean_deviance <- mean(as.vector(samples[[1]]))
# n <- datalist$N
# k <- 6+J+L
# aic_approx <- post_mean_deviance + 2*k
# bic_approx <- post_mean_deviance + k * log(n)

# results_gen <- run.jags(model = univ,
#                         params_to_monitor <- c("beta0","beta1","beta2","sigma","nu"),
#                         data = datalist,
#                         n.chains = 3,
#                         burnin= 20000,
#                         sample = 10000,
#                         adapt = 10000)
# 
# # summary(results_naive)
# # summary(results_beta)
# summary(results_gen)