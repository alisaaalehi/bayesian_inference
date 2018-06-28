

#By: 
#	Ali Salehi 
#	Md Shiplu Hawlader


rm(list = ls()) # remove all variables from workspace
setwd("C:/Users/ali/Desktop/Google Drive/Courses/Bayesian Inference/project")

# Input data
Dose <- c(1.69, 1.72, 1.76, 1.78, 1.81, 1.84, 1.86, 1.88)
N <- c(60, 62, 63, 60, 64, 60, 62, 64)
A <- c(6, 13, 20, 30, 53, 55, 61, 62)
frequency <- A/N


# Convert to Bernoulli
d <- c() # empty vector
y <- c()
for (i in 1:length(Dose))
{
	# create dose vector
	tempD <- c(rep(Dose[i], N[i]))
	d <- c(d , tempD)

	# creatse y vector
	temp0 <- c(rep(0, N[i]-A[i]))
	temp1 <- c(rep(1, A[i]))
	y <- c(y, temp0, temp1)
}

print(d)
print(y)

N <- length(y) # Number of samples
N1  <- sum(y)  # Number of successes
N0  <- N - N1  # Number of failures

# Form the X matrix
x <- cbind(d, d^2) # In this project we want to use d^2 also
D <- 3 # Number of predictors
X <- matrix(c(rep(1, length(d)), x), ncol = D) # X is complete data matrix

# MLE: fit the model to the data
fit <- glm(y ~ x , family = binomial(link = probit))
summary(fit)

mle_beta <- fit$coefficients


###################################################### Gibbs sampling 

require(mvtnorm) # for sampling from Multivariate Normal distribution
require(truncnorm)# for sampling from Truncated Normal distribution

# Initialize parameters
beta <- rep(0, D)
z <- rep(0, N) # Latent variables

theta = 0.000001 # I use theta = 1/ (sigma^2)
prior_galpha <- 0.5 #0.1 # alpha of the prior distribution for (1/sigma^2 or thera)
prior_gbeta <- 0.005 #0.5 # beta of the prior distribution for (1/sigma^2 or thera)

# Gibbs sampler Parameters
N_sim <- 10000 # Number of simulations
burn_in <- 5000 # Burn in period
beta_chain <- matrix(0, nrow = N_sim, ncol = D) # Matrix storing samples of the \beta parameter

# Compute posterior variance of beta
X_trans_X_inv <- solve(crossprod(X, X)) #inverse of (X-transpose times X)
V <- X_trans_X_inv
X_trans_X_inv_X_trans <- tcrossprod (X_trans_X_inv, X) # Transpose of second parameter times the first param

post_galpha = N/2 - prior_galpha

method_flag = 2 # 0: non informative, 1: prior on sigma^2, 2: informative on beta
if (method_flag == 0){
	print("Non informative!")
  	for (t in 2:N_sim) 
	{
	  # Update Mean of z
	  mu_z <- X %*% beta

	  # Draw latent variable z from its full conditional: z | beta, y, X, sigma^2
	  z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1 , a = -Inf, b = 0)
	  z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1 , a = 0, b = Inf)
	  
	  # Update beta
	  M <- X_trans_X_inv_X_trans %*% z  # posterior mean of beta
	  V <- X_trans_X_inv # posterior covariance of beta
	  beta <- c(rmvnorm(1, M, V))# Draw variable beta from its full conditional: \beta | z, X, sigma^2

	  # Store the \beta draws
	  beta_chain[t, ] <- beta
	}
 curve_file <- "cr_nonInformative.pdf"
 hist_file <- "hist_nonInformative.pdf"

} else if (method_flag == 1){
	print("Prior on 1/Sima^2!")
	for (t in 2:N_sim) {
	  # Update Mean of z
	  mu_z <- X %*% beta

	  # Draw latent variable z from its full conditional: z | beta, y, X, sigma^2
	  z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1 , a = -Inf, b = 0)
	  z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1 , a = 0, b = Inf)
	  
	  # Update beta
	  M <- X_trans_X_inv_X_trans %*% z  # posterior mean of beta
	  V <- X_trans_X_inv * (1/theta) # posterior covariance of beta
	  beta <- c(rmvnorm(1, M, V))# Draw variable beta from its full conditional: \beta | z, X, sigma^2
		
	  # update theta (1 over sigma squared)
	  yxbeta <- y - (X %*% beta)
	  post_gbeta <- prior_gbeta + ((crossprod(yxbeta,yxbeta))/2)
	  theta <- rgamma(1, shape=post_galpha , rate = post_gbeta)
	  #print (1/theta)

	  # Store the \beta draws
	  beta_chain[t, ] <- beta
	}
 curve_file <- "cr_sigmaprior.pdf"
 hist_file <- "hist_sigmaprior.pdf"

} else {  
	#informative prior on beta  
	print("Informative prior on beta!")
	Q_0 <- diag(10, D)
	prec_0 <- solve(Q_0)
	beta_0 <- rep(0, D)
	V <- solve(prec_0 + crossprod(X, X))
	for (t in 2:N_sim) {
		  # Update Mean of z
		  mu_z <- X %*% beta
		  # Draw latent variable z from its full conditional: z | \beta, y, X
		  z[y == 0] <- rtruncnorm(N0, mean = mu_z[y == 0], sd = 1, a = -Inf, b = 0)
		  z[y == 1] <- rtruncnorm(N1, mean = mu_z[y == 1], sd = 1, a = 0, b = Inf)
  
		  # Compute posterior mean of theta
		  M <- V %*% (prec_0 %*% beta_0 + crossprod(X, z))
		  # Draw variable \beta from its full conditional: \theta | z, X
		  beta <- c(rmvnorm(1, M, V))
  
		  # Store the \beta draws
		  beta_chain[t, ] <- beta
	}
 curve_file <- "cr_informative.pdf"
 hist_file <- "hist_informative.pdf"
}

# Get posterior mean of beta
post_beta <- colMeans(beta_chain[-(1:burn_in), ])


# Plot covariates x versus observations y

# Show the fitted function using the posterior mean estimates
xt <- seq(from = min(Dose), to = max(Dose ), by = 0.0001)
Xt <- matrix(c(rep(1, length(xt )), xt, xt^2), ncol = D)

dev.new()
pdf(file = curve_file)
plot(Dose, frequency, ylab="p")
lines(x = xt , y = pnorm(Xt %*% mle_beta), col = "red3", lwd = 2)
lines(x = xt, y = pnorm(Xt %*% post_beta), col = "blue3", lwd = 2)
legend("bottomright", legend=c("MLE","Posterior Mean"), col=c("red3","blue3"), 
       bty = 'n', lwd = 2, inset = c(0.02, 0.08), lty = 1, cex = 0.9)
dev.off()


# Plot histograms and Trace plots
dev.new()
pdf(file = hist_file)
par(mfrow = c(2,3))

hist(beta_chain[-(1:burn_in),1],breaks=30, main="Posterior of Beta 1", 
     xlab=paste('Mean B1',toString(round(post_beta[1],1)), sep=" = "),ylab="", col="cornflowerblue")
abline(v = post_beta[1], col="goldenrod2", lwd=3)

hist(beta_chain[-(1:burn_in),2], breaks=30, main="Posterior of Beta 2", 
     xlab=paste('Mean B2',toString(round(post_beta[2],1)), sep=" = "),ylab="", col="cornflowerblue")
abline(v = post_beta[2], col="goldenrod2", lwd=3)

hist(beta_chain[-(1:burn_in),3], breaks=30, main="Posterior of Beta 3", 
     xlab=paste('Mean B3',toString(round(post_beta[3],1)), sep=" = "),ylab="", col="cornflowerblue")
abline(v = post_beta[3], col="goldenrod2", lwd=3)

legend("topright", c("Mean"), lty=1, lwd=2,
       col=c("goldenrod2"), bty='n', cex=.95)


plot(beta_chain[, 1], type = "l", xlab="Iteration" , ylab="", 
     main = "Chain values of beta1")
lines(cumsum(beta_chain[, 1])/(1:N_sim), col="red", lwd=2)

plot(beta_chain[, 2], type = "l", xlab="Iteration" , ylab="",  
     main = "Chain values of beta2")
lines(cumsum(beta_chain[, 2])/(1:N_sim), col="red", lwd=2)

plot(beta_chain[, 3], type = "l", xlab="Iteration" , ylab="", 
     main = "Chain values of beta3")
lines(cumsum(beta_chain[, 3])/(1:N_sim), col="red", lwd=2)
legend("bottomright", c("Mean chain"), lty=1, lwd=2,
       col=c("red"), bty='n', cex=.95)

dev.off()

print (matrix(round(mu_z,2), nrow = 10, byrow = FALSE))
print (post_beta)
print (mle_beta)


