### define Z(t) as function of t

q = 4
n = 100 ## number of sample components
T = seq(1, 4, length = 10) 
m = length(T)  ##number of time points
Z = function(t) {
  return(matrix(c(12-0.5*t^2, -2*t, log(3+1/t), exp(-t)), ncol = 1))
}

sigma_e = 0.05

### define paramters for Theta_i
library(MASS)

theta = c(4, 3, 5, 5)
Sigma = matrix(c(0.75, -.27, .092, -.21, -.27, .86, 
                 .15, .049, .092, .15, .75, -.071, 
                 -.21, .049, -0.071, .24), nrow = q)



##### Simulate

Simu = function(n, m, q, theta, Sigma, sigma_e) {
  Datmat = matrix(nrow = m, ncol = n)
  for(i in 1:n) {
    Theta = mvrnorm(1, theta, Sigma)
    for(j in 1:m){
      Datmat[j, i] = t(Z(T[j])) %*% Theta + rnorm(1, mean = 0, sd = sigma_e)
    }
  }
  mylist = list("data" = Datmat)
  return(mylist)
}

### declare cutoff
eta = 20

#### Estimation of R-hat(t)

Estim = function(n, m, q, S, Z, T) {
  Zmat = matrix(nrow = m, ncol = q)
  for(i in 1:m){
    Zmat[i, ] = Z(T[i])
  }
  
  Ybar = 0
  for(i in 1:n){
    Ybar = Ybar + S[, i]
  }
  Ybar = Ybar/n
  
  #### hat(theta)_n
  
  theta_hat = solve(t(Zmat) %*% Zmat) %*% t(Zmat) %*% Ybar
  
  #### hat(sigma)_e^2
  
  sigma_e2_hat = 0
  for(i in 1:n) {
    sigma_e2_hat = sigma_e2_hat + t(S[,i]) %*% S[,i] - t(S[,i]) %*% Zmat %*% solve(t(Zmat)%*%Zmat) %*% t(Zmat)%*%S[,i]
  }
  
  sigma_e2_hat = sigma_e2_hat/(n*(m-q))
  
  #### hat(Sigma_Theta)
  
  Sigma_Theta = matrix(rep(0, q^2), nrow = q)
  for(i in 1:n) {
    X = solve(t(Zmat) %*% Zmat) %*% t(Zmat) %*% (S[,i]-Ybar) %*% t(S[,i]-Ybar) %*% Zmat %*% solve(t(Zmat) %*% Zmat)
    Sigma_Theta = Sigma_Theta + X
  }
  
  Sigma_Theta = Sigma_Theta/n - as.numeric(sigma_e2_hat) * solve(t(Zmat) %*% Zmat)
  
  #### hat(s(t))
  
  st = c()
  for(i in 1:m) {
    x = sqrt(t(Zmat[i, ]) %*% Sigma_Theta %*% Zmat[i, ])
    st = c(st, x)
  }
  
  ##### hat(R(t))
  
  Rt = c()
  for(i in 1:m) {
    x = (t(Zmat[i, ]) %*% theta_hat - eta) / st[i]
    Rt = c(Rt, pnorm(x))
  }
  
  mylist = list("theta_hat" = theta_hat, "sigma_e2_hat" = sigma_e2_hat, "Sigma_Theta" = Sigma_Theta, "Rt" = Rt)
  return(mylist)
  
}


########RUN
set.seed(1)
S = Simu(n,m,q,theta, Sigma, sigma_e)$data
#result = Estim(n,m,q,S,Z,T)
#result
#is.positive.definite(result$Sigma_Theta)

#### Now that we are done with Simulation and Estimation, lets jump into Jackknife and Bootstrap!


###### JACKKNIFE ############
Jmat = matrix(nrow = n, ncol = m)
for(i in 1:n) {
  S_new = S[, -i]
  result = Estim((n-1), m, q, S_new, Z, T)
  Jmat[i, ] = result$Rt
}

v_jack = ((n-1)^2/n) * var(Jmat) ####TAA_DAA




###### BOOTSTRAP ############
B = 100
Bmat = matrix(nrow = B, ncol = m)
for(b in 1:B) {
  S_new = S[ , sample(n, n, replace = TRUE)]
  result = Estim(n, m, q, S_new, Z, T)
  Bmat[b, ] = result$Rt
}

v_boot = var(Bmat) ####TAA-DAA
