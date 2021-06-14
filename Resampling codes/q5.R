##### Defining global things
### define Z(t) as function of t

q = 4
n = 100 ## number of sample components
T = seq(1, 4, length = 10) 
m = length(T)  ##number of time points
Z = function(t) {
  return(matrix(c(12-0.5*t^2, -2*t, log(3+1/t), exp(-t)), ncol = 1))
}

sigma_e = 0.05

### declare cutoff
eta = 20

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
  
  mylist = list("theta_hat" = theta_hat, "sigma_e2_hat" = sigma_e2_hat,
                "Sigma_Theta_hat" = Sigma_Theta, "Rt" = Rt)
  return(mylist)
  
}


#### Now that we are done with Simulation and Estimation, lets jump into Jackknife and Bootstrap!

vjack = function(n, m, S) {
  Jmat = matrix(nrow = n, ncol = m)
  for(i in 1:n) {
    S_new = S[, -i]
    result = Estim((n-1), m, q, S_new, Z, T)
    Jmat[i, ] = result$Rt
  }
  return(((n-1)^2/n) * diag(var(Jmat)))
}

vboot = function(B = 200, m, S) {
  Bmat = matrix(nrow = B, ncol = m)
  for(b in 1:B) {
    S_new = S[ , sample(n, n, replace = TRUE)]
    result = Estim(n, m, q, S_new, Z, T)
    Bmat[b, ] = result$Rt
  }
  return(((B-1)/B) * diag(var(Bmat)))
}

### repeating Jackknife and Bootstrap multiple times
N = 100; v_jack = c(); v_boot = c(); B = 2*n
for(i in 1:N) {
  set.seed(i)
  S = Simu(n,m,q,theta, Sigma, sigma_e)$data
  v_jack = rbind(v_jack, vjack(n, m, S))
  v_boot = rbind(v_boot, vboot(B, m, S))
}


##### Now trying to find the tedious one -- the linear estimate

gradg = function(x1, x2, x3, i) {
  d = x2 - x1^2 - x3*(t(Zmat[i, ])%*%(solve(t(Zmat)%*%Zmat))%*%Zmat[i,])
  g1 = dnorm((x1 - eta)/sqrt(d)) * (sqrt(d) + (x1 - eta)*x1/sqrt(d))/d
  g2 = dnorm((x1 - eta)/sqrt(d)) * (-0.5*(x1 - eta)/sqrt(d))/d
  g3 = dnorm((x1 - eta)/sqrt(d)) * 
    (0.5*(x1 - eta)*(t(Zmat[i, ])%*%(solve(t(Zmat)%*%Zmat))%*%Zmat[i,])/sqrt(d))/d
  
  return(c(g1, g2, g3))
}

N = 100; v_L = c()
Zmat = matrix(nrow = m, ncol = q)
for(i in 1:m){
  Zmat[i, ] = Z(T[i])
}
for(i in 1:N) {
  set.seed(i)
  S = Simu(n,m,q,theta, Sigma, sigma_e)$data
  M=solve(t(Zmat)%*%Zmat)
  K=Zmat%*%M%*%t(Zmat)
  v_linear = c()
  for(t in 1:m) {
    P=t(Zmat[t, ])%*%M%*%t(Zmat)
    f1<-function(v)
    {
      return(P%*%v)
    }
    x1=apply(S, 2, f1)
    x2=x1^2
    f2<-function(v)
    {
      return(t(v)%*%v-t(v)%*%K%*%v)
    }
    x3=apply(S, 2, f2)/(m-q)
    X=cbind(x1, x2, x3)
    Sigma_hat_X=var(X)
    u=colMeans(X)
    
    g=gradg(u[1],u[2],u[3], t)
    v=(t(g)%*%Sigma_hat_X%*%g)/n
    v_linear = c(v_linear, v)
  }
  v_L = rbind(v_L, v_linear)
}

report_mean_matrix = matrix(nrow = m, ncol = 3)
report_sd_matrix = matrix(nrow = m, ncol = 3)

for(j in 1:m) {
  report_mean_matrix[j, 1] = mean(sqrt(v_jack[, j]))
  report_mean_matrix[j, 2] = mean(sqrt(v_boot[, j]))
  report_mean_matrix[j, 3] = mean(sqrt(v_L[, j]))
  report_sd_matrix[j, 1] = sd(sqrt(v_jack[, j]))
  report_sd_matrix[j, 2] = sd(sqrt(v_boot[, j]))
  report_sd_matrix[j, 3] = sd(sqrt(v_L[, j]))
}

report_mean_matrix
report_sd_matrix

library(beepr)
beep(4)