Twotrim <- function (X, Z, y, Kn = NULL, c1 = 5, HDIC_Type = "HDBIC", c2 = 2, 
                     c3 = 2.01, c_n, intercept = TRUE) 
{
###############################1.First step Ohit################################
  
  #############################1.0.Setup and Check##############################
  if (!is.vector(y)) 
    stop("y should be a vector")
  if (!is.matrix(X)) 
    stop("X should be a matrix")
  n = nrow(X)
  pL = ncol(X)
  if (n != length(y)) 
    stop("the number of observations in y is not equal 
         to the number of rows of X")
  if (n == 1) 
    stop("the sample size should be greater than 1")
  
  ##################################1.1.OGA_L###################################
  if (is.null(Kn)) 
    K = max(1, min(floor(c1 * sqrt(n/log(pL))), pL))
  else {
    if ((Kn < 1) | (Kn > pL)) 
      stop(paste("Kn should between 1 and ", pL, sep = ""))
    if ((Kn - floor(Kn)) != 0) 
      stop("Kn should be a positive integer")
    K = Kn
  }
  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))
  Jhat = sigma2hat = rep(0, K)
  XJhat = matrix(0, n, K)
  RL = as.matrix(dy)
  xnorms = sqrt(colSums((dX)^2))
  aSSE = (abs(t(RL) %*% dX)/xnorms)
  Jhat[1] = which.max(aSSE)
  XJhat[, 1] = (dX[, Jhat[1]]/sqrt(sum((dX[, Jhat[1]])^2)))
  RL = RL - XJhat[, 1] %*% t(XJhat[, 1]) %*% RL
  sigma2hat[1] = mean(RL^2)
  if (K > 1) {
    for (k in 2:K) {
      aSSE = (abs(t(RL) %*% dX)/xnorms)
      aSSE[Jhat[1:(k - 1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dX[, Jhat[k]] - XJhat[, 1:(k - 1)] %*%
        t(XJhat[, 1:(k - 1)]) %*% dX[, Jhat[k]]
      XJhat[, k] = (rq/sqrt(sum((rq)^2)))
      RL = RL - XJhat[, k] %*% t(XJhat[, k]) %*% RL
      sigma2hat[k] = mean(RL^2)
    }
  }
  
  ##################################1.2.HDIC_L##################################
  #need sigma2hat and Jhat from step 1
  if ((HDIC_Type != "HDAIC") & (HDIC_Type != "HDBIC") & (HDIC_Type != 
                                                         "HDHQ")) 
    stop("HDIC_Type should be \"HDAIC\", \"HDBIC\" or \"HDHQ\"")
  if (HDIC_Type == "HDAIC") 
    omega_n = c2
  if (HDIC_Type == "HDBIC") 
    omega_n = log(n)
  if (HDIC_Type == "HDHQ") 
    omega_n = c3 * log(log(n))
  hdic = (n * log(sigma2hat)) + ((1:K) * omega_n * (log(pL)))
  kn_hat = which.min(hdic)
  benchmark = hdic[kn_hat]
  # Outpue of HDIC
  J_HDIC = sort(Jhat[1:kn_hat])
  
  ##################################1.3.Trim_L##################################
  # Initialize output of Trim
  J_Trim = Jhat[1:kn_hat]
  trim_pos = rep(0, kn_hat)
  if (kn_hat > 1) {
    for (l in 1:(kn_hat - 1)) {
      JDrop1 = J_Trim[-l]
      fit = lm(dy ~ . - 1, data = data.frame(dX[, JDrop1]))
      RLDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(RLDrop1^2))) + ((kn_hat - 
                                                  1) * omega_n * (log(pL)))
      if (HDICDrop1 > benchmark) 
        trim_pos[l] = 1
    }
    trim_pos[kn_hat] = 1
    J_Trim = J_Trim[which(trim_pos == 1)]
  }
  # Outpue of Trim
  J_Trim = sort(J_Trim)
  
  # Prepare the results
  X_HDIC = as.data.frame(as.matrix(X[, J_HDIC]))
  X_Trim = as.data.frame(as.matrix(X[, J_Trim]))
  X = data.frame(X)
  colnames(X_HDIC) = names(X)[J_HDIC]
  colnames(X_Trim) = names(X)[J_Trim]
  if (intercept == TRUE) {
    fit_HDIC = lm(y ~ ., data = X_HDIC)
    fit_Trim = lm(y ~ ., data = X_Trim)
  }
  else {
    fit_HDIC = lm(y ~ . - 1, data = X_HDIC)
    fit_Trim = lm(y ~ . - 1, data = X_Trim)
  }
  betahat_HDIC = summary(fit_HDIC)
  betahat_Trim = summary(fit_Trim)

  #############################1.4.Obtain residuals#############################
  # projection matrix H_{L,\hat{N}_{L,n}}
  if (intercept == TRUE) {
    A <- cbind(matrix(1,n,1), dX[, result[["J_Trim"]]])
  }
  else {
    A <- dX[, result[["J_Trim"]]]
  }
  proj <- A %*% solve(t(A) %*% A) %*% t(A)
  eta_hat2 <- ((diag(1,n) - proj) %*% dy)^2
  eta_tilde2 <- apply(rbind(dy, rep(c_n,length(dy))), 2, function(x) max(x))
  r_t <- log(eta_tilde2)
  dr <- r_t - mean(r_t)
  
##############################2.Second step Ohit################################
  
  #############################2.0.Setup and Check##############################
  if (!is.matrix(Z)) 
    stop("Z should be a matrix")
  n = nrow(Z)
  pD = ncol(Z)
  if (n != length(y)) 
    stop("the number of observations in y is not equal to the number of rows of Z")
  if (n == 1) 
    stop("the sample size should be greater than 1")
  
  ##################################2.1.OGA_D###################################
  if (is.null(Kn)) 
    K = max(1, min(floor(c1 * sqrt(n/log(pD))), pD))
  else {
    if ((Kn < 1) | (Kn > pD)) 
      stop(paste("Kn should between 1 and ", pD, sep = ""))
    if ((Kn - floor(Kn)) != 0) 
      stop("Kn should be a positive integer")
    K = Kn
  }
  dy = y - mean(y)
  dZ = apply(Z, 2, function(x) x - mean(x))
  Jhat = sigma2hat = rep(0, K)
  ZJhat = matrix(0, n, K)
  RD = as.matrix(dy)
  znorms = sqrt(colSums((dZ)^2))
  aSSE = (abs(t(RD) %*% dZ)/znorms)
  Jhat[1] = which.max(aSSE)
  ZJhat[, 1] = (dZ[, Jhat[1]]/sqrt(sum((dZ[, Jhat[1]])^2)))
  RD = RD - ZJhat[, 1] %*% t(ZJhat[, 1]) %*% RD
  sigma2hat[1] = mean(RD^2)
  if (K > 1) {
    for (k in 2:K) {
      aSSE = (abs(t(RD) %*% dZ)/znorms)
      aSSE[Jhat[1:(k - 1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dZ[, Jhat[k]] - ZJhat[, 1:(k - 1)] %*%
        t(ZJhat[, 1:(k - 1)]) %*% dZ[, Jhat[k]]
      ZJhat[, k] = (rq/sqrt(sum((rq)^2)))
      RD = RD - ZJhat[, k] %*% t(ZJhat[, k]) %*% RD
      sigma2hat[k] = mean(RD^2)
    }
  }
  
  ##################################2.2.HDIC_D##################################
  #need sigma2hat and Jhat from step 1
  if ((HDIC_Type != "HDAIC") & (HDIC_Type != "HDBIC") & (HDIC_Type != 
                                                         "HDHQ")) 
    stop("HDIC_Type should be \"HDAIC\", \"HDBIC\" or \"HDHQ\"")
  if (HDIC_Type == "HDAIC") 
    omega_n = c2
  if (HDIC_Type == "HDBIC") 
    omega_n = log(n)
  if (HDIC_Type == "HDHQ") 
    omega_n = c3 * log(log(n))
  hdic = (n * log(sigma2hat)) + ((1:K) * omega_n * (log(pD)))
  kn_hat = which.min(hdic)
  benchmark = hdic[kn_hat]
  # Outpue of HDIC
  J_HDIC = sort(Jhat[1:kn_hat])
  
  ##################################2.3.Trim_D##################################
  # Initialize output of Trim
  J_Trim = Jhat[1:kn_hat]
  trim_pos = rep(0, kn_hat)
  if (kn_hat > 1) {
    for (l in 1:(kn_hat - 1)) {
      JDrop1 = J_Trim[-l]
      fit = lm(dy ~ . - 1, data = data.frame(dZ[, JDrop1]))
      RDDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(RDDrop1^2))) + ((kn_hat - 
                                                  1) * omega_n * (log(pD)))
      if (HDICDrop1 > benchmark) 
        trim_pos[l] = 1
    }
    trim_pos[kn_hat] = 1
    J_Trim = J_Trim[which(trim_pos == 1)]
  }
  # Outpue of Trim
  J_Trim = sort(J_Trim)
  
  # Prepare the results
  Z_HDIC = as.data.frame(as.matrix(Z[, J_HDIC]))
  Z_Trim = as.data.frame(as.matrix(Z[, J_Trim]))
  Z = data.frame(Z)
  colnames(Z_HDIC) = names(Z)[J_HDIC]
  colnames(Z_Trim) = names(Z)[J_Trim]
  if (intercept == TRUE) {
    fit_HDIC = lm(y ~ ., data = Z_HDIC)
    fit_Trim = lm(y ~ ., data = Z_Trim)
  }
  else {
    fit_HDIC = lm(y ~ . - 1, data = Z_HDIC)
    fit_Trim = lm(y ~ . - 1, data = Z_Trim)
  }
  betahat_HDIC = summary(fit_HDIC)
  betahat_Trim = summary(fit_Trim)
  
  
  
  return(list(n = n, p = p, Kn = K, J_OGA = Jhat, HDIC = hdic, 
              J_HDIC = J_HDIC, J_Trim = J_Trim, betahat_HDIC = betahat_HDIC, 
              betahat_Trim = betahat_Trim))
}
