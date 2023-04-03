Twotrim <- function (X, Z, y, Kn = NULL, c1 = 5, HDIC_Type = "HDBIC", c2 = 2, 
                     c3 = 2.01, intercept = TRUE) 
{
  ###############################1.First step Ohit################################
  if (!is.vector(y)) 
    stop("y should be a vector")
  if (!is.matrix(X)) 
    stop("X should be a matrix")
  n = nrow(X)
  p = ncol(X)
  if (n != length(y)) 
    stop("the number of observations in y is not equal 
         to the number of rows of X")
  if (n == 1) 
    stop("the sample size should be greater than 1")
  
  ##################################1.1.OGA_L###################################
  if (is.null(Kn)) 
    K = max(1, min(floor(c1 * sqrt(n/log(p))), p))
  else {
    if ((Kn < 1) | (Kn > p)) 
      stop(paste("Kn should between 1 and ", p, sep = ""))
    if ((Kn - floor(Kn)) != 0) 
      stop("Kn should be a positive integer")
    K = Kn
  }
  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))
  Jhat = sigma2hat = rep(0, K)
  XJhat = matrix(0, n, K)
  u = as.matrix(dy)
  xnorms = sqrt(colSums((dX)^2))
  aSSE = (abs(t(u) %*% dX)/xnorms)
  Jhat[1] = which.max(aSSE)
  XJhat[, 1] = (dX[, Jhat[1]]/sqrt(sum((dX[, Jhat[1]])^2)))
  u = u - XJhat[, 1] %*% t(XJhat[, 1]) %*% u
  sigma2hat[1] = mean(u^2)
  if (K > 1) {
    for (k in 2:K) {
      aSSE = (abs(t(u) %*% dX)/xnorms)
      aSSE[Jhat[1:(k - 1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dX[, Jhat[k]] - XJhat[, 1:(k - 1)] %*%
        t(XJhat[, 1:(k - 1)]) %*% dX[, Jhat[k]]
      XJhat[, k] = (rq/sqrt(sum((rq)^2)))
      u = u - XJhat[, k] %*% t(XJhat[, k]) %*% u
      sigma2hat[k] = mean(u^2)
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
  hdic = (n * log(sigma2hat)) + ((1:K) * omega_n * (log(p)))
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
      uDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(uDrop1^2))) + ((kn_hat - 
                                                  1) * omega_n * (log(p)))
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
  
  ##############################2.Second step Ohit################################
  if (!is.vector(y)) 
    stop("y should be a vector")
  if (!is.matrix(X)) 
    stop("X should be a matrix")
  n = nrow(X)
  p = ncol(X)
  if (n != length(y)) 
    stop("the number of observations in y is not equal to the number of rows of X")
  if (n == 1) 
    stop("the sample size should be greater than 1")
  
  ##################################2.1.OGA_D###################################
  if (is.null(Kn)) 
    K = max(1, min(floor(c1 * sqrt(n/log(p))), p))
  else {
    if ((Kn < 1) | (Kn > p)) 
      stop(paste("Kn should between 1 and ", p, sep = ""))
    if ((Kn - floor(Kn)) != 0) 
      stop("Kn should be a positive integer")
    K = Kn
  }
  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))
  Jhat = sigma2hat = rep(0, K)
  XJhat = matrix(0, n, K)
  u = as.matrix(dy)
  xnorms = sqrt(colSums((dX)^2))
  aSSE = (abs(t(u) %*% dX)/xnorms)
  Jhat[1] = which.max(aSSE)
  XJhat[, 1] = (dX[, Jhat[1]]/sqrt(sum((dX[, Jhat[1]])^2)))
  u = u - XJhat[, 1] %*% t(XJhat[, 1]) %*% u
  sigma2hat[1] = mean(u^2)
  if (K > 1) {
    for (k in 2:K) {
      aSSE = (abs(t(u) %*% dX)/xnorms)
      aSSE[Jhat[1:(k - 1)]] = 0
      Jhat[k] = which.max(aSSE)
      rq = dX[, Jhat[k]] - XJhat[, 1:(k - 1)] %*%
        t(XJhat[, 1:(k - 1)]) %*% dX[, Jhat[k]]
      XJhat[, k] = (rq/sqrt(sum((rq)^2)))
      u = u - XJhat[, k] %*% t(XJhat[, k]) %*% u
      sigma2hat[k] = mean(u^2)
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
  hdic = (n * log(sigma2hat)) + ((1:K) * omega_n * (log(p)))
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
      fit = lm(dy ~ . - 1, data = data.frame(dX[, JDrop1]))
      uDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(uDrop1^2))) + ((kn_hat - 
                                                  1) * omega_n * (log(p)))
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
  
  
  
  return(list(n = n, p = p, Kn = K, J_OGA = Jhat, HDIC = hdic, 
              J_HDIC = J_HDIC, J_Trim = J_Trim, betahat_HDIC = betahat_HDIC, 
              betahat_Trim = betahat_Trim))
}
