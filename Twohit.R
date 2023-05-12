Twotrim <- function (X, Z, y, KLn = NULL, KDn = NULL, c1 = 5,
                     HDIC_Type_L = "HDBIC", HDIC_Type_D = "HDBIC",
                     c2 = 2, c3 = 2.01, intercept = TRUE) 
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
  if (is.null(KLn)) 
    KLn = max(1, min(floor(c1 * sqrt(n/log(pL))), pL))
  else {
    if ((KLn < 1) | (KLn > pL)) 
      stop(paste("KLn should between 1 and ", pL, sep = ""))
    if ((KLn - floor(KLn)) != 0) 
      stop("KLn should be a positive integer")
    # K = KLn
  }
  dy = y - mean(y)
  dX = apply(X, 2, function(x) x - mean(x))
  Jhat_L = sigma2hat = rep(0, KLn)
  XJhat_L = matrix(0, n, KLn)
  RL = as.matrix(dy)
  xnorms = sqrt(colSums((dX)^2))
  aSSE = (abs(t(RL) %*% dX)/xnorms)
  Jhat_L[1] = which.max(aSSE)
  XJhat_L[, 1] = (dX[, Jhat_L[1]]/sqrt(sum((dX[, Jhat_L[1]])^2)))
  RL = RL - XJhat_L[, 1] %*% t(XJhat_L[, 1]) %*% RL
  sigma2hat[1] = mean(RL^2)
  if (KLn > 1) {
    for (k in 2:KLn) {
      aSSE = (abs(t(RL) %*% dX)/xnorms)
      aSSE[Jhat_L[1:(k - 1)]] = 0
      Jhat_L[k] = which.max(aSSE)
      rq = dX[, Jhat_L[k]] - XJhat_L[, 1:(k - 1)] %*%
        t(XJhat_L[, 1:(k - 1)]) %*% dX[, Jhat_L[k]]
      XJhat_L[, k] = (rq/sqrt(sum((rq)^2)))
      RL = RL - XJhat_L[, k] %*% t(XJhat_L[, k]) %*% RL
      sigma2hat[k] = mean(RL^2)
    }
  }
  
  ##################################1.2.HDIC_L##################################
  #need sigma2hat and Jhat from step 1
  if ((HDIC_Type_L != "HDAIC") & (HDIC_Type_L != "HDBIC") & (HDIC_Type_L != 
                                                         "HDHQ")) 
    stop("HDIC_Type_L should be \"HDAIC\", \"HDBIC\" or \"HDHQ\"")
  if (HDIC_Type_L == "HDAIC") 
    omega_n = c2
  if (HDIC_Type_L == "HDBIC") 
    omega_n = log(n)
  if (HDIC_Type_L == "HDHQ") 
    omega_n = c3 * log(log(n))
  hdic_L = (n * log(sigma2hat)) + ((1:KLn) * omega_n * (log(pL)))
  kln_hat = which.min(hdic_L)
  benchmark = hdic_L[kln_hat]
  # Outpue of HDIC
  J_HDIC_L = sort(Jhat_L[1:kln_hat])
  
  ##################################1.3.Trim_L##################################
  # Initialize output of Trim
  J_Trim_L = Jhat_L[1:kln_hat]
  trim_pos = rep(0, kln_hat)
  if (kln_hat > 1) {
    for (l in 1:(kln_hat - 1)) {
      JDrop1 = J_Trim_L[-l]
      fit = lm(dy ~ . - 1, data = data.frame(dX[, JDrop1]))
      RLDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(RLDrop1^2))) + ((kln_hat - 
                                                  1) * omega_n * (log(pL)))
      if (HDICDrop1 > benchmark) 
        trim_pos[l] = 1
    }
    trim_pos[kln_hat] = 1
    J_Trim_L = J_Trim_L[which(trim_pos == 1)]
  }
  # Outpue of Trim
  J_Trim_L = sort(J_Trim_L)
  
  # Prepare the results
  X_HDIC = as.data.frame(as.matrix(X[, J_HDIC_L]))
  X_Trim = as.data.frame(as.matrix(X[, J_Trim_L]))
  X = data.frame(X)
  colnames(X_HDIC) = names(X)[J_HDIC_L]
  colnames(X_Trim) = names(X)[J_Trim_L]
  if (intercept == TRUE) {
    fit_HDIC = lm(y ~ ., data = X_HDIC)
    fit_Trim = lm(y ~ ., data = X_Trim)
  }
  else {
    fit_HDIC = lm(y ~ . - 1, data = X_HDIC)
    fit_Trim = lm(y ~ . - 1, data = X_Trim)
  }
  betahat_HDIC_L = summary(fit_HDIC)
  betahat_Trim_L = summary(fit_Trim)

  #############################1.4.Obtain residuals#############################
  # projection matrix H_{L,\hat{N}_{L,n}}
  if (intercept == TRUE) {
    A <- cbind(matrix(1,n,1), dX[, result[["J_Trim_L"]]])
  }
  else {
    A <- dX[, result[["J_Trim_L"]]]
  }
  proj <- A %*% solve(t(A) %*% A) %*% t(A)
  eta_hat2 <- ((diag(1,n) - proj) %*% dy)^2
  eta_tilde2 <- apply(rbind(dy, rep(c_n,length(dy))), 2, function(x) max(x))
  r <- log(eta_tilde2)
  
##############################2.Second step Ohit################################
  
  #############################2.0.Setup and Check##############################
  if (!is.matrix(Z)) 
    stop("Z should be a matrix")
  n = nrow(Z)
  pD = ncol(Z)
  if (n != length(r)) 
    stop("the number of observations in r is not equal to the number of rows of Z")
  if (n == 1) 
    stop("the sample size should be greater than 1")
  
  ##################################2.1.OGA_D###################################
  if (is.null(KDn)) 
    KDn = max(1, min(floor(c1 * sqrt(n/log(pD))), pD))
  else {
    if ((KDn < 1) | (KDn > pD)) 
      stop(paste("KDn should between 1 and ", pD, sep = ""))
    if ((KDn - floor(KDn)) != 0) 
      stop("KDn should be a positive integer")
    # K = KDn
  }
  dr = r - mean(r)
  dZ = apply(Z, 2, function(x) x - mean(x))
  Jhat_D = sigma2hat = rep(0, KDn)
  ZJhat_D = matrix(0, n, KDn)
  RD = as.matrix(dr)
  znorms = sqrt(colSums((dZ)^2))
  aSSE = (abs(t(RD) %*% dZ)/znorms)
  Jhat_D[1] = which.max(aSSE)
  ZJhat_D[, 1] = (dZ[, Jhat_D[1]]/sqrt(sum((dZ[, Jhat_D[1]])^2)))
  RD = RD - ZJhat_D[, 1] %*% t(ZJhat_D[, 1]) %*% RD
  sigma2hat[1] = mean(RD^2)
  if (KDn > 1) {
    for (k in 2:KDn) {
      aSSE = (abs(t(RD) %*% dZ)/znorms)
      aSSE[Jhat_D[1:(k - 1)]] = 0
      Jhat_D[k] = which.max(aSSE)
      rq = dZ[, Jhat_D[k]] - ZJhat_D[, 1:(k - 1)] %*%
        t(ZJhat_D[, 1:(k - 1)]) %*% dZ[, Jhat_D[k]]
      ZJhat_D[, k] = (rq/sqrt(sum((rq)^2)))
      RD = RD - ZJhat_D[, k] %*% t(ZJhat_D[, k]) %*% RD
      sigma2hat[k] = mean(RD^2)
    }
  }
  
  ##################################2.2.HDIC_D##################################
  #need sigma2hat and Jhat from step 1
  if ((HDIC_Type_D != "HDAIC") & (HDIC_Type_D != "HDBIC") & (HDIC_Type_D != 
                                                         "HDHQ")) 
    stop("HDIC_Type_D should be \"HDAIC\", \"HDBIC\" or \"HDHQ\"")
  if (HDIC_Type_D == "HDAIC") 
    omega_n = c2
  if (HDIC_Type_D == "HDBIC") 
    omega_n = log(n)
  if (HDIC_Type_D == "HDHQ") 
    omega_n = c3 * log(log(n))
  hdic_D = (n * log(sigma2hat)) + ((1:KDn) * omega_n * (log(pD)))
  kdn_hat = which.min(hdic_D)
  benchmark = hdic_D[kdn_hat]
  # Output of HDIC
  J_HDIC_D = sort(Jhat_D[1:kdn_hat])
  
  ##################################2.3.Trim_D##################################
  # Initialize output of Trim
  J_Trim_D = Jhat_D[1:kdn_hat]
  trim_pos = rep(0, kdn_hat)
  if (kdn_hat > 1) {
    for (l in 1:(kdn_hat - 1)) {
      JDrop1 = J_Trim_D[-l]
      fit = lm(dr ~ . - 1, data = data.frame(dZ[, JDrop1]))
      RDDrop1 = fit$residuals
      HDICDrop1 = (n * log(mean(RDDrop1^2))) + ((kdn_hat - 
                                                  1) * omega_n * (log(pD)))
      if (HDICDrop1 > benchmark) 
        trim_pos[l] = 1
    }
    trim_pos[kdn_hat] = 1
    J_Trim_D = J_Trim_D[which(trim_pos == 1)]
  }
  # Outpue of Trim
  J_Trim_D = sort(J_Trim_D)
  
  # Prepare the results
  Z_HDIC = as.data.frame(as.matrix(Z[, J_HDIC_D]))
  Z_Trim = as.data.frame(as.matrix(Z[, J_Trim_D]))
  Z = data.frame(Z)
  colnames(Z_HDIC) = names(Z)[J_HDIC_D]
  colnames(Z_Trim) = names(Z)[J_Trim_D]
  if (intercept == TRUE) {
    fit_HDIC = lm(y ~ ., data = Z_HDIC)
    fit_Trim = lm(y ~ ., data = Z_Trim)
  }
  else {
    fit_HDIC = lm(y ~ . - 1, data = Z_HDIC)
    fit_Trim = lm(y ~ . - 1, data = Z_Trim)
  }
  betahat_HDIC_D = summary(fit_HDIC)
  betahat_Trim_D = summary(fit_Trim)
  
  return(list(n = n, pL = pL, pD = pD, KLn = KLn, KDn = KDn, J_OGA_L = Jhat_L,
              J_OGA_D = Jhat_D, HDIC_L = hdic_L, HDIC_D = hdic_D,
              J_HDIC_L = J_HDIC_L, J_HDIC_D = J_HDIC_D,
              J_Trim_L = J_Trim_L, J_Trim_D = J_Trim_D,
              betahat_HDIC_L = betahat_HDIC_L, betahat_HDIC_D = betahat_HDIC_D,
              betahat_Trim_L = betahat_Trim_L, betahat_Trim_D = betahat_Trim_D))
}
