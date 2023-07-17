##Compute centerized theta.

centerized_theta(dat,burnin){
  Theta_ij <- dat[["theta_ij"]][(burnin+1):dim(dat[["theta_ij"]])[2]]
  Theta_i <- dat[["theta_i"]][(burnin+1):dim(dat[["theta_i"]])[2]]
  TF_nondup_id <- !duplicated(dat[["TF_names"]])

  sigmaS <- dat[["sigmaS"]][(burnin+1):dim(dat[["sigmaS"]])[2]]
  tau0S <- dat[["tau0S"]][(burnin+1):length(dat[["tau0S"]])]
  tau1S <- dat[["tau1S"]][(burnin+1):length(dat[["tau1S"]])]

  centerized_theta_ij <- data.frame(matrix(ncol=2,nrow=dim(Theta_ij)[1]))


}
