###For the simulation

simu_data_generator<-function(M=25,Mc=25,N=10000,n=20,mu=-2,tau2=0.5,tau2_1=5,sigma_i=1){
    simu_data_list<-list()
    for(i in 1:M){
      simu_data_list[[i]]<-list()
      simu_data_list[[i]][["theta_i"]]=rnorm(1,mu,sqrt(tau2))
      simu_data_list[[i]][["sigma_i"]]=sigma_i
      simu_data_list[[i]][["theta_ij"]]=rnorm(n,simu_data_list[[i]][["theta_i"]],sqrt(simu_data_list[[i]][["sigma_i"]]))
      x_ij_vec<-c()
      for(k in 1:n){
        x_ij_vec<-c(x_ij_vec,rbinom(1,N,dlogit(simu_data_list[[i]][["theta_ij"]][k])))
      }
      simu_data_list[[i]][["x_ij"]]<-x_ij_vec
      simu_data_list[[i]][["n_ij"]]<-rep(N,length(x_ij_vec))
    }
    for(j in (M+1):(M+Mc)){
      simu_data_list[[j]]<-list()
      simu_data_list[[j]][["theta_i"]]=rnorm(1,mu,sqrt(tau2))
      simu_data_list[[j]][["theta_ij"]]=rnorm(1,simu_data_list[[j]][["theta_i"]],sqrt(tau2_1))
      simu_data_list[[j]][["x_ij"]]<-rbinom(1,N,dlogit(simu_data_list[[j]][["theta_ij"]]))
      simu_data_list[[i]][["n_ij"]]<-c(N)
    }
    return(simu_data_list)
}

simu_data_generator_from_vec<-function(TF_quantity_vec,N=10000,mu=0,tau2=0.5,tau2_1=5,sigma_i=1){
  simu_data_list<-c()
  for(i in 1:length(TF_quantity_vec)){
    simu_data_list[[i]]<-list()
    simu_data_list[[i]][["theta_i"]]=rnorm(1,mu,sqrt(tau2))
    simu_data_list[[i]][["sigma_i"]]=sigma_i
    if(TF_quantity_vec[i]==1){
      simu_data_list[[i]][["theta_ij"]]=rnorm(TF_quantity_vec[i],simu_data_list[[i]][["theta_i"]],
                                              sqrt(tau2_1))
    }else{
      simu_data_list[[i]][["theta_ij"]]=rnorm(TF_quantity_vec[i],simu_data_list[[i]][["theta_i"]],
                                              sqrt(simu_data_list[[i]][["sigma_i"]]))
    }
    x_ij_vec<-c()
    for(k in 1:TF_quantity_vec[i]){
      x_ij_vec<-c(x_ij_vec,rbinom(1,N,dlogit(simu_data_list[[i]][["theta_ij"]][k])))
    }
    simu_data_list[[i]][["x_ij"]]<-x_ij_vec
    simu_data_list[[i]][["n_ij"]]<-rep(N,length(x_ij_vec))
  }
  return(simu_data_list)
}


dlogit<-function(theta_ij){
  return(exp(theta_ij)/(1+exp(theta_ij)))
}

simu_data_constructor<-function(rdata){
  xij<-c()
  nij<-c()
  label_vec<-c()
  theta_i<-c()
  theta_ij<-c()
  for(i in 1:length(rdata)){
    xij<-c(xij,rdata[[i]]$x_ij)
    nij<-c(nij,rdata[[i]]$n_ij)
    theta_i<-c(theta_i,rdata[[i]]$theta_i)
    theta_ij<-c(theta_ij,rdata[[i]]$theta_ij)
    label_vec<-c(label_vec,rep(i,length(rdata[[i]]$x_ij)))
  }
  return(list(xij=xij,nij=nij,label_vec=label_vec,theta_i=theta_i,theta_ij=theta_ij))
}
