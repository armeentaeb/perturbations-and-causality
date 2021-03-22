# libraries
library(MASS)
library(pcalg)
library(xlsx)





# Setup model parameters ----------------------------------------------
# connectivity matrix
p <- 10
sarr <- c(200,5,5,5,5,5,5)
B <- generate_DAG(0.9)

# 1 latent variable -- create the latent effect matrix
A <- svd(matrix(rnorm(p^2,mean = 0,sd = 1), ncol=p))
condition <- 0
Sigma_eps <- 0.5*diag(p)
while (!condition){
  A<- t(runif(p,0,1))
  A[which(A<0.1)]<- 0
  Sigma_latent <- 0.2*t(A)%*%(A)
  M<- t(diag(p)-B)%*%solve(Sigma_eps+Sigma_latent)%*%t(diag(p)-B)
  latent_effect <- t(diag(p)-B)%*%solve(Sigma_eps)%*%t(diag(p)-B) - M;
  if (max(latent_effect[which((t(diag(p)-B)%*%solve(Sigma_eps)%*%t(diag(p)-B))==0)]) > 0){
    condition <- 1
  }
}


## 2 latent variable -- create the latent effect matrix
A <- svd(matrix(rnorm(p^2,mean = 0,sd = 1), ncol=p))
condition <- 0
Sigma_eps <- 0.5*diag(p)
while (!condition){
  A1<- t(runif(p,0,1))
  A1[which(A1<0.1)]<- 0
  A2<- t(runif(p,0,1))
  A2[which(!A1==0)]<-0
  A2[which(A2<0.1)]<- 0
  Sigma_latent2 <- 0.2/sqrt(2)*t(rbind(A1,A2))%*%(rbind(A1,A2))
  M<- t(diag(p)-B)%*%solve(Sigma_eps+Sigma_latent2)%*%t(diag(p)-B)
  latent_effect <- t(diag(p)-B)%*%solve(Sigma_eps)%*%t(diag(p)-B) - M;
  if (max(latent_effect[which((t(diag(p)-B)%*%solve(Sigma_eps)%*%t(diag(p)-B))==0)]) > 0){
    condition <- 1
  }
}


#################################
## Setting 1: perturbation 5, 1 hiddden variable
root <- dirname(getwd())
link <- paste0(root,"/files/mismatch_nonlinear/")
dir.create(link)
str1<-link
# observed perturbation magnitudes
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[,1]<-0
# latent perturbation magnitudes
latent_perturbed <- c(t(1+runif(length(sarr))/2))


# setup covariance matrices
covariance_data<-matrix(list(), length(sarr), 1)
var_lat <- 1*c(1)
Gamma_latent <-  sqrt(0.3)*t(rbind(A1))%*%sqrt(diag((as.numeric(var_lat))))
xi <- 0.2





for (iter_out in c(1,2,4,16,64)){
  sarr <- c(300,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
  str1<-link
  str1<-paste0(str1,as.character(iter_out*5),"_samples")
  dir.create(str1)
  for (iter in 1:10){
   
    Latent_data <- Gamma_latent%*%t(mvrnorm(n =sarr[1], matrix(0,1,1), diag(1), tol = 1e-6, empirical = FALSE))
    epsilon <- 1/sqrt(2)*t(mvrnorm(n =sarr[1], mu,diag(p), tol = 1e-6, empirical = FALSE))
    
    order<- topological_sort(t(adjacency))
    data <- matrix(0,p,sarr[1])
    for (i in 1:p){
      data[as.numeric(order[i]),] <-  epsilon[as.numeric(order[i]),] + (Latent_data[as.numeric(order[i]),] + B[as.numeric(order[i]),]%*%data)+xi*(Latent_data[as.numeric(order[i]),] + B[as.numeric(order[i]),]%*%data)^2
    }
    data <- t(data)
    for (j in 2:as.numeric(length(sarr))){
      Latent_data <- sqrt(1+latent_perturbed[j])*Gamma_latent%*%t(mvrnorm(n =sarr[j], matrix(0,1,1), diag(1), tol = 1e-6, empirical = FALSE))
      epsilon <- 1/sqrt(2)*t(mvrnorm(n =sarr[j], matrix(0,p,1),diag(p), tol = 1e-6, empirical = FALSE))
      perturbed_data <- sqrt(obs_perturbed[j])*t(mvrnorm(n =sarr[j], matrix(0,p,1),diag(p), tol = 1e-6, empirical = FALSE))
      
      data_n <- matrix(0,p,as.numeric(sarr[j]))
      for (i in 1:p){
        data_n[as.numeric(order[i]),] <-  epsilon[as.numeric(order[i]),] + perturbed_data[as.numeric(order[i]),]+(Latent_data[as.numeric(order[i]),] + B[as.numeric(order[i]),]%*%data_n)+xi*(Latent_data[as.numeric(order[i]),] + B[as.numeric(order[i]),]%*%data_n)^2
      }
      data <- rbind((data), t(data_n))
    }
    
    write.table(data, file =  paste0(str1,"/overall_data",as.character(iter)),append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
    score <- new("GaussL0penIntScore", data,lambda = 1*log(nrow(data)))
    ges.fit <- ges(score)
    A <- ges.fit$essgraph
    M <- matrix(0,p,p)
    for (i in 1:p){M[i,A[['.in.edges']][[names(A$.in.edges)[i]]]]<-1}
    res1<- pdag2allDags(M)
    amat1 <- matrix(res1$dags[1,],p,p, byrow = TRUE)
    if (dim(res1$dags)[1]>1){
      for (j in 2:dim(res1$dags)[1]){
        amat1 <- rbind(amat1,matrix(res1$dags[j,],p,p, byrow = TRUE))
      }
    }
    write.table(amat1, file =  paste0(str1,"/equivalent_DAGS_GES",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  }
}






