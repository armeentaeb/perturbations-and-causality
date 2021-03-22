

# libraries
library(MASS)
library(pcalg)
library(xlsx)




B <- rbind(c(0,0.7,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,-0.7,0,0,0),c(0,0,0,0,0,0,0,0,-0.7,0),c(0,0,0,0,0,-0.7,0,0,0,0),c(0,0,0,-0.7,-0.7,0,0,0,0,0),c(0,0,-0.7,0,0,0,0,0,0,0),c(0,0,0,-0.7,0,0,0,0,0,0))



# Setup model parameters ----------------------------------------------
# connectivity matrix
p <- 10
sarr <- c(300,5,5,5,5,5,5)

link <- paste0(root,"/comparison_to_two_stage_non_dense/")

# 1 latent variable -- create the latent effect matrix
set.seed(1)
Sigma_eps <- 0.5*diag(p)

#################################
## non-deseness, three hidden variables
p<-10
B <- rbind(c(0,0.7,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,-0.7,0,0,0),c(0,0,0,0,0,0,0,0,-0.7,0),c(0,0,0,0,0,-0.7,0,0,0,0),c(0,0,0,-0.7,-0.7,0,0,0,0,0),c(0,0,-0.7,0,0,0,0,0,0,0),c(0,0,0,-0.7,0,0,0,0,0,0))
Sigma_eps <- 0.5*diag(p)

root <- dirname(getwd())
link <- paste0(root,"/files/comparison_to_two_stage_non_dense/")
dir.create(link)
str1<-link
sarr <- c(10000,5,5,5,5,5,5)
set.seed(1)

# observed perturbation magnitudes
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 2+runif(length(sarr))}# latent perturbation magnitudes

obs_perturbed[,1] <- 0
#obs_perturbed = c(0,2+runif(1),2+runif(1),2+runif(1),2+runif(1),2+runif(1),2+runif(1))
# latent perturbation magnitudes
latent_perturbed <- runif(length(sarr))/2
latent_perturbed[1]<-0

A1 <- as.matrix(c(0,0,0,0,0,1,0,0,0,0))
A2 <- as.matrix(c(0,0,0,0,1,0,0,0,0,0))
A3<- t(runif(p,0,1))
A3[which(A3<0.5)]<- 0
Sigma_latent <- (A1)%*%t(A1)+(A2)%*%t(A2)+t(A3)%*%(A3)

# setup covariance matrices
covariance_data<-matrix(list(), length(sarr), 1)
Mj <- solve(diag(p)-B)%*%(Sigma_eps+Sigma_latent)%*%(t(solve(diag(p)-B)))
covariance_data[[1,1]] <- Mj
for (j in 1:as.numeric(length(sarr)-1)){
  Mj <- solve(diag(p)-B)%*%(Sigma_latent*as.numeric(1+latent_perturbed[j+1])+Sigma_eps+diag(obs_perturbed[,j+1]))%*%(t(solve(diag(p)-B)))
  covariance_data[[j+1,1]] <- Mj
}

num_obs <- 1000
for (iter_out in c(200)){
  sarr <- c(num_obs,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
  str1<-link
  str1<-paste0(str1,as.character(iter_out*5),"_samples")
  dir.create(str1)
  for (iter in 1:10){
    mu <- matrix(0,p,1)
    obs_data <- (mvrnorm(n = as.numeric(sarr[1]), mu, covariance_data[[1,1]], tol = 1e-6, empirical = FALSE))
    data <- (mvrnorm(n = as.numeric(sarr[1]), mu, covariance_data[[1,1]], tol = 1e-6, empirical = FALSE))
    
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
    write.table(amat1, file =  paste0(str1,"/equivalent_DAGS_GES_obs",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
    
    
    for (j in 2:as.numeric(length(sarr))){
      data1 <- (mvrnorm(n = as.numeric(sarr[j]), mu, covariance_data[[j,1]], tol = 1e-6, empirical = FALSE))
      data <- rbind(data,data1)
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
  }
  
}
