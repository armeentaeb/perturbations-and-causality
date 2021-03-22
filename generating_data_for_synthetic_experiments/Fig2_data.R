# libraries
library(MASS)
library(pcalg)
library(xlsx)




B <- rbind(c(0,0.7,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,-0.7,0,0,0),c(0,0,0,0,0,0,0,0,-0.7,0),c(0,0,0,0,0,-0.7,0,0,0,0),c(0,0,0,-0.7,-0.7,0,0,0,0,0),c(0,0,-0.7,0,0,0,0,0,0,0),c(0,0,0,-0.7,0,0,0,0,0,0))



# Setup model parameters ----------------------------------------------
# connectivity matrix
p <- 10
sarr <- c(300,5,5,5,5,5,5)

# 1 latent variable -- create the latent effect matrix
set.seed(1)
A <- svd(matrix(rnorm(p^2,mean = 0,sd = 1), ncol=p))
Sigma_eps <- 0.5*diag(p)
  A<- t(runif(p,0,1))
  A[which(A<0.5)]<- 0
  Sigma_latent <- 0.3*t(A)%*%(A)
 


## 2 latent variable -- create the latent effect matrix
A <- svd(matrix(rnorm(p^2,mean = 0,sd = 1), ncol=p))
condition <- 0
Sigma_eps <- 0.5*diag(p)
A1<- t(runif(p,0,1))
A1[which(A1<0.5)]<- 0
A2<- t(runif(p,0,1))
A2[which(A2<0.5)]<- 0
Sigma_latent2 <- 0.3/sqrt(2)*t(rbind(A1,A2))%*%(rbind(A1,A2))
  

#################################
## Setting 1: perturbation 5, 1 hiddden variable
root <- dirname(getwd())
link <- paste0(root,"/files/Fig2_zeta5_h1_psi1_latent/")
sarr <- c(300,5,5,5,5,5,5)
if (dir.exists(link)) {
  #Delete file if it exists
  unlink(link, recursive=TRUE)
}
dir.create(link)
str1<-link
# observed perturbation magnitudes
set.seed(1)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
obs_perturbed[k,] = 5+runif(length(sarr))}
obs_perturbed[,1]<-0
# latent perturbation magnitudes
latent_perturbed = 0.5+0.5*runif(length(sarr))
latent_perturbed[1] <- 0

# setup covariance matrices
covariance_data<-matrix(list(), length(sarr), 1)
Mj <- solve(diag(p)-B)%*%(Sigma_eps+Sigma_latent)%*%(t(solve(diag(p)-B)))
covariance_data[[1,1]] <- Mj
for (j in 1:as.numeric(length(sarr)-1)){
  Mj <- solve(diag(p)-B)%*%(Sigma_latent*as.numeric(1+latent_perturbed[j+1])+Sigma_eps+diag(obs_perturbed[,j+1]))%*%(t(solve(diag(p)-B)))
  covariance_data[[j+1,1]] <- Mj
}


for (iter_out in c(1,2,4,16,64)){
  sarr <- c(300,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
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
    write.table(amat1, file =  paste0(str1,"/equivalent_DAGS_GES",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  }
}

link <- paste0(root,"/files/Fig2_zeta5_h2_psi1_latent/")
sarr <- c(300,5,5,5,5,5,5)
if (dir.exists(link)) {
  #Delete file if it exists
  unlink(link, recursive=TRUE)
}
dir.create(link)
str1<-link
# observed perturbation magnitudes
set.seed(1)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+runif(length(sarr))}
obs_perturbed[,1]<-0
# latent perturbation magnitudes
latent_perturbed = 0.5+0.5*runif(length(sarr))
latent_perturbed[1] <- 0


# setup covariance matrices
covariance_data<-matrix(list(), length(sarr), 1)
Mj <- solve(diag(p)-B)%*%(Sigma_eps+Sigma_latent)%*%(t(solve(diag(p)-B)))
covariance_data[[1,1]] <- Mj
for (j in 1:as.numeric(length(sarr)-1)){
  Mj <- solve(diag(p)-B)%*%(Sigma_latent*as.numeric(1+latent_perturbed[j+1])+Sigma_eps+diag(obs_perturbed[,j+1]))%*%(t(solve(diag(p)-B)))
  covariance_data[[j+1,1]] <- Mj
}


for (iter_out in c(1,2,4,16,64)){
  sarr <- c(300,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
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
    write.table(amat1, file =  paste0(str1,"/equivalent_DAGS_GES",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  }
}




#################################
## Setting 2: perturbation 2, 1 hiddden variable
link <- paste0(root,"/files/Fig2_zeta2_h1_psi1_latent/")
if (dir.exists(link)) {
  #Delete file if it exists
  unlink(link, recursive=TRUE)
}
dir.create(link)
str1<-link
# observed perturbation magnitudes
for (k in 1:p){
  obs_perturbed[k,] = 2+runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[,1]<-0
latent_perturbed = 0.5+0.5*runif(length(sarr))
latent_perturbed[1] <- 0


# setup covariance matrices
covariance_data<-matrix(list(), length(sarr), 1)
Mj <- solve(diag(p)-B)%*%(Sigma_eps+Sigma_latent)%*%(t(solve(diag(p)-B)))
covariance_data[[1,1]] <- Mj
for (j in 1:as.numeric(length(sarr)-1)){
  Mj <- solve(diag(p)-B)%*%(Sigma_latent*as.numeric(1+latent_perturbed[j+1])+Sigma_eps+diag(obs_perturbed[,j+1]))%*%(t(solve(diag(p)-B)))
  covariance_data[[j+1,1]] <- Mj
}


for (iter_out in c(1,2,4,16,64)){
  sarr <- c(300,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
  str1<-link
  str1<-paste0(str1,as.character(iter_out*5),"_samples/")
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
    write.table(amat1, file =  paste0(str1,"/equivalent_DAGS_GES",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  }
}





#################################
## Setting 5: perturbation 5,1 hiddden variable,do_interventions
link <- paste0(root,"/files/Fig2_zeta5_h1_psi1_latent_do/")
if (dir.exists(link)) {
  #Delete file if it exists
  unlink(link, recursive=TRUE)
}
dir.create(link)
str1<-link
# observed perturbation magnitudes
for (k in 1:p){
  obs_perturbed[k,] = 5+runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[,1]<-0
latent_perturbed = 0.5+0.5*runif(length(sarr))

do_magnitude <- 5


# setup covariance matrices
covariance_data<-matrix(list(), length(sarr), 1)
Mj <- solve(diag(p)-B)%*%(Sigma_eps+Sigma_latent)%*%(t(solve(diag(p)-B)))
covariance_data[[1,1]] <- Mj
#for (j in 1:as.numeric(length(sarr)-2)){
j <- 1
Mj <- solve(diag(p)-B)%*%(Sigma_latent*as.numeric(1+latent_perturbed[j+1])+Sigma_eps+diag(obs_perturbed[,j+1]))%*%(t(solve(diag(p)-B)))
covariance_data[[j+1,1]] <- Mj
#}
# the do-intervention environments
j <- 3
ind_do1<- sample.int(p, 2, replace = FALSE)
do_perturb_mat <- matrix(0,p,p)
do_perturb_mat[ind_do1[1],ind_do1[1]]<- do_magnitude
do_perturb_mat[ind_do1[2],ind_do1[2]]<- do_magnitude
Proj_offdo <- diag(p)
Proj_offdo[ind_do1[1],] <- 0
Proj_offdo[ind_do1[2],] <- 0
N <- Sigma_eps+Sigma_latent
N <- Proj_offdo%*%N%*%Proj_offdo+do_perturb_mat
Mj <- solve(diag(p)-Proj_offdo%*%B)%*%N%*%(t(solve(diag(p)-Proj_offdo%*%B)))
covariance_data[[j,1]] <- Mj
ind_do <- ind_do1

for (j in 4:as.numeric(length(sarr))){
ind_do1<- sample.int(p, 2, replace = FALSE)
do_perturb_mat <- matrix(0,p,p)
do_perturb_mat[ind_do1[1],ind_do1[1]]<- do_magnitude
do_perturb_mat[ind_do1[2],ind_do1[2]]<- do_magnitude
Proj_offdo <- diag(p)
Proj_offdo[ind_do1[1],] <- 0
Proj_offdo[ind_do1[2],] <- 0
N <- Sigma_eps+Sigma_latent
N <- Proj_offdo%*%N%*%Proj_offdo+do_perturb_mat
Mj <- solve(diag(p)-Proj_offdo%*%B)%*%N%*%(t(solve(diag(p)-Proj_offdo%*%B)))
covariance_data[[j,1]] <- Mj
ind_do <- rbind(ind_do,ind_do1)
}
write.table(ind_do, file =  paste0(str1,"/do_locs"), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)


for (iter_out in c(1,2,4,16,64)){
  sarr <- c(300,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
  str1<-link
  str1<-paste0(str1,as.character(iter_out*5),"_samples/")
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
    write.table(amat1, file =  paste0(str1,"/equivalent_DAGS_GES",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
  }
}

