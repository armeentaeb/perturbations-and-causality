


library(MASS)
library(pcalg)
library(xlsx)





#Functions ---------------------------------------

other_methods <- function(link){
  if (dir.exists(link)) {
    #Delete file if it exists
    unlink(link, recursive=TRUE)
    }
  dir.create(link)
  # setup covariance matrices
  s<-0
  covariance_data<-matrix(list(), length(sarr), 1)
  Mj <- solve(diag(p)-B)%*%(Sigma_eps+Sigma_latent)%*%(t(solve(diag(p)-B)))
  covariance_data[[1,1]] <- Mj
  for (j in 1:as.numeric(length(sarr)-1)){
    Mj <- solve(diag(p)-B)%*%(Sigma_latent*as.numeric(1+latent_perturbed[j+1])+Sigma_eps+diag(obs_perturbed[,j+1]))%*%(t(solve(diag(p)-B)))
    covariance_data[[j+1,1]] <- Mj
  }
  M = matrix(0,p,p)
  obs_vec<- c(60,200)
  FD_cd<-matrix(0,1,length(obs_vec))
  TD_cd<-matrix(0,1,length(obs_vec))
  FD_bs<-matrix(0,1,length(obs_vec))
  TD_bs<-matrix(0,1,length(obs_vec))
  
  l <- 1
  for (iter_out in obs_vec){
    str1<-link
    
    num_obs <- iter_out*5
    str1<-paste0(str1,as.character(iter_out*5),"_samples")
    if (dir.exists(str1)) {
      #Delete file if it exists
      unlink(str1,recursive=TRUE)
    }
    dir.create(str1)
    
    for (iter in 1:10){
      
      sarr <- c(num_obs,iter_out*5,iter_out*5,iter_out*5,iter_out*5)
      
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
        t<-matrix(0,p,1)
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
      
      
      environment <- c(rep(1,num_obs),rep(2,iter_out*5),rep(3,iter_out*5),rep(4,iter_out*5),rep(5,iter_out*5))
      # causal dantzig
      Y <- data[,10]
      X <- data[,c(1,2,3,4,5,6,7,8,9)]
      result <- hiddenICP(X, Y, environment)
      N <- matrix(0,p,p)
      
      
      
      FD_cd[l] <- FD_cd[l] + length(setdiff(which(result$pvalues < 10^(-2)),which(B[10,]<0)))
      TD_cd[l] <- TD_cd[l] + length(intersect(which(result$pvalues < 10^(-2)),which(B[10,]<0)))
      
      
      
      
      Latent_data <- t(Gamma_latent%*%t(mvrnorm(n = sarr[1], matrix(0,1,1), diag(1), tol = 1e-6, empirical = FALSE)))
      for (j in 2:length(sarr)){
        Latent_data <- rbind(Latent_data,t(Gamma_latent%*%t(mvrnorm(n = sarr[1], matrix(0,1,1), diag(1), tol = 1e-6, empirical = FALSE)))*sqrt((1+latent_perturbed[j])))
      }
      
      n <- sum(sarr)
      # the noise variance across all the environments
      epsilon <-  sqrt(0.5)*matrix(rnorm(n * p), nrow = n)
      
      ### here creating the perturbation noise, the first is set to zero since observational
      perturb <-  0*matrix(rnorm(sarr[1] * p), nrow = sarr[1])
      for (j in 2:as.numeric(length(sarr))){
        multiplier <- sqrt(obs_perturbed[,j])
        perturb_n <- sweep(matrix(rnorm(sarr[j] *p), ncol = p), 2, multiplier, FUN = "*") ## Christina: intervention variances must not always be the same
        perturb <- rbind((perturb), (perturb_n))
      }
      ### here creating the data
      data <- (perturb+epsilon+Latent_data)%*%(solve(diag(p)-B)) ## Christina: no transpose here
      

      
      #backShift
      network <- backShift(data, environment, ev = 1,threshold =0.75)  
      M <- network$AhatAdjacency
      write.table(M, file =  paste0(str1,"/backShiftSolution",as.character(iter)), append = FALSE, sep = " ", dec = ".",row.names = FALSE, col.names = FALSE)
      
      FD_bs[l] <- FD_bs[l] + length(setdiff(which(M[10,] > 0),c(3,4)))
      TD_bs[l] <- TD_bs[l] + length(intersect(which(M[10,] > 0),c(3,4)))
      
      # 
    }
    l <- l+1
  }
  return(rbind(FD_cd,TD_cd,FD_bs,TD_bs))
  }




# 1 latent variable -- create the latent effect matrix
A <- svd(matrix(rnorm(p^2,mean = 0,sd = 1), ncol=p))
condition <- 0
p <- 10
Sigma_eps <- 0.5*diag(p)
set.seed(2)
A<- t(runif(p,0,1))
A[which(A<0)]<- 0
var_lat <- 1*c(1)
Gamma_latent <-  sqrt(0.3)*t(rbind(A))%*%sqrt(diag((as.numeric(var_lat))))

Sigma_latent <- 0.3*t(A)%*%(A)

B <- rbind(c(0,0.7,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,-0.7),c(0,0,0,0,0,0,0,0,0,0),c(0,0,0,0,0,0,-0.7,0,0,0),c(0,0,0,0,0,0,0,0,-0.7,0),c(0,0,0,0,0,-0.7,0,0,0,0),c(0,0,0,-0.7,-0.7,0,0,0,0,0),c(0,0,-0.7,0,0,0,0,0,0,0),c(0,0,0,-0.7,0,0,0,0,0,0))
B[10,3] <- -0.7
B[,10] <- 0
B[7:9,10] <- -0.7
sarr <- c(1000,1000,1000,1000,1000)
root <-dirname(getwd())

link <- paste0(root,"/files/comparison_no_intervention_DAG/")
# observed perturbation magnitudes
set.seed(2)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+1*runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[10,] <- 0+0*runif(length(sarr))
obs_perturbed[,1] <- 0
latent_perturbed <- 0+0*runif(length(sarr))
latent_perturbed[1]<-0

result_NoInt<- other_methods(link)



link <- paste0(root,"/files/comparison_no_intervention_DAG/")
# observed perturbation magnitudes
set.seed(2)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+1*runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[10,] <- 5+1*runif(length(sarr))
obs_perturbed[,1] <- 0
latent_perturbed <- 0+0*runif(length(sarr))
latent_perturbed[1]<-0

result_TargetInt<- other_methods(link)

link <- paste0(root,"/files/comparison_no_intervention_DAG/")
# observed perturbation magnitudes
set.seed(2)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+1*runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[10,] <- 0+0*runif(length(sarr))
obs_perturbed[,1] <- 0
latent_perturbed <- 1+1*runif(length(sarr))
latent_perturbed[1]<-0
result_LatentInt<- other_methods(link)

link <- paste0(root,"/files/comparison_no_intervention_DAG/")
# observed perturbation magnitudes
set.seed(2)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+1*runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[10,] <- 5+1*runif(length(sarr))
obs_perturbed[,1] <- 0
latent_perturbed <- 1+1*runif(length(sarr))
latent_perturbed[1]<-0
result_AllInt<- other_methods(link)




set.seed(2)
obs_perturbed <- matrix(0,p,length(sarr))
for (k in 1:p){
  obs_perturbed[k,] = 5+5*runif(length(sarr))}# latent perturbation magnitudes
obs_perturbed[10,] <- 5+5*runif(length(sarr))
obs_perturbed[,1] <- 0
latent_perturbed <- 0*runif(length(sarr))
latent_perturbed[1]<-0
result_AllInt_Strong<- other_methods(link)








