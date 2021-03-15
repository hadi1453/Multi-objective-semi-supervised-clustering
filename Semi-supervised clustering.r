library(ggplot2) 
library(reshape2)
library(gplots)
library(RColorBrewer)
library(boot)
library(DAAG)
library(gridExtra)
library(cluster)
library(Gmedian)
library(fclust)
library(tree)
library(MASS)
library(e1071)
set.seed(1234)

########################################################################################################################
### SGD for k-median and Regression non-dominant different K

# auxiliary function to allocate data points into clusters
Allocation_Distance <- function(Z, Xr){
  Dis <- matrix(0,nrow=nrow(Z), ncol=nrow(Xr))
  Amat <-  matrix(0,nrow=nrow(Z), ncol=nrow(Xr))
  for (j in 1:nrow(Xr)) {
    Dis[,j] <-  matrix(apply(Z, 1, FUN = function(x) {sum(abs(x-Xr[j,]))}),,1)
  }
  min.dist <- Dis==apply(Dis, 1, min) # for each point find the cluster with the minimum distance
  for (n in 1:nrow(min.dist)) {
    if(length(which(min.dist[n,]))>1){
      min.dist[n,sample(which(min.dist[n,]),1)] <- FALSE
    }
  }
  Amat[min.dist] <- 1 # assign each point to the cluster with the highest probability
  Amat[!min.dist] <- 0 # remove points from clusters with lower probabilites
  return(list(Distance=Dis, Allocation=Amat))
} 

# auxiliary function to calculate a cost function for clustering
Cost_func <- function (Z, Xr){
  a <- Allocation_Distance(Z,Xr)
  A <- sum(colSums(a$Allocation*a$Distance)/apply(a$Allocation*a$Distance,2,FUN = function(x){length(x[x!=0])}),na.rm=T)
  return(A)
}

# auxiliary function to calculate an error function for regression 
error_func <- function (X, Y, m){
  df <- data.frame(X=X,Y=Y)
  fit <- glm(Y~as.factor(X), data=df, family = "gaussian")
  library(boot)
  set.seed(1)
  CV_fit <- cv.glm(df, fit, cost = function(r,p) sqrt(mean((r - p)^2)), m)
  return(CV_fit$delta[1])
}


# INITIALIZATION
X <- as.matrix(Service_L1_Cost_W1[,15:33]) 
Y <- Service_L1_Cost_W1$Cost_3Y


# Set the parameters:
N <- nrow(X) # number of samples
K <- 11    # number of clusters
D <- ncol(X)    # number of dimensions
C <- 100      # number of initial population

tau.max <- 100 # maximum number of iterations
eta <- 5000 # learning rate
#epsilon <- 5000 # a threshold on the cost (to terminate the process)

error <- data.frame('tau'=rep(1:tau.max,each=C), 'center'=rep(1:C, tau.max))  # to be used to trace the test and training errors in each iteration

tau <- 1 # iteration counter
terminate <- FALSE

Z <- X
YZ <- Y
for (c in 1:(C/10)) {
  for (k in 2:K) {
    assign(paste("Xr",(c-1)*10+(k-1),sep = ""), X[sample(N,k),]) # randomly choose K samples as cluster centers
  }
}

# Main loop
while(!terminate){
  print(tau)
  # check termination criteria:
  terminate <- tau >= tau.max #| Cost_func(Z, Xr)<=epsilon
  
  for (n in 1:C) {
    a <- Allocation_Distance(Z, get(paste("Xr",n,sep = "")))
    Cl <- as.data.frame(a$Allocation)
    for (i in 1:ncol(Cl)) {
      colnames(Cl)[i] <- paste("Cl",i,sep = "")
    }
    assign(paste("Cl",n,sep = ""), Cl)
    assign(paste("Dis",n,sep = ""), as.data.frame(a$Distance))
    
    p.comp <- prcomp(Z)
    X.comp <- as.data.frame(p.comp$x[,1:2])
    X_comp_Cl <- cbind(X.comp,Cl)
    X_comp_Cl.m <- data.frame(X_comp_Cl[,1:2], cluster=NA)
    for (ti in 1:nrow(X_comp_Cl)) {
      for (tj in 3:ncol(X_comp_Cl)) {
        if(X_comp_Cl[ti,tj]==1){X_comp_Cl.m[ti,3] <- colnames(X_comp_Cl)[tj]}
      }
    }
    if (tau %% 100 == 1){
      print(ggplot(data=X_comp_Cl.m, aes(x=PC1, y=PC2, col=cluster)) + 
              geom_point() +
              ggtitle (paste('Clustering (tau=', tau,', C=', n,')')) + theme_minimal())
    }
    
    error[error$tau==tau & error$center==n, 'Cluster'] <- Cost_func(Z, get(paste("Xr",n,sep = "")))
    error[error$tau==tau & error$center==n, 'Regression'] <- error_func(X_comp_Cl.m$cluster, YZ, 10)
  }
  # for each population
  ND=0
  indx <- NULL
  for (m in 1:C) {
    Dom=0
    for (mm in setdiff(1:C, m)) {
      if(error[error$tau==tau & error$center==m,3] > error[error$tau==tau & error$center==mm,3] & error[error$tau==tau & error$center==m,4] > error[error$tau==tau & error$center==mm,4]){Dom = Dom + 1}}
    if(Dom > 0){
      Xc <- matrix(0,nrow(get(paste("Xr",m,sep = ""))),D)
      # for each center:
      for (j in 1:nrow(get(paste("Xr",m,sep = "")))){
        # update the coefficient:
        Zq <- Z[get(paste("Cl",m,sep = ""))[,j]==1,][sample(nrow(Z[get(paste("Cl",m,sep = ""))[,j]==1,]),1),]
        Xc[j,] <- if(colSums(get(paste("Cl",m,sep = "")))[j]>0){get(paste("Xr",m,sep = ""))[j,] + eta/((1+colSums(get(paste("Cl",m,sep = "")))[j])^0.75) * (Zq-get(paste("Xr",m,sep = ""))[j,]) / sum(abs(Zq-get(paste("Xr",m,sep = ""))[j,]))}else{get(paste("Xr",m,sep = ""))[j,]}
        if(any(is.na(Xc[j,]))){
          print(paste("m=",m,sep = ""))
          print(colSums(get(paste("Cl",m,sep = "")))[j])
          print(get(paste("Dis",m,sep = ""))[get(paste("Cl",m,sep = ""))[,j]==1,][1,j])
          Xc[j,] <- get(paste("Xr",m,sep = ""))[j,]
        }
      }
      assign(paste("Xr",m,sep = ""), Xc)
    }else{
      ND = ND + 1
      indx[ND]=m
    }
  }
  print(ggplot(data=error[error$tau==tau,], aes(x=Cluster, y=Regression)) + geom_point(col="blue") +
          geom_point(data=error[error$tau==tau,][indx,], colour="red") + ggtitle(paste('Population pool (tau=', tau,')')) + theme_minimal())
  # update the counter:
  tau <- tau + 1
 
}
