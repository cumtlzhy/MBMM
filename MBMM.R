customKmeans<-function(dataset,k){
  if(is.na(dataset) || is.na(k)){
    stop("You must input valid parameters!!")
  }
  
  Eudist<-function(x,y){
    distance<-sqrt(sum((x-y)^2))
    return (distance)
  }
  rows.dataset<-nrow(dataset)
  cols.dataset<-ncol(dataset)
  continue.change=TRUE
  set.seed(0)
  random_number <- sample.int(rows.dataset,size = k)
  initPoint<-dataset[random_number,]
  
  formerPoint<-initPoint
  iterPoint<-matrix(0,nrow = k,ncol = ncol(dataset))
  
  error.matrix<-matrix(0,nrow=rows.dataset,ncol=k)
  while(continue.change){
    
    cluster.matrix<-matrix(0,nrow=rows.dataset,ncol=k)
    for(i in 1:rows.dataset){
      for(j in 1:k){
        error.matrix[i,j]<-Eudist(dataset[i,],formerPoint[j,])
      }
      
    }
    
    
    for(i in 1:rows.dataset){
      cluster.matrix[i,which.min(error.matrix[i,])]<-1
      
    }
   
    
    for(i in 1:k){
      iterPoint[i,]<-apply(dataset[which(cluster.matrix[,i] == 1),],2,"mean")
      
    }
    all.true<-c()
    
    for(i in 1:k){
      if(all(formerPoint[i,] == iterPoint[i,]) == T){
        all.true[i]<-TRUE
      }
      
    }
    formerPoint = iterPoint
    continue.change=ifelse(all(all.true) == T,F,T)
    
  }
  colnames(iterPoint)<-colnames(dataset)
  
  cluster_lable <- rep(1,rows.dataset)
  for(i in 1:rows.dataset){
    
    cluster_lable[i]<-which(cluster.matrix[i,] == 1)
    
  }
  
  yita_m <- (table(cluster_lable) + 1) / (length(cluster_lable) + k) 
  out=list()
  out[["centers"]] <- iterPoint
  out[["distance"]] <- error.matrix
  out[["cluster_lable"]] <- cluster_lable   
  
  out[["yita_m"]] <- yita_m
  
  out[["random_number"]] <- random_number
  return(out)
  
}


moments_estimate <- function(x,y,k){
  J <-ncol(x)
  
  a = matrix(0, k,J) 
  b = matrix(0, k,J)
  means <- matrix(0,k,J)
  vars <- matrix(0,k,J)
  for (j in 1:J) {
    means[,j] <- tapply(x[,j], y, mean)
    vars[,j] <- tapply(x[,j], y, var)
    
  }
  com <- means * (1 - means) / vars - 1 
  a <- means * com
  b <- (1 - means) * com
  out=list()
  out[["alpha"]] <- a
  out[["beta"]] <- b
  return(out)
}


expect_z_estimate <- function(x,k,yita,a,b){
  J <- ncol(x)
  n <- nrow(x)
  
  z <- matrix(0,n,k)
  for (i in 1:n) {
    zik_numerator <- matrix(0,1,k)  
    
    for (m in 1:k){
      D_beta <- matrix(0,J,1) 
      for (j in 1:J){
        # browser()
        
        D_beta[j] <- dbeta(x[i,j],a[m,j],b[m,j]) 
        
      }
      D_beta[D_beta == 0] <- .Machine$double.xmin
      zik_numerator[m] <- yita[m]*prod(D_beta)  
    } 
    z[i,] <- zik_numerator/sum(zik_numerator)
    
  }
  y <- max.col(z)  
  out=list()
  out[["expect_z"]] <- z
  out[["lable"]] <- y
  return(out)
  
}


log_observed_likelihood <- function(x,yita,a,b,k){
  J <- ncol(x)
  n <- nrow(x)
  likelihood_of_each_i_for_all_m <- matrix(0,n,k)
  for (i in 1:n) {
    likelihood_i_m <- matrix(0,1,k)  
    
    for (m in 1:k){
      D_beta <- matrix(0,J,1) 
      for (j in 1:J){
        
        D_beta[j] <-  dbeta(x[i,j],a[m,j],b[m,j]) 
      }
      likelihood_i_m[m] <- yita[m]*prod(D_beta)  
      
    }
    likelihood_i_m[likelihood_i_m == 0] <- .Machine$double.xmin  
    likelihood_of_each_i_for_all_m[i,] <- likelihood_i_m  
    
  }
  log_observed_likelihood <- sum(log(apply(likelihood_of_each_i_for_all_m,1,sum)))
  
  return(log_observed_likelihood)
  
}


ICL_BIC <- function(expect_z,n,log_likelihood,J,k){
  l <- J*2*k+k-1
 
  Bayesian_IC <- -2*log_likelihood + l*log(n)
  expect_z[expect_z == 0] <- .Machine$double.xmin  
  
  Entropy_z <- -sum(mapply(function(x,y) x*y, expect_z,log(expect_z)))
  ICL_bayesian_infor_criterion <- Bayesian_IC +2*Entropy_z
  return(ICL_bayesian_infor_criterion)
}


MBMM <- function(x,eps = 1e-5,k){
  
  J <- ncol(x)
  n <- nrow(x)
  index <- 1
  
  iterations <- 1 
  
  observed_likelihood <- 0 
  result_kmeans <- customKmeans(dataset = x,k = k)
  y_0 <- result_kmeans$cluster_lable
  
  yita_0 <- result_kmeans$yita_m
  random_number <- result_kmeans$random_number
  
  result_moment_0 <- moments_estimate(x = x, y = y_0, k = k)
  a_0 <- result_moment_0$alpha
  b_0 <- result_moment_0$beta
  z_y_estimate_1 <- expect_z_estimate(x = x,k = k,yita = yita_0, a = a_0, b = b_0)
  update_z_1 <- z_y_estimate_1$expect_z
  update_y_1 <- z_y_estimate_1$lable
  
  update_yita_1 <- colSums(update_z_1)/n
  
  result_moment_1 <- moments_estimate(x = x, y = update_y_1, k = k)
  update_a_1 <- result_moment_1$alpha
  update_b_1 <- result_moment_1$beta
  
  log_observed_likelihood_1 <- log_observed_likelihood(x = x ,yita = update_yita_1,a = update_a_1, b = update_b_1,k = k)
  
  likelihood_pre <- log_observed_likelihood_1
  observed_likelihood[index] <- log_observed_likelihood_1
  
  repeat{
    
    
    # E-step：
    
    z_y_estimate <- expect_z_estimate(x = x,k = k,yita = update_yita_1,a = update_a_1, b = update_b_1)
    update_z <- z_y_estimate$expect_z
    update_y <- z_y_estimate$lable
    # M-step：
    
    update_yita <- colSums(update_z)/n 
    
    result_moment <- moments_estimate(x = x, y = update_y, k = k)
    update_a <- result_moment$alpha
    update_b <- result_moment$beta
    
    update_log_observed_likelihood <- log_observed_likelihood(x = x ,yita = update_yita,a = update_a, b = update_b,k = k)
    
    observed_likelihood[index+1] <- update_log_observed_likelihood
    
    iterations[index+1] <- index + 1
    
    f.e <- update_log_observed_likelihood - likelihood_pre
    
    if (f.e < eps ){
      
      compute_ICL_BIC <- ICL_BIC(expect_z = update_z_1, n = n, log_likelihood = likelihood_pre, J = J, k = k)
      
      break}
    
    update_yita_1 <- update_yita
    update_a_1 <- update_a
    update_b_1 <-  update_b
    update_y_1 <- update_y
    update_z_1 <- update_z
    
    likelihood_pre <- update_log_observed_likelihood
    
    index <- index+1}
  
  out = list()
  assign(paste("iterations", k, sep="_"),iterations)
  out[[paste("iterations", k, sep="_")]] <- get(paste("iterations", k, sep="_"))
  out[["observed_likelihood"]] <- observed_likelihood
  out[["alpha"]] <- update_a_1
  out[["beta"]] <- update_b_1
  out[["it_max"]] <- index + 1
  out[["expect_z"]] <- update_z_1
  out[["y_lable"]] <- update_y_1
  out[["yita_k"]] <- update_yita_1
  out[["random_number"]] <-  random_number
  out[["ICL_BIC"]] <- compute_ICL_BIC
  
  return(out)
}
model_select <- function(cmax,x,eps = 1e-5){
  
  criterion <- 0
  out = list()
  for (c in 2:cmax) {
    result_MBMM <- MBMM(x = x,eps = 1e-5, k = c)
    assign(paste("model", c, sep="_"),result_MBMM)
    out[[paste("model", c, sep="_")]] <- get(paste("model", c, sep="_"))
    
    criterion[c-1] <- result_MBMM$ICL_BIC
    correct_model <- which.min(criterion) + 1
   
  }
  
  out[["ICL_BIC_OF_EACH_MODEL"]] <- criterion
  out[["correct_model"]] <- correct_model
  
  return(out)
}
