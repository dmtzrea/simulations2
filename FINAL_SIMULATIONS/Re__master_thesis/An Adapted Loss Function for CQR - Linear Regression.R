##########################################################################
################### An Adapted Loss Function for #########################
################### Censored Quantile Regression #########################
###################  R Code - Linear Regression  #########################
##########################################################################

# Date: 8/12/2017
# Structure: 
#     - 1/ Literature Functions
#     - 2/ Function for NEW procedure
#     - 3/ Example on one simulated data set
#     - 4/ Example of simulation study (using parallelization on servers)


rm(list=ls(all=TRUE))

  # Folder Adress
    path <- "/home/micdebac/Quantile Regression/New Loss Function/Simple Linear/Summary Project"
    setwd(path)
 
    # --------------------------------------------- #
    # ----------- 1. LOAD LITERATURE -------------- #
    # --------------------------------------------- #

        library("quantreg")
        library("survival")
    
      # Load Functions Literature 
        source("Functions Literature CQR - Linear Regression.R")

    # --------------------------------------------- #
    # ----------- 2. MM algorithm NEW  ------------ #
    # --------------------------------------------- #

      # MM algorithm for QR with complete data
        MMLQR <- function(T.obs,X,tau,beta,toler=1e-8,maxit=5000){
  iteration <- 0
  n.obs     <- length(T.obs)
  df.mat    <- cbind(rep(1,n.obs),X)
  
  # Calculation of epsilon
  tn        <- toler/n.obs
  e0        <- -tn/log(tn)
  eps       <- (e0-tn)/(1+log(e0))
  
  # Initialization of condition for break
  cond <- T
  while(cond)
  { 
    beta.prev <- beta
    
    r.vec     <- T.obs-df.mat%*%as.matrix(beta)
    A.entries <- c(1/(eps+abs(r.vec)))/2
    A.mat     <- diag(A.entries)
    
    B.vec     <- as.matrix(rep(tau-.5,n.obs))
    
    beta      <- solve(t(df.mat)%*%A.mat%*%df.mat)%*%(t(df.mat)%*%(A.mat%*%as.matrix(T.obs)+B.vec))
    cond      <- max(abs(beta-beta.prev)) > toler
    iteration <- iteration + 1
    if(iteration > maxit){warning("WARNING: Algorithm did not converge"); break} 
  }
  return(list("beta"=c(beta),"IterN"=iteration))
}

      # MM algorithm for CQR (univariate)
        MMLQR.cens <- function(Y,delta,X,tau,beta,h,toler=1e-10,maxit=5000){
         
         iteration <- 0
         n.obs     <- length(Y)
         df.mat    <- cbind(1,X)
         
         # Calculation of epsilon
         tn        <- toler/n.obs
         e0        <- -tn/log(tn)
         eps       <- (e0-tn)/(1+log(e0))
         
         # Conditional Kaplan-Meier for G_C(.|X): same code as Wang and Wang
         Bnk.func    <- function(x0, x, h){
           xx             <- (x-x0)/h  
           xx[abs(xx)>=1] <- 1
           w              <- 15*(1-xx^2)^2/16 
           w              <- w/ifelse(sum(w)==0,1,sum(w))
           return(w)
         } 
         tauhat.func <- function(y0, x0, z, x, delta,h){
           
           n<-length(z)
           ###kernel weights#########################
           Bn = Bnk.func(x0, x, h)
           if (y0<max(z))
           {
             # sort the data z, and the delta, Bn correspondingly to the order of sorted y
             z2 = sort(z)
             Order = order(z) # so z[Order] = z2
             Bn2 = Bn[Order]
             delta2 = delta[Order]
             eta = which(delta2==1 & z2<=y0) # the index of those observations satisfying delta2==1 & z2<=y0
             Bn3 = Bn2[n:1]                  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
             tmp = 1- Bn2 /cumsum(Bn3)[n:1]  
             out = 1-prod(tmp[eta], na.rm=T) # na.rm=T, as some of those tmp=NA as the denom =0
           } 
           else out<-1 
           return(out)
         }
         Beran.fun   <- function(a,h){ # rewritten to use apply on rows of matrix of (beta X,X) in MM afterwards
           tauhat.func(a[1],a[2], Y, X, 1-delta,h)
         }
         
         # Initialization of condition for break
         cond <- T
         while(cond)
         { 
           beta.prev <- beta
           
           r.vec     <- Y-df.mat%*%as.matrix(beta)
           A.entries <- c(1/(eps+abs(r.vec)))/2
           A.mat     <- diag(A.entries)
           
           B.mat     <- as.matrix(rep(tau-.5,n.obs))
           G.mat     <- as.matrix((1-tau)*apply(cbind(df.mat%*%as.matrix(beta.prev),X),1,FUN="Beran.fun",h=h))
           
           beta      <- solve(t(df.mat)%*%A.mat%*%df.mat)%*%(t(df.mat)%*%(A.mat%*%as.matrix(Y)+B.mat+G.mat))
           cond      <- max(abs(beta-beta.prev)) > toler
           iteration <- iteration + 1
           
           if(iteration > maxit){warning("WARNING: Algorithm did not converge"); break} 
         }
         return(list("beta"=c(beta),"IterN"=iteration))
       }
       
      # CV for MM algo with Beran (univariate)
        MMLQR.cv <- function(Y,X,delta, nfold=5, h, tau,initial){ 
          # cross validation for selecting the bandwidth
          n <- length(Y)
          m <- n/nfold
          
          pp  <- 1:n
          pp1 <- sample(pp)   ##a random permutation of 1:n
          
          fpred.5cv = NULL
          for (i in 1:nfold)
          {
            #validation data
            ind    <- pp1[seq((i-1)*m+1,i*m,1)]
            xv     <- X[ind]
            yv     <- Y[ind]
            deltav <- delta[ind]        
            #training data
            xt     <- X[-ind]
            yt     <- Y[-ind]
            deltat <- delta[-ind]
            
            beta <- suppressWarnings(MMLQR.cens(yt,deltat,xt,tau,initial,h=h,toler=1e-6)$beta) # bigger tolerance than for real estimation to save time
            pred <- beta[1]+xv*beta[-1]
            
            ind2      <- which(deltav==1)
            tmp.error <- sum((yv[ind2]-pred[ind2])*(tau-as.numeric(yv[ind2]<=pred[ind2])))  #absolute prediction error based on the ith validation sample
            fpred.5cv <- c(fpred.5cv, tmp.error)
          }
          result<- mean(fpred.5cv)
          result
        } 
        
      # MM algorithm for CQR (multivariate)
        MMLQR.cens.multi <- function(Y,delta,X,tau,beta,h,toler=1e-10,maxit=5000){
          require(survival)
          
          iteration <- 0
          n.obs     <- length(Y)
          df.mat    <- cbind(rep(1,n.obs),X)
          
          # Calculation of epsilon
          tn        <- toler/n.obs
          e0        <- -tn/log(tn)
          eps       <- (e0-tn)/(1+log(e0))
          
          
          # Conditional Kaplan-Meier for G_C(.|X) (same code as Wang et al.)
          tauhat.func  <- function(y0, x0, y, x, delta, h, kernel.type="4th"){
            # estimate the cond. censoring probability P(T<y0|x0)
            # x0: k-dimensional covariate vector
            # y: observed survival time = T^C, n-vector
            # x: k-dimensional covariate, n*k matrix
            
            Univ.Kernel.func <- function(x0, x, h, kernel.type){
              # the kernel weight function Bnk(x0, x), where x0 is a scalar, and x is a vector
              # returns a vector
              # h is the bandwidth
              xx<-(x-x0)/h 
              if(kernel.type=="normal")
              {
                w=dnorm(xx)
              }
              if(kernel.type=="4th")
              { 
                xx[abs(xx)>=1]<-1
                w<-105*(1-5*xx^2+7*xx^4-3*xx^6)/64  #biquadratic kernel 
              }
              w<-w/ifelse(sum(w)==0,1,sum(w))
              return(w)
            }
            Kernel.func      <- function(U, U0, h, kernel.type){
              # U: n*k matrix
              # U0: k-vector
              # return: K((U-U0)/h)
              n = nrow(U)
              if(kernel.type=="4th")
              {
                tt = rbind(U, U0)
                tmp = apply(tt, 2, function(x) 
                {
                  Univ.Kernel.func(x[1:n],x[n+1],h,kernel.type)
                })
                tmp = apply(tmp, 1, prod)
                tmp = tmp/sum(tmp)
              }
              if(kernel.type=="normal") # use multivariate normal density as the kernel function
              {
                k = length(U0)
                tt = t(t(U)-U0)/h
                tmp = dmvnorm(tt, mean=rep(0,k), sigma=diag(k), log=FALSE)
                tmp = tmp/sum(tmp)
              }
              return(tmp)    
            }
            
            # this is the same for model1 with interaction, and model2: without interaction
            n<-length(y)
            w<-rep(0,n)  
            
            # modified by Huixia on 04/15/2008
            p = qr(x)$rank
            
            if(p>1)  Bn = Kernel.func(x, x0, h, kernel.type)
            if(p==1) Bn = Univ.Kernel.func(x, x0, h, kernel.type)
            
            if (y0<max(y))
            {
              # sort the data y, and the delta, Bn correspondingly to the order of sorted y
              y2     = sort(y)
              Order  = order(y) # so y[Order] = z2
              Bn2    = Bn[Order]
              delta2 = delta[Order]
              eta    = which(delta2==1 & y2<=y0) # the index of those observations satisfying delta2==1 & y2<=y0
              Bn3    = Bn2[n:1]  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
              tmp    = 1- Bn2 /cumsum(Bn3)[n:1]  
              out    = 1-prod(tmp[eta], na.rm=T) # na.rm=T, as some of those tmp=NA as the denom =0
            } 
            else out=1
            out
          }
          Beran.fun    <- function(a,h){ # rewritten to use apply on rows of matrix of (beta X,X) in MM afterwards
            tauhat.func(a[1],a[2:(dim(X)[2]+1)], Y, X, 1-delta,h)
          }
          
          # Initialization of condition for break
          cond <- T
          while(cond)
          { 
            beta.prev <- beta
            
            r.vec     <- Y-df.mat%*%as.matrix(beta)
            A.entries <- c(1/(eps+abs(r.vec)))/2
            A.mat     <- diag(A.entries)
            
            B.mat     <- as.matrix(rep(tau-.5,n.obs))
            G.mat     <- as.matrix((1-tau)*apply(cbind(df.mat%*%as.matrix(beta.prev),X),1,FUN="Beran.fun",h=h))
            
            beta      <- solve(t(df.mat)%*%A.mat%*%df.mat)%*%(t(df.mat)%*%(A.mat%*%as.matrix(Y)+B.mat+G.mat))
            cond      <- max(abs(beta-beta.prev)) > toler
            iteration <- iteration + 1
           
            if(iteration > maxit){warning("WARNING: Algorithm did not converge"); break} 
          }
          return(list("beta"=c(beta),"IterN"=iteration))
        }
      
      # CV for MM algo with Beran (multivariate)
        MM.multi.cv <- function(Y,X,delta, nfold=5, h, tau,initial){ 
          # cross validation for selecting the bandwidth
          n <- length(Y)
          m <- n/nfold
          
          pp  <- 1:n
          pp1 <- sample(pp)   ##a random permutation of 1:n
          
          fpred.5cv = NULL
          for (i in 1:nfold)
          {
            #validation data
            ind    <- pp1[seq((i-1)*m+1,i*m,1)]
            xv     <- X[ind,]
            yv     <- Y[ind]
            deltav <- delta[ind]        
            #training data
            xt     <- X[-ind,]
            yt     <- Y[-ind]
            deltat <- delta[-ind]
            
            beta <- suppressWarnings(MMLQR.cens.Gcond(yt,deltat,xt,tau,initial,h=h,toler=1e-6)$beta) # bigger tolerance than for real estimation to save time
            pred <- beta[1]+xv*beta[-1]
            
            ind2      <- which(deltav==1)
            tmp.error <- sum((yv[ind2]-pred[ind2])*(tau-as.numeric(yv[ind2]<=pred[ind2])))  #absolute prediction error based on the ith validation sample
            fpred.5cv <- c(fpred.5cv, tmp.error)
          }
          result<- mean(fpred.5cv)
          result
        } 
    
        
    # --------------------------------------------- #
    # ---------- 3. ONE SAMPLE EXAMPLE  ----------- #
    # --------------------------------------------- #   
       
        # Generate Data for DGP1 in our article
          DataGen.DGP1 <- function(n,tau=.5,prop.cens=.4){
            b0  <- 3
            b1  <- 5 
            X   <- runif(n)
            eta <- rnorm(n)
            eps <- eta - qnorm(tau)
            
            T <- b0+b1*X+eps
            
            censoring <- c(.15,.4)
            max.C     <- c(36,14)[which(censoring==prop.cens)]
            C         <- runif(n,min=0,max=max.C)
            Y         <- pmin(T,C)
            delta     <- T <= C
            
            result     <- matrix(nrow= n, ncol=5)
            result[,1] <- T
            result[,2] <- Y
            result[,3] <- delta
            result[,4] <- X
            result[,5] <- C
            
            return(result)
          }
        
        # Setting
          tau       <- .5
          prop.cens <- .4
          n.obs     <- 200
          beta      <- c(3,5)
          
        # Simulate data
          set.seed(2015)
        
          data.sim  <- DataGen.DGP1(n.obs,tau,prop.cens)
          T.sim     <- data.sim[,1]
          Y.sim     <- data.sim[,2]
          delta.sim <- data.sim[,3]
          X.sim     <- data.sim[,4]
          
        # Estimation
        
          # Fully observed data
            (m.est.full <- rq(T.sim~X.sim,tau)$coef)
          
          # Censored data 
          
            # 1/ Wang and Wang
              (m.est.WW <- WW.cens(Y.sim, X.sim, delta.sim, tau=tau, h=.05)$coeff)          
            
            # 2/ Bang and Tsiatis
              (m.est.BT <- BT.cens(Y.sim, delta.sim, X.sim, tau=tau))
            
            # 3/ New estimator (could be fastened up with C version of Beran estimator)
              (m.est.NEW <- MMLQR.cens(Y.sim,delta.sim,X.sim,tau,m.est.BT+1e-2,h=.05))
          
    # --------------------------------------------- #
    # ----------- 4. SIMULATION STUDY  ------------ #
    # --------------------------------------------- #             

      # Setting
        IterN     <- 500
        n.obs     <- 100
        prop.cens <- .15
        
      # Storage
        T.sim.mat     <- matrix(nrow= IterN, ncol=n.obs)
        Y.sim.mat     <- matrix(nrow= IterN, ncol=n.obs)
        delta.sim.mat <- matrix(nrow= IterN, ncol=n.obs)
        X.sim.mat     <- matrix(nrow= IterN, ncol=n.obs)
        
      # 0/ Data Generation - Store every data set for i in 1:IterN
        set.seed(12345678)
        for(i in 1:IterN){    
          
          data.sim <- DataGen.DGP1(n.obs,tau,prop.cens)
          
          T.sim.mat[i,]     <- data.sim[,1]
          Y.sim.mat[i,]     <- data.sim[,2]
          X.sim.mat[i,]     <- data.sim[,4]
          delta.sim.mat[i,] <- data.sim[,3]
          
        }
        plot(density(apply(delta.sim.mat,1,mean)),main="Average censoring proportion")
        
      # Settings parralelization
        library(doParallel)
        library(foreach)
        
      # Timer
        p2 <- Sys.time()
      
      # 1/ Full Observations
        do.full <- TRUE
        if(do.full){
          cl <- makeCluster(6, outfile="")
          registerDoParallel(cl)
          
          M.est.full <- foreach(i=1:IterN, .combine=rbind) %dopar% {    
        
            T.sim     <- T.sim.mat[i,]
            Y.sim     <- Y.sim.mat[i,]
            delta.sim <- delta.sim.mat[i,]
            X.sim     <- X.sim.mat[i,]
          
            require(quantreg)
            rq(T.sim~X.sim,tau)$coef
        }
      
          Sys.sleep(.15)
          stopCluster(cl) 
        }
        
      # 2/ Wang and Wang estimator
        do.WW <- TRUE
        if(do.WW){
          cl <- makeCluster(6, outfile="")
          registerDoParallel(cl)
        
         M.est.WW <- foreach(i=1:IterN, .combine=rbind) %dopar% {    
          require(quantreg)
        
          T.sim     <- T.sim.mat[i,]
          Y.sim     <- Y.sim.mat[i,]
          delta.sim <- delta.sim.mat[i,]
          X.sim     <- X.sim.mat[i,]
        
          # Bandwidth CV
            h.vect.WW <- seq(0.05, .5, length=15)
            cv.WW     <- NULL
            for (j in 1:length(h.vect.WW)){
              h.temp <- h.vect.WW[j]
              tmp.cv <- WW.cv(Y.sim,X.sim, delta.sim, nfold=5, h.temp, tau)
              cv.WW  <- c(cv.WW, tmp.cv)
            }
            h.WW  <- h.vect.WW[which.min(cv.WW)]
        
            cat("\r  WW estimator",
                "- EXECUTION:", i,"of",IterN,
                "- TIME",round(difftime(Sys.time(), p2, units="mins"),2),"min",collapse="")
            
          WW.cens(Y.sim, X.sim, delta.sim, tau=tau, h=h.WW)$coeff          
      }
  
         Sys.sleep(.15)
         stopCluster(cl) 
        }
        
      # 3/ MM algorithm for NEW procedure
        do.DBEGVK <- TRUE
        if(do.DBEGVK){
          cl <- makeCluster(8, outfile="")
          registerDoParallel(cl)
        
          M.est.MM.G <- foreach(i=1:IterN, .combine=rbind) %dopar% {
        
            T.sim     <- T.sim.mat[i,]
            Y.sim     <- Y.sim.mat[i,]
            delta.sim <- delta.sim.mat[i,]
            X.sim     <- X.sim.mat[i,]

            # Cross-validation for bandwidth selection (takes long time without C version of Beran)            
             # h.vect.MM <- seq(0.05, .5, length=10)
             # cv.MM     <- NULL
             # for (j in 1:length(h.vect.MM)){
             #   h.temp <- h.vect.MM[j]
             #   tmp.cv <- MMLQR.cv(Y.sim,X.sim, delta.sim, nfold=5, h.temp, tau,initial=M.est.WW[i,]+1e-4)
             #   cv.MM  <- c(cv.MM, tmp.cv)
             # }
            
            h.MM  <- .1#h.vect.MM[which.min(cv.MM)]
            
            cat("\r NEW estimator (No CV to fasten up)",
                "- EXECUTION:", i,"of",IterN,
                "- TIME",round(difftime(Sys.time(), p2, units="mins"),2),"min",collapse="")
        
            MMLQR.cens(Y.sim,delta.sim,X.sim,tau,M.est.WW[i,]+1e-1,h=h.MM,toler=1e-8)$beta
        
      }
      
          Sys.sleep(.15)
          stopCluster(cl) 
        }
        
      # 4/ Bang and Tsiatis estimator
        do.Qout <- TRUE
        if(do.Qout){
          cl <- makeCluster(6, outfile="")
          registerDoParallel(cl)
        
          M.est.Qout <- foreach(i=1:IterN, .combine=rbind) %dopar% {    
          
            T.sim     <- T.sim.mat[i,]
            Y.sim     <- Y.sim.mat[i,]
            delta.sim <- delta.sim.mat[i,]
            X.sim     <- X.sim.mat[i,]
            
            require(quantreg)
  
            cat("\r BT estimator",
                "- EXECUTION:", i,"of",IterN,
                "- TIME",round(difftime(Sys.time(), p2, units="mins"),2),"min",collapse="")
            
          rq(Y.sim~X.sim,tau,weights=delta.sim/(1-KM.str.fun(Y.sim,Y.sim,1-delta.sim)))$coef
        }
          
          Sys.sleep(.15)
          stopCluster(cl) 
        }
  
        
        #### Graphs  ####
        #################
  
          a1 <- 3
          a2 <- 5
  
          Beta0.mat <- cbind(Omni=M.est.full[,1],WW=M.est.WW[,1],NEW=M.est.MM.G[,1],BT=M.est.Qout[,1])
          Beta1.mat <- cbind(Omni=M.est.full[,2],WW=M.est.WW[,2],NEW=M.est.MM.G[,2],BT=M.est.Qout[,2])
        
          par(mfrow=c(2,1))
          boxplot(Beta0.mat, main="Beta0");abline(a=a1,b=0)
          boxplot(Beta1.mat, main="Beta1");abline(a=a2,b=0)    
  
        #### Results  ####
        ##################
  
          # BIAS
            (beta0.BIAS <- round(apply(Beta0.mat,2,mean,na.rm=TRUE)-a1,3))
            (beta1.BIAS <- round(apply(Beta1.mat,2,mean,na.rm=TRUE)-a2,3))
          
          # RMSE
            (beta0.RMSE <- round(sqrt(apply((Beta0.mat-a1)^2,2,mean,na.rm=TRUE)),3))
            (beta1.RMSE <- round(sqrt(apply((Beta1.mat-a2)^2,2,mean,na.rm=TRUE)),3))
          
          # MAE
            (beta0.MAE <- round(apply(abs(Beta0.mat-a1),2,median,na.rm=TRUE),3))
            (beta1.MAE <- round(apply(abs(Beta1.mat-a2),2,median,na.rm=TRUE),3))
         
          end.time    <- Sys.time()
          (time.taken <- end.time - p2)      
          
###########################################################################################
###########################################################################################
###########################################################################################