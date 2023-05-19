##########################################################################
################### An Adapted Loss Function for #########################
################### Censored Quantile Regression #########################
#################  R Code - Functions Literature  ########################
##########################################################################

# Date: 8/12/2017
# This file contains the basic functions needed for the estimation of a
# censored quantile regression model using several estimators in the literature
# These include: 
#         - Bang and Tsiatis (2002)
#         - Portnoy (2003)
#         - Peng and Huang (2008)
#         - Wang and Wang (2009)
#         - Wang et al. (2013)
#         - Wey et al. (2014)

library("quantreg")
library("survival")

    # ------------------------------------------- #
    # ----------- 1. BASIC FUNCTIONS ------------ #
    # ------------------------------------------- #
    
      # Indicator function #
        ind            <- function(x,a,b)  as.numeric(x >= a & x <= b)
        ind.str        <- function(x,a,b)  as.numeric(x >= a & x < b)
        
      # Handmade Kaplan-Meier for distribution (not elegantly written)
        fonc.G         <- function(ordre,Y,delta){
          n<-length(Y)
          indice<-order(Y)
          prod<-1
          if(ordre>1){ordre<-ordre-1}
          for(j in 1:ordre){prod<-prod*((n-j)/(n-j+1))^(1-delta[indice[j]])}
          return(1-prod)
}
        poids.KM       <- function(ordre,Y,delta){ # Attention divis? par (n+1)!!
          pds<-delta[ordre]/((length(Y)+1)*(1-fonc.G(ordre,Y,delta)))
          return(pds)
}
        KM.fun         <- function(t,Y,delta){ #Attention K-M de fonction de distribution!!!!
          indice<-order(Y)
          Y1<-sort(Y)
          delta1<-c()
          for(j in 1:length(Y)){delta1[j]<-delta[indice[j]]}
          m<-length(Y);sum<-0
          for(i in 1:m){sum<-sum+poids.KM(i,Y1,delta1)*ind(Y1[i],-Inf,t)}
          return(sum)
}
        KM.str.fun     <- function(t,Y,delta){ #Attention K-M de fonction de distribution!!!! INEGALITE STRICTE DANS INDICATRICE
          indice<-order(Y)
          Y1<-sort(Y)
          delta1<-c()
          for(j in 1:length(Y)){delta1[j]<-delta[indice[j]]}
          m<-length(Y);sum<-0
          for(i in 1:m){sum<-sum+poids.KM(i,Y1,delta1)*ind.str(Y1[i],-Inf,t)}
          return(sum)
}
  
    # --------------------------------------------- #
    # -------------- 2. LITERATURE ---------------- #
    # --------------------------------------------- #

      # Bang and Tsiatis (2002)
        BT.cens <- function(y,delta,tau,X){
          rq(y~X,tau,weights=delta/(1-KM.str.fun(y,y,1-delta)))$coef
        }
        
      # Portnoy (2003)
        Port.cens <- function(y,delta,tau,X){
  portnoy.cqr  <- crq(Surv(y,delta) ~ X, tau=tau, method = "Portnoy")
  coef.port    <- summary(portnoy.cqr, c(tau,tau),R=10)[[1]]$coef
  beta.port    <- as.numeric(coef.port[,1])
  return(beta.port)
}
      
      # Peng and Huang (2008)
        PH.cens <- function(y,delta,tau,X){
  portnoy.cqr  <- crq(Surv(y,delta) ~ X, tau=tau, method = "PengHuang")
  coef.port    <- summary(portnoy.cqr, c(tau,tau),R=10)[[1]]$coef
  beta.port    <- as.numeric(coef.port[,1])
  return(beta.port)
}
      
      # Wang and Wang (2009) (Code from webpage Huixia Wang)
        Bnk.func    <- function(x0, x, h){
  # the kernel weight function Bnk(x0, x), where x0 is a scalar, and x is a vector
  # returns a vector
  # h is the bandwidth
  xx<-(x-x0)/h  
  xx[abs(xx)>=1]<-1
  w<-15*(1-xx^2)^2/16  #biquadratic kernel 
  w<-w/sum(w)
  return(w)
} 
        tauhat.func <- function(y0, x0, z, x, delta,h){
  # tau0(y0, x0) = F(T<y0|x0); so y0 is the C_i, and x0 is the xi in the paper
  # z is observed vector of response variable
  # x is the observed covariate
  # delta is the censoring indicator function
  # h is the bandwidth
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
    Bn3 = Bn2[n:1]  # change the order of Bn2, make the first obs of Bn2 to be the last of Bn3
    tmp = 1- Bn2 /cumsum(Bn3)[n:1]  
    out = 1-prod(tmp[eta], na.rm=T) # na.rm=T, as some of those tmp=NA as the denom =0
  } 
  else out<-1 
  return(out)
}
        WW.cens     <- function(y, x, delta, tau, h){
  # x is one dimensional, not including the intercept yet
  # y is the observed survival time = min(T, C)
  # delta is the censoring indicator function with 1 standing for uncensored, and 0 censored
  # tau is the quantile level of interest
  # h is the handwidth used in calculating tauhat for each censored observation
  n = length(y)
  ind = which(delta==0)
  w = rep(1, n)  # the weight
  if(length(ind)>=1)
  {
    for(i in 1:length(ind))
    {
      x0 = x[ind[i]]
      y0 = y[ind[i]]
      tau.star = tauhat.func(y0,x0, y, x, delta,h=h)
      if (tau>tau.star) w[ind[i]] = (tau-tau.star)/(1-tau.star)
    }
    # pseudo observations
    ind2 = which(w!=1)
    y.pse = rep(max(y)+100, length(ind2))
    x.pse = x[ind2]
    yy = c(y, y.pse)
    xx = c(x, x.pse)
    ww = c(w, 1-w[ind2])
  }
  else
  {
    yy=y; xx=x; ww=w
  }
  rq1 = rq(yy~xx, weights=ww, tau=tau)
  result<- rq1$coeff
  return(list(coeff=result,weights=w))    
}
        
        # Cross-validation code
          WW.cv <- function(y,x,delta, nfold=5, h, tau){ 
            # cross validation for selecting the bandwidth
            n <- length(y)
            m <- n/nfold
            
            pp  <- 1:n
            pp1 <- sample(pp)   ##a random permutation of 1:n
            
            fpred.5cv = NULL
            for (i in 1:nfold)
            {
              #validation data
              ind    <- pp1[seq((i-1)*m+1,i*m,1)]
              xv     <- x[ind]
              yv     <- y[ind]
              deltav <- delta[ind]        
              #training data
              xt     <- x[-ind]
              yt     <- y[-ind]
              deltat <- delta[-ind]
              
              beta <- WW.cens(yt, xt, deltat, tau=tau, h)$coeff
              pred <- beta[1]+xv*beta[-1]
              
              ind2      <- which(deltav==1)
              tmp.error <- sum((yv[ind2]-pred[ind2])*(tau-as.numeric(yv[ind2]<=pred[ind2])))  # prediction error based on the ith validation sample
              fpred.5cv <- c(fpred.5cv, tmp.error)
            }
            result<- mean(fpred.5cv)
            result
          }
        
        # MM version of Wang and Wang's estimator to check if MM code OK
          MMLQR.WW    <- function(y,X,delta,tau,beta,h,toler=1e-8,maxit=5000){
          # MM version of Wang and Wang's estimator to check if MM code OK
          iteration <- 0
          n.obs     <- length(y)
          df.mat    <- cbind(rep(1,n.obs),X)
          
          # Calculation of epsilon
          tn        <- toler/n.obs
          e0        <- -tn/log(tn)
          eps       <- (e0-tn)/(1+log(e0))
          
          # Calculation of weights \hat{w}
          ind <- which(delta==0)
          w   <- rep(1, n.obs)  # the weight
          
          for(i in 1:length(ind))
          {
            x0 = X[ind[i]]
            y0 = y[ind[i]]
            tau.star = tauhat.func(y0,x0, y, X, delta,h=h)
            if (tau>tau.star) w[ind[i]] = (tau-tau.star)/(1-tau.star)
          }
          
          # pseudo observations
          ind2       <- which(w!=1)
          y.pse      <- rep(max(y)+100, length(ind2))
          x.pse      <- X[ind2]
          df.mat.pse <- cbind(rep(1,length(ind2)),x.pse)
          
          yy    <- c(y, y.pse)
          xx    <- c(X, x.pse)
          ww    <- c(w, 1-w[ind2])
          
          # Initialization of condition for break
          cond <- T
          while(cond)
          { 
            beta.prev <- beta
            
            r.vec     <- y-df.mat%*%as.matrix(beta)
            r.vec.pse <- y.pse-df.mat.pse%*%as.matrix(beta)
            
            A.entries     <- c(w/(eps+abs(r.vec)))/2
            A.mat         <- diag(A.entries)
            A.entries.pse <- c((1-w[ind2])/(eps+abs(r.vec.pse)))/2
            A.mat.pse     <- diag(A.entries.pse)
            
            B.vec     <- as.matrix(rep(tau-.5,n.obs))
            
            beta      <- solve(t(df.mat)%*%A.mat%*%df.mat+t(df.mat.pse)%*%A.mat.pse%*%df.mat.pse)%*%(t(df.mat)%*%(A.mat%*%as.matrix(y)+B.vec)+t(df.mat.pse)%*%A.mat.pse%*%as.matrix(y.pse))
            cond      <- max(abs(beta-beta.prev)) > toler
            iteration <- iteration + 1
            if(iteration > maxit){warning("WARNING: Algorithm did not converge"); break} 
          }
          return(list("beta"=c(beta),"IterN"=iteration))
        } 

      # Wang et al. (2013) (i.e. multivariate version Wang and Wang)
        tauhat.multi.func <- function(y0, x0, y, x, delta, h, kernel.type="4th"){
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
            w<-w/sum(w)
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
        Pseudo.func       <- function(y, x, delta, tau, h, kernel.type="4th"){
          # estimate the local K-M weights for reweighting (based on cond. distribution of y given x)
          # return the pseudo weights, Y, and X
          
          n   <- length(y)
          ind <- which(delta==0)
          w   <- rep(1, n)  # the weight
          
          for(i in 1:length(ind))
          {
            x0       <- x[ind[i],]
            y0       <- y[ind[i]]
            tau.star <- tauhat.multi.func(y0, x0 , y, x, delta, h=h, kernel.type)
            if (tau>tau.star) w[ind[i]] = (tau-tau.star)/(1-tau.star)
            else w[ind[i]]=0
            
          }
          # pseudo observations
          ind2  <- which(w!=1)
          y.pse <- rep(max(y)+100, length(ind2))
          yy    <- c(y, y.pse)   
          xx    <- rbind(x, x[ind2,])
          ww    <- c(w, 1-w[ind2])
          
          out = cbind(yy, ww, xx)   # xx does not include the intercept
          return(out)
        }
        WW.multi.cens     <- function(y, x, delta, tau, h, kernel.type="4th"){
          # Locally weighted censored quantile regression method (without any penalization)
          # DIR: indices
          # local weights are obtained by the cond. distribution of T given DIR
          # x: does not include the first column of 1's corresponding to the intercept
          tmp = Pseudo.func(y, x, delta, tau, h, kernel.type)
          yy = tmp[,1]
          ww = tmp[,2]
          xx = tmp[,-(1:2)]
          rq1 = rq(yy~xx, weights=ww, tau=tau)
          result<-rq1$coeff
          result    
        }
        
        # Cross-validation code
        WW.multi.cv       <- function(y,x,delta, nfold=5, h, tau, kernel.type="4th"){ 
          # cross validation for selecting the bandwidth
          n <- length(y)
          m <- n/nfold
          
          pp  <- 1:n
          pp1 <- sample(pp)   ##a random permutation of 1:n
          
          fpred.5cv = NULL
          for (i in 1:nfold)
          {
            #validation data
            ind    <- pp1[seq((i-1)*m+1,i*m,1)]
            xv     <- x[ind,]
            yv     <- y[ind]
            deltav <- delta[ind]        
            #training data
            xt     <- x[-ind,]
            yt     <- y[-ind]
            deltat <- delta[-ind]
            
            beta <- WW.multi.cens(yt, xt, deltat, tau=tau, h, kernel.type)
            pred <- beta[1]+xv%*%beta[-1]
            
            ind2      <- which(deltav==1)
            tmp.error <- sum(abs(yv[ind2]-pred[ind2]))  #absolute prediction error based on the ith validation sample
            fpred.5cv <- c(fpred.5cv, tmp.error)
          }
          result<- mean(fpred.5cv)
          result
        }

      # Wey et al. (2014)
        Wey.RPcrq <- function(y, x, delta, tau, minatrisk = 15, pruneZ = 0, track = FALSE, bagN = 20){
           require("RPcrq")
          first <- order(y)
          Nnodes <- NULL
          X <- as.matrix(x)[first, , drop = FALSE]
          y <- y[first]
          delta <- delta[first]
          zdat <- 1:length(y)
          subgrp <- rep(TRUE, length(y))
          nObs <- length(y)
          p <- dim(X)[2]
          ux <- jx <- NULL
          for (j in 1:p) {
            u <- sort(unique(X[subgrp, j]))
            if (length(u) > 1) {
              ux <- c(ux, u[-1] - diff(u)/2)
              jx <- c(jx, rep(j, length(u) - 1))
            }
          }
          xsplits.all <- matrix(F, nObs, length(ux))
          for (j in 1:p) {
            jp <- jx == j
            if (any(jp)) {
              xsplits.all[, jp] <- matrix(X[, j] >= rep(ux[jp], 
                                                        each = nObs), nObs)
            }
          }
          fEst <- rep(0, nObs)
          for (i in 1:bagN) {
            if (i == 1) {
              ind <- zdat
            }
            if (i > 1) {
              ind <- sample(zdat, length(y), replace = TRUE)
            }
            y1 <- y[ind]
            X1 <- X[ind, , drop = FALSE]
            delta1 <- delta[ind]
            ord <- order(y1)
            X.tmp <- X1[ord, , drop = FALSE]
            y.tmp <- y1[ord]
            event <- delta1[ord]
            orgInd <- ind
            ordInd <- ind[ord]
            xsplits <- xsplits.all[ordInd, ]
            xsplits <- 1 * xsplits
            Ctree <- .Call("nodefxn", y.tmp, event, 1 * subgrp, xsplits, 
                           c(0, 0, 1, 1), c(0, 1, 0, 1), c(minatrisk, tau))
            nodes <- length(unlist(Ctree, recursive = F))/length(Ctree[[1]])
            Crslt <- NULL
            for (i in 1:nodes) {
              Crslt <- rbind(Crslt, as.numeric(Ctree[[i]][c(1:5, 
                                                            9:13)]))
            }
            ind <- order(Crslt[, 1])
            for (nde in 1:nodes) {
              Ctree[[nde]][13] <- F
              Ctree[[nde]][12] <- Ctree[[nde]][12] == 1
            }
            done <- F
            while (!done) {
              nbrPrune <- 0
              for (nde in ind) {
                if (!Ctree[[nde]][[12]]) {
                  if (abs(Ctree[[nde]][[5]]) <= pruneZ & Ctree[[ind[Ctree[[nde]][[9]]]]][[12]] & 
                      Ctree[[ind[Ctree[[nde]][[10]]]]][[12]]) {
                    nbrPrune <- nbrPrune + 1
                    Ctree[[nde]][12] <- T
                    Ctree[[ind[Ctree[[nde]][[9]]]]][13] <- T
                    Ctree[[ind[Ctree[[nde]][[10]]]]][13] <- T
                  }
                }
              }
              if (nbrPrune == 0) 
                done <- T
            }
            back.track <- function(crslt, j) {
              cur <- j
              stp <- TRUE
              boo.all <- boo <- rep(1, length(y))
              while (stp) {
                parent.n <- Crslt[cur, 8]
                if (Crslt[parent.n, 7] == cur) {
                  blah <- 1 - xsplits[, Crslt[parent.n, 4]]
                  blah.all <- 1 - xsplits.all[, Crslt[parent.n, 
                                                      4]]
                }
                if (Crslt[parent.n, 6] == cur) {
                  blah <- xsplits[, Crslt[parent.n, 4]]
                  blah.all <- xsplits.all[, Crslt[parent.n, 4]]
                }
                boo <- boo * blah
                boo.all <- boo.all * blah.all
                if (parent.n == 1) {
                  stp <- FALSE
                }
                if (parent.n > 1) {
                  cur <- parent.n
                }
              }
              return(cbind(boo, boo.all))
            }
            sEst <- rep(1, length(y))
            Crslt <- CStree <- CsimMemtree <- allMemMat <- allMemMat.tst <- NULL
            for (i in 1:nodes) {
              if (nodes == 1) {
                sEst <- sapply(y, function(x) {
                  if (x < min(y.tmp)) {
                    return(1)
                  }
                  if (x >= min(y.tmp)) {
                    return(unlist(min(unlist(Ctree[[1]][7])[y.tmp <= 
                                                              x], na.rm = TRUE)))
                  }
                })
              }
              if (nodes > 1) {
                Crslt <- rbind(Crslt, as.numeric(Ctree[[ind[i]]][c(1:5, 
                                                                   9:13)]))
                CStree <- unlist(Ctree[[ind[i]]][7])
                if (Crslt[i, 9] == 1) {
                  node.rslt <- back.track(Crslt, i)
                  sEst[node.rslt[, 2] == 1] <- sapply(y[node.rslt[, 
                                                                  2] == 1], function(x) {
                                                                    if (x < min(y.tmp)) {
                                                                      return(1)
                                                                    }
                                                                    if (x >= min(y.tmp)) {
                                                                      return(unlist(min(CStree[y.tmp <= x], na.rm = TRUE)))
                                                                    }
                                                                  })
                }
              }
            }
            fEst <- fEst + (1 - sEst)/bagN
            Nnodes <- c(Nnodes, sum(Crslt[, 9], na.rm = TRUE))
          }
          n = length(y)
          ind = which(delta == 0)
          wR = rep(1, n)
          if (length(ind) >= 1) {
            boo <- delta == 0 & fEst < tau
            wR[boo] <- (tau - fEst[boo])/(1 - fEst[boo])
            ind2 = which(wR != 1)
            y.pse = rep(10000 * (max(y) + 1), length(ind2))
            x.pse = X[ind2, , drop = FALSE]
            yy = c(y, y.pse)
            xx = rbind(X, x.pse)
            wwR = c(wR, 1 - wR[ind2])
          }
          else {
            yy = y
            xx = X
            wwR = wR
          }
          rq1 = rq(yy ~ xx, weights = wwR, tau = tau)
          result1 <- rq1$coeff
          if (track) {
            return(list(rslt = result1, weights = matrix(wR, nrow = 1), 
                        nodeN = Nnodes))
          }
          else {
            return(result1)
          }
}
        
        
        