###### MATRIX OF SETTINGS  ##########

matrix = matrix(c(0, 7, 0.5, 200, 500,0, 7, 0.5, 500, 500,0, 7, 0.7, 200, 500,0, 7, 0.7, 500, 500),
                nrow = 4, ncol= 5, byrow=TRUE)
colnames(matrix) = c("min", "max", "tau", "n", "N")
h = c(1,2,3,4,5,6,7,8,9,10)
type = c("Legendre")
m=2
m_tilde=2

datasets = list(array(dim = c(200,6,500)) ,array(dim = c(500, 6, 500)), array(dim = c(200, 6, 500)),
                array(dim = c(500, 6, 500)))


set.seed(1234)


for(k in 1:nrow(matrix)){
  for (i in 1:dim(datasets[[k]])[3]){

    tau  <- matrix[k, 3]
    min = matrix[k,"min"]
    max = matrix[k,"max"]
    n = matrix[k,"n"]
    N = matrix[k,"N"]
    beta <- c(2,1)

    epsilon_raw=rnorm(n, mean=0, sd=1) - qnorm(tau, mean= 0 , sd=1)
    x     <- runif(n, min = -1, max = 1)
    X     <- matrix(0,ncol=3,nrow=n)
    X[,1] <- 1
    X[,2] <- x
    X_s = as.matrix(X[,2])
    #epsilon <- (ifelse(x_f == 1, cos(X[,2]), 0.5*cos(1.5*X[,2])))*epsilon_raw
    epsilon <- (4*cos(1.5*(X[,2])))*epsilon_raw
    #epsilon <- ((3.2+3*cos(8*(X[,2]))))*epsilon_raw

    T     <- X[,1:2]%*%beta+epsilon

    C     <- runif(n,min=min,max=max)

    Y<- pmin(T,C)

    Delta <- as.numeric(T==Y)

    datasets[[k]][,,i] = cbind(X, Y, Delta, T)
    rm(list=c("tau","min","max","n","N","beta","epsilon_raw","x","X","X_s","epsilon","T",
              "C","Y","Delta"))
  }
}

rm(list=c("i","k"))


#### array of seeds
seeds = array(1:(length(h)*dim(matrix)[1]*10000), dim=c(length(h), dim(matrix)[1],
                                                        10000))

###### FUNCTION TO AGGREGATE THE DIFFERENT ARRAYS ######
###### arguments x and y are arrays here
cube = function(...){
  return(abind(..., along=3))
}




