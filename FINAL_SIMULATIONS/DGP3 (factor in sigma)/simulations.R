## Load Libraries ----
library("abind")
library(foreach)
library(doParallel)
library(dplyr)
library(Laguerre)
library(pracma)
library(tidyr)
library(ggplot2)
library(kableExtra)

## SET WD  ----
setwd(dir = "C:/Users/dzr/Documents/KU Leuven PC/Documents/Master Thesis/FINAL_SIMULATIONS/DGP3 (factor in sigma)")

## Find the number of cores in your system ----
clno <- detectCores()
cl   <- makeCluster(clno,outfile="test2")
registerDoParallel(cl)

## LOAD LITERATURE AND DATASETS ----
source(file = "Loading Literature.R")
source(file = "DGP3.R")

## Identity link ----

id = function(x){return(x)}
idd = function(x){return(x^0)}
link = list(id, idd)

## absolute value link ----

id = function(x){return(abs(x))}
idd = function(x){return(ifelse(x>=0, 1, -1))}
link2 = list(id, idd)

# Initialize list ----
h_list = vector(mode = "list", length = length(h))

## Iterative loop ----

for (H in h){
  out7 =
    foreach(k = 1:(dim(matrix)[1]), .combine = 'cube', .packages = 'abind', .multicombine = TRUE)%:%
    foreach(i=1:(dim(datasets[[k]])[3]-498),.packages=c('nloptr','SphericalCubature', 'EQL','orthopolynom',
                                                  'quantreg', 'survival', 'Laguerre'),
            .combine=rbind) %dopar% {

              ### setting random seed
              set.seed(seeds[which(h==H,arr.ind = TRUE), k, i])


              cat("Step ",i," of ",matrix[k,"N"]," from simulation ",k, " ", "h = ", H, "\n")

              X = datasets[[k]][,1:2,i]
              Y = datasets[[k]][,4,i]
              Delta = datasets[[k]][,5,i]
              T = datasets[[k]][,6,i]
              X_s = as.matrix(datasets[[k]][,2:3, i])
              tau  <- matrix[k, 3]

              # Bandwidth CV
              h.vect.WW <- seq(0.05, .5, length=15)
              cv.WW     <- NULL
              for (j in 1:length(h.vect.WW)){
                h.temp <- h.vect.WW[j]
                tmp.cv <- WW.cv(Y,X[,2], Delta, nfold=5, h.temp, tau)
                cv.WW  <- c(cv.WW, tmp.cv)
              }
              h.WW  <- h.vect.WW[which.min(cv.WW)]

              omni = rq(T~X[,2], tau = tau)
              crq = crq(Surv(Y,Delta, type='right')~X[,2], tau=tau, method = "Portnoy")
              
              estexp = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="exp"))
              while(estexp$objective %in% c(-Inf, Inf)){
                estexp = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="exp"))
              }
              
              estquad = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="quad"))
              while(estquad$objective %in% c(-Inf, Inf)){
                estquad = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link="quad"))
              }
              
              est_id = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link))
              while(est_id$objective %in% c(-Inf, Inf)){
                est_id = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link))
              }
              
              est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link2))
              while(est_abs$objective %in% c(-Inf, Inf)){
                est_abs = try(laguerre_estimator_het(m,m_tilde,H,X,X_s,type, Y, Delta, tau,trials=32, verbose=0,link=link2))
              }
              beta_h = laguerre_estimator_het(m,m_tilde,0, X=X,X_s,type, Y=Y, Delta=Delta, tau=tau,trials=32, verbose = 0)$beta
              adapted = MMLQR.cens(Y,Delta,X[,2],tau,h=0.5, beta=c(4,5))
              W.W = WW.cens(Y, X[,2], Delta, tau, 0.1)
              W.W_cv = WW.cens(Y, X[,2], Delta, tau, h.WW)
              PHuang = PH.cens(Y, Delta, tau, X[,2])

              ### Collecting the results
              c(omni$coefficients,
                PHuang,
                W.W$coeff,
                W.W_cv$coeff,
                crq$sol[2:3,which.min(abs(tau - crq$sol["tau",]))],
                adapted$beta,
                beta_h$beta,
                estexp$beta,estquad$beta, est_id$beta, est_abs$beta,
                estexp$H, estquad$H, est_id$H, est_abs$H,
                estexp$theta, estquad$theta, est_id$theta, est_abs$theta,
                estexp$theta_tilde, estquad$theta_tilde, est_id$theta_tilde, est_abs$theta_tilde)
            }
  h_list[[which(h==H,arr.ind = TRUE)]] = out7
}

## Save results ----
save(list=c("h_list"), file="results.RData")


## Compute statistics ----

true_beta = c(2,1) # CHANGE IF YOU CHANGE THE TRUE MODEL?

#BIAS COMPUTATION

bias = vector(mode = "list", length = length(h))

for(i in 1:10){
  bias[[i]] = array(dim = c(11, length(true_beta), nrow(matrix)))
  for(k in 1:nrow(matrix)){
colnames(h_list[[i]]) = NULL
rownames(h_list[[i]]) = NULL
bias[[i]][,,k] = (h_list[[i]][,1:(11*length(true_beta)),k] %>% colMeans() %>% matrix(ncol = length(true_beta), byrow = TRUE)) -
  matrix(true_beta, ncol = length(true_beta), nrow = 11, byrow = TRUE)
  }
}

# MSE COMPUTATION

MSE = vector(mode = "list", length = length(h))

for(i in 1:10){
  MSE[[i]] = array(dim = c(11, length(true_beta), nrow(matrix)))
  for(k in 1:nrow(matrix)){
    colnames(h_list[[i]]) = NULL
    rownames(h_list[[i]]) = NULL
    MSE[[i]][,,k] = ((h_list[[i]][,1:(11*length(true_beta)),k] - matrix(true_beta, nrow = dim(h_list[[i]])[1], ncol = (11*length(true_beta)), byrow = TRUE))^2 %>%
                        colMeans() %>% matrix(ncol = length(true_beta), byrow = TRUE))
  }
}


# SIGMAS COMPUTATION ----
SIGMA_BIG = vector(mode = "list", length = nrow(matrix))
SIGMA = vector(mode = "list", length = length(h))

N = 10 #CHANGE THIS TO 500 FOR REAL SIMULATION
links = list("exp", "quad", link, link2)

for(k in 1:nrow(matrix)){
for(i in 1:10){
  SIGMA[[i]] = array(dim = c(dim(datasets[[k]])[1], 5, length(1:N)))
  for(n in 1:N){

    colnames(h_list[[i]]) = NULL
    rownames(h_list[[i]]) = NULL
    X_s = as.matrix(datasets[[k]][,2,n])
    SIGMA[[i]][,1, n] = X_s

    # Compute sigma ----
    start = (11*length(true_beta)) + 1 #To fetch the H coefficients.
    start_theta = (11*length(true_beta)) + 4*(i + 1) + 1
    start_theta_tilde = (11*length(true_beta)) + 4*(i + 1) + 8 + 1

    # LOOP OVER THE SIGMA ESTIMATORS
    for(l in 1:4){
    link_temp = links[[l]]
    H = h_list[[i]][n, start:(start + i), k]
    start = start + i + 1

    theta = h_list[[i]][n, start_theta:(start_theta + 1), k]
    start_theta = start_theta + 2


    theta_tilde = h_list[[i]][n, start_theta_tilde:(start_theta_tilde + 1), k]
    start_theta_tilde = start_theta_tilde + 2



    Her = Her(X_s, deg=i, type=type)
    if(link_temp == "exp"){
      sigma = exp(Her%*%H)/exp(1)
    }

    if(link_temp=="quad"){
      sigma = (Her%*%H)^2
    }

    if (link_temp!="exp" & link_temp!="quad"){
      sigma = as.vector(unlist(lapply(link_temp, function(f) f(Her%*%H))[1]))
      dsigma = as.vector(unlist(lapply(link_temp, function(f) f(Her%*%H))[2]))
    }
    sigma = sigma*sqrt(laguerre_var(theta, theta_tilde, matrix[k, 'tau']))

    SIGMA[[i]][,l + 1, n] = sigma

    }




}
}

  SIGMA_BIG[[k]] = SIGMA
}

# Arrange sigmas in a dataframe
x_s  = as.matrix(datasets[[1]][,2,1])
sigmas = cbind(x_s, (0.2+2*(x_s-0.5)^2)) %>% as.data.frame() %>% # TRUE SIGMA
  mutate(type = "true sigma", iter = NA, degree = NA, dataset = NA) %>%
  rename(c("x" = "V1", "sigma" = "V2"))

for(k in 1:nrow(matrix)){
for(i in 1:10){
    for(n in 1:N){
    sigmas_temp = SIGMA_BIG[[k]][[i]][,,n] %>% as.data.frame() %>%
    rename(c("x" = "V1", "exp" = "V2", "quad" = "V3", "id" = "V4", "abs" = "V5")) %>%
    pivot_longer(cols = exp:abs, names_to = "type", values_to = "sigma") %>%
    mutate(iter = n, degree = i, dataset = k) %>%
    arrange(type, x)

  sigmas = rbind(sigmas, sigmas_temp)
}
}
}


# PLOT AND SAVE IMAGES
# New facet label names for degree variable
labs <- c("Hermite degree 2", "Hermite degree 3", "Hermite degree 4",
               "Hermite degree 5", "Hermite degree 6", "Hermite degree 7",
               "Hermite degree 8", "Hermite degree 9")
names(labs) <- c("2", "3", "4", "5", "6", "7", "8", "9")


for(link in c("exp", "quad", "id", "abs")){
  for(k in 1:4){
  ggsave(plot = sigmas %>%
  filter(type %in% c("true sigma", link), dataset == k, degree %in% c(2, 3, 4, 5, 6, 7, 8, 9)) %>%
  ggplot(aes(x = x, y = sigma, group = iter)) +
  geom_line(color = 'gray') +
  facet_wrap(degree~., ncol = 2, labeller = labeller(degree = labs)) +
  geom_line(data = sigmas %>% filter(type == 'true sigma') %>%
              select(-degree), aes(x = x, y = sigma), color = 'black') +
  ylim(c(0,30)) +
  ylab(label = expression(paste(theta, "(", sigma, ")"))) +
  theme_bw() + theme(panel.border = element_blank(), panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank(), axis.line = element_line(colour = "black")),
  path = paste0("Results Wang/PLOTS/"),
  filename = paste0('Link function_', link, '_', 'quantile_', matrix[k, 'tau'], '_',
                   'Sample_size_', matrix[k, 'n'], ".png"),
  width = 5,
  height = 7

  )
  }
}

# GENERATE BIAS TABLES FOR LATEX FILE ----

est_names = c("Omni", "P & W", "W & W", "W & W (CV)",
              "Portnoy", "DB et. al.", "Laguerre", "Laguerre H. (exp link)",
              "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

bias_tables = c()
for(i in 1:10){

  bias_temp = c()
  for(k in 1:nrow(matrix)){
    bias_temp2 = bias[[i]][,,k] %>% round(digits = 4) %>% as.data.frame() %>%
      mutate(Estimator = est_names, tau = matrix[k, "tau"], n = matrix[k, "n"], degree = i,
             dataset = k) %>%
      rename(c("beta_0" = "V1", "beta_1" = "V2"))

    bias_temp = rbind(bias_temp, bias_temp2)
  }
  bias_tables = rbind(bias_tables, bias_temp)
}

# TABLES LITERATURE
bias_tables_literature = bias_tables %>%
  filter(degree == 1, Estimator %in% c("Omni", "P & W", "W & W", "W & W (CV)",
                                       "Portnoy", "DB et. al.", "Laguerre"))

# apply styling to tables
table1_styled <- kable(bias_tables_literature %>% filter(dataset == 1) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table2_styled <- kable(bias_tables_literature %>% filter(dataset == 2) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table3_styled <- kable(bias_tables_literature %>% filter(dataset == 3) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table4_styled <- kable(bias_tables_literature %>% filter(dataset == 4) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

# combine tables into a grid
table_grid <- cbind(table1_styled, table2_styled, table3_styled, table4_styled)


# TABLES FOR PROPOSED ESTIMATOR

est_proposed = c("Laguerre H. (exp link)",
              "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

bias_tables_temp = bias_tables %>%
  filter(Estimator == 'Laguerre') %>%
  mutate(degree = ifelse(Estimator == 'Laguerre', 0, degree)) %>%
  arrange(dataset) %>%
  group_by(dataset) %>%
  slice_head(n = 4) %>%
  mutate(Estimator = est_proposed) %>%
  ungroup()

bias_tables_proposed = bias_tables %>%
  filter(Estimator %in% est_proposed) %>%
  bind_rows(bias_tables_temp) %>%
  arrange(dataset, Estimator, degree) %>%
  pivot_wider(names_from = Estimator, values_from = c(beta_0, beta_1))

bias_tables_proposed = bias_tables_proposed  %>%
  select(tau, n, degree, dataset,
         ends_with("(Id link)"),
         ends_with("(abs link)"),
         ends_with("(exp link)"),
         ends_with("(quad link)"))

table1 = kable(bias_tables_proposed %>% filter(dataset == 1) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)
table2 = kable(bias_tables_proposed %>% filter(dataset == 2) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)
table3 = kable(bias_tables_proposed %>% filter(dataset == 3) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)
table4 = kable(bias_tables_proposed %>% filter(dataset == 4) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)


# GENERATE MSE TABLES FOR LATEX FILE ----

est_names = c("Omni", "P & W", "W & W", "W & W (CV)",
              "Portnoy", "DB et. al.", "Laguerre", "Laguerre H. (exp link)",
              "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

mse_tables = c()
for(i in 1:10){

  mse_temp = c()
  for(k in 1:nrow(matrix)){
    mse_temp2 = MSE[[i]][,,k] %>% round(digits = 4) %>% as.data.frame() %>%
      mutate(Estimator = est_names, tau = matrix[k, "tau"], n = matrix[k, "n"], degree = i,
             dataset = k) %>%
      rename(c("beta_0" = "V1", "beta_1" = "V2"))

    mse_temp = rbind(mse_temp, mse_temp2)
  }
  mse_tables = rbind(mse_tables, mse_temp)
}

# TABLES LITERATURE
mse_tables_literature = mse_tables %>%
  filter(degree == 1, Estimator %in% c("Omni", "P & W", "W & W", "W & W (CV)",
                                       "Portnoy", "DB et. al.", "Laguerre"))

# apply styling to tables
table1_styled <- kable(mse_tables_literature %>% filter(dataset == 1) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table2_styled <- kable(mse_tables_literature %>% filter(dataset == 2) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table3_styled <- kable(mse_tables_literature %>% filter(dataset == 3) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

table4_styled <- kable(mse_tables_literature %>% filter(dataset == 4) %>%
                         select(-dataset, -degree), "latex", booktabs = T) %>%
  kable_styling(bootstrap_options = c("striped", "hover")) %>%
  column_spec(1, width = "50px")

# combine tables into a grid
table_grid <- cbind(table1_styled, table2_styled, table3_styled, table4_styled)


# TABLES FOR PROPOSED ESTIMATOR

est_proposed = c("Laguerre H. (exp link)",
                 "Laguerre H. (quad link)", "Laguerre H. (Id link)", "Laguerre H. (abs link)")

mse_tables_temp = mse_tables %>%
  filter(Estimator == 'Laguerre') %>%
  mutate(degree = ifelse(Estimator == 'Laguerre', 0, degree)) %>%
  arrange(dataset) %>%
  group_by(dataset) %>%
  slice_head(n = 4) %>%
  mutate(Estimator = est_proposed) %>%
  ungroup()

mse_tables_proposed = mse_tables %>%
  filter(Estimator %in% est_proposed) %>%
  bind_rows(mse_tables_temp) %>%
  arrange(dataset, Estimator, degree) %>%
  pivot_wider(names_from = Estimator, values_from = c(beta_0, beta_1))

mse_tables_proposed = mse_tables_proposed  %>%
  select(tau, n, degree, dataset,
         ends_with("(Id link)"),
         ends_with("(abs link)"),
         ends_with("(exp link)"),
         ends_with("(quad link)"))

table1 = kable(mse_tables_proposed %>% filter(dataset == 1) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)
table2 = kable(mse_tables_proposed %>% filter(dataset == 2) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)
table3 = kable(mse_tables_proposed %>% filter(dataset == 3) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)
table4 = kable(mse_tables_proposed %>% filter(dataset == 4) %>%
                 select(-dataset, -tau, -n), "latex", booktabs = T)





