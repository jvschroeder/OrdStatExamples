#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

data <- read.table(textConnection("0.60 5 70 0.15 0.704 0.426 0.524 0.216 0.249 0.570 1.500
0.60 5 80 0.15 0.795 0.584 0.681 0.308 0.396 0.622 1.500
0.60 5 90 0.15 0.860 0.738 0.797 0.409 0.538 0.673 1.467
0.60 5 100 0.15 0.906 0.867 0.875 0.518 0.657 0.721 1.380
0.60 20 50 0.15 0.664 0.265 0.307 0.042 0.028 0.607 1.500
0.60 20 60 0.15 0.781 0.610 0.630 0.139 0.157 0.695 1.500
0.60 20 70 0.15 0.859 0.894 0.877 0.320 0.378 0.765 1.343
0.60 20 80 0.15 0.911 0.991 0.961 0.563 0.599 0.819 1.225
0.60 60 40 0.15 0.706 0.280 0.315 0.005 0.002 0.666 1.500
0.60 60 50 0.15 0.823 0.903 0.891 0.088 0.099 0.771 1.360
0.60 60 60 0.15 0.895 1.000 0.995 0.453 0.492 0.842 1.183
0.60 100 30 0.15 0.632 0.037 0.031 0.000 0.000 0.609 1.500
0.60 100 40 0.15 0.788 0.786 0.789 0.011 0.006 0.750 1.450
0.60 100 50 0.15 0.879 1.000 0.999 0.278 0.270 0.838 1.200
0.80 5 40 0.15 0.697 0.417 0.496 0.212 0.252 0.566 1.500
0.80 5 50 0.15 0.844 0.696 0.784 0.380 0.502 0.659 1.500
0.80 5 60 0.15 0.924 0.915 0.895 0.574 0.704 0.742 1.333
0.80 20 30 0.15 0.692 0.330 0.352 0.057 0.036 0.626 1.500
0.80 20 40 0.15 0.859 0.892 0.893 0.319 0.357 0.764 1.350
0.80 60 20 0.15 0.616 0.062 0.060 0.001 0.000 0.590 1.500
0.80 60 30 0.15 0.845 0.963 0.952 0.148 0.147 0.791 1.333
0.80 100 20 0.15 0.717 0.286 0.284 0.001 0.003 0.684 1.500
0.80 100 30 0.15 0.896 1.000 1.000 0.453 0.506 0.855 1.167
1.00 5 30 0.15 0.790 0.574 0.680 0.304 0.392 0.617 1.500
1.00 20 20 0.15 0.699 0.350 0.431 0.064 0.045 0.630 1.500
1.00 20 30 0.15 0.913 0.992 0.966 0.579 0.614 0.821 1.200
1.00 60 20 0.15 0.853 0.978 0.959 0.183 0.225 0.799 1.300
1.00 100 20 0.15 0.904 1.000 0.999 0.549 0.586 0.863 1.150"), sep = "" , header = F,
           na.strings ="", stringsAsFactors= F)

data <- dplyr::select(data,c(V1:V3,V9))

i <- 7
lambda <- 0.01
m <- 200

n <- data$V3[i]
r <- data$V2[i]/200
alpha <- 0.01
theta <- data$V1[i]

df <- 2*n-2
ncp <- sqrt(n/2) * theta

t <- alpha * (1:m)/m
pr <- 1-r
Fn <- function(t) pt(qt(t/2,df,0,F),df,ncp,F)-pt(-qt(t/2,df,0,F),df,ncp)
print(system.time(res <- fd_fm_lambda_pwr(1-rev(t),1-rev(Fn(t)),m,pr,lambda)))
print(system.time(res1 <- fd_fm_lambda_pwr_steck(1-rev(t),1-rev(Fn(t)),m,pr,lambda)))
print(correctDigits(res,res1))

data$res <- 1:nrow(data)
data$cDigits <- data$res
data$cDigits1 <- data$res
data$our_sim <- 1:nrow(data)

for(i in 1:nrow(data)) {
  print(i)
  lambda <- 0.9
  m <- 200
  
  n <- data$V3[i]
  r <- data$V2[i]/200
  alpha <- 0.15
  theta <- data$V1[i]
  
  df <- 2*n-2
  ncp <- sqrt(n/2) * theta
  
  t <- alpha * (1:m)/m
  pr <- 1-r
  Fn <- function(t) pt(qt(t/2,df,0,F),df,ncp,F)-pt(-qt(t/2,df,0,F),df,ncp)
  
  data$res[i] <- fd_fm_lambda_pwr(1-rev(t),1-rev(Fn(t)),m,pr,lambda)
  res_double <- fd_fm_lambda_pwr_double(1-rev(t),1-rev(Fn(t)),m,pr,lambda)
  data$cDigits[i] <- correctDigits(data$res[i],res_double)
  res_steck <- fd_fm_lambda_pwr_steck(1-rev(t),1-rev(Fn(t)),m,pr,lambda)
  data$cDigits1[i] <- correctDigits(data$res[i],res_steck)
  
  set.seed(0);
  data$our_sim[i] <- var(replicate(3e2,{
    fn <- function(df,ncp){
      m0 <- rbinom(1,m,1-r)
      p <- 2*pt(abs(c(rt(m0,df,0),rt(m-m0,df,ncp))),df,0,F)
      sum(p.adjust(p,"BH")[-(1:m0)]<=alpha)/(m-m0)
    }
  
    res <- replicate(1e3,fn(df,ncp))
    sum(res>=0.9,na.rm=T)/1e3
  }))
  print(data$our_sim[i])
}

colnames(data) <- c("Eff Sz. $\\theta$","$\\mathbb{E}(M_m)$","n","est. $\\lambda_{90}$-pwr","$\\lambda_{90}$-pwr", "Diff in std")

data[,6] <- abs(data[,4]-data[,5])/sqrt(data[,6])

tbl <- xtable::xtable(data,digits=5,auto=T,caption="Exact calculation of the $\\lambda_{90}-$power cf. \\cite[Table 3]{izmirlian2018average}.",label="table:lambda_power")
#print(tbl,file="../tables/lambda_power.tex",sanitize.text.function=function(x){x})