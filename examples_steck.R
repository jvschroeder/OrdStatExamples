res <- outer(1:50,1:50)

for(n1 in 1:50) {
  for(n2 in 1:50) {
    print(c(n1,n2))
    res[n1,n2] <- steck2_slow_opcount(n1,n2)
  }
}

coeffs <- c()
for(i in 1:50) {
  df <- data.frame(n=1:50,k=res[,i])
  fit <- lm(k ~ poly(n,3,raw=T),df)
  coeffs <- rbind(coeffs,fit$coefficients)
}

for(i in 1:4) {
  df <- data.frame(n=1:50,coeff=coeffs[,i])
  cf <- lm(coeff ~ poly(n,3,raw=T),df)$coefficients
  print(cf[abs(cf)>1e-6])
}


#kfn <- function(n1,n2) 2-5*n2+6.5*n2^2+0.5*n2^3 + n1 * (-5-5*n2+9.5*n2^2+0.5*n2^3) + n1^2 *(6.5+9.5*n2+3*n2^2) + n1^3*(0.5+0.5*n2)

kfn <- function(n1,n2) 0.5 * (n1^3*(n2+1)+n2^3*(n1+1)) + 3 * n1^2*n2^2 + 9.5 * (n1^2*n2 + n2^2*n1) + 6.5 * (n1^2+n2^2) - 5*(n1*n2+n1+n2) + 2

sum(abs(outer(1:50,1:50,kfn)-res))

res1 <- outer(1:50,1:50)
for(n1 in 1:50) {
  for(n2 in 1:50) {
    print(c(n1,n2))
    res1[n1,n2] <- steck2_opcount(n1,n2)
  }
}

hist((res1-res)/res)

hist(outer(1:500,1:500,function(n1,n2) kfn1(n1,n2)/kfn(n1,n2)),xlim=c(0,10),breaks=100)