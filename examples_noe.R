fn <- function(i) bolshev2_exact_p(0.05*(1:(2*i))/(2*i),0.05*(1:(2*i))/(2*i),i,i)
fn1 <- function(i) noe2_faithful_p(0.05*(1:(2*i))/(2*i),0.05*(1:(2*i))/(2*i),i,i)
fn2 <- function(i) steck2_exact_p(0.05*(1:(2*i))/(2*i),0.05*(1:(2*i))/(2*i),i,i)

res1 <- data.frame(x=1:50,fn="Bolshev Rational",lq=0,median=0,uq=0)
res2 <- data.frame(x=1:50,fn="Noe Faithful",lq=0,median=0,uq=0)
res3 <- data.frame(x=1:50,fn="Steck Rational",lq=0,median=0,uq=0)
for(n in 1:50) {
  print(n)
  m <- microbenchmark::microbenchmark(fn(n),fn1(n),fn2(n), times=10)
  r <- summary(m,"ms")[,c("lq","median","uq")]
  res1[n,c("lq","median","uq")] <- r[1,]
  res2[n,c("lq","median","uq")] <- r[2,]
  res3[n,c("lq","median","uq")] <- r[3,]
  print(n)
}

res <- tidyr::gather(rbind(res1,res2,res3),key,value,-x,-lq,-uq,-median,-fn)

plot1 <- ggplot2::ggplot(res,ggplot2::aes(fill=fn)) + ggplot2::geom_line(ggplot2::aes(x=x,y=median,linetype=fn)) +
  ggplot2::xlab("$\\ell$") + ggplot2::ylab("Execution time in ms") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::guides(linetype=ggplot2::guide_legend(title="Algorithm"))


complexity.noe <- function(n1,n2) 0.5*n2+1.5*n2^2+n2^3+(1.5+1.5*n2+4.25*n2^2+1.25*n2^3)*n1+(1.5+4.25*n2+3*n2^2+0.25*n2^3)*n1^2+(1+1.25*n2+0.25*n2^2)*n1^3
complexity.bolshev <- function(n1,n2) 2 + 3 * (n1 + n1^2 + n2 + n2^2) + 7.5 * n1 * n2 + 4.5 * (n1^2 * n2 + n1 * n2^2) + 1.5 * n1^2 * n2^2
complexity.steck <- function(n1,n2) 0.5 * (n1^3*(n2+1)+n2^3*(n1+1)) + 3 * n1^2*n2^2 + 9.5 * (n1^2*n2 + n2^2*n1) + 6.5 * (n1^2+n2^2) - 5*(n1*n2+n1+n2) + 2

res$complexity <- res$median
comp <- complexity.bolshev
fn <- "Bolshev Rational"
res$complexity[res$fn == fn] <- res$complexity[res$fn == fn] / comp(res$x[res$fn == fn],res$x[res$fn == fn])
comp <- complexity.steck
fn <- "Steck Rational"
res$complexity[res$fn == fn] <- res$complexity[res$fn == fn] / comp(res$x[res$fn == fn],res$x[res$fn == fn])
comp <- complexity.noe
fn <- "Noe Faithful"
res$complexity[res$fn == fn] <- res$complexity[res$fn == fn] / comp(res$x[res$fn == fn],res$x[res$fn == fn])

plot2 <- ggplot2::ggplot(res,ggplot2::aes(fill=fn)) + ggplot2::geom_line(ggplot2::aes(x=x,y=log10(complexity),linetype=fn)) +
  ggplot2::xlab("$\\ell$") + ggplot2::ylab("$\\log_{10}$ of ms/arithmetic operation") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  #ggplot2::theme(aspect.ratio = 1) +
  ggplot2::guides(linetype=ggplot2::guide_legend(title="Algorithm"))

plot3 <- ggplot2::ggplot(res[res$fn=="Noe Faithful",],ggplot2::aes(fill=fn)) + ggplot2::geom_line(ggplot2::aes(x=x,y=median,linetype=fn)) +
  ggplot2::xlab("$\\ell$") + ggplot2::ylab("Execution time in ms") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::guides(linetype=ggplot2::guide_legend(title="Algorithm"))

tikzDevice::tikz(file = "../graphics/comp_bolshev_noe.tex", width = 10, height = 15)
print(gridExtra::grid.arrange(plot1,plot3,plot2,layout_matrix=rbind(c(1,2),c(3,3))))
dev.off()
fix_tikz("../graphics/comp_bolshev_noe.tex")


#Calculation of k(n_1,n_2) and the actual complexity (as a function of n_1 and n_2)

t <- 0.5 * (1:800)/800
print(system.time(res <- noe_faithful_k(t,t,200,200)))

coeffs <- c()
for(i in 2:201) {
  df <- data.frame(n=1:200,k=res[2:201,i])
  fit <- lm(k ~ poly(n,1,raw=T),df)
  coeffs <- rbind(coeffs,fit$coefficients)
}

for(i in 1:2) {
  df <- data.frame(n=1:200,coeff=coeffs[1:200,i])
  print(round(lm(coeff ~ poly(n,1,raw=T),df)$coefficients,2))
}

kfn <- function(n1,n2) -7+8*(n1+n2)+n1*n2
#outer(0:200,0:200,kfn)-res
sum(abs(outer(2:200,2:200,kfn)-res[3:201,3:201]))



res <- outer(1:50,1:50)

for(n1 in 1:50) {
  for(n2 in 1:50) {
    print(c(n1,n2))
    res[n1,n2] <- noe2_opcount(n1,n2)
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
  cf <-  lm(coeff ~ poly(n,3,raw=T),df)$coefficients
  print(cf[abs(cf)>1e-6])
}

kfn <- function(n1,n2) 0.5*n2+1.5*n2^2+n2^3+(1.5+1.5*n2+4.25*n2^2+1.25*n2^3)*n1+(1.5+4.25*n2+3*n2^2+0.25*n2^3)*n1^2+(1+1.25*n2+0.25*n2^2)*n1^3

sum(abs(outer(1:50,1:50,kfn)-res))