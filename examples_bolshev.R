#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

fn <- function(i) bolshev_exact(1-rev(0.05*(1:i)/i))

res1 <- data.frame(x=1:150,fn="Bolshev Rational",lq=0,median=0,uq=0)
for(n in 1:150) {
  m <- microbenchmark::microbenchmark(fn(n), times=200)
  r <- summary(m,"ms")[,c("lq","median","uq")]
  res1[n,c("lq","median","uq")] <- r[1,]
  print(n)
}

res <- tidyr::gather(rbind(res1),key,value,-x,-lq,-uq,-median,-fn)

tikzDevice::tikz(file = "../graphics/bolshev_single_execution_time.tex", width = 5, height = 5)
plot <- ggplot2::ggplot(res,ggplot2::aes(fill=fn)) + ggplot2::geom_ribbon(ggplot2::aes(x=x,ymin=lq,ymax=uq),alpha=0.3,show.legend=F) + ggplot2::geom_line(ggplot2::aes(x=x,y=median,colour=fn),show.legend=F) +
  ggplot2::xlab("$n$") + ggplot2::ylab("Execution time in ms") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) + ggplot2::scale_color_grey() + ggplot2::scale_fill_grey()
print(plot)
dev.off()
fix_tikz("../graphics/bolshev_single_execution_time.tex")

res <- outer(1:100,1:100)

for(n1 in 1:100) {
  for(n2 in 1:100) {
    print(c(n1,n2))
    res[n1,n2] <- bolshev2_opcount(n1,n2)
  }
}

coeffs <- c()
for(i in 1:100) {
  df <- data.frame(n=1:100,k=res[,i])
  fit <- lm(k ~ poly(n,2,raw=T),df)
  coeffs <- rbind(coeffs,fit$coefficients)
}

for(i in 1:3) {
  df <- data.frame(n=1:100,coeff=coeffs[,i])
  print(lm(coeff ~ poly(n,2,raw=T),df)$coefficients)
}

#kfn <- function(n1,n2) 2+3*n2+3*n2^2+(3+7.5*n2+4.5*n2^2)*n1+(3+4.5*n2+1.5*n2^2)*n1^2
kfn <- function(n1,n2) 2 + 3 * (n1 + n1^2 + n2 + n2^2) + 7.5 * n1 * n2 + 4.5 * (n1^2 * n2 + n1 * n2^2) + 1.5 * n1^2 * n2^2

sum(abs(outer(1:100,1:100,kfn)-res))


fn <- function(i) {
  i <- max(i,10)
  bolshev2_exact(0.05*(1:(2*i))/(2*i),0.05*(1:(2*i))/(2*i),i,i)
} 
fn3 <- function(i) {
  i <- max(i,10);
  bolshev2_double(0.05*(1:(2*i))/(2*i),0.05*(1:(2*i))/(2*i),i,i)
}

res1 <- data.frame(x=1:120,fn="Bolshev Rational",lq=0,median=0,uq=0)
res4 <- data.frame(x=1:120,fn="Bolshev (double)",lq=0,median=0,uq=0)

for(n in 1:50) {
  print(n)
  m <- microbenchmark::microbenchmark(fn(n),fn3(n), times=10)
  r <- summary(m,"ms")[,c("lq","median","uq")]
  res1[n,c("lq","median","uq")] <- r[1,]
  res4[n,c("lq","median","uq")] <- r[2,]
  print(n)
  save(res1,res4,file='D:/test.RData')
}

res <- tidyr::gather(rbind(res1,res2,res3),key,value,-x,-lq,-uq,-median,-fn)
plot(res1$median[10:50] / kfn(10:50,10:50))
plot(res4$median[10:50] / kfn(10:50,10:50))

res4 <- res4[10:50,]
res4$fn <- kfn(res4$x,res4$x)

tikzDevice::tikz(file = "../graphics/double_bolshev.tex", width = 5, height = 5)
plot <- ggplot2::ggplot(res4,ggplot2::aes(fill=fn)) + ggplot2::geom_point(ggplot2::aes(x=x,y=median/fn)) + ggplot2::guides(fill=FALSE) +
  ggplot2::xlab("$\\ell$") + ggplot2::ylab("Execution time in ms / \\#ops") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1)
print(plot)
dev.off()
fix_tikz("../graphics/double_bolshev.tex")

res1 <- res1[10:50,]
res1$fn <- kfn(res1$x,res1$x)

tikzDevice::tikz(file = "../graphics/rational_bolshev.tex", width = 5, height = 5)
plot <- ggplot2::ggplot(res1,ggplot2::aes(fill=fn)) + ggplot2::geom_point(ggplot2::aes(x=x,y=median/fn)) + ggplot2::guides(fill=FALSE) +
  ggplot2::xlab("$\\ell$") + ggplot2::ylab("Execution time in ms / \\#ops") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1)
print(plot)
dev.off()
fix_tikz("../graphics/rational_bolshev.tex")