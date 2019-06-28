#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

t <- 0.05 * (1:5)/5
m0 <- 3
m <- 5
N <- 5
#(24) in Glueck (2008)
Fn <- function(t) 1 + pnorm(qnorm(t/2)-sqrt(N)) - pnorm(qnorm(1-t/2)-sqrt(N))

n <- 150

system.time({res <- matrix(NA,n-1,n)
for(m in 2:n) {
  t <- 0.05 * (1:m)/m
  res[m-1,1:m] <- fd_fm_pwr(1-rev(t),1-rev(Fn(t)),m)
}
})

system.time({res1 <- matrix(NA,n-1,n)
for(m in 2:n) {
  t <- 0.05 * (1:m)/m
  res1[m-1,1:m] <- fd_fm_pwr_bolshev(1-rev(t),1-rev(Fn(t)),m)
}
})

system.time({res2 <- matrix(NA,n-1,n)
for(m in 2:n) {
  t <- 0.05 * (1:m)/m
  res2[m-1,1:m] <- fd_fm_pwr_bolshev_exact(1-rev(t),1-rev(Fn(t)),m)
}
})

m <- 90
t <- sort(runif(100))#0.05 * (1:m)/(m*cumsum(1/(1:m)))
order_stat <- steck2_lower(m+1,m+1,1-rev(t),1-rev(Fn(t)))
print(res <- fd_fm_pwr_order_stat(1-rev(t),1-rev(Fn(t)),m,order_stat))

res- fd_fm_pwr(1-rev(t),1-rev(Fn(t)),m)

plot(log10(abs(res2[,1]-res1[,1])/res2[,1]))
plot(log10(abs(res2[,1]-res[,1])/res2[,1]))

rownames(res) <- 2:50
colnames(res) <- 0:49
tbl <- xtable::xtable(res[,1:6],digits=5,auto=T)
print(tbl,file="../tables/average_power.tex")

fn <- function(m){
  t <- 0.05 * (1:m)/m
  fd_fm_pwr(1-rev(t),1-rev(Fn(t)),m)
}

res1 <- data.frame(x=1:150,fn="C++ Implementation",lq=0,median=0,uq=0)
for(n in 1:150) {
  m <- microbenchmark::microbenchmark(fn(n), times=10)
  r <- summary(m,"ms")[,c("lq","median","uq")]
  res1[n,c("lq","median","uq")] <- r[1,]
  print(n)
}

res <- tidyr::gather(rbind(res1),key,value,-x,-lq,-uq,-median,-fn)
res$op_complexity <- res$median/complexity.noe(res$x,res$x)

tikzDevice::tikz(file = "../graphics/average_power_time.tex", width = 10, height = 10)
plot1 <- ggplot2::ggplot(res,ggplot2::aes(fill=fn)) + ggplot2::geom_ribbon(ggplot2::aes(x=x,ymin=lq,ymax=uq),alpha=0.3,show.legend=F) + ggplot2::geom_line(ggplot2::aes(x=x,y=median,colour=fn),show.legend=F) +
  ggplot2::xlab("$m$") + ggplot2::ylab("Execution time in ms") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) + ggplot2::scale_color_grey() + ggplot2::scale_fill_grey()
plot2 <- ggplot2::ggplot(res,ggplot2::aes(fill=fn)) + ggplot2::geom_line(ggplot2::aes(x=x,y=log10(op_complexity),colour=fn),show.legend=F) +
  ggplot2::xlab("$m$") + ggplot2::ylab("$\\log_{10}$ of ms/arithmetic operation") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) + ggplot2::scale_color_grey() + ggplot2::scale_fill_grey()
print(gridExtra::grid.arrange(plot1,plot2,layout_matrix=rbind(c(1,2))))
dev.off()
fix_tikz("../graphics/average_power_time.tex")