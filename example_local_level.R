local_level <- function(n,alpha) {
  interval <- c(alpha/n,alpha)
  i <- 1:n
  
  expectation <- function(v) 1-noe2_faithful_p(1-rev(qbeta(v,i,n-i+1)),1:n,n,0)[n+1]
  
  for(j in 1:1000) {
    e <- expectation(mean(interval))
    if(e > alpha) {
      interval[2] <- mean(interval)
    } else {
      interval[1] <- mean(interval)
    }
    if( diff(interval) < .Machine$double.eps) break
  }
  print(expectation(mean(interval)))
  return(mean(interval))
}

n <- 1000
a <- local_level(n,0.025)
i <- 1:n
t <- qbeta(a,i,n-i+1)
Fn <- function(x) pnorm(qnorm(x),mean=-1,sd=1.2)
t1 <- 1-rev(t)
t2 <- 1-rev(Fn(t))

#rowIndices + nrow(valueMatrix) * (colIndices - 1)

#k <- 65
#print(system.time(m <- 1-noe2_faithful_p(t1,t2,n,k)))
#pow <- m[cbind((n+1):(n+1-k),1:(k+1))]
#plot((0:k)/n,pow,type='l')

#print(system.time(m <- 1-noe2_p(t1,t2,n,k)))
#pow1 <- m[cbind((n+1):(n+1-k),1:(k+1))]
#plot((0:k)/n,pow1,col='red')

k <- 65
n <- 100
a <- local_level(n,0.025)
i <- 1:n
t <- qbeta(a,i,n-i+1)

# pows <- c()
# for(d in seq(-0.1,-1,length.out=300)){
#   Fn <- function(x) pnorm(qnorm(x),mean=d,sd=1.2)
#   t1 <- 1-rev(t)
#   t2 <- 1-rev(Fn(t))
#   m <- 1-noe2_p(t1,t2,n,k)
#   pows <- c(pows,m[n+1-k,k+1])
#   print(d)
# }
# 
# plot(seq(-0.1,-1,length.out=300),pows,type='l')

pows <- c()
for(d in seq(0,-1,length.out=300)){
  Fn <- function(x) pnorm(qnorm(x),mean=d)
  t1 <- 1-rev(t)
  t2 <- 1-rev(Fn(t))
  m <- 1-noe2_faithful_p(t1,t2,n,k)
  pows <- c(pows,m[n+1-k,k+1])
  print(d)
}

res <- data.frame(x=seq(0,-1,length.out=300),y=pows)

tikzDevice::tikz(file = "../graphics/power_local_mu.tex", width = 5, height = 5)
plot <- ggplot2::ggplot(res) + 
  ggplot2::geom_line(ggplot2::aes(x=x,y=y)) +
  ggplot2::xlab("$\\mu_1$") + ggplot2::ylab("Power") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::theme(legend.position = "none")
print(plot)
dev.off()
fix_tikz("../graphics/power_local_mu.tex")

pows <- c()
for(k in 1:(n-1)){
  Fn <- function(x) pnorm(qnorm(x),mean=-1)
  t1 <- 1-rev(t)
  t2 <- 1-rev(Fn(t))
  m <- 1-noe2_faithful_p(t1,t2,n,k)
  pows <- c(pows,m[n+1-k,k+1])
  print(k)
}

res <- data.frame(x=1:(n-1),y=pows)

breaks = c(1,seq(25, 100, by=25), 31)
# and labels
labels = as.character(breaks)

breaksy = seq(0,1,by=.2)
labelsy = as.character(breaksy)

tikzDevice::tikz(file = "../graphics/power_local_num_cont.tex", width = 5, height = 5)
plot <- ggplot2::ggplot(res) + 
  ggplot2::geom_vline(xintercept=min(which(pows>=0.8)),colour = "grey") +
  ggplot2::geom_hline(yintercept=0.8,colour = "grey") +
  ggplot2::geom_point(ggplot2::aes(x=x,y=y)) +
  ggplot2::xlab("$k$") + ggplot2::ylab("Power") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1) +
  ggplot2::theme(legend.position = "none") +
  ggplot2::scale_x_continuous(limits = c(1, 100), breaks = breaks, labels = labels) +
  ggplot2::scale_y_continuous(limits = c(0, 1), breaks = breaksy, labels = labelsy)
print(plot)
dev.off()
fix_tikz("../graphics/power_local_num_cont.tex")