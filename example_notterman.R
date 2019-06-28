ncps <- seq(2.8,6.7,0.1)
do_calc <- function(n=125,R=95) {
  t <- (1:n)/n*0.05
  k <- 89
  res <- matrix(0,length(ncps),k)
  for(ncp_i in seq_along(ncps)) {
    ncp <- ncps[ncp_i]
    print(c(n,R,ncp_i))
    gamma_star <- max(fdp.conf.upper(0.05,alpha=0.025,df=17,ncp=ncp,n))
    m0.upper <- floor(n-(1-gamma_star)*R)
    for(df in 1:k) {
      Fn <- function(t) pt(qt(t/2,df,0,F),df,ncp,F)-pt(-qt(t/2,df,0,F),df,ncp)
      print(res[ncp_i,df] <- min(fd_fm_pwr_lim(1-rev(t),1-rev(Fn(t)),n,m0.upper+1)[1:(m0.upper+1)]))
    }
  }
  res
}

res1 <- do_calc(79, 56)
save(res1,file='D:/results/res1_samp.RData')

plot(res[1,],ylim=c(0,1))
sampSize1 <- 1:length(ncps)
for(i in 1:length(ncps)) {
  points(res1[i,])
  print(sampSize1[i] <- min(which(res1[i,]>=0.8)))
}
abline(h=0.8)

res3 <- do_calc(126, 96)
save(res3,file='D:/results/res3_samp.RData')
sampSize3 <- 1:length(ncps)
for(i in 1:length(ncps)) {
  print(sampSize3[i] <- min(which(res3[i,]>=0.8)))
}

res2 <- do_calc(225, 190)
save(res2,file='D:/results/res2_samp.RData')
sampSize2 <- 1:length(ncps)
for(i in 1:length(ncps)) {
  print(sampSize2[i] <- min(which(res2[i,]>=0.8)))
}


n=126
R=96

t <- (1:n)/n*0.05
k <- 89
df <- 17
ncp <- 3.0
m0 <- 60
Fn <- function(t) pt(qt(t/2,df,0,F),df,ncp,F)-pt(-qt(t/2,df,0,F),df,ncp)
res <- fd_fm_pwr(1-rev(t),1-rev(Fn(t)),n)

tikzDevice::tikz(file = "../graphics/example_power.tex", width = 5, height = 5)
ggplot2::ggplot(data.frame(res,m0=0:(length(res)-1))) + ggplot2::geom_point(size=1,ggplot2::aes(x=m0,y=res)) +
  ggplot2::xlab("$m_0$") + ggplot2::ylab("$\\textrm{Pow}_{\\textrm{avg}}$") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1,legend.title = ggplot2::element_blank())
dev.off()
fix_tikz("../graphics/example_power.tex")


tikzDevice::tikz(file = "../graphics/example_sampsize.tex", width = 5, height = 5)
sampSize1[!is.finite(sampSize1)] <- NA
sampSize2[!is.finite(sampSize2)] <- NA
sampSize3[!is.finite(sampSize3)] <- NA
res <- tidyr::gather(data.frame(ncps,sampSize1,sampSize3,sampSize2),key,sampSize,-ncps)
plot(ggplot2::ggplot(res,ggplot2::aes(fill=key)) + ggplot2::geom_point(size=2,ggplot2::aes(x=ncps,y=sampSize)) + ggplot2::geom_line(size=1,ggplot2::aes(x=ncps,y=sampSize,linetype=key,colour=key)) +
  ggplot2::xlab("Effect Size $\\mu$") + ggplot2::ylab("Sample Size $n$") +
  ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
  ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::theme(aspect.ratio = 1,legend.title = ggplot2::element_blank()))
dev.off()
fix_tikz("../graphics/example_sampsize.tex")

tikzDevice::tikz(file = "../graphics/example_sampsize1.tex", width = 5, height = 5)
sampSize3[!is.finite(sampSize3)] <- NA
res <- tidyr::gather(data.frame(ncps,sampSize3),key,sampSize,-ncps)
plot(ggplot2::ggplot(res,ggplot2::aes(fill=key)) + ggplot2::geom_point(size=2,ggplot2::aes(x=ncps,y=sampSize)) + ggplot2::geom_line(size=1,ggplot2::aes(x=ncps,y=sampSize,linetype=key)) +
       ggplot2::xlab("Effect Size $\\mu$") + ggplot2::ylab("Sample Size $n$") +
       ggplot2::theme(axis.title.x = ggplot2::element_text(size = 15, vjust=-.2)) +
       ggplot2::theme(axis.title.y = ggplot2::element_text(size = 15, vjust=0.3)) +
       ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
       ggplot2::theme(aspect.ratio = 1,legend.title = ggplot2::element_blank()))
dev.off()
fix_tikz("../graphics/example_sampsize1.tex")