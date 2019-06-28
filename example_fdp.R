#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

fdp_dist_r <- function(order_stat,t1,t2,m,m0) {
  dist <- fd_fm_order_stat(t1,t2,m0,m,order_stat)
  
  # vals <- matrix(0,m0+1,m+1)
  # for(k in 1:(m+1)) {
  #   for(j in 1:min(m0+1,k)) {
  #     vals[j,k] <- (j-1)/max(k-1,1)
  #   }
  # }
  
  vals <- outer(1:(m0+1),1:(m+1),function(j,k) (j-1)/pmax(k-1,1) * (j<=k))
  
  supp <- sort(unique(as.vector(vals)))
  p <- supp
  
  for(i in seq_along(supp)) {
    p[i] <- sum(dist[vals==supp[i]],na.rm=T)
  }
  
  list(supp=supp,p=p)
}

fdp.conf.upper <- function(conf,alpha=0.05,df=17,ncp=2.5,m) {
  Fn <- function(t) pt(qt(t/2,df,0,F),df,ncp,F)-pt(-qt(t/2,df,0,F),df,ncp)
  t <- alpha * (1:m)/m
  t1 <- 1-rev(t)
  t2 <- 1-rev(Fn(t))
  order_stat <- noe2_p(t1,t2,m,m)
  res <- 0:m
  for(m0 in 0:m) {
    pdist <- fdp_dist(order_stat,t1,t2,m,m0)
    res[m0+1] <- min(pdist$supp[cumsum(pdist$p)>=1-conf])
  }
  res
}

N <- 50
Fn <- function(t) 1 + pnorm(qnorm(t/2)-sqrt(N)) - pnorm(qnorm(1-t/2)-sqrt(N))
m <- 50
m0 <- 5

plts <- list()
l <- 1


correctDigits <- function(p1,p2) {
  res <- floor(log2(abs(p1-p2)))-floor(log2(p1))
  res[p1-p2 == 0] <- -53
  res[p1 == 0] <- ifelse(p2[p1 == 0]==0,-53,floor(log2(abs(p2))))
  res[res>0] <- 0
  -res
}

cDigits1 <- c()
cDigits2 <- c()

for(alpha in c(0.15,0.2,0.3,0.4)) {
  t <- alpha * (1:m)/m
  dist <- fd_fm(1-rev(t),1-rev(Fn(t)),m0,m)
  dist1 <- fd_fm_order_stat(1-rev(t),1-rev(Fn(t)),m0,m,steck2_double_p(1-rev(t),1-rev(Fn(t)),m0,m))
  dist2 <- fd_fm_double(1-rev(t),1-rev(Fn(t)),m0,m)
  probs <- matrix(0,m0+1,m+1)
  probs1 <- probs
  probs2 <- probs
  vals <- matrix(0,m0+1,m+1)
  for(k in 0:m) {
    for(j in 0:min(m0,k)) {
      vals[j+1,k+1] <- j/max(k,1)
      probs[j+1,k+1] <- dist[j+1,k+1]
      probs1[j+1,k+1] <- dist1[j+1,k+1]
      probs2[j+1,k+1] <- dist2[j+1,k+1]
    }
  }
  supp <- sort(unique(as.vector(vals)))
  supp <- supp[supp<=0.2]
  p <- supp
  p1 <- supp
  p2 <- supp
  for(i in seq_along(supp)) {
    p[i] <- sum(probs[vals==supp[i]])
    p1[i] <- sum(probs1[vals==supp[i]])
    p2[i] <- sum(probs2[vals==supp[i]])
  }
  cDigits1 <- c(cDigits1,correctDigits(p,p1))
  cDigits2 <- c(cDigits2,correctDigits(p,p2))
  df <- data.frame(supp=supp,p=p)
  plts[[l]] <- (ggplot2::ggplot(df,ggplot2::aes(x=supp,y=p)) + ggplot2::geom_bar(stat="identity",width=0.0025, fill='lightgrey') + ggplot2::geom_vline(xintercept=sum(vals * probs),col="black",linetype='dotted') +
        ggplot2::xlab("x") + ggplot2::ylab("$\\mathbb{P}(\\text{FDP}\\leq x)$") +
        ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
        ggplot2::ggtitle(paste0("$\\alpha=",alpha,"$")))
  l <- l+1
}

library(tikzDevice)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}\n"))

tikzDevice::tikz(file = "../graphics/fdp_dist_digits.tex", width = 5, height = 5)
tbl <- table(cDigits1)
ggplot2::ggplot(data.frame(x=as.numeric(names(tbl)),y=as.numeric(tbl)),
                ggplot2::aes(x=x, y=y)) +
  ggplot2::geom_bar(stat="identity") +
  ggplot2::xlab("number of correct digits (in base 2)") + ggplot2::ylab("frequency") +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank())
dev.off()
fix_tikz("../graphics/fdp_dist_digits.tex")

tbl1 <- table(cDigits2)
plt2 <- ggplot2::ggplot(data.frame(x=as.numeric(names(tbl1)),y=as.numeric(tbl1)),
                        ggplot2::aes(x=x, y=y)) +
  ggplot2::geom_bar(stat="identity") +
  ggplot2::xlab("number of correct digits (in base 2)") + ggplot2::ylab("frequency") +
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank())
dev.off()

tikzDevice::tikz(file = "../graphics/fdp_dist.tex", width = 5, height = 5)
print(do.call(gridExtra::grid.arrange,plts))
dev.off()
fix_tikz("../graphics/fdp_dist.tex")
