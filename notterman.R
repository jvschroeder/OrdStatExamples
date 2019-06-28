library(dplyr)
df <- readxl::read_xls("./data/CarcinomaNormaldatasetCancerResearch.xls",range="Cancer!A9:AN7465",col_names=F)
header <- unlist(readxl::read_xls("./data/CarcinomaNormaldatasetCancerResearch.xls",range="Cancer!A1:AN1",col_names=F), use.names = F)
colnames(df) <- header
df <- df %>% dplyr::select(-Sample)
df.tumor <- df[,3:20]
df.normal <- df[,21:38]

D <- df.tumor-df.normal
n1 <- ncol(df.tumor)
df$Stat <- sqrt(n1) * apply(D,1,mean) / apply(D,1,sd)
df$p <- 2*pt(abs(df$Stat),n1-1,lower.tail = F)

df$tumor.mean <- apply(df.tumor,1,mean)
df$normal.mean <- apply(df.normal,1,mean)
df$over <- df$tumor.mean/pmax(df$normal.mean,10)
df$under <- df$normal.mean/pmax(df$tumor.mean,10)
df$over_under <- pmax(df$over,df$under)

df <- df[(df$tumor.mean > 10 | df$normal.mean > 10),]

rejections <- function(t,alpha) {
  c(sum(df$over_under>=t),
    sum(p.adjust(df$p[df$over_under>=t],'BH')<=alpha))
}

alpha <- 0.025
rejections(3,alpha)
rejections(4,alpha)
rejections(5,alpha)