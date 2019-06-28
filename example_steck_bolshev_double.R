#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

p <- c(rep(2^-10,10),rep(2^(-1),1))
df <- data.frame(exact=bolshev_exact(p))
df$steck <- steck_double(p)
df$steck_rel <- abs(df$exact-df$steck)/df$exact
df$bolshev <- bolshev_double(p)
df$bolshev_rel <- abs(df$exact-df$bolshev)/df$exact

stopifnot(all(df$exact-noe2_p(p,p,11,0)==0))

colnames(df) <- c("Exact Probability","Steck","Rel. Err. (Steck)","Bolshev","Rel. Err. (Bolshev)")
tbl <- xtable::xtable(df[-1,],digits=-5,auto=T,label='tbl:steck_bolshev',caption="test")
print(tbl,file="../tables/bolshev_steck.tex")