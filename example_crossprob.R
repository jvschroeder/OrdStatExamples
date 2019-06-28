#setwd(dirname(rstudioapi::getSourceEditorContext()$path))

fix_tikz <- function(fname) {
  # remove all lines that invisibly mess up the bounding box
  lines <- readLines(con=fname)
  lines <- lines[-which(grepl("\\path\\[clip\\]*", lines,perl=F))]
  lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines,perl=F))]
  writeLines(lines,con=fname)
}

tikzDevice::tikz(file = "../graphics/relative_numerical_error.tex", width = 5, height = 5)
result <- jsonlite::fromJSON("./test.json")
result <- matrix(unlist(result, use.names = F), ncol=4, byrow=T)
result <- result[1:200,]
df <- data.frame(result[,1],result[,2],result[,3],result[,4])
names(df) <- c("n","Numerical Integration","Poisson O($n^2$)","Poisson O($n^2\\log(n)$)")
res <- tidyr::gather(df,Algorithm,value,-n)
res$Algorithm <- factor(res$Algorithm)
lvls <- levels(res$Algorithm)
lvls[1:2] <- lvls[2:1]
res$Algorithm <- factor(res$Algorithm,levels=lvls)
plot <- ggplot2::ggplot(res) + ggplot2::geom_line(ggplot2::aes(x=n,y=value,linetype=Algorithm,color=Algorithm),size=1) + ggplot2::ylab("Magnitude of the relative error") +
  ggplot2::xlab("$n$") + ggplot2::coord_equal() + 
  ggplot2::theme(axis.line = ggplot2::element_line(colour = "grey"), panel.background = ggplot2::element_blank()) +
  ggplot2::scale_color_grey()
print(plot)
dev.off()
fix_tikz("../graphics/relative_numerical_error.tex")