#Note - Compiling test.cc requires catch.hpp (https://github.com/catchorg/Catch2) to be in your include/ directory
path <- 'D:/Seafile/Publications'#insert the path of your libgmp here
Sys.setenv(	"PKG_LIBS" 		= paste0("-lgmp -lgmpxx -L'",path,"/libs/lib'"),
			"PKG_CPPFLAGS" 	= paste0("-mfma -I'",path,"/libs/include'"))
bin_path <- paste0(path,"/libs/bin;")
bin_path <- gsub("/","\\\\",bin_path)
if(!grepl(bin_path,Sys.getenv("path"))) Sys.setenv(path = paste0(bin_path, Sys.getenv("path")))
loadCpp <- (function(path) function(f){
  require(Rcpp)
  Rcpp::sourceCpp(paste0(path,'/',f))
})(dirname(rstudioapi::getSourceEditorContext()$path))

loadCpp('./cpp/test.cc')
loadCpp('./cpp/test_k.cc')

library(tikzDevice)
options(tikzLatexPackages 
        =c(getOption( "tikzLatexPackages" ),"\\usepackage{amssymb}\n"))

fix_tikz <- function(fname) {
  # remove all lines that invisibly mess up the bounding box
  lines <- readLines(con=fname)
  lines <- lines[-which(grepl("\\path\\[clip\\]*", lines,perl=F))]
  lines <- lines[-which(grepl("\\path\\[use as bounding box*", lines,perl=F))]
  writeLines(lines,con=fname)
}