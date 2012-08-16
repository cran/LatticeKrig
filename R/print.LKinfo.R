print.LKinfo<- function(x,...){
  LKinfo<- x
  L<-  LKinfo$nlevel
  cat( "Number of levels:", L, fill=TRUE)
  cat("grid sizes and number of basis functions", fill=TRUE)
   print( cbind( LKinfo$mx, LKinfo$my, diff(LKinfo$offset) ))
  cat("total number of basis functions:", LKinfo$m, fill=TRUE)
  cat("grid info", fill=TRUE)
  print( LKinfo$grid.info)
}
  
   
