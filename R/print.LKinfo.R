print.LKinfo <- function(x, ...) {
    LKinfo <- x
    L <- LKinfo$nlevel
    cat("Number of levels:", L, fill = TRUE)
    cat(" ", fill = TRUE)
    cat("alpha:", unlist(x$alpha), fill = TRUE)
    if (!is.null(x$nu)) {
        cat("based on nu: ", x$nu, fill = TRUE)
    }
    cat("a.wght", unlist(x$a.wght), fill = TRUE)
       cat(" ", fill = TRUE)
    cat("Ranges of spatial domain (excluding the buffer regions)", 
        fill = TRUE)
       temp <- unlist(LKinfo$grid.info)
       names(temp) <- names(LKinfo$grid.info)
       print(temp)
    cat("Distance type: ", LKinfo$distance.type, fill = TRUE)
        cat(" ", fill = TRUE)
   
# Details on basis functions at ecah level
        temp <- cbind(LKinfo$mx - 2 * LKinfo$NC.buffer.x, LKinfo$my - 
            2 * LKinfo$NC.buffer.y, LKinfo$m.domain)
        temp<- cbind( temp, LKinfo$delta)
        dimnames(temp) <- list(paste("level", 1:LKinfo$nlevel), 
            c("mx", "my", "total", "grid spacing"))
        cat("Total number of basis functions is ", LKinfo$m,"using buffer of", LKinfo$NC.buffer, "points.",
            fill = TRUE)
        if( LKinfo$normalize){
        cat("Basis functions will be normalized to approximate stationarity", fill=TRUE)
        }
          cat(" ", fill = TRUE)
        cat("Total number of basis functions _within_ domain: ", 
                 sum(temp[, 3]), fill = TRUE)
           cat(" ", fill = TRUE)
        cat("Grid sizes and number of basis functions _within_ spatial domain", 
            fill = TRUE)
        print(temp)   
}



