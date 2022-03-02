data_process_death <- function (x) {
  
  n <- 5
  
  x[x < 0] <- NA
  
  x <- x[complete.cases(x), ]
  
  x <- x[x[,2] < quantile(x[,2],prob=1-n/100),]
}
