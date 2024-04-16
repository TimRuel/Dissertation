cl <- makeCluster(6, type="SOCK")
registerDoSNOW(cl)
m <- matrix(rnorm(90000), 30000, 3)
foreach(i=1:nrow(m), .combine=rbind) %dopar% (m[i,] / mean(m[i,]))
stopCluster(cl)




library(doSNOW)
library(itertools)
cl <- makeCluster(6, type="SOCK")
registerDoSNOW(cl)
m <- matrix(rnorm(90000), 30000, 3)

foreach(s=isplitRows(m, chunks=getDoParWorkers()), .combine='rbind') %dopar% {
  t(apply(s, 1, function(x) x / mean(x)))
}