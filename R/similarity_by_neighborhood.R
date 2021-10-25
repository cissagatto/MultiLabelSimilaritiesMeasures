
# df exemplo
#a1 <- c(1,0,1,1,1,0,0,0,1,1)
#a2 <- c(0,1,1,1,0,1,0,1,1,0)
#a3 <- c(1,1,1,0,0,0,1,1,1,0)
#a4 <- c(0,1,1,1,0,0,0,1,1,0)
#a5 <- c(0,1,1,1,0,0,0,1,1,0)
#c1 <- c(1,0,0,1,0,0,0,1,1,0)
#c2 <- c(1,0,0,1,0,0,0,1,1,0)
#c3 <- c(1,0,0,0,1,1,1,0,1,1)
#df <- data.frame(a1,a2,a3,a4,a5,c1,c2,c3)
#firstLabel <- 6


# df yeast
setwd("/Users/mauriferrandin/Dropbox/UFSC/ComCerri/part_correlation")
df = read.csv("yeast.csv")
firstLabel <- 104

k = 100

lastLabel <- ncol(df)
lastAtt <- firstLabel -1

classes <- seq(firstLabel,lastLabel)

p <- data.frame()
for (i in 1:length(classes)){
  pi = colSums(df[which(df[,classes[i]] == 1),1:lastAtt])/sum(df[,classes[i]])
  p <- rbind(p,pi)
}
colnames(p) <- colnames(df)[1:lastAtt]


findknn <- function(df, pi, k){
  
  dists <- sqrt(rowSums(sweep(as.matrix(df[,1:lastAtt]),2, as.numeric(pi),"-")^2))
  
  knn <- data.frame(rownames(df),dists);
  colnames(knn) <- c("nearest", "dist")
  
  knn <- knn[order(knn$dist),]
  knn <- knn[1:k,]
  findknn <- knn
}



sim = matrix(ncol = length(classes), nrow = length(classes))
rownames(sim) = colnames(df)[firstLabel:lastLabel]
colnames(sim) = colnames(df)[firstLabel:lastLabel]
for (i in 1:length(classes)){
  knndist <- findknn(df, p[i,], k)
  knnids <- knndist[1]
  dfknn <- df[as.list(knnids)[[1]],]
  for (j in 1:length(classes)){
    
    pj = colSums(dfknn[which(dfknn[,classes[j]] == 1),1:lastAtt])/k
    
    sim[i,j] = sqrt(sum((pi - pj)^2))
    #print(pi - pj)
  }
}

################################################################################
# any errors, please, contact me: elainececiliagatto@gmail.com                 #
################################################################################