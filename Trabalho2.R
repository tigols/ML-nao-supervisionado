# Trabalho 2 ML nao supervisionado

library(mclust)
library(dbscan)
library(ggplot2)

rm(list = ls())
set.seed(21)

load("snp50k.rda")

plot(alleleA50k[1,], alleleB50k[1,], xlab = "Alelo A", ylab = "Alelo B")


########################################### FUNCOES ###########################################


ObjectiveKmeansTracew <- function(veccenters, K, X, A){
  centers <- matrix(veccenters, nrow = K)
  W <- crossprod(X - A%*%centers)
  return(log(det(W))) # optim waits for numeric output
}
GradientKmeansTracew <- function(veccenters, K, X, A){
  centers <- matrix(veccenters, nrow = K)
  dcenters <- 0*centers # Mesmas dimensoes, matriz de zeros
  W <- crossprod(X - A%*%centers)
  W.inv <- solve(W) # Inversa, mas é pequena
  for(j in 1:K){
    if(sum(A[,j]) == 0) next() # Gradient set to zero (arbitrarily) if cluster is empty
    temp <- sum(A[,j])*centers[j, ] - apply(x[as.logical(A[,j]), ], 1, sum) # elements of j-th iteration
    for(w in 1:ncol(centers)){
      dO.dc <- 0*W # Mesmas dimensoes, matriz de zeros
      dO.dc[w,] <- temp
      dO.dc[,w] <- temp
      dO.dc[w,w] <- 2*temp[w] # changes main cell
      dcenters[j,w] <- sum(diag(W.inv%*%dO.dc))
    }
  }
  return(as.numeric(dcenters)) # optim waits for vector output
}
sqrtM.inv <- function(W){
  Eig <- eigen(W, symmetric = TRUE)
  return(Eig$vectors %*% diag(1/sqrt(Eig$values)) %*% t(Eig$vectors))
}
# A função que vou declarar, kmeansTracew, copia comportamento de kmeans()
# se centers é um inteiro (ou quase), então é o número de clusteres
# caso contrário é uma matriz de coordenadas pros centros
#
# A minha implementação não está ótima principalmente porque eu corro
# o risco de ter clusteres vazios, mas deve funcionar se K for pequeno
kmeansTracew <- function(X, centers, maxit = 10, tol = 0.01){
  X <- as.matrix(X)
  if(!is.matrix(centers)){
    if(centers > 1 & centers < nrow(X)) centers <- X[sample.int(nrow(X), size = round(centers)), ]
  } else if(is.matrix(centers)) {
    if(ncol(centers) != ncol(X)) stop("centers need to be either an integer or a K times p matrix of centers")
  } else {
    stop("centers need to be either an integer or a K times p matrix of centers")
  }
  K <- nrow(centers); n <- nrow(X); p <- ncol(X); it <- 1; W.05 <- diag(p)
  div <- function(x, y) (x-y)^2
  repeat{
    A <- matrix(0, ncol = K, nrow = n); D <- matrix(0, ncol = K, nrow = n); newcenters <- centers
    for(w in 1:p){
      D <- D + outer((X%*%W.05)[,w], (newcenters%*%W.05)[,w], FUN = "div")
    }
    closestCentroid <- apply(D, 1, which.min)
    for(w in seq_along(closestCentroid)) A[w, closestCentroid[w]] <- 1
    veccenters <- as.numeric(newcenters)
    W <- crossprod(X - A%*%centers); W.05 <- sqrtM.inv(W)
    newveccenters <- optim(veccenters, fn = ObjectiveKmeansTracew, gr = GradientKmeansTracew, K = K, X = X, A = A)$par
    centers <- matrix(newveccenters, nrow = K)
    if(max(abs(centers - newcenters))/max(abs(centers)) <= tol | it >= maxit) break()
    it <- it + 1
  }
  return(list(labels = closestCentroid, centers = centers))
}


################################################################################

snp = sample(nrow(alleleA50k), 3, replace = FALSE)



par(ask =TRUE)
for( i in snp[3]){
  message("Estamos no SNP ",i)
  titulo = paste0("SNP ",i)
  MS_tmp = NULL
  
  M = log2(alleleA50k[i,]/alleleB50k[i,])
  S = (log2(alleleA50k[i,]) + log2(alleleB50k[i,]))/2
  MS = cbind(S,M)
  plot(MS,main = titulo)
  
  ########   Kmeans   #############
  cluster = hclust(dist(MS), method = "average")
  
  N = nrow(MS)
  h = seq(min(cluster$height), max(cluster$height), by = 0.01)
  
  totalVar = numeric(length(h))
  
  K =  numeric(length(h))
  
  for(i in seq_along(h)){
    groups <- factor(cutree(cluster, h = h[i]))
    K[i] <- length(levels(groups))
    B <- summary(manova(MS ~ groups))$SS$groups
    W <- summary(manova(MS ~ groups))$SS$Residuals
    totalVar[i] <- det(B)/(det(B+W))
  }
  
  plot(K, totalVar, xlim = c(0,100), main = "Elbolw")
  axis(1, at = 1:25, labels = 1:25)
  abline(v = 3, col = "red", lty = 2)
  
  k.medias = kmeans(MS, 21)
  
  MS_tmp = cbind(MS, k.medias$cluster)
  
  plot(MS, col = k.medias$cluster, main = titulo)
  
  
  ##########    Kmeans Optimization    #############
  
  k.opt = kmeansTracew(MS, centers = 3, maxit = 10, tol = 0.01)
  k = cbind(MS,k.opt$labels)
  
  k = as.data.frame(k)
  colnames(k) = c("S","M","cluster")
  k$cluster = as.factor(k$cluster)
  
  
  print(ggplot(k, aes(x = S, y = M , color = cluster)) + geom_point() +
          labs(title = paste0("Kmeans Optimization - ",titulo))) 
  
  # o modelo nao ajustou bem (talvez erro nos chutes iniciais)
  

  
}
par(ask = FALSE)

#########     DBSCAN    #################
# snp1
M_tmp = log2(alleleA50k[snp[1],]/alleleB50k[snp[1],])
S_tmp = (log2(alleleA50k[snp[1],]) + log2(alleleB50k[snp[1],]))/2
MS_tmp = cbind(S_tmp,M_tmp)

dens = densityMclust(MS_tmp)
plot(dens, what = "density", type = "persp", main = titulo)

R <- 0.3 # por serem elipses em sua maioria achatadas estimei por olho o menor "raio" da elipse

model <- dbscan(MS_tmp, eps = R, minPts = 5)
hullplot(MS_tmp, model,main = titulo)

############
# snp2
M_tmp = log2(alleleA50k[snp[2],]/alleleB50k[snp[2],])
S_tmp = (log2(alleleA50k[snp[2],]) + log2(alleleB50k[snp[2],]))/2
MS_tmp = cbind(S_tmp,M_tmp)

dens = densityMclust(MS_tmp)
plot(dens, what = "density", type = "persp", main = titulo)

R <- 0.3 # por serem elipses em sua maioria achatadas estimei por olho o menor "raio" da elipse

model <- dbscan(MS_tmp, eps = R, minPts = 50)
hullplot(MS_tmp, model,main = "SNP 12089")

#########
# snp3
M_tmp = log2(alleleA50k[snp[3],]/alleleB50k[snp[3],])
S_tmp = (log2(alleleA50k[snp[3],]) + log2(alleleB50k[snp[3],]))/2
MS_tmp = cbind(S_tmp,M_tmp)

plot(MS_tmp)

dens = densityMclust(MS_tmp)
plot(dens, what = "density", type = "persp", main = titulo)

R <- 0.3 # por serem elipses em sua maioria achatadas estimei por olho o menor "raio" da elipse

model <- dbscan(MS_tmp, eps = R, minPts = 25)
hullplot(MS_tmp, model,main = titulo)





# 
# 
# M_tmp = log2(alleleA50k[1,]/alleleB50k[1,])
# S_tmp = (log2(alleleA50k[1,]) + log2(alleleB50k[1,]))/2
# MS_tmp = cbind(S_tmp,M_tmp)
# plot(MS_tmp, xlab = "S", ylab = "M")
# 
# smoothScatter(MS_tmp)
# 
# 
# 
# 
# # KMEANS 
# 
# # Determinando numero de clusters
# 
# cluster = hclust(dist(MS_tmp), method = "average")
# 
# plot(cluster, labels = FALSE)
# 
# N = nrow(MS_tmp)
# h = seq(min(cluster$height), max(cluster$height), by = 0.01)
# 
# totalVar = numeric(length(h))
# 
# K =  numeric(length(h))
# 
# for(i in seq_along(h)){
#   groups <- factor(cutree(cluster, h = h[i]))
#   K[i] <- length(levels(groups))
#   B <- summary(manova(MS_tmp ~ groups))$SS$groups
#   W <- summary(manova(MS_tmp ~ groups))$SS$Residuals
#   totalVar[i] <- det(B)/(det(B+W))
# }
# 
# plot(K, totalVar, xlim = c(0,25))
# axis(1, at = 1:25, labels = 1:25)
# abline(v = 3, col = "red", lty = 2)
# # numero ideal igual a 3
# 
# 
# for(i in 1:50){
#   k.medias = kmeans(MS_tmp, 3)
#   
#   MS_tmp = cbind(MS_tmp, k.medias$cluster)
#   
#   plot(MS_tmp, col = k.medias$cluster)
#   
# }
# # Kmeans ta com problema de otimizacao dos chutes iniciais dos centroides
# 
# 
# 
# for(i in 1:50){
#   kmeans = kmeansTracew(MS_tmp, centers = n_cluster, maxit = 10, tol = 0.01)
#   k = cbind(MS_tmp,kmeans$labels)
#   
#   k = as.data.frame(k)
#   colnames(k) = c("S","M","cluster")
#   k$cluster = as.factor(k$cluster)
#   
#   
#   print(ggplot(k, aes(x = S, y = M , color = cluster)) + geom_point() +
#           labs(title = "Kmeans Optimization")) 
# }
# 




