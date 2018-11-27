# NonParametric CoDa PCA
#
# R functions associated to "Representation Learning of Compositional Data", NIPS18
# Methods, Algorithms and Python implementation developed by R Nock, C S Ong and K Sun, Data61, CSIRO, Australia
# R version coded by B Xu, SISTM, INRIA and INSERM U1219, France, verified by M Avalos   
# 27/10/2018
#
# NOTE: to get better performance, use "control" option in optim.

#clr_transform:  centered log-ratio transformation of each row
clr_transform = function(X, eps=1e-8 ){
  X = log( X + eps )
  X =  X - rowMeans(X, na.rm = TRUE, dims = 1)
  return(X)
}

# the loss function to be minimized
loss_func <- function(x, n_components, N, dim, checkX) {
  U = matrix(x[1:((n_components+1)*dim)], nrow = n_components+1, byrow = TRUE)
  U = U - rowMeans(U, na.rm = TRUE, dims = 1)
  
  B = matrix(x[((n_components+1)*dim)+1:(N*dim)], nrow = N, ncol = n_components, byrow = TRUE )
  B = cbind(B, rep(1,N))
  
  Y = B %*% U
  
  loss = sum(exp( Y ) - checkX * Y)
  return(loss)
}

# gradient of the loss function
grad_func <- function(x, n_components, N, dim, checkX) {
  U = matrix(x[1:((n_components+1)*dim)], nrow = n_components+1, byrow = TRUE)
  U = U - rowMeans(U, na.rm = TRUE, dims = 1)
  
  B = matrix(x[((n_components+1)*dim)+1:(N*dim)], nrow = N, ncol = n_components, byrow = TRUE )
  B = cbind(B, rep(1,N))
  
  Y = B %*% U
  
  diff  = exp( Y ) - checkX
  
  gradB = diff %*% t(U)
  gradB = gradB[,1:dim(gradB)[2]-1]
  
  gradU = t(B) %*% diff
  gradU = gradU - rowMeans(gradU, na.rm = TRUE, dims = 1)
  
  grad = c(as.vector(t(gradU)), as.vector(t(gradB)))

  return(grad)
}


SCodaPCA = function(x0, loss_func, grad_func,n_components,N,dim,checkX) {
  result <- optim(x0, fn = loss_func, gr = grad_func, method="L-BFGS-B",
                  n_components=n_components, N=N, dim=dim, checkX=checkX)
  
  pcadata     = list()
  pcadata$U   = matrix(result$par[1:((n_components+1)*dim)], nrow = n_components+1, byrow = TRUE)
  pcadata$pca = matrix(result$par[((n_components+1)*dim)+1:(N*dim)], nrow = N, ncol = n_components, byrow = TRUE )

  return(pcadata)
}

options(stringsAsFactors=F)
X     = read.table('atlas.csv',        header = TRUE, sep = ',')
group = read.table('atlasfactors.csv', header = TRUE, sep = ',')
DNA_extraction_method = group$DNA_extraction_method
DNA_extraction_method[is.na(DNA_extraction_method)]<-"NA"
clist = match(DNA_extraction_method, unique(DNA_extraction_method))

X   = as.matrix(X)
N   = dim(X)[1]
dim = dim(X)[2]

n_components = 5;

checkX = exp( clr_transform( X ) )
x0 = rnorm( (n_components+1)*dim + N*n_components )

result.scoda   = SCodaPCA(x0, loss_func, grad_func,n_components,N,dim,checkX)
result.clrpca  = prcomp(clr_transform(X), center = TRUE, scale = TRUE)
result.pca     = prcomp(X, center = TRUE, scale = TRUE)

par(mfrow = c(1, 3))
plot(result.scoda$pca, col=clist, xlab = 'PC1', ylab = 'PC2',main = 'NonParametricCodaPCA')
plot(result.clrpca$x,  col=clist, xlab = 'PC1', ylab = 'PC2',main = 'Atlas by CLR-PCA')
plot(result.pca$x,     col=clist, xlab = 'PC1', ylab = 'PC2',main = 'Atlas by PCA')