# CoDa PCA
# R functions associated to "Representation Learning of Compositional Data", NIPS18
# Methods, Algorithms and Python implementation developed by R Nock, C S Ong and K Sun, Data61, CSIRO, Australia
# R version coded by J Rouar, verified by M Avalos and B Xu, SISTM, INRIA and INSERM U1219, France  
# 20/10/2018

source("CoDa-PCA.R")

# .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
# Demo: CoDa-PCA on simplex arms or take a n x d data matrix.
n    = 100
d    = 10
data = buildSimplexArms(n = n, d = d)

# Parameters for CoDa-PCA
loss         = 'coda'
dimension    = 2
method       = 'Adam'
cv.criterion = 10 ^ (-7)
Beta1        = 0.9
Beta2        = 0.999
eps          = 10 ^ (-3)
epochs       = 50
eps_adam     = 10 ^ (-3)
ratioBatches = 1 / 3
max_iter     = 600
cycles       = 2
lrate        = 0.005

# CoDa-PCA running
result.codapca = gPCA(
  data,
  epochs       = epochs,
  ratioBatches = ratioBatches,
  method       = 'Adam',
  lrate        = lrate,
  dimension    = dimension,
  cv.criterion = cv.criterion,
  cycles       = cycles,
  max_iter     = max_iter,
  loss         = loss,
  Beta1        = Beta1,
  Beta2        = Beta2,
  eps          = eps,
  eps_adam     = eps_adam
)

result.clrpca  = prcomp(clr_transform(data), center = TRUE, scale = TRUE)
result.pca     = prcomp(data, center = TRUE, scale = TRUE)

# Plotting .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.
par(mfrow = c(1, 3)) 

clist = sort(rep(1:d, n))

plot(result.codapca$A[,1], result.codapca$A[,2], col=clist, xlab = 'PC1', ylab = 'PC2',main = '10-arm simplex by CoDa-PCA')
plot(result.clrpca$x[,1],  result.clrpca$x[,2],  col=clist, xlab = 'PC1', ylab = 'PC2',main = '10-arm simplex by CLR-PCA')
plot(result.pca$x[,1],     result.pca$x[,2],     col=clist, xlab = 'PC1', ylab = 'PC2',main = '10-arm simplex by PCA')



