library(compositions)
library(tensr)

#gPCA: main fonction, compute a generalized PCA for the coda loss by using a 
# alternating minimization process which uses all the parameters.
#
#gPCA return a list with:
# A:   the score matrix in the PCA space, 
# V:   the matrix basis of the new space,
# Y:   the matrix projections in the usual space
# dim: the number of dimensions.
# 
#parameters are explained in the following functions (SGDAlternatePRocess_mat 
# and SGD_mat)
gPCA = function(data,
                epochs = 50,
                ratioBatches = 1 / 3,
                method = 'Adam',
                lrate = 0.005,
                dimension = 2,
                cv.criterion = 10 ^ (-7),
                cycles = 2,
                max_iter = 600,
                loss = 'coda',
                Beta1 = 0.9,
                Beta2 = 0.999,
                eps = 10 ^ (-3),
                eps_adam = 10 ^ (-3)) {
  
  resp = SGDAlternateProcess_mat(
    data,
    epochs = epochs,
    ratioBatches = ratioBatches,
    method = 'Adam',
    lrate  = lrate,
    dimension = dimension,
    cv.criterion = cv.criterion,
    cycles   = cycles,
    max_iter = max_iter,
    loss  = loss,
    Beta1 = Beta1,
    Beta2 = Beta2,
    eps      = eps,
    eps_adam = eps_adam
  )
  
  temp = orthoLQGPCA(resp$A, resp$V)
  
  objpca   = list()
  objpca$A = temp$A
  objpca$V = temp$V
  objpca$Y = resp$A %*% resp$V
  objpca$dim = dimension
  
  return(objpca)
}

#SCGAlternateProcess_mat: compute an alternating minimization process by using
# two Stochastic Gradient Minimization (SGD_mat) up to the coda loss converges 
#
#data: n x d matrix, n is the length of the sample and d number of variables
#ratioBatches: between 0 and 1, it's the fraction of each batch.
#loss: it's the name of loss function, only 'coda' is available for now.!!!!!!!!
#cv.criterion and max_iter: two ways to stop the number of SGD_mat, max_iter is 
# the easiest to use.
#cycles: the number of times we initialize V, not need to be a high value.
#dimension: number of principal components to find.
#eps: the estimated value of zeros input data.
#
SGDAlternateProcess_mat = function(data,
                                   epochs = 50,
                                   ratioBatches = 1 / 3,
                                   method = 'Adam',
                                   lrate  = 0.005,
                                   dimension = 2,
                                   cv.criterion = 10 ^ (-5),
                                   cycles   = 2,
                                   max_iter = 600,
                                   loss  = 'coda',
                                   Beta1 = 0.9,
                                   Beta2 = 0.999,
                                   eps   = 10 ^ (-3),
                                   eps_adam = 10 ^ (-3)) {
  f  = F
  gC = F
  data_ = data
  data_nozero = data
  data_nozero[which(data == 0, arr.ind = T)] = eps
  data_nozero = CoDa(data_nozero)
  if (loss == 'coda') {
    f  = codaLoss
    gC = T
    g  = exp
    gradfa = gradCodaA_mat
    gradfv = gradCodaV_mat
  }
  
  if (gC) {
    data_ = unclass(exp(clr(data_nozero)))
  }
  
  d  = ncol(data_)
  m  = apply(data_, 2, mean)
  n  = nrow(data_)
  m2 = var(data_)
  prev_cost = 0
  count = 0
  cost_ = matrix(Inf, 1, max_iter)
  A = matrix(0, n, dimension)
  V = matrix(0, dimension, d)
  list_resp = list()
  
  for (cycle in 1:cycles) {
    V[, ] = t(svd(data)$v[, 1:dimension])
    cv    = 1
    count = 1
    neg   = 0
    while ((cv > cv.criterion) & (count < max_iter)) {
      batches = list(1, d)
      A = SGD_mat(
        init = A,
        f    = f,
        gradf  = gradfa,
        method = method,
        dataOption = list(x = data_, a = A, v = V),
        epochs  = epochs,
        batches = batches,
        gamma = gamma,
        lrate = lrate,
        cv.criterion = cv.criterion,
        Beta1 = Beta1,
        Beta2 = Beta2,
        eps_adam = eps_adam
      )
      batches = list(max(1, round(1 / ratioBatches)), n)
      V = SGD_mat(
        init = V,
        f = f,
        gradf  = gradfv,
        method = method,
        dataOption = list(
          x = data_,
          a = A,
          v = V,
          ind = F,
          cp  = F,
          var = F
        ),
        epochs  = epochs,
        batches = batches,
        gamma = gamma,
        lrate = lrate,
        cv.criterion = cv.criterion,
        Beta1 = Beta1,
        Beta2 = Beta2,
        eps_adam = eps_adam
      )
      cost = f(x = data_, a = A, v = V)
      # print(cost)
      cost_[count + 1] = cost
      neg = max(neg + sign(cost - prev_cost), 0)
      cv  = abs(prev_cost - cost) / abs(prev_cost) * (neg < 2)
      prev_cost = cost
      count = count + 1
      list_resp[[count]] = list(A = A, V = V)
    }
  }
  return(list_resp[[which.min(cost_)]])
}

#SGD_mat: compute a Stochastic Gradient Descent with mini-batches for a matrix
#
#init: initialization matrix, it's the first iteration of A or V.
#dataOption: parameters and matrix needed to compute the gradient.
#method: only 'Adam' is available up to now.
#f, gradf : loss function and his gradient.
#batches: a list with the number of batches and the length of sample.
#lrate: learning rate of gradient descent.
#Beta1, Beta2, eps_adam: parameters of Adam gradient descent, #eps_adam is used 
# to replace zeros.
#
SGD_mat = function(init,
                   dataOption = list(x = F, a = F, v = F),
                   method = 'Adam',
                   f,
                   gradf,
                   epochs = 100,
                   batches = list(3, nrow(init)),
                   lrate = 0.005,
                   Beta1 = 0.9,
                   Beta2 = 0.999,
                   eps_adam = 10 ^ (-3),
                   ...) {
  a   = dataOption$a
  v   = dataOption$v
  ind = dataOption$ind
  var = dataOption$var
  cp  = dataOption$cp
  x   = dataOption$x
  
  resp = init
  resp_prev = resp
  list_resp = list()
  mm = list()
  vv = list()
  pp = list()
  qq = list()
  ss = list()
  alpha = list()
  list_grad = list()
  yy = list()
  list_grad[[1]] = 0
  list_resp[[1]] = init
  
  for (e in 1:epochs) {
    Batches = buildBatches(n = batches[[2]], batches[[1]])
    for (b in 1:batches[[1]]) {
      bb = Batches[[b]]
      resp_prev = resp
      
      if (method == 'Adam') {
        if (e == 1) {
          mm[[b]] = matrix(0, nrow(init), ncol(init))
          vv[[b]] = mm[[b]]
        }
        tempGrad = clipGradNorm(gradf(
          x  = x,
          a  = a,
          bb = bb,
          v  = v,
          value = resp
        ))
        mm[[b]] = Beta1 * mm[[b]] + (1 - Beta1) * tempGrad
        vv[[b]] = Beta2 * vv[[b]] + (1 - Beta2) * tempGrad ^ 2
        mm[[b]] = mm[[b]] / (1 - Beta1 ^ e)
        vv[[b]] = vv[[b]] / (1 - Beta2 ^ e)
        resp = resp - lrate * mm[[b]] / (sqrt(vv[[b]]) + eps_adam)
      }
    }
  }
  return(resp)
}

#orthoLQGPCA: return A and V with V orthogonal, it's not needed for this coda 
# algorithm but it's free.
orthoLQGPCA = function(A, V) {
  decompo = lq(V)
  V_ = t(decompo$Q)
  A_ = A %*% decompo$L
  return(list(A = A_, V = V_))
}

#buildBatches: return a list with the id of the batches of a n-sample.
buildBatches = function(n, numberBatches) {
  resp = list()
  rand = sample(1:n)
  start = 1
  for (i in 1:numberBatches) {
    end = round(i * n / numberBatches)
    cut = rand[start:end]
    resp[[i]] = cut
    start = end + 1
  }
  return(resp)
}

#clipGradNorm: if the L2 norm of gradient is higher than max, it returns the 
# gradient with L2 norm equal to max
clipGradNorm = function(grad, max = 10) {
  norm = sqrt(sum(grad ^ 2))
  if (norm > max) {
    grad = grad / norm * max
  }
  return(grad)
}

#buildSimplexArms: build the artifical data set with (n*d) rows and d variables 
# which is the d simplex arms
buildSimplexArms = function(n = 100, d = 5) {
  resp = matrix(0, n * d, d)
  resp[1, ] = rep(1 / d, d)
  for (i in 1:(n - 1)) {
    resp[i+1,] = resp[i,] + c(1/(n-1) * (d-1)/d, -rep(1/(n-1)/d, d - 1))
  }
  for (j in 2:d) {
    resp[((j - 1) * n + 1):(j * n), ] = resp[1:n, ]
    temp = resp[((j - 1) * n + 1):(j * n), j]
    resp[((j - 1) * n + 1):(j * n), j] = resp[((j - 1) * n + 1):(j * n), 1]
    resp[((j - 1) * n + 1):(j * n), 1] = temp
  }
  resp[which(resp < 10 ^ (-5), arr.ind = T)] = 0
  return(resp)
}

#codaLoss: compute the coda loss with x geo-centered and y=av
codaLoss = function(x_, a, v, ...) {
  y_   = scale(a %*% v, scale = F)
  resp = exp(y_) - x_ * (y_ - log(x_)) - x_
  resp = sum(mean(resp))
  return(resp)
}

#gradCodaA_mat: compute the gradient of codaloss with v fixed.
#bb are the id of the batche.
gradCodaA_mat = function(x, a, v, value, bb = F, ...) {
  p = ncol(x)
  if (sum(bb == F)) {
    bb = (1:ncol(x))
  }
  a[, ] = value
  resp  = exp(a %*% v[, bb]) %*% t(v[, bb]) - x[, bb] %*% t(v[, bb])
  return(resp)
}

#gradCodaV_mat: compute the gradient of codaloss with a fixed.
#bb are the id of the batche.
gradCodaV_mat = function(x, a, v, value, bb = F, ...) {
  p = ncol(x)
  if (sum(bb == F)) {
    bb = (1:nrow(x))
  }
  v[, ] = value
  resp = t(a[bb, ]) %*% exp(a[bb, ] %*% v[, ]) - t(a[bb, ]) %*% x[bb, ]
  return(resp)
}

#CoDa: normalize the data in coda.
CoDa = function(data) {
  data = as.matrix(data)
  if (ncol(data) == 1) {
    data = t(data)
  }
  resp = data / apply(data, 1, sum)
  return(resp)
}



# .-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-.-
# Example: construct the the simplex arms or take a n x d data matrix.
n = 100
d = 5
data = buildSimplexArms(n = n, d = d)
#data=read.csv("name_file.csv")

#parameters used for the buildSimplexArms(n=100,d=5).
loss = 'coda'
dimension = 2
method = 'Adam'
cv.criterion = 10 ^ (-7)
Beta1 = 0.9
Beta2 = 0.999
eps   = 10 ^ (-3)
epochs = 50
eps_adam = 10 ^ (-3)
ratioBatches = 1 / 3
max_iter = 600
cycles = 2
lrate  = 0.005

#run the generalized PCA for the coda loss.
resp = gPCA(
  data,
  epochs = epochs,
  ratioBatches = ratioBatches,
  method = 'Adam',
  lrate = lrate,
  dimension = dimension,
  cv.criterion = cv.criterion,
  cycles = cycles,
  max_iter = max_iter,
  loss  = loss,
  Beta1 = Beta1,
  Beta2 = Beta2,
  eps   = eps,
  eps_adam = eps_adam
)

par(mfrow = c(1, 2)) 

#plot the score matrix of gPCA with one color by arm.
plot(
  resp$A,
  col = sort(rep(1:5, n)),
  xlab = 'PC 1',
  ylab = 'PC 2',
  xlim = c(min(resp$A), max(resp$A)),
  ylim = c(min(resp$A), max(resp$A)),
  main = 'gPCA'
)

data.pca <- prcomp(data, center = TRUE,scale. = TRUE)

plot(
  data.pca$x[,1:2],
  col = sort(rep(1:5, n)),
  xlab = 'PC 1',
  ylab = 'PC 2',
  xlim = c(min(data.pca$x[,1:2]), max(data.pca$x[,1:2])),
  ylim = c(min(data.pca$x[,1:2]), max(data.pca$x[,1:2])),
  main = 'PCA'
)


