
library(stops)
data(Pendigits500)

library(ProjectionBasedClustering)
dis <- as.matrix(dist(Pendigits500[,1:16]))
?CCA
train.len <- 100
lambda0 <- 1
proj <- CCA(dis,Epochs=train.len,alpha0=0.5,lambda0=lambda0) #default lamda 3*max(apply(dis,2,sd))
proj2 <- CCA(dis,Epochs=100,alpha0=0.5)
PlotProjectedPoints(proj$ProjectedPoints,Pendigits500[,17]+1)
PlotProjectedPoints(proj2$ProjectedPoints,Pendigits500[,17]+1)

potency_curve = function(v0, vn, l) return(v0 * (vn/v0)^((0:(l - 1))/(l - 1)))
potency_curve(lambda0,0.01,train.len)

#need to change ProjectedPoints to conf and calculate a stress (there it is error but that one is always 0)

library(Rdimtools)
Xd <- Pendigits500[,1:16]
X <- as.matrix(Xd)
cca1 <- do.crca(X,lambda=130,alpha=0.5,maxiter=100)
cca1 <- do.crca(X,lambda=1000000000,alpha=0.5,maxiter=100)

library(smacof)
mds1 <- mds(dis)

par(mfrow=c(1,3))
plot(cca1$Y,pch=19,col=Pendigits500[,17]+1)
plot(proj$ProjectedPoints,pch=19,col=Pendigits500[,17]+1)
plot(mds1$conf,pch=19,col=Pendigits500[,17]+1)
par(mfrow=c(1,1))
plot(smacof::Procrustes(cca1$Y,mds1$conf))
plot(smacof::Procrustes(proj$ProjectedPoints,mds1$conf))

#so I tested for very high lambda: shoudl look like MDS because they are all with weight 1. Seems do.crca does this for a constant lambda. CCA also does it if we set epochs to 1. If there are more epochs it dynamically changes lambda (making it always smaller) so eventually no values get a positive lambda and the error is 0. The latter seems strange, so we need to work with epochs=1.  


#####
load("../data/Swissroll.rda")
dis <- as.matrix(dist(Swissroll[,1:3]))

library(ProjectionBasedClustering)
dis <- as.matrix(dist(Swissroll[,1:3]))
train.len <- 100
lambda0 <- 10
proj <- CCA(dis,Epochs=train.len,alpha0=0.5,lambda0=lambda0,PlotIt=TRUE) #default lamda 3*max(apply(dis,2,sd))
plot(proj$ProjectedPoints,col=Swissroll[,4])

potency_curve = function(v0, vn, l) return(v0 * (vn/v0)^((0:(l - 1))/(l - 1)))
potency_curve(lambda0,0.01,train.len)

#need to change ProjectedPoints to conf and calculate a stress (there it is error but that one is always 0)

library(Rdimtools)
Xd <- Swissroll[,1:3]
X <- as.matrix(Xd)
cca1 <- do.crca(X,lambda=10,alpha=0.5,maxiter=100)
cca1 <- do.crca(X,lambda=1000000000,alpha=0.5,maxiter=100)

library(smacof)
mds1 <- mds(dis)

par(mfrow=c(1,3))
plot(cca1$Y,pch=19,col=Swissroll[,4])
plot(proj$ProjectedPoints,pch=19,col=Swissroll[,4])
plot(mds1$conf,pch=19,col=Swissroll[,4])
par(mfrow=c(1,1))
plot(smacof::Procrustes(cca1$Y,mds1$conf))
plot(smacof::Procrustes(proj$ProjectedPoints,mds1$conf))



library(stops)
data(Pendigits500)
dis <- as.matrix(dist(Swissroll[,1:3]))
dis <- as.matrix(dist(Pendigits500[,1:16]))
X <- as.matrix(Pendigits500[,1:16])

epo <- 100
lambda0 <- 10
alpha0 <- 0.5
cca1 <- do.crca(X,lambda=lambda0,alpha=alpha0,maxiter=epo)
proj <- my.cca(dis,Epochs=epo,alpha0=alpha0,lambda0=lambda0) #default lamda 3*max(apply(dis,2,sd))
projo <- ProjectionBasedClustering::CCA(dis,Epochs=epo,alpha0=alpha0,lambda0=lambda0)
par(mfrow=c(1,3))
plot(cca1$Y,pch=19,col=Swissroll[,4])
plot(proj$ProjectedPoints,pch=19,col=Swissroll[,4])
plot(projo$ProjectedPoints,pch=19,col=Swissroll[,4])

###HAck CCA

DataOrDistances <- dis
Epochs <- 1
OutputDimension <- 2
method <- "euclidean"
alpha0 <- 0.5
PlotIt <- FALSE
Cls <- Swissroll[,4]
lambda0 <- 10

OrigCCA <- function (DataOrDistances, Epochs, OutputDimension = 2, method = "euclidean", 
    alpha0 = 0.5, lambda0, PlotIt = FALSE, Cls) 
{
#    if (missing(DataOrDistances)) 
#        stop("No DataOrDistances given")
#    DataOrDistances
#    if (!is.matrix(DataOrDistances)) 
#        stop("DataOrDistances has to be a matrix, maybe use as.matrix()")
#    if (missing(Epochs)) {
#        warning("scalar value for number of eppochs is missing. Setting epochs=20 which may not be prerable in order to continue the algorithm. There is no default setting for this parameter!")
#        Epochs = 20
#    }
#    else {
        epochs = Epochs
#    }
#    if (missing(lambda0)) 
#        lambda0 = NULL
#    if (isSymmetric(unname(DataOrDistances))) {
        Mdist = DataOrDistances
        AnzVar = ncol(DataOrDistances)
        AnzData = nrow(DataOrDistances)
#    }
#    else {
#        AnzVar = ncol(DataOrDistances)
#        AnzData = nrow(DataOrDistances)
#        Mdist = NULL
#    }
    squareform = function(X) {
        requireNamespace("pracma")
        if (is.matrix(X)) {
            if (isSymmetric(X)) {
                return(X[upper.tri(X)])
            }
            else {
                return(pracma::squareform(X))
            }
        }
        else {
            return(pracma::squareform(X))
        }
    }
    D = DataOrDistances
    P = OutputDimension
    potency_curve = function(v0, vn, l) return(v0 * (vn/v0)^((0:(l - 
        1))/(l - 1)))
    cca_error = function(P, Mdist, lambda) {
        noc = nrow(P)
        odim = ncol(P)
        noc_x_1 = rep(1, noc)
        odim_x_1 = rep(1, odim)
        error = 0
        for (i in 1:noc) {
            known = which(!is.nan(Mdist[, i]))
            if (length(known) > 0) {
                y = t(as.matrix(P[i, ]))
                Dy = P[known, ] - y[noc_x_1[known], ]
                dy = sqrt(rowSums(Dy^2))
                fy = exp(-dy/lambda)
                error = error + sum(((Mdist[known, i] - dy)^2) * 
                  fy)
            }
        }
        error = error/2
        return(error)
    }
    noc = nrow(D)
    dim = ncol(D)
    noc_x_1 = rep(1, noc)
    me = matrix(rep(0, dim), nrow = 1)
    st = matrix(rep(0, dim), nrow = 1)
    for (i in 1:dim) {
        me[i] = mean(D[which(is.finite(D[, i]) == TRUE), i])
        st[i] = sd(D[which(is.finite(D[, i]) == TRUE), i])
    }
    P  <-  P1  <-  ((2 * matrix(runif(noc * P), noc) - 1) * st[noc_x_1, 
        1:P]) + me[noc_x_1, 1:P]
    dummy = nrow(P)
    odim = ncol(P)
    odim_x_1 = matrix(1, odim, 1)
    train_len = epochs * noc
    sample_inds = ceiling(runif(train_len, 0, noc))
  #  if (is.null(Mdist)) 
  #      Mdist = as.matrix(dist(D, diag = TRUE, upper = TRUE))
  #  else {
  #      if (nrow(Mdist) == 1) 
  #          Mdist = squareform(Mdist)
  #      if (nrow(Mdist) != noc) 
  #          stop("Mutual distance matrix size and data set size do not match")
  #  }
    alpha = potency_curve(alpha0, alpha0/100, train_len)
  #  if (is.null(lambda0)) 
  #      lambda0 = max(st) * 3
    lambda = potency_curve(lambda0, 0.01, train_len)
    k = 0
    print(sprintf("iterating: %d / %d epochs", k, epochs))
    for (i in 1:train_len) {
        ind = sample_inds[i]
        dx = Mdist[, ind]
        known = which(!is.nan(dx))
        if (length(known) > 0) {
            y = t(as.matrix(P[ind, ]))
            Dy = P[known, ] - y[noc_x_1[known], ]
            dy = sqrt(rowSums(Dy^2))
            dy[which(dy == 0)] = 1
            fy = as.matrix(exp(-dy/lambda[i]) * (dx[known]/dy - 
                1))
            P[known, ] = P[known, ] + alpha[i] * fy[, odim_x_1] * 
                Dy
        }
        if (i%%noc == 0) {
            k = k + 1
            print(sprintf("iterating: %d / %d epochs", k, epochs))
        }
    }
    error = cca_error(P, Mdist, lambda[train_len])
    print(sprintf("%d iterations, error %f", epochs, error))
    unknown = which(sum(t(is.nan(D))) == dim)
    P[unknown, ] = NaN
    ProjectedPoints = P
   # if (PlotIt) {
   #     if (missing(Cls)) {
   #         AnzData = nrow(ProjectedPoints)
   #         Cls = rep(1, AnzData)
   #     }
   #     string = paste0("CCA with error ", round(error, 4), " and epochs ", 
   #         Epochs)
   #     PlotProjectedPoints(ProjectedPoints, Cls, main = string)
   # }
    return(list(ProjectedPoints = ProjectedPoints, Error = error))
}



my.cca <- function (DataOrDistances, Epochs, OutputDimension = 2, method = "euclidean", 
    alpha0 = 0.5, lambda0, PlotIt = FALSE, Cls) 
{
#    if (missing(DataOrDistances)) 
#        stop("No DataOrDistances given")
#    DataOrDistances
#    if (!is.matrix(DataOrDistances)) 
#        stop("DataOrDistances has to be a matrix, maybe use as.matrix()")
#    if (missing(Epochs)) {
#        warning("scalar value for number of eppochs is missing. Setting epochs=20 which may not be prerable in order to continue the algorithm. There is no default setting for this parameter!")
#        Epochs = 20
#    }
#    else {
        epochs = Epochs
#    }
#    if (missing(lambda0)) 
#        lambda0 = NULL
#    if (isSymmetric(unname(DataOrDistances))) {
        Mdist = DataOrDistances
        AnzVar = ncol(DataOrDistances)
        AnzData = nrow(DataOrDistances)
#    }
#    else {
#        AnzVar = ncol(DataOrDistances)
#        AnzData = nrow(DataOrDistances)
#        Mdist = NULL
#    }
    squareform = function(X) {
        requireNamespace("pracma")
        if (is.matrix(X)) {
            if (isSymmetric(X)) {
                return(X[upper.tri(X)])
            }
            else {
                return(pracma::squareform(X))
            }
        }
        else {
            return(pracma::squareform(X))
        }
    }
    D = DataOrDistances
    P = OutputDimension
    potency_curve = function(v0, vn, l) return(v0 * (vn/v0)^((0:(l - 
        1))/(l - 1)))
    cca_error = function(P, Mdist, lambda) {
        noc = nrow(P)
        odim = ncol(P)
        noc_x_1 = rep(1, noc)
        odim_x_1 = rep(1, odim)
        error = 0
        for (i in 1:noc) {
            known = which(!is.nan(Mdist[, i]))
            if (length(known) > 0) {
                y = t(as.matrix(P[i, ]))
                Dy = P[known, ] - y[noc_x_1[known], ]
                dy = sqrt(rowSums(Dy^2))
                fy = exp(-dy/lambda)
                error = error + sum(((Mdist[known, i] - dy)^2) * 
                  fy)
            }
        }
        error = error/2
        return(error)
    }
    noc = nrow(D)
    dim = ncol(D)
    noc_x_1 = rep(1, noc)
    me = matrix(rep(0, dim), nrow = 1)
    st = matrix(rep(0, dim), nrow = 1)
    for (i in 1:dim) {
        me[i] = mean(D[which(is.finite(D[, i]) == TRUE), i])
        st[i] = sd(D[which(is.finite(D[, i]) == TRUE), i])
    }
    P  <-  P1  <-  ((2 * matrix(runif(noc * P), noc) - 1) * st[noc_x_1, 
        1:P]) + me[noc_x_1, 1:P]
    dummy = nrow(P)
    odim = ncol(P)
    odim_x_1 = matrix(1, odim, 1)
    train_len = epochs * noc
    sample_inds = ceiling(runif(train_len, 0, noc))
  #  if (is.null(Mdist)) 
  #      Mdist = as.matrix(dist(D, diag = TRUE, upper = TRUE))
  #  else {
  #      if (nrow(Mdist) == 1) 
  #          Mdist = squareform(Mdist)
  #      if (nrow(Mdist) != noc) 
  #          stop("Mutual distance matrix size and data set size do not match")
  #  }
  #  alpha = potency_curve(alpha0, alpha0/100, train_len)
   alpha = rep(alpha0,train_len) 
  #  if (is.null(lambda0)) 
  #      lambda0 = max(st) * 3
                                        #  lambda = potency_curve(lambda0, 0.01, train_len)
    lambda <- rep(lambda0,train_len)  
    k = 0
    print(sprintf("iterating: %d / %d epochs", k, epochs))
    for (i in 1:train_len) {
        ind = sample_inds[i]
        dx = Mdist[, ind]
        known = which(!is.nan(dx))
        if (length(known) > 0) {
            y = t(as.matrix(P[ind, ]))
            Dy = P[known, ] - y[noc_x_1[known], ]
            dy = sqrt(rowSums(Dy^2))
            dy[which(dy == 0)] = 1
            fy = as.matrix(exp(-dy/lambda[i]) * (dx[known]/dy - 
                1))
            P[known, ] = P[known, ] + alpha[i] * fy[, odim_x_1] * 
                Dy
        }
        if (i%%noc == 0) {
            k = k + 1
            print(sprintf("iterating: %d / %d epochs", k, epochs))
        }
    }
    error = cca_error(P, Mdist, lambda[train_len])
    print(sprintf("%d iterations, error %f", epochs, error))
    unknown = which(sum(t(is.nan(D))) == dim)
    P[unknown, ] = NaN
    ProjectedPoints = P
   # if (PlotIt) {
   #     if (missing(Cls)) {
   #         AnzData = nrow(ProjectedPoints)
   #         Cls = rep(1, AnzData)
   #     }
   #     string = paste0("CCA with error ", round(error, 4), " and epochs ", 
   #         Epochs)
   #     PlotProjectedPoints(ProjectedPoints, Cls, main = string)
   # }
    return(list(ProjectedPoints = ProjectedPoints, Error = error))
}


