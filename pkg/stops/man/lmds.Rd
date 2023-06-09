% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/lmds.R
\name{lmds}
\alias{lmds}
\title{An MDS version for local MDS (Chen & Buja 2006)}
\usage{
lmds(delta, init = NULL, ndim = 3, k = 2, tau = 1, itmax = 5000, verbose = 0)
}
\arguments{
\item{delta}{dissimilarity or distance matrix}

\item{init}{initial configuration. If NULL a classical scaling solution is used.}

\item{ndim}{the dimension of the configuration}

\item{k}{the k neighbourhood parameter}

\item{tau}{the penalty parameter (suggested to be in [0,1])}

\item{itmax}{number of optimizing iterations, defaults to 5000.}

\item{verbose}{prints progress if > 4.}
}
\description{
An MDS version for local MDS (Chen & Buja 2006)
}
\details{
Note that k and tau are not independent. It is possible for normalized stress to become negative if the tau and k combination is so that the absolute repulsion for the found configuration dominates the local stress substantially less than the repulsion term does for the solution of D(X)=Delta, so that the local stress difference between the found solution and perfect solution is nullified. This can typically be avoided if tau is between 0 and 1. If not, set k and or tau to a smaller value.
}
\examples{
dis<-smacof::kinshipdelta
res<- lmds(as.matrix(dis),k=2,tau=0.1)
res
summary(res)
plot(res)

}
\author{
Lisha Chen & Thomas Rusch
}
