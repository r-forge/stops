% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stop_isomap.R
\name{stop_isomap1}
\alias{stop_isomap1}
\title{STOPS version of isomap to optimize over integer k.}
\usage{
stop_isomap1(
  dis,
  theta = 3,
  type = "ratio",
  weightmat = NULL,
  ndim = 2,
  init = NULL,
  stressweight = 1,
  structures = c("cclusteredness", "clinearity", "cdependence", "cmanifoldness",
    "cassociation", "cnonmonotonicity", "cfunctionality", "ccomplexity", "cfaithfulness",
    "cregularity", "chierarchy", "cconvexity", "cstriatedness", "coutlying",
    "cskinniness", "csparsity", "cstringiness", "cclumpiness", "cinequality"),
  strucweight = rep(1/length(structures), length(structures)),
  strucpars,
  verbose = 0,
  stoptype = c("additive", "multiplicative"),
  itmaxi = NULL,
  registry = struc_reg
)
}
\arguments{
\item{dis}{numeric matrix or dist object of a matrix of proximities}

\item{theta}{the number of shortest dissimilarities retained for a point (nearest neighbours), the isomap parameter. Must be a numeric scalar. Defaults to 3.}

\item{type}{MDS type. Is "ratio".}

\item{weightmat}{(optional) a matrix of nonnegative weights}

\item{ndim}{number of dimensions of the target space}

\item{init}{(optional) initial configuration}

\item{stressweight}{weight to be used for the fit measure; defaults to 1}

\item{structures}{which structuredness indices to be included in the loss}

\item{strucweight}{weight to be used for the structuredness indices; ; defaults to 1/#number of structures}

\item{strucpars}{the parameters for the structuredness indices}

\item{verbose}{numeric value hat prints information on the fitting process; >2 is extremely verbose}

\item{stoptype}{How to construct the target function for the multi objective optimization? Either 'additive' (default) or 'multiplicative'}

\item{itmaxi}{placeholder for compatibility in stops call; not used}

\item{registry}{registry object with c-structuredness indices.}
}
\value{
A list with the components
   \itemize{
        \item{stress:} Not really stress but 1-GOF[2] where GOF is the second element returned from smacofx::cmdscale (the sum of the first ndim eigenvalues divided by the sum of all absolute eigenvalues).
        \item{stress.m:} default normalized stress (sqrt explicitly normalized stress; really the stress this time)
        \item{stoploss:} the weighted loss value
        \item{indices:} the values of the structuredness indices
        \item{parameters:} the parameters used for fitting 
        \item{fit:} the returned object of the fitting procedure
         \item{stopobj:} the stopobj object
}
}
\description{
Free parameter is k.
}
\details{
Currently this version is a bit less flexible than the vegan one, as the only allowed parameter for isomap is the theta (k in isomap, no epsilon) and the shortest path is always estimated with argument "shortest". Also note that fragmentedOK is always set to TRUE which means that for theta that is too small only the largest conected group will be analyzed. If that's not wanted just set the theta higher.
}
\keyword{multivariate}
