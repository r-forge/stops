% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structurednessindices.R
\name{c_faithfulness}
\alias{c_faithfulness}
\title{c-faithfulness 
calculates the c-faithfulness based on the index by Chen and Buja 2013 (M_adj) with equal input neigbourhoods}
\usage{
c_faithfulness(confs, obsdiss, k = 3, ...)
}
\arguments{
\item{confs}{a numeric matrix or a dist object}

\item{obsdiss}{a symmetric numeric matrix or a dist object. Must be supplied.}

\item{k}{the number of nearest neighbours to be looked at}

\item{...}{additional arguments passed to dist()}
}
\value{
a numeric value; faithfulness
}
\description{
c-faithfulness 
calculates the c-faithfulness based on the index by Chen and Buja 2013 (M_adj) with equal input neigbourhoods
}
\examples{
delts<-smacof::kinshipdelta
dis<-smacofSym(delts)$confdist
c_faithfulness(dis,obsdiss=delts,k=3)
}
