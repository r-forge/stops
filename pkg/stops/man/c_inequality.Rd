% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/structurednessindices.R
\name{c_inequality}
\alias{c_inequality}
\title{c-inequality
Calculates c-inequality (as in an economic measure of inequality) as Pearsons coefficient of variation of the fitted distance matrix. This can help with avoiding degenerate solutions.   
This is one of few c-structuredness indices not between 0 and 1, but 0 and infinity.}
\usage{
c_inequality(confs, ...)
}
\arguments{
\item{confs}{a numeric matrix or data frame}

\item{...}{additional arguments (don't do anything)}
}
\value{
a numeric value; inequality (Pearsons coefficient of variation of the fitted distance matrix)
}
\description{
c-inequality
Calculates c-inequality (as in an economic measure of inequality) as Pearsons coefficient of variation of the fitted distance matrix. This can help with avoiding degenerate solutions.   
This is one of few c-structuredness indices not between 0 and 1, but 0 and infinity.
}
\examples{
x<-1:10
y<-2+3*x+rnorm(10)
z<- sin(y-x)
confs<-cbind(z,y,x)
c_inequality(confs)
}
