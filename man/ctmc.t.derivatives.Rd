% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{ctmc.t.derivatives}
\alias{ctmc.t.derivatives}
\title{CTMC transition probability derivatives with respect to time interval length}
\usage{
ctmc.t.derivatives(t, subst_mod, max_order)
}
\arguments{
\item{t}{A nonnegative numeric scalar representing the CTMC time interval 
length.}

\item{subst_mod}{An S3 object of class \code{"substitution.model"}.}

\item{max_order}{An integer specifying the maximum derivative order of 
interest.}
}
\value{
A three-dimensional array that stores the transition probability 
  derivative matrices, where the third dimension indexes the different 
  derivative orders.  Specifically, if \code{outp} denotes the output array,
  then \code{outp[, , k]} represents the transition probability derivative 
  matrix of order \code{k - 1} for all integers \code{k} between \code{1} 
  and \code{max_order + 1}.
}
\description{
Computes the CTMC transition probability derivative matrices with respect to
the time interval length for all orders less than or equal to the 
user-specified maximum order of interest.
}
\details{
The zeroth-order transition probability derivative matrix is defined to be 
the transition probability matrix.
}
\references{
Kenney T and Gu H (2012) \dQuote{Hessian Calculation for 
  Phylogenetic Likelihood based on the Pruning Algorithm and its 
  Applications}, \emph{Statistical Applications in Genetics and Molecular 
  Biology}, 11(4).
  
  OUR PAPER!!
}
\seealso{
\code{\link{ctmc.Q.derivatives}}, \code{\link{phylo.t.derivatives}}
}

