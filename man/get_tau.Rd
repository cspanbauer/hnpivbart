\name{get_tau}
\title{Heterogenous Linear instrumental variables using BART}
\alias{get_tau}

\description{ Get terminal node prior variances }

\usage{
get_tau(sigmaf1,sigmaf2,k,m1,m2,rT,rY)
}

\arguments{
   \item{sigmaf1}{ User defined (if given) value for setting prior
   variances in stage 1. }
   \item{sigmaf2}{ User defined (if given) value for setting prior
   variances in stage 2. }
   \item{k}{ Regularization hyperparameter. Higher values pull the
   predictions toward the prior mean (usually the sample mean) }
   \item{m1}{ Number of trees in stage 1. }
   \item{m2}{ Number of trees in stage 2. }
   \item{rT}{ Range of exposure \code{T}. }
   \item{rY}{ Range of outcome \code{Y}. }
}

\author{
Charles Spanbauer: \email{spanb008@umn.edu},\cr
}
