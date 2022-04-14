\name{check_build_training_data}
\title{Heterogenous Linear instrumental variables using BART}
\alias{check_build_training_data}

\description{ Checks inputs \code{Z}, \code{T}, and
  \code{X} (if desired) and processes them to be in the expected
  format for the C++ algorithm. Returns \code{X} data. }

\usage{
check_build_training_data(Z,T,X,n,method)
}

\arguments{
   \item{Z}{ Training data for the instruments. }
   \item{T}{ Training data for the exposure. }
   \item{X}{ Training data for the other covariates. }
   \item{n}{ \code{length(Y)}}
   \item{method}{ String giving the name of the function. }
}

\author{
Charles Spanbauer: \email{spanb008@umn.edu},\cr
}
