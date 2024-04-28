\name{check_build_test_data}
\title{Heterogenous Linear instrumental variables using BART}
\alias{check_build_test_data}

\description{ Checks inputs \code{Z.test}, \code{T.test}, and
  \code{X.test} (if desired) and processes them to be in the expected
  format for the C++ algorithm. Returns a list of the processed
  test data}

\usage{
check_build_test_data(Z.test,T.test,X.test,method)
}

\arguments{
   \item{Z.test}{ Test data for the instruments. }
   \item{T.test}{ Test data for the exposure. }
   \item{X.test}{ Test data for the other covariates. }
   \item{method}{ String giving the name of the function. }
}

\author{
Charles Spanbauer: \email{spanb008@umn.edu},\cr
}
