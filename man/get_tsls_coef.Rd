\name{get_tsls_coef}
\title{Heterogenous Linear instrumental variables using BART}
\alias{get_tsls_coef}

\description{ Computes two-stage-least square estimates from
  the training data and returns inputs exposure coefficient and
  the estimated Cholesky decomposition of \eqn{Sigma} which are
  used for starting values. }

\usage{
get_tsls_coef(T,Z,X,Y)
}

\arguments{
   \item{T}{ Training data for the exposure. }
   \item{Z}{ Training data for the instruments. }
   \item{X}{ Training data for the other covariates. }
   \item{Y}{ Training data for the outcomes. }
}

\author{
Charles Spanbauer: \email{spanb008@umn.edu},\cr
}
