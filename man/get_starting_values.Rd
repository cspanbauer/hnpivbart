\name{get_starting_values}
\title{Heterogenous Linear instrumental variables using BART}
\alias{get_starting_values}

\description{ Gets starting values (if not given) from two-stage least
  squares results. }

\usage{
get_starting_values(betas,sTs,sYs,gammas,L,null2)
}

\arguments{
   \item{betas}{ Training data for the exposure. }
   \item{sTs}{ Training data for the instruments. }
   \item{sYs}{ Training data for the other covariates. }
   \item{gammas}{ Training data for the outcomes. }
   \item{L}{ Cholesky of \eqn{Sigma}{Sigma}}
   \item{null2}{ Boolean describing whether we fit the second stage. }
}

\author{
Charles Spanbauer: \email{spanb008@umn.edu},\cr
}
