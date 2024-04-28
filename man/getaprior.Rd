\name{getaprior}
\title{Heterogenous Linear instrumental variables using BART}
\alias{getaprior}

\description{ Returns an \eqn{\alpha}{alpha} concentration parameter prior.
}

\usage{
getaprior(n, Imin, Imax, psi, gs)
}

\arguments{
   \item{n}{ The number of observations.}
   \item{Imin}{ The minimum number of clusters; typically, 1. }
   \item{Imax}{ The maximum number of clusters. }
   \item{psi}{ Hyperparameter for \eqn{\alpha}{alpha}, DPM sparsity parameter}
   \item{gs}{ Number of gridpoints for hyperparameter \eqn{\alpha}{alpha}}
}

\author{
Robert McCulloch: \email{robert.e.mcculloch@gmail.com},\cr
Rodney Sparapani: \email{rsparapa@mcw.edu}.
}
