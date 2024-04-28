\name{hetivbart}
\title{Heterogenous Linear instrumental variables using BART}
\alias{hetivbart}

\description{ hetivBART is a Bayesian \dQuote{sum-of-trees} model.\cr For a
  numeric treatment \eqn{T}, we have the IV first stage: \eqn{T =
    f[1](Z) + eT}{T = f_1(Z) + eT} where \eqn{f[1]}{f_1} follows the BART prior.  For a numeric outcome
  \eqn{Y}, we have the IV second stage: \eqn{Y = f[21](X) + f[22](X) T + eY}{Y
  = f_21(X)T + f_22(X) + eY} where
  \eqn{f[21](X)T+f[22](X)}{f_21(X)T+f_22(X)} follows a semiparametric
  BART prior.
    Errors \eqn{eT} and \eqn{eY} are correlated and their 
  joint distribution is estimated parametrically using a bivariate
  normal specification or nonparametrically using Dirichlet Process Mixtures (DPM).}

\usage{
hetivbart(Z, T, Y, X=NULL,
                  Z.test=NULL,
                  T.test=NULL, 
                  X.test=NULL,
                  burn=1000, nd=2000, burnf1=1000, burnf2=1000,
                  keepevery=10,
                  m1=200, m2=200, nc=100,
                  power=2, base=0.95,
                  k=2, sigmaf1=NA, sigmaf2=NA,
                  doDP=TRUE, do.tsls=TRUE,
                  v=0.17, nu=2.004, a=0.016,
                  Imin=1, Imax=floor(0.1*n)+1, psi=0.5, gs=100,
                  centermeans=TRUE,
                  n=length(Y),
                  f1s = rep(0, n), f2s = rep(0, n),
                  mTs=0, mYs=0, sTs=NA, sYs=NA, gammas=NA,
                  null.f1=FALSE, null.f2=FALSE,
                  printevery=100,quiet=TRUE
                  )
}
\arguments{
   \item{Z}{ A vector of a single instrument or a matrix of multiple instruments. }
   \item{T}{ A vector of the continuous treatment. }
   \item{Y}{ A vector of the continuous outcome. }
   \item{X}{ A vector of a single confounder or a matrix of multiple
     confounders.}
   \item{Z.test}{ A vector or matrix of the instrument to obtain
     \eqn{f_1}{f_1} predictions. }
   \item{T.test}{ A vector or matrix of the instrument to obtain
     \eqn{f_{21}}{f_21} predictions. }
   \item{X.test}{ A vector or matrix of the instrument to obtain
     \eqn{f_{22}}{f_22} predictions. }
   \item{type1}{ Set to \code{1} for a continuous outcome in
     \eqn{f_1}{f_1}, set to \code{2} for a binary outcome w/ probit link
     in \eqn{f_1}{f_1}. }
   \item{type2}{ Set to \code{1} for a continuous outcome in
     \eqn{f_2}{f_2}, set to \code{2} for a binary outcome w/ probit link
     in \eqn{f_2}{f_2}. }
   \item{burn}{ The number of \eqn{\beta}{beta} draws to burn-in and discard. }
   \item{nd}{ The number of \eqn{\beta}{beta} draws to return. }
   \item{burnf1}{ The number of \eqn{f_1} draws to burn-in before the initial value. }
   \item{burnf2}{ The number of \eqn{f_2} draws to burn-in before the initial value. }
   \item{keepevery}{ Every \code{keepevery} draw of \eqn{\beta}{beta} is returned. }
   \item{m1}{ The number of trees for the BART prior on \eqn{f_1}{f_1}. }
   \item{m2}{ The number of trees for the BART prior on \eqn{f[21](X)T+f[22](X)}{f_21(X)T+f_22(X)}. }
   \item{nc}{ The number of cut-points for the BART prior. }
   \item{power}{ Power parameter for trees in the BART prior. }
   \item{base}{ Base parameter for tree in the BART prior. }
   \item{k}{ The number of prior standard deviations \eqn{E(Y|x) = f(x)}
     is away from +/-0.5. }
   \item{doDP}{ The DPM error setting is used when this is \code{TRUE}.}
   \item{do.tsls}{ Do two-stage least squares to get starting values.}
   \item{sigmaf1}{ The SD of \eqn{f}. }
   \item{sigmaf2}{ The SD of \eqn{h}. }
   \item{v}{ Hyperparameter for DPM. }
   \item{nu}{ Hyperparameters for DPM. }
   \item{a}{ Hyperparameters for DPM. }
      \item{sparse1}{ Sparsity parameter for \eqn{f_1}. Use \code{1} for
     traditional BART where every predictor dimension in (z,x) is selected with
     probability 1/p. Use \code{2} for DART where this probability has a
     Dirichlet prior for adaptivity.}
   \item{sparse2}{ Sparsity parameter for \eqn{f_2}. Use \code{1} for
     traditional BART where every predictor dimension in (t,x) is selected with
     probability 1/p. Use \code{2} for DART where this probability has a
     Dirichlet prior for adaptivity.}
   \item{n}{ The number of observations as calculated from \code{length(Y)}
     for convenience. Specifying this parameter on the command line is NOT recommended. }
   \item{Imin}{ The minimum number of DPM clusters. }
   \item{Imax}{ The maximum number of DPM clusters. }
   \item{psi}{ Hyperparameter for \eqn{\alpha}{alpha}, DPM sparsity parameter}
   \item{gs}{ Number of gridpoints for hyperparameter \eqn{\alpha}{alpha}}
   \item{centermeans}{ By default, center \code{T} and \code{Y}.  Set this parameter
     to \code{FALSE} to refrain from centering. }
   \item{in_sample}{ Boolean giving whether to return the in-sample
     predictions. This is useful when the sample size is large and
     memory is a concern if only test predictions are required. }
   \item{aggregate_sigma}{ Currently not used. }
   \item{f1s}{ Starting values for \eqn{f_{1}}{f_1}. }
   \item{f2s}{ Starting values for \eqn{f_{2}}{f_2}. }
   \item{mTs}{ The starting value for \eqn{\mu_{Ts}}{muTs}: TSLS is used
     when not provided. }
   \item{mYs}{ The starting value for \eqn{\mu_{Ys}}{muYs}: TSLS is used
     when not provided. }
   \item{sTs}{ The starting value for \eqn{\sigma_{Ts}}{sdTs}: TSLS is used
     when not provided. }
   \item{sYs}{ The starting value for \eqn{\sigma_{Ys}}{sdYs}: TSLS is used
     when not provided. }
   \item{gammas}{ The starting value for \eqn{\gamma}{gamma}: TSLS is used
     when not provided. }
   \item{null.f1}{ Boolean determining if \eqn{f_1=0}{f_1=0} is fixed. }
   \item{null.f2}{ Boolean determining if \eqn{f_2=0}{f_2=0} is fixed. }
   \item{printevery}{ During MCMC, a message is printed every
     \code{printevery} draws. }
   \item{quiet}{ When \code{TRUE}, output (including \code{printevery}
     updates) will be suppressed during BART
     algorithm run. }
}
\details{
   BART is an Bayesian MCMC method.
   At each MCMC interation, we produce a draw from the joint posterior
   \eqn{(f,\sigma) | (x,y)}{(f,sigma) \| (x,y)} in the numeric \eqn{y} case.
   %and just \eqn{f} in the binary \eqn{y} case.

   Thus, unlike a lot of other modelling methods in R, we do not produce a single model object
   from which fits and summaries may be extracted.  The output consists of values
   \eqn{f^*(x)}{f*(x)} (and \eqn{\sigma^*}{sigma*} in the numeric case) where * denotes a particular draw.
   The \eqn{x} is either a row from the training data (x.train) or the test data (x.test).
}
\value{
  \code{hetivbart} returns a list containing the following items.
  
   \item{dnpart}{ The number of DPM clusters. }
   \item{dalpha}{ Draws of the DPM concentration parameters. }
   \item{dsigma1}{ Draws of the SD from the DPM with a Normal prior for \eqn{eT}. }
   \item{dsigma2}{ Draws of the SD from the DPM with a Normal prior for
     \eqn{eY}. }
   \item{dcov}{ Draws of the covariance from the DPM with a Normal prior for \eqn{eY}. }
   \item{df1}{ Draws of \eqn{f_1}{f_1}. }
   \item{df22}{ Draws of \eqn{f_{22}}{f_22}. }
   \item{df2}{ Draws of \eqn{f_{2}}{f_2}. }
   \item{df1.test}{ Draws of \eqn{f_1}{f_1} using \code{Z.test}}
   \item{df22.test}{ Draws of \eqn{f_{22}}{f_22} using \code{X.test}}
   \item{df2.test}{ Draws of \eqn{f_{2}}{f_2} using \code{T.test} and
     \code{X.test}}
   \item{df1.varcount}{ Counts of splitting rules using each predictor
     for \eqn{f_1}{f_1} at each MCMC iteration. }
   \item{df1.varcount}{ Counts of splitting rules using each predictor
     for \eqn{f_1}{f_1} at each MCMC iteration. } 
}
\references{
  Spanbauer, C., and Pan, W. (20XX)
  Flexible Instrumental Variable Models Using Bayesian Additive Regression Tree Ensembles
  \emph{Submitted}
  
Chipman, H., George, E., and McCulloch R. (2010)
   Bayesian Additive Regression Trees.
   \emph{The Annals of Applied Statistics}, \bold{4,1}, 266-298 <doi:10.1214/09-AOAS285>.
}
\author{
Charles Spanbauer: \email{spanb008@umn.edu}.
}
\examples{

\dontrun{
## Test installation
set.seed(32687)
n <- 250
pZ <- 10
cov <- 0.7
h <- function(v,x) {v*ifelse(x<0,0,1)}
sig. <- matrix(c(1,cov,cov,1),nrow=2)
chl. <- chol(sig.)
Z. <- matrix(rnorm(n*pZ,mean=0,sd=0.5),ncol=pZ)
ZP. <- matrix(rnorm((10e4)*pZ,mean=0,sd=0.5),ncol=pZ)
f1. <- sin(pi*Z.[,1]*Z.[,2])+Z.[,3]^2+2*Z.[,4]+Z.[,5]
f1P. <- sin(pi*ZP.[,1]*ZP.[,2])+ZP.[,3]^2+2*ZP.[,4]+ZP.[,5] 
stdNorm <- matrix(rnorm(n*2),ncol=2)
T. <- (f1. - mean(f1P.))/sd(f1P.) + chl.[1,1]*stdNorm[,1]
X. <- matrix(runif(n*10,min=-1,max=1),nrow=n)
f2. <- h(T.,X.[,1])
Y. <- f2. + chl.[1,2]*stdNorm[,1] + chl.[2,2]*stdNorm[,2]                      
dat <- data.frame(Y=Y.,T=T.,Z=Z.,X=X.)
colnames(dat) <- c('Y','T',paste0('Z',1:pZ),paste0('X',1:10))
Tc <- dat$T
Yc <- dat$Y-mean(dat$Y)
Zc <- dat[,3:12]
Xc <- BART3::bartModelMatrix(dat[,13:ncol(dat)])
Tgrd <- seq(-2.5,2.5,length.out=5)
Xt <- cbind(Tc,Xc)[rep(1:length(Tc),each=10),]
Xt[,1] <- rep(rep(Tgrd,2),length(Tc))
Xt[,2] <- rep(rep(c(-.5,.5),each=5),length(Tc))
post = hetivbart(Z=Zc, T=Tc, Y=Yc, X=Xc,
                 burn=1, burnf1=1, burnf2=1, keepevery=1,
                 nd=1, sigmaf1=1, sigmaf2=1,
                 printevery=1,quiet=FALSE)
print(names(post))


    post = hetivbart(Z=Zc, T=Tc, Y=Yc, X=Xc.
                     nd=1000,burnf1=1000,keepevery=5,burnf2=500,burn=500)                  

                }
              
}
