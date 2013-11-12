\name{PEMM_fun}
\alias{PEMM_fun}
\title{
A penalized EM algorithm incorporating missing-data mechanism for 
multivariate parameter estimation
}
\description{
The PEMM function utilizes a penalized EM algorithm to estimate 
the mean and covariance of multivariate Gaussian data with 
ignorable or abundence-dependent missing-data mechanism. For 
example, in proteomics data, the smaller the abundance value of a protein is, the more likely the protein cannot be detected in the experiment.
The PEMM function incorporates the known or estimated abundance-dependent 
missing-data mechanism in the penalized likelihood, and obtains the 
maximum penalized likelihood estimates for multvariate mean and 
covariance via the PEMM algorithm.
}
\usage{
  PEMM_fun(X, phi, lambda = NULL, K = NULL, pos = NULL, tol = 0.001, 
          maxIter = 100)
}
\arguments{
  \item{X}{ 
  a n-by-p matrix, the value of the multivariate incomplete data. 
Each rown is a sample and each column is one feature (e.g. protein). 
}
  \item{phi}{ 
  the parameter in the missing-data mechanism. It is a non-negative 
number. We model the abundence-dependent missing-data mechanism as  
          P(missing Xi)=constant*exp(-phi*Xi) 
if Xi is strictly positive; or 
          P(missing Xi)=constant*exp(-phi*Xi^2) 
if Xi can be positive and negative. When phi is set to 0, 
the data will be treated as missing-at-random (MAR) and a penalized EM 
algorithm will be performed without the consideration of 
missing-data mechanism. When phi is specified as being greater 
than zero, the PEMM algorithm will be utilized.   
}
  \item{lambda}{ 
  the tuning parameter in the Inverse-Wishart penalty function.
Appropriate value of lambda helps to bound the eigen-values of the covariance matrix away from zero 
and thus assures the positive-definiteness of the estimated covariance .
}
  \item{K}{ 
  the second tuning parameter in the Inverse-Wishart penalty function.
Appropriate choice of K helps to stabalize the estimation of the covariance matrix. Here we suggest K being a small 
positive integer (=5), which seems work well when dimensionality and 
sample size are both moderate (<100).
}
  \item{pos}{ 
  a logical value. If pos=TRUE, all X are one-sided and the
abundance-dependent missing-data mechanism is given by
             P(missing Xi)=constant*exp(-phi*Xi).
If pos=FALSE, X can be two-sided. The 
abundance-dependent missing-data mechanism is given by
             P(missing Xi)=constant*exp(-phi*Xi^2).
}
  \item{tol}{ 
  the tolerance to assess the convergence of the iterative algorithm. We set tol=0.001 as default, 
i.e., when the change in log-likelihood of the current iteration 
versus the last iteration is less than 0.1\%, the algorithm stops.
}
  \item{maxIter}{ 
  the maximum number of iterations allowed. We set maxIter=100 as default,
i.e., if after 100 iterations the algorithm still do not
converge, we force the algorithm to stop and give the users a warning.
}
}
\details{
The algorithm is motivated by data characteristics in proteomics data 
with substantial abundance-dependent missing values. The algorithm 
calculates the maximum penalized likelihood estimates 
of mean and covariance of multivariate incomplete data, with abundence-dependent missing data mechanism incorporated.
}
\value{
The algorithm will return a list of estimated mean and covariance, and
the imputed missing-values from the last iteration.
\item{mu}{ the estimated mean of the multivariate incomplete data}
\item{Sigma}{ the estimated covariance matrix of the multivariate 
incomplete data}
\item{Xhat}{ the "imputed" missing-values from the last iteration.}
}
\references{
Lin S. Chen, Ross Prentice and Pei Wang. (2014) A penalized EM algorithm
incorporating missing data mechanism for Gaussian parameter estimation. 
Biometrics, in revision.
}
\author{
Lin S. Chen and Pei Wang

Maintainer: Lin S. Chen <lchen@health.bsd.uchicago.edu>
}
\note{
Motivated by characteristics and the abundance-dependent missing-data 
mechanism in proteomics data, the current function takes the known or
estimated abundance-dependent missing-data parameter as input and 
estimates the mean and covariance of multivariate incomplete data.

When the missing-data mechansim is MAR or MCAR, one could specify phi 
being zero so that a penalized EM algorithm will be performed without modelling
of the missing-data mechanism.

When the missing data are non-ignorable and the missing-data mechanism
is unknown, one would first examine if the data follow abundance-dependent
pattern and if so, one can estimate the abundance-dependent parameters and perform
the analysis based on the PEMM algorithm. Otherwise, if the missing data are not
abundence-dependent, one would need alternative approaches.
}
\examples{

set.seed(111)
library(PEMM)
data(sim_dat)
X.mar=sim_dat
X.mar[sample(1:length(X.mar),round(length(X.mar)*0.2))]<-NA

## If data are MAR or MCAR, by only specifying phi=0, 
## a penalized EM algorithm will be performed at default.
PEM.result = PEMM_fun(X.mar, phi=0) 

## By specifying phi=0, lambda=0, K=0, an EM algorithm will be performed. 
## Although when n is small, EM may not converge.
EM.result = PEMM_fun(X.mar, phi=0, lambda=0, K=0) 

## Generate data with non-ignorable missingness -- observations with 
## lower absolute values are more likely to be missing
phi=1
prob <- 0.5*exp(-phi*(sim_dat)^2)
X.mnar=sim_dat
X.mnar[which(rbinom(length(X.mnar),1,prob)==1)] <- NA
mean(is.na(X.mnar))  ## proportion of data being missing 

## Geting the estimated results
PEMM.result.mnar = PEMM_fun(X.mnar, phi=1) 
PEM.result.mnar = PEMM_fun(X.mnar, phi=0)  ## ignoring MNAR mechanism
EM.result.mnar = PEMM_fun(X.mnar, phi=0, lambda=0, K=0)  ## ignoring MNAR

## Compare the mean estimates for data with MNAR from different methods
## complete data
colMeans(sim_dat)

## EM results ignoring MNAR mechanism
EM.result.mnar$mu

## PEMM estimates
PEMM.result.mnar$mu

cor(colMeans(sim_dat),PEMM.result.mnar$mu)
cor(colMeans(sim_dat),EM.result.mnar$mu)

}
\keyword{ abundence-dependent missing-data mechanism }
\keyword{ a penalized EM algorithm }
\keyword{ non-ignorable missing data }
