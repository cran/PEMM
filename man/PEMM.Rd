\name{PEMM}
\alias{PEMM}
\docType{package}
\title{
A penalized EM algorithm incorporating missing-data mechanism for multivariate parameter estimation
}
\description{
For many modern high-throughput technologies, missing values arise at high
rates and the missingness probabilities may depend on the values
to be measured. In mass spectrometry based proteomics experiments, the smaller the abundance value of a protein is, the harder the protein can be detected.
That is, the probability of values being missing depends on the values to be 
measured.

  Motivated by data characteristics in mass spectrometry based proteomics studies, 
we consider the problem of estimating mean and covariance of multivariate data 
with ignorable and non-ignorable missingness. The current R package will provide
functions to perform a penalized Expectation-Maximization (EM) algorithm in which
abundance-dependent missing-data mechanisms if present will be incorporated. The 
package is tailored for but not limited to proteomics data, in which sample sizes 
are typically small, and a large proportion of the data are missing-not-at-random. 
The package can be used to jointly estimate the mean abundance and covariance 
structure of multiple functionally-related proteins.
}
\details{
\tabular{ll}{
Package: \tab PEMM \cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2013-11-12\cr
License: \tab GPL \cr
LazyLoad: \tab yes\cr
}
The package contains a PEMM function, which utilizes a penalized EM algorithm 
to estimate the mean and covariance of multivariate Gaussian data with ignorable 
or abundence-dependent missing-data mechanisms. The PEMM function incorporate the 
abundance-dependent missing-data mechanism in the penalized likelihood, 
and obtain the maximum penalized likelihood estimates for multvariate mean and covariance
via the PEMM algorithm.
}
\author{
Lin S. Chen and Pei Wang

Maintainer: Lin S. Chen <lchen@health.bsd.uchicago.edu>
}
\references{
Lin S. Chen, Ross Prentice and Pei Wang. (2014) A penalized EM algorithm
incorporating missing data mechanism for Gaussian parameter estimation. Biometrics, in revision.
}
\keyword{package}
\seealso{
\code{\link{PEMM_fun}} 
}