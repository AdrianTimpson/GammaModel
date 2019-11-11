\name{GammaModel} 
\alias{GammaModel} 
\title{GammaModel} 
\docType{package}
\description{
	Tools to model count data in discrete classes using the 	Gamma distribution
	}
\details{
	Package  to analyse age-at-death count data across 	discrete age classes.

	Inferential tools where there are large uncertainties 	from a) small sample sizes, and b) multi-class 	assignments. Therefore typically applicable to 	Archaeological data.

	Calculates the exact probability of observed data given 	either a Gamma model or an age-class model. Therefore 	unsuitable for data with large sample sizes in multi-	class assignments, since all possible combinations are 	computed.

	MLE (under either model) found using DEoptimR search 	algorithm.

	Likelihood distribution (under Gamma model) found using 	MCMC Metropolis-Hastings algorithm.

	Performs Goodness of Fit tests, and model comparison 	using AIC and BIC.
	}
\references{
	'Modelling caprine age-at-death profiles using the Gamma 	distribution' by Adrian Timpson, Rosalind E. Gillis, 	Katie Manning and Mark G. Thomas, published in Journal of 	Archaeological Science 2018
	}
