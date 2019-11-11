\name{mcmc} 
\alias{mcmc} 
\title{GammaModel::mcmc} 
\usage{
	mcmc(counts, class.ages = NULL, N = 30000, burn = 2000, thin = 5, prop = 0.4, plot.chain = TRUE)
	}
\arguments{
  \item{counts}{
	A one-row data frame comprising integer counts, with 	column names in capital letters.  Formatting requirements 	detailed in \code{\link{checkCountsFormat}}.
	}
  \item{class.ages}{
	An optional one-row data frame specifying the starting 	age of each age-class.  Age-classes as column names.
	}
  \item{N}{
	An optional integer specifying how many samples in the 	chain.
	}
  \item{burn}{
	An optional integer specifying how many initial samples 	in the chain should be discarded for burn-in.
	}
  \item{thin}{
	An optional integer specifying the proportion of samples 	to discard. I.e. 5 = every 5th sample in the chain is 	retained.
	}
  \item{prop}{
	An optional numerical value controlling the average jump 	size in the proposal function.
	}  
\item{plot.chain}{
	An optional logical value indicating if the chain should 	be plotted after completion.
	}
}
\value{
	Data frame of two columns, giving the Gamma 	parameters ('shape' and 'mean') sampled from the chain, 	after burn-in and thinning.
	}
\references{
	Hastings, W.K., 1970. Monte Carlo sampling methods using 	Markov chains and their applications, Biometrika 57, 	97-109.
	}
\description{
	Generates a single Markov Chain Monte Carlo (MCMC) chain 	of Gamma distribution parameter samples, using the 	Metropolis-Hastings algorithm.
	}
\details{
	Prints progress every 1000th sample, and the final 	acceptance rate. For an efficient chain this should 	be around 0.3 to 0.6. If the rate is too low, efficiency 	can be improved by lowering the value of the prop 	argument.

	The chain can be inspected to ensure good mixing using 	the default plot.chain = TRUE.
	}
\examples{

    # Sampling from the joint likelihood distribution for TRA1
    # using the default Payne class ages.
    data(Neolithic)
    counts <- Neolithic['TRA1',]
    pars <- mcmc(counts)

    # plot the final samples
    par(mfrow = c(1, 1)  )
    plot(pars,  xlab = 'shape',  ylab = 'mean' )

	}
