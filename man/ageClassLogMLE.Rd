\name{ageClassLogMLE} 
\alias{ageClassLogMLE} 
\title{GammaModel::ageClassLogMLE} 
\usage{
	ageClassLogMLE(counts, N = 200, I = 40, C = 2)
	}
\arguments{
  \item{counts}{
	A one-row data frame comprising integer counts, with 	column names in capital letters.  Formatting requirements 	detailed in \code{\link{checkCountsFormat}}.
	}
  \item{N}{
	An integer specifying the population size of particles in each iteration.
	}
  \item{I}{
	An integer specifying the number of iterations.
	}
  \item{C}{
	A numeric value greater than 1, average parameter jump 	size.
	}
}
\value{
	A single numeric, giving the Maximum Likelihood Estimate.
	}
\seealso{
	\code{\link{ageClassMLparameters}}
}
\references{
	Rastrigin, L.A. (1963). The convergence of the random 	search method in the extremal control of a many parameter 	system. Automation and Remote Control. 24 (10): 1337-1342
	}
\description{
	Finds the Maximum Likelihood of an age-class model, given 	some 	age-class count data.
	}
\details{
A search for the MLE is constrained by the fact we are searching for a vector of model probabilities which must by definition sum to 1. Therefore our general approach is to sample from the n-simplex using the Dirichlet distribution. We use a Random Search Algorithm which samples a population of N particles from the Simplex, and retains the particle with the greatest likelihood. At each subsequent iteration, a new population of N particles is sampled from a region of the simplex close to the previous best particle (also including the previous best particle). At each iteration, the average distance to the previous best particle is reduced, the rate of which is controlled by C.
	}
\examples{
    # Maximum Likelihood Age-class parameters of 'TRA1' data
    # using the default Payne caprine age classes.
    data(Neolithic)
    counts <- Neolithic['TRA1',]
    ageClassLogMLE(counts)
	}
