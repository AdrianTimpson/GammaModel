\name{gammaMLparameters} 
\alias{gammaMLparameters} 
\title{GammaModel::gammaMLparameters} 
\usage{
	gammaMLparameters(counts, class.ages = NULL, trace = FALSE)
	}
\arguments{
  \item{counts}{
	A one-row data frame comprising integer counts, with 	column names in capital letters.  Formatting requirements 	detailed in \code{\link{checkCountsFormat}}.
	}
  \item{class.ages}{
	An optional one-row data frame specifying the starting 	age of each age-class.  Age-classes as column names.}
  	\item{trace}{An optional logical value indicating if the 	search progress should be printed.
	}
}
\value{
	Data frame giving the estimated Gamma parameters 	('shape' and 'mean')
	}
\seealso{
	\code{\link[DEoptimR]{JDEoptim}}
}
\references{
	Kill-off Patterns in Sheep and Goats: The Mandibles from 	Asvan Kale. Payne 1973.
	}
\description{
	Finds the Maximum Likelihood Gamma parameters given some 	age-class count data
	}
\details{
	class.ages must be a one-row data frame giving the 	starting age of each age class or age multi-class.  These 	must be sequential, such that the start of the (n+1)th 	class equals the end of the nth class.  Since all 	possible ages must be accounted for, the first starting 	age must be zero, whilst the final age class will 	automatically include ages above this final age.

	If the default class.ages = NULL is used, Payne's caprine 	age classes are used as follows, where units are in 	years:

	class.ages <- data.frame(A=0, B=1/6, C=1/2, D=1, E=2, F=	3, G=4, H=6, I=8)

	The search for the Maximum Likelihood parameters is 	acheived using a bespoke jDE variant of the Differential 	Evolution stochastic algorithm for global optimization, 	as implemented by \code{\link[DEoptimR]{JDEoptim}} in 	Eduardo L T Conceicao's DEoptimR package. 
	}
\examples{

    # Maximum Likelihood Gamma Parameters of the 'TRA1' data
    # using the default Payne caprine age classes.
    data(Neolithic)
    counts <- Neolithic['TRA1',]
    par <- gammaMLparameters(counts)

    # Plot the Gamma distribution defined by these parameters
    # Note, dgamma() requires the shape and rate

    x <- seq(0,10,length.out=1000)
    y <- dgamma(x, shape = par$shape, rate = par$shape/par$mean )
    plot(x, y, type='l', xlab='age', ylab='PD', main='Gamma')

	}
