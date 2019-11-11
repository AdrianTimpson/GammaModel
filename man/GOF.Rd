\name{GOF} 
\alias{GOF} 
\title{GammaModel::GOF} 
\usage{GOF(counts, model, N = 1e+04)}
\arguments{
  \item{counts}{
	A one-row data frame comprising integer counts, with 	column names in capital letters.  Formatting requirements 	detailed in \code{\link{checkCountsFormat}}.
	}
  \item{model}{
	A one-row data frame specifying the age-class model 	probabilities.  Age-classes as column names.
	}
  \item{N}{
	Positive integer specifying the number of random 	arrangements to generate.
	}
}
\value{
	P-value in scientific notation, with two significant 	figures.
	}
\references{
	Pearson, K., 1900. On the criterion that a given 	system 	of deviations from the probable in the case of a 	correlated system of variables is such that it can be 	reasonably supposed to have arisen from random sampling, 	The London, Edinburgh, and Dublin Philosophical Magazine 	and Journal of Science 50, 157-175
	}
\description{
	Estimates Goodness of Fit of the counts to the age-class 	model, using the Chi-squared test.
	}
\details{
	The general approach to calculating a p-value for n 	observed counts given the model probabilities is to 	calculate the probability of each possible arrangement of 	the n counts, then sum the probabilities that are smaller 	or equal to the probability of the observed arrangement.

	This approach can be easily achieved if each count was 	assigned to only one class, but becomes computationally 	expensive when dealing with the additional problem of 	multi-class assignments.

	Instead, the p-value can be estimated by generating 	random arrangements of the observed data (this naturally 	takes care of the different frequencies of each 	arrangement), calculating a p-value for each using a chi-	squared test, then calculating the mean p-value.  With a 	large number of random arrangements this approaches the 	exact probability of the observed data being as or more 	extreme.

	This should not be confused methods of combining p-values 	such as Fisher's, which solves the different problem of 	combining p-values from independent trials.
	}
\examples{
    # How well does the 'FON1' data fit Payne's meat model 
    data(Neolithic)
    data(models.payne)
    counts <- Neolithic['FON1',] 
    model <- models.payne['meat',]
    GOF(counts, model)
	}
