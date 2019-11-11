\name{allArrangements} 
\alias{allArrangements} 
\title{GammaModel::allArrangements} 
\usage{allArrangements(counts)}
\arguments{
  \item{counts}{
	A one-row data frame comprising integer counts, with 	column names in capital letters. Formatting requirements 	detailed in \code{\link{checkCountsFormat}}.
	}
}
\value{
	Data frame of all possible arrangements of the counts in 	rows, and age classes in columns.
	}
\description{
	Generates a data frame of all possible arrangements of 	the counts in rows.  Age classes that always share counts 	remain as multi-classes.
	}
\details{
	Calculating the probability of the observed counts under 	some model probabilities requires calculating the sum of 	the probabilities across all possible arrangements (each 	row). 

	Memory demands from the total number of arrangements 	rapidly becomes huge if there are many counts across many 	multi-classes, rendering the calculation of exact 	probabilities inpractical.
	}
\examples{
    # All 30 possible arrangements for site 'PPI'
    # Notice classes E & F (also H & I) always share counts 
    data(Neolithic)
    counts <- Neolithic['PPI',]
    allArrangements(counts)

    # Number of possible arrangements can be very large
    counts <- Neolithic['TRA2',]
    x <- allArrangements(counts)
    nrow(x)
	}