\name{checkCountsFormat} 
\alias{checkCountsFormat} 
\title{GammaModel::checkCountsFormat} 
\usage{checkCountsFormat(counts)}
\arguments{
  \item{counts}{
	A one-row data.frame comprising integer counts, with column 	names in capital letters
	}
}
\value{
	The string 'OK' is returned if no errors are found. 	Otherwise a string stating the number of errors.
	}
\references{
	Kill-off Patterns in Sheep and Goats: The Mandibles from 	Asvan Kale. Payne 1973. 
	}
\description{
	Helper function to check count data are in the required 	format
	}
\details{
	Performs a few basic checks to ensure the raw count data are 	formatted as required: 

	A one-row data frame comprising column names of age classes 	(or multi-classes) in capital letters from 'A'. For example, 	Payne's 1973 caprine age classes are 'A', 'B', 'C', 'D', 	'E', 'F', 'G', 'H', and 'I'. However, it is more likely that 	the raw counts will include some multi-class assignments, 	for example 'ABCD', 'A', 'AB', 'B', 'BC' etc. Order does not 	matter.

	Row names are unconstrained. In the 'Neolithic' dataset 	provided, these specify the site codes.

	Raw counts must by definition be integers.

	Deviations from these formatting requirements are output as 	warnings.
	}
\examples{
    # load the raw data counts for all ten Neolithic sites
    data(Neolithic)

    # select the 'TRA1' subset
    counts <- Neolithic['TRA1',]

    # eyeball data
    print(counts)

    # check format
    checkCountsFormat(counts)

    # change to a non-integer count
    mistake <- counts
    mistake[3] <- 2.5
    checkCountsFormat(mistake)
	}
