% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/steigerTest.R
\name{steigerTest}
\alias{steigerTest}
\title{Comparison of two dependent correlation coefficients}
\usage{
steigerTest(r, n, type, alternative = "greater")
}
\arguments{
\item{r}{numeric vector or matrix containing the observed correlations. If type 1 , then r must be a vector containing the 
correlation coefficients in the following order r(a,b), r(a,c), r(b,c). If type 2, then r must be a correlation matrix.}

\item{n}{numeric scalar indicating the sample size.}

\item{type}{numeric scalar indicating which type of test. Possible values are 1 and 2 (see description).}

\item{alternative}{character string indicating the kind of alternative hypothesis. The possible values are "greater" (the default), "less" and "two-tailed".}
}
\value{
a list of class "htest" containing the results of the hypothesis test. See the help file for \code{\link{htest.object}} for details.
}
\description{
This function allows to compare dependent correlation coefficients. It tests whether two correlation coefficients do not differ from each other. 
The first possible test  is to compare the correlation between variable a and variable b with the correlation between variable a and variable b (type 1)
\deqn{H0: \rho(a,b) = \rho(a,c).} The second possible test is to compare the correlation between variable a and variable b with the correlation between 
variable c and variable d (type 2)\deqn{H0: \rho(a,b) = \rho(c,d).}
}
\examples{
 r <- matrix(c(1, 0.5, 0.8, 0.5, 0.5, 1, 0.5, 0.7, 0.8, 0.5, 1, 0.6, 0.5, 0.7, 0.6, 1), ncol = 4, nrow = 4)
 results <- steigerTest(r, 103, 2)
}
\author{
J�rn Alexander Quent, \email{alexander.quent@rub.de}
}
\references{
Bortz, J., & Schuster, C. (2010). Statistik f�r Human- und Sozialwissenschaftler. Berlin, Heidelberg: Springer-Verlag Berlin Heidelberg. Retrieved from http://dx.doi.org/10.1007/978-3-642-12770-0

Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87(2), 245-251. http://doi.org/10.1037/0033-2909.87.2.245
}
\keyword{Steiger}
\keyword{coefficients}
\keyword{comparison}
\keyword{correlation}
\keyword{of}
\keyword{test}

