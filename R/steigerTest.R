#' Comparison of two dependent correlation coefficients
#'
#' This function allows to compare dependent correlation coefficients. It tests whether two correlation coefficients do not differ from each other. 
#' The first possible test  is to compare the correlation between variable a and variable b with the correlation between variable a and variable b (type 1)
#' \deqn{H0: \rho(a,b) = \rho(a,c).} The second possible test is to compare the correlation between variable a and variable b with the correlation between 
#' variable c and variable d (type 2)\deqn{H0: \rho(a,b) = \rho(c,d).}
#' @author Jörn Alexander Quent, \email{alexander.quent@rub.de}
#' @keywords comparison of correlation coefficients
#' @keywords Steiger test
#' @param r numeric vector or matrix containing the observed correlations. If type 1 , then r must be a vector containing the 
#' correlation coefficients in the following order r(a,b), r(a,c), r(b,c). If type 2, then r must be a correlation matrix. 
#' @param n numeric scalar indicating the sample size.
#' @param type numeric scalar indicating which type of test. Possible values are 1 and 2 (see description).
#' @param alternative character string indicating the kind of alternative hypothesis. The possible values are "greater" (the default), "less" and "two-tailed".
#' @return a list of class "htest" containing the results of the hypothesis test. See the help file for \code{\link{htest.object}} for details.
#' @export
#' @references Bortz, J., & Schuster, C. (2010). Statistik für Human- und Sozialwissenschaftler. Berlin, Heidelberg: Springer-Verlag Berlin Heidelberg. Retrieved from http://dx.doi.org/10.1007/978-3-642-12770-0
#' @references Steiger, J. H. (1980). Tests for comparing elements of a correlation matrix. Psychological Bulletin, 87(2), 245-251. http://doi.org/10.1037/0033-2909.87.2.245
#' @examples
#'  r <- matrix(c(1, 0.5, 0.8, 0.5, 0.5, 1, 0.5, 0.7, 0.8, 0.5, 1, 0.6, 0.5, 0.7, 0.6, 1), ncol = 4, nrow = 4)
#'  results <- steigerTest(r, 103, 2)

steigerTest <- function(r, n, type, alternative = 'greater'){
  
  results <- c()
  if(type == 1){
    Z1     <- 0.5 * log((1 + r[1])/(1 - r[1]))
    Z2     <- 0.5 * log((1 + r[2])/(1 - r[2]))
    ra.    <- (r[1] + r[2])/2
    CV     <- (1/(1 - ra.**2))*(r[3]*(1 - 2 * ra. ** 2) - 0.5 * ra.**2 * (1 - 2 * ra.**2 - r[3]**2))
    resukts$alternative <- alternative
    results$null.value  <- 'rho(a,c)'
    names(results$null.value) <- 'rho(a,b)'
  } else if ( type == 2){
    Z1     <- 0.5 * log((1 + r[1,2])/(1 - r[1,2]))
    Z2     <- 0.5 * log((1 + r[3,4])/(1 - r[3,4]))
    rab.cd <- (r[1,2] + r[3,4])/2
    Zae    <- 0.5 * ((r[1,3] - r[1,2] * r[2,3]) * (r[2,4] - r[2,3] * r[3,4]) + (r[1,4] - r[1,3] * r[3,4]) * (r[2,3] - r[1,2] * r[1,3]) + (r[1,3] - r[1,4] * r[3,4]) * (r[2,4] - r[1,2] * r[1,4]) + (r[1,4] - r[1,2] * r[2,4]) * (r[2,3] - r[2,4] * r[3,4]))
    CV     <- Zae/(1 - rab.cd**2)**2
    results$alternative <- alternative
    results$null.value  <- 'rho(a,b)'
    names(results$null.value) <- 'rho(c,d)'
  } else {
    stop('Choose one of the possible tests (1 or 2).')
  }
  spaces <- '\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t'
  results$method       <- paste('Steiger test for comparing \n', spaces, 'two dependent correlation\n', spaces ,'coefficients (Steiger, 1980)')
  results$data.name    <- deparse(substitute(r))
  results$sample.size  <- n
  results$statistic    <- (sqrt(n - 3) * (Z1 - Z2))/sqrt(2 - 2 * CV)
  names(results$statistic) <- 'Z'
  if (alternative == 'greater'){
    results$p.value      <- 1 - pnorm(results$statistic)
  } else if (alternative == 'less'){
    results$p.value      <- pnorm(results$statistic)
  } else if (alternative == 'two-tailed'){
    results$p.value      <- (1 - pnorm(abs(results$statistic)))*2
  } else{
    stop('Choose valid alternative: less, greater or two-tailed')
  }
  oldClass(results)    <- 'htest'
  return(results)
}