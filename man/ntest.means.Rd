\name{ntest.means}
\alias{ntest.means}
\title{
Normality Check for Interval-Censored Data with Repeated Measurements - Means
}
\description{
The function checks whether interval-censored data with two repeated measurements meet the normality assumption for subjects' means. This is a prerequisite for the random effects model used in ICC calculation.
}
\usage{
ntest.means(r1, r2, predefined.classes=FALSE, classes, c.limits, optim.method=1, bins=10)
}
\arguments{
  \item{r1}{
argument passed to \code{\link{intervalICC}}; see documentation for that function.
}
  \item{r2}{
argument passed to \code{\link{intervalICC}}; see documentation for that function.
}
  \item{predefined.classes}{
argument passed to \code{\link{intervalICC}}; see documentation for that function.
}
  \item{classes}{
argument passed to \code{\link{intervalICC}}; see documentation for that function.
}
  \item{c.limits}{
argument passed to \code{\link{intervalICC}}; see documentation for that function.
}
  \item{optim.method}{
argument passed to \code{\link{intervalICC}}; see documentation for that function.
}
  \item{bins}{
number of categories in chi-square test; see details below (default is 10).
}
}
\details{
For ICC estimation the random effects data model
\deqn{Y_{ij} = \mu + b_i + e_{ij},} 
is used, where \eqn{b_i} and \eqn{e_{ij}} are normally distributed with means 0 and variances \eqn{\sigma^2_b} and \eqn{\sigma^2}, respectively. 
This function assesses the assumption that the subjects' means \eqn{0.5 (Y_{i1}+Y_{i2})} are normally distributed with mean \eqn{\mu} and variance \eqn{\sigma^2_b + 0.5 \sigma^2}, as is expected under the specified model.

To test normality, chi-square goodness-of-fit test with \code{bins} subsequent data categories is used (call to \code{chisq.test} from package \code{stats}). 
The categories (bins) are determined using the equidistant quantiles of expected normal distribution, with corresponding maximum likelihood 
parameters. Maximum likelihood estimates for parameters \eqn{\mu}, \eqn{\sigma^2_b} and \eqn{\sigma^2} are obtained by calling the function \code{\link{intervalICC}}.
The probability corresponding to each bin is \code{1/bins} (expected relative frequencies; this corresponds to \code{p = rep(1/bins,bins)} in \code{chisq.test} function). 
Since means are interval-censored and censoring intervals overlap, the observed relative frequencies are calculated in the following way. If one of the original intervals representing subjects mean
spans multiple bins, each bin receives a share of votes from the original interval. This share is calculated using the
expected normal density function and it is proportional to the probability of data falling within the intersection of the
original interval and bin.
}
\value{
An object of class "ntestMeans". The object is a list with the components:
  \item{statistic }{value of chi-squared statistic; \code{statistic} in the output of \code{chisq.test}}
  \item{parameter }{number of degrees of freedom for chi-squared distribution; \code{parameter} in the output of \code{chisq.test}}
  \item{p.value }{p-value of test; \code{p.value} in the output of \code{chisq.test}}
  \item{data }{character string with value "means"}
  \item{mu }{mean of the expected normal distribution for subjects' means; equal to maximum likelihood estimate for \eqn{\mu} from \code{intervalICC}}
  \item{var }{variance of the expected normal distribution for subjects' means; equal to maximum likelihood estimate for \eqn{\sigma^2_b + 0.5 \sigma^2} from \code{intervalICC}}
  \item{bins }{number of categories in chi-square test}
}
\references{
Kovacic J, Varnai VM. Intraclass correlation coefficient for grouped data. \emph{Epidemiology} 2014;25(5):769--770.
}
\author{
Jelena Kovacic \email{jkovacic@imi.hr}
}
\note{
This function was designed as a help in assessing goodness of model fit. However, it has not been tested in simulations nor in any other way.
It is the responsibility of the user to provide appropriate number of bins; the function checks only if \code{bins} is a positive integer.
Testing normality with low number of bins is unreliable. On the other hand, if the number of bins is too large, \code{chisq.test} will complain since the expected frequencies will be too low.
}
\seealso{
\code{\link{summary.ntestMeans}}, \code{\link{intervalICC}}, \code{\link[stats]{chisq.test}}
}
\examples{
# Example with 6 predefined classes (grouped data)
classes <- 1:6
class.limits <- cbind(classes-0.5,classes+0.5)
r1 <- sample(classes,30,replace=TRUE)
r2 <- sample(classes,30,replace=TRUE)
ntest.means(r1,r2,predefined.classes=TRUE,classes,class.limits,bins=10)
}
