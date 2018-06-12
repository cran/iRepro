\name{summary.ntestRes}
\alias{summary.ntestRes}
\alias{print.summary.ntestRes}
\title{
Summary for ntestRes Objects
}
\description{
The function summarizes the results of normality check for residuals.
}
\usage{
\method{summary}{ntestRes}(object, ...)
}
\arguments{
  \item{object}{
object of the class ntestRes (output of the \code{\link{ntest.res}} function)
}
\item{... }{additional arguments passed to the function (they do not affect the summary produced)
}
}
\details{
For more details about normality check and output values shortly described below, please refer to the documentation for \code{\link{ntest.res}}.
}
\value{
An object of class "summary.ntestRes". The object is a list with the components:
  \item{test.res }{a data frame containing the chi-squared statistics and p-values for normality tests}
  \item{mu }{mean of the expected normal distribution for residuals}
  \item{stdev }{standard deviation of the expected normal distribution for residuals}
  \item{bins }{number of categories in chi-squared normality test}
  \item{df }{number of degrees of freedom in chi-squared normality test}
}
\references{
Kovacic J, Varnai VM. Intraclass correlation coefficient for grouped data. \emph{Epidemiology} 2014;25(5):769--770.
}
\author{
Jelena Kovacic \email{jkovacic@imi.hr}
}
\seealso{
\code{\link{ntest.res}}
}
\examples{
# Example with 6 predefined classes (grouped data)
classes <- 1:6
class.limits <- cbind(classes-0.5,classes+0.5)
r1 <- sample(classes,30,replace=TRUE)
r2 <- sample(classes,30,replace=TRUE)
nr <- ntest.res(r1,r2,predefined.classes=TRUE,classes,class.limits,bins=10)
summary(nr)
}
