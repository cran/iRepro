\name{summary.ntestMeans}
\alias{summary.ntestMeans}
\alias{print.summary.ntestMeans}
\title{
Summary for ntestMeans Objects
}
\description{
The function summarizes the results of normality check for means.
}
\usage{
\method{summary}{ntestMeans}(object, ...)
}
\arguments{
  \item{object }{
object of the class ntestMeans (output of the \code{\link{ntest.means}} function)
}
\item{... }{additional arguments passed to the function (they do not affect the summary produced)
}
}
\details{
For more details about normality check and output values shortly described below, please refer to the documentation for \code{\link{ntest.means}}.
}
\value{
An object of class "summary.ntestMeans". The object is a list with the components:
  \item{test.res }{a data frame containing the chi-squared statistic and p-value for normality test}
  \item{mu }{mean of the expected normal distribution for means}
  \item{stdev }{standard deviation of the expected normal distribution for means}
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
\code{\link{ntest.means}}
}
\examples{
# Example with 6 predefined classes (grouped data)
classes <- 1:6
class.limits <- cbind(classes-0.5,classes+0.5)
r1 <- sample(classes,30,replace=TRUE)
r2 <- sample(classes,30,replace=TRUE)
nm <- ntest.means(r1,r2,predefined.classes=TRUE,classes,class.limits,bins=10)
summary(nm)
}
