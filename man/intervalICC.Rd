\name{intervalICC}
\alias{intervalICC}
\alias{intervalICC.default}
\alias{print.ICCfit}
\title{
Intraclass Correlation Coefficient for Interval-Censored Data
}
\description{
The function calculates intraclass correlation coefficient (ICC) for interval-censored data with two repeated measurements. ICC is estimated by maximum likelihood from model with one fixed and one random effect (both intercepts).}
\usage{
intervalICC(r1, r2, predefined.classes=FALSE, classes, c.limits, optim.method=1)
}
\arguments{
  \item{r1}{
data corresponding to the first measurement. If \code{predefined.classes=TRUE} (appropriate for grouped data), this is a vector of length \eqn{n}, where each observation is one of the labels given in \code{classes}. 
Otherwise, if \code{predefined.classes=FALSE}, \code{r1} is a matrix or a data frame with \eqn{n} rows and 2 columns, with columns representing lower and upper bounds of censoring intervals (e.g., if \eqn{i}-th observation
lies in the interval [\eqn{a}, \eqn{b}], then \code{r1[i,]=c(a,b)}).
}
  \item{r2}{
data corresponding to the second measurement. If \code{predefined.classes=TRUE} (appropriate for grouped data), this is a vector of length \eqn{n}, where each observation is one of the labels given in \code{classes}. 
Otherwise, if \code{predefined.classes=FALSE}, \code{r2} is a matrix or a data frame with \eqn{n} rows and 2 columns, with columns representing lower and upper bounds of censoring intervals (e.g., if \eqn{i}-th observation
lies in the interval [\eqn{a}, \eqn{b}], then \code{r2[i,]=c(a,b)}).
}
  \item{predefined.classes}{
logical, indicating whether observations belong to predefined classes (e.g. grouped data in questionnaires) or each observation has its own lower and upper limit (default; \code{FALSE}).
}
  \item{classes}{
a vector with unique labels for the \emph{k} predefined classes. Required if \code{predefined.classes=TRUE}.
}
  \item{c.limits}{
a matrix or a data frame with \emph{k} rows and 2 columns, corresponding to lower and upper bounds of censoring intervals for classes. Required if \code{predefined.classes=TRUE}.
}
  \item{optim.method}{
an integer (1 or 2) specifying the optimization method to be used in maximum likelihood estimation (default is 1). Details are given below.
}
}
\details{
ICC is estimated by maximum likelihood from random effects model
\deqn{Y_{ij} = \mu + b_i + e_{ij},} 
where \eqn{b_i} and \eqn{e_{ij}} are independent and normally distributed with means 0 and variances \eqn{\sigma^2_b} and \eqn{\sigma^2}, respectively. If data were uncensored, this would be analogous to 

\code{lme(ratings~1, random=~1|id, method="ML", data=observed)}

in \code{nlme} package, where 

\code{observed=as.data.frame(rbind(cbind(r1,1:n), cbind(r2,1:n)))}

and \code{colnames(observed)=c("ratings","id")}.
To maximize log-likelihood, \code{constrOptim} from \code{stats} package is used (\code{method=BFGS}).

Two available optimization methods, specified by \code{optim.method}, correspond to two mathematically equivalent expressions for log-likelihood. The option \code{optim.method=1} resulted in slightly more accurate estimates in 
simulations with grouped data, but \code{optim.method=2} was more numerically stable. See the reference for more details.
}
\value{
An object of class "ICCfit". The object is a list with the components:
  \item{icc }{maximum likelihood estimate (MLE) of ICC}
  \item{sigma2.b }{MLE of between-class variance \eqn{\sigma^2_b}}
  \item{sigma2.w }{MLE of within-class variance \eqn{\sigma^2}}
  \item{mu }{MLE of mean \eqn{\mu}}
  \item{loglikelihood }{log-likelihood evaluated at MLE parameters}
}
\references{
Kovacic J, Varnai VM. Reproducibility of grouped data assessed by intraclass correlation coefficient. Accepted for publication in \emph{Epidemiology}.
}
\author{
Jelena Kovacic \email{jkovacic@imi.hr}
}
\note{
If there are many observations with same values (i.e. with the same lower and upper bounds), it is advisable to group all observations into \code{classes} and use option \code{predefined.classes=TRUE}; this will reduce 
computation time.

Subjects with only one measurement are omitted from ICC calculation. 
}
\seealso{
\code{\link{print.ICCfit}}, \code{\link{summary.ICCfit}}
}
\examples{
# Example with 6 predefined classes (grouped data)
classes <- 1:6
class.limits <- cbind(classes-0.5,classes+0.5)
r1 <- sample(classes,30,replace=TRUE)
r2 <- sample(classes,30,replace=TRUE)

intervalICC(r1,r2,predefined.classes=TRUE,classes,class.limits)

# The same result can be obtained with predefined.classes=FALSE option, 
# although with slower computation time
rtg1 <- matrix(nrow=30,ncol=2)
rtg2 <- matrix(nrow=30,ncol=2)
# when predefined.classes=FALSE, ratings must be given with lower and upper bounds 
# for each observation:
for(i in 1:length(classes)){
  rtg1[r1==classes[i],1] <- class.limits[i,1]
  rtg1[r1==classes[i],2] <- class.limits[i,2]
  rtg2[r2==classes[i],1] <- class.limits[i,1]
  rtg2[r2==classes[i],2] <- class.limits[i,2]
}

intervalICC(rtg1,rtg2,predefined.classes=FALSE)
}