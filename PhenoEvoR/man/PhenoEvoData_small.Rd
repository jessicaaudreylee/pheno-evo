\name{PhenoEvoData_small}
\alias{PhenoEvoData_small}
\docType{data}
\title{
A small PhenoEvoDataset for example use - a few runs, last 10 timepoints.
}
\description{
This dataset is a small sample dataset, a subset of the main PhenoEvoData example dataset available with our tutorial online. It contains only the last 10 timepoints of data from only run numbers 2, 4, and 6. This is useful for testing functions such as extract.endpoint() and summarize.endpoint().
}
\usage{data("PhenoEvoData_small")}
\format{
  A data frame with 30 observations on the following 12 variables.
  \describe{
    \item{\code{run.number}}{a numeric vector}
    \item{\code{toxin.conc}}{a numeric vector}
    \item{\code{env.noise}}{a numeric vector}
    \item{\code{step}}{a numeric vector}
    \item{\code{count.turtles}}{a numeric vector}
    \item{\code{mean.toxin}}{a numeric vector}
    \item{\code{degrade.rate}}{a character vector}
    \item{\code{switch.rate}}{a character vector}
    \item{\code{response.error}}{a character vector}
    \item{\code{barcode}}{a character vector}
    \item{\code{generation}}{a character vector}
    \item{\code{xydr}}{a character vector}
  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
}
\references{
%%  ~~ possibly secondary sources and usages ~~
}
\examples{
data(PhenoEvoData_small)
PE.ends_small<-extract.endpoint(PhenoEvoData_small, c(2, 4, 6))
pop.level.traits<-c('degrade.rate','switch.rate','response.error','generation')
PE.ends_small<-summarize.endpoint(PE.ends_small, traitlist=pop.level.traits)
}
\keyword{datasets}
