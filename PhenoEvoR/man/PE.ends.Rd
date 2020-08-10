\name{PE.ends}
\alias{PE.ends}
\docType{data}
\title{
A sample final-timepoint dataset from a Pheno-Evo experiment
}
\description{
This is an example Pheno-Evo dataset containing the final timepoints from each of 9 runs in a BehaviorSpace parameter sweep. It is useful for testing any functions for analyzing or visualizing data from the final timepoints. It was generated from the main example dataset PhenoEvoData (available online with our tutorial) with the following commands:
PE.ends<-extract.endpoint(PhenoEvoData, run.numbers)
pop.level.traits<-c('degrade.rate','switch.rate','response.error','generation')
PE.ends<-summarize.endpoint(PE.ends, pop.level.traits)
dr.fundiv<-calc.fun.div(PE.ends, 'degrade.rate')
PE.ends<-merge(PE.ends, dr.fundiv, by='run.number')
}
\usage{data("PE.ends")}
\format{
  A data frame with 9 observations on the following 20 variables.
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
    \item{\code{mean.degrade.rate}}{a numeric vector}
    \item{\code{mean.switch.rate}}{a numeric vector}
    \item{\code{mean.response.error}}{a numeric vector}
    \item{\code{mean.generation}}{a numeric vector}
    \item{\code{nbsp}}{a numeric vector}
    \item{\code{FRic}}{a numeric vector}
    \item{\code{FEve}}{a numeric vector}
    \item{\code{FDis}}{a numeric vector}
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
data(PE.ends)
generation.heatmap(PE.ends, toxin.conc, env.noise)
phenotype.histogram(PE.ends, xvar=env.noise, yvar=toxin.conc)
}
\keyword{datasets}
