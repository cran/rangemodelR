% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/range_shuffle.R
\name{range_shuffle}
\alias{range_shuffle}
\title{Range shuffle models for for range extents reccorded along gradients}
\usage{
range_shuffle(x, boundary, var, interval, sites, reps, degen, lowest = NA,
  highest = NA)
}
\arguments{
\item{x}{Input data for elevational extents of species.data.dataframe
with names 'genus_species','min','max','range','mid',
'num_zones'. See Wang et.al (2012) for details. Species names
column is optional.}

\item{boundary}{nature of boundaries at the extremes of the gradien ie.
 either 'hard boundaries' that species cannot cross, or
 'soft boundaries' that species can move across. Can be
one of the following choises "hh", "sh", "hs", or "ss"}

\item{var}{Predictor variable for constraining randomizations.
dataframe with columns 'mid' and 'weights'. Where 'wights' provide
relative chance of selecting a range location, typically based on
environmental predictors}

\item{interval}{Numeric. Interval between}

\item{sites}{Numeric. Locations on domain for calculating species richness}

\item{reps}{number of iterations}

\item{degen}{logical. If TRUE save each randomized distribution.}

\item{lowest}{Numeric. If analysis is only for a subset of the sampled
gradient then the minimum point within the subset}

\item{highest}{Numeric. If analysis is only for a subset of the sampled
gradient then the maximum point within the subset}
}
\value{
If degen is FALSE, a data frame with four colums for mean, SD and
       confidence intervals of expected richness

\itemize{
 \item{"mod.rich"}{ mean richness of each site}
 \item{"mod.sd"}{ standard deviation of species richness}
 \item{"q2.5"}{ lower limit of the confidence interval}
 \item{"q97.5"}{ upper limit of the confidence interval}
}
       If degen is TRUE, then a list containing above data frame and a list
       of all the randomized matrices
}
\description{
This function was used by Wang et. al.(2012) to test geometric
             constraints, on elevational gradients. The function randomizes
             the range extents and range location for each species and
             returns expected pattern in species richness under geometric
             constraints.
}
\details{
range_shuffle impliments simulations described by Wang et.al (2012)
         to estimate effect of geometric  and environmental constraints on
         pattern in species richness across spatial graidients. It calculates
         a vector of all possible range locations for each given range extent
         based on the conditions for geometric constraints given by
         'boundary'. Range locations are randomized by sampling from this
         vector.
}
\examples{
data(rangedata)
range_shuffle_rnd <- range_shuffle(x=rangedata,boundary = 'hh',
                                   interval = 200, var = NULL,
                                   sites = seq(600,2400,200),
                                   reps =10,degen = FALSE)
range_shuffle_rnd
plot(range_shuffle_rnd$mean,ylab = "Species Richness",pch = 19)
}
\references{
Wang, X., and J. Fang. 2012. Constraining null models with
            environmental gradients: a new method for evaluating the effects
            of environmental factors and geometric constraints on geographic
            diversity patterns. Ecography 35:1147-1159.
}
