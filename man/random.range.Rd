% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomrange.R
\name{random.range}
\alias{random.range}
\title{random.range - supporting function for other rangemodel fucntions}
\usage{
random.range(uid, nb, range.size, var, first)
}
\arguments{
\item{uid}{a vector of unique ids for selection}

\item{nb}{a neighbour object similar to generated from 'shp2nb'}

\item{range.size}{a vector of number of sites occupied by each species}

\item{var}{an optional vector of variables for constraining the randomization}

\item{first}{If true, var is used while choosing the first occurence as well.
if var is null, first is always set FALSE}
}
\value{
a numeric vector specifying selected possitions in 'uid'
}
\description{
random.range is used within other rangemodel functions to
             radomly place given number of species occurences
}
\details{
this function is not intended for any direct use but is called
          within other functions of 'rangemodelR'.
}
