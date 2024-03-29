% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CALF-package.R
\encoding{UTF-8}
\name{CALF-package}
\alias{CALF-package}
\title{Coarse Approximation Linear Function}
\description{
Forward selection linear regression greedy algorithm.
}
\details{
The Coarse Approximation Linear Function (CALF) algorithm is a type of forward selection
linear regression greedy algorithm. Nonzero weights are restricted to the values +1 and -1 and
their number limited by an input parameter.  CALF operates similarly on two different types of samples,
binary and nonbinary, with some notable distinctions between the two.
All sample data is provided to CALF as a data matrix. A binary sample must contain a distinguished first
column with at least one 0 entries (e.g. controls) and at least one 1 entry (e.g. cases); at least one
other column contains predictor values of some type.  A nonbinary sample is similar but must contain a 
first column with real dependent (target) values. Columns containing values other that 0 or 1 must be 
normalized, e.g. as z-scores.
As its score of differentiation, CALF uses either the Welch t-statistic p-value or AUC for binary samples
and the Pearson correlation for non-binary samples, selected by input parameter.  When initiated CALF
selects from all predictors (markers) (first in the case of a tie) the one that yields the best score.
CALF then checks if the number of selected markers is equal to the limit provided and terminates if so.
Otherwise, CALF seeks a second marker, if any, that best improves the score of the sum function generated
by adding the newly selected marker to the previous markers with weight +1 or weight -1.
The process continues until the limit is reached or until no additional marker can be included in the sum
to improve the score.
By default, for binary samples, CALF assumes control data is designated with a 0 and case data with a 1.
It is allowable to use the opposite convention, however the weights in the final sum may need to be reversed.
}
\author{
{ Stephanie Lane [aut, cre],\cr
   John Ford [aut],\cr
   Clark Jeffries [aut],\cr
   Diana Perkins [aut]
}
Maintainer: John Ford \email{JoRuFo@gmail.com}
}
\keyword{calf}
