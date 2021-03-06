% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/annoyance.R
\name{lannoyed}
\alias{lannoyed}
\title{Inverse of the Miedema-Oudshoorn curves.}
\usage{
lannoyed(p, type = "road", annoyance = "HA", Lden_max = 80, C = 72)
}
\arguments{
\item{p}{percentage of people annoyed A, highly HA or little LA annoyed.}

\item{type}{available source type:
'road', 'aircraft', 'railway', 'wtn'.}

\item{annoyance}{annoyance category. Use the default 'HA' for the percentage of highly annoyed people, 'A' for annoyed and 'LA' for little annoyed.}

\item{Lden_max}{maximum admitted Lden value for the specific source.
If exceeded the function returns \code{NA}.}

\item{C}{cut off value for HA, A, LA definitions. Modify only to use different thresholds in the statistical model. If modified the '\code{annoyance}' parameter is ignored.}
}
\description{
The function '\code{lannoyed}' computes the Lden noise level
corresponding to a given percentage of annoyed people '\code{p}'.
It is the inverse function of '\code{\link{pannoyed}}', and
can be computed for annoyed 'A', highly annoyed 'HA' or little
annoyed 'LA' percentages of people exposed to noise and
for different noise sources.
}
\examples{

## noise levels at equal percentage \%HA=20
lannoyed(p = 20, type = "wtn", annoyance = "HA")
lannoyed(p = 20, type = "road", annoyance = "HA")
lannoyed(p = 20, type = "aircraft", annoyance = "HA")
lannoyed(p = 20, type = "railway", annoyance = "HA")

## plotting inverse Miedema-Oudshoorn curves
p_HA <- seq(3,30)
Lden <- lannoyed(p_HA, type = "wtn", annoyance = "HA")

plot(Lden~p_HA)

}
\references{

}
\seealso{
\code{\link{pannoyed}}, \code{\link{pannoyed_conf}}
}
