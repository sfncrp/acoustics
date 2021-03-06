% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic.R
\name{eqloudness}
\alias{eqloudness}
\title{Normal equal loudness level contours ISO 226:2003}
\usage{
eqloudness(phon = 40)
}
\arguments{
\item{phon}{loudness level, the default value is 40 dB(A)}
}
\description{
\code{eqloudness} return the normal equal loudness level contours
defined in the standard ISO 226:2003 for the corresponding loudness value.
The validity of the curves is for loudness values from a lower limit of 20 dB(A) to the following upper limits:
20 Hz to 4000 Hz 90 phon,
5000 Hz to 12500 Hz 80 phon.
}
\examples{

eqloudness(40)
plot(eqloudness(40), ylim =c(20,100), type = "l")
for(l in seq(20,80,5) )
    lines(eqloudness(l), lty=2)

}
