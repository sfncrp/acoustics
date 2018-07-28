#' Equivalent Continuous Sound Pressure Level
#'
#' The function '\code{leq}' return the Equivalent Sound Level of a vector of noise levels measured in the same base time period.
#'
#' @param levels a vector of noise levels
#' @param na.rm a logical value indicating whether ‘NA’ values should be
#'        ignored before the computation proceeds.
#'
#' @examples
#' levels <- c(60, 60, 63, 66)
#' leq(levels)
#'
#' levels <- c(60, 60, 63, NA, NA, 66)
#' leq(levels, na.rm=TRUE)
#'
#'
#' @export
leq<-function(levels,na.rm=FALSE){
    if(na.rm==FALSE)
        {
            l=10*log10(sum(10^(levels/10))/length(levels))
        } else{
            l=10*log10(sum(10^(na.omit(levels)/10))/length(na.omit(levels))
                )
        }
  return(l)
}

#' Sum of sound pressure levels
#'
#' \code{lsum} return the energetic sum of a vector of sound pressure levels
#'
#' @inheritParams leq
#'
#' @examples
#'
#' lsum(c(40,40))
#'
#' @export
lsum<-function(levels,na.rm=FALSE){
    if(na.rm==FALSE)
        {
            l<-10*log10(sum(10^(levels/10)))
        } else{
            l<-10*log10(sum(10^(na.omit(levels)/10)))
        }
    return(l)
  }

#' Difference between two sound pressure levels
#'
#' \code{ldiff} return the difference of two sound pressure levels
#'
#' @param level1 sound pressure level
#' @param level2 sound pressure level to subtract from level1
#' @export
#'
#' @examples
#'
#' ldiff(63,60)
ldiff<-function(level1,level2)
{
    diff<-10*log10(10^(level1/10)-10^(level2/10))
    return(diff)
}


#' Normal equal loudness level contours ISO 226:2003
#'
#' \code{eqloudness} return the normal equal loudness level contours
#' defined in the standard ISO 226:2003 for the corresponding loudness value.
#' The validity of the curves is for loudness values from a lower limit of 20 dB(A) to the following upper limits:
#' 20 Hz to 4000 Hz 90 phon,
#' 5000 Hz to 12500 Hz 80 phon.
#'
#' @param phon loudness level, the default value is 40 dB(A)
#'
#' @examples
#'
#' eqloudness(40)
#' plot(eqloudness(40), ylim =c(20,100), type = "l")
#' for(l in seq(20,80,5) )
#'     lines(eqloudness(l), lty=2)
#'
#' @export

eqloudness<-function(phon=40)
  {
## fromISO 226
## Tf<- hearing threshold

## af exponent for loudness perception

## Lu magnitude of the linear transfer function normalized at 1000 Hz

## this applies for values between 20 phon and

## 20 Hz to 4000 Hz 90 phon
## 5000 Hz to 12500 Hz 80 phon

freq<-c(  20.0,    25.0,    31.5,    40.0 ,   50.0,    63.0 ,   80.0,   100.0 ,  125.0, 160.0  , 200.0 ,  250.0 ,  315.0 ,  400.0 ,  500.0  , 630.0  , 800.0,  1000.0,  1250.0 , 1600.0 , 2000.0 , 2500.0 , 3150.0,  4000.0 , 5000.0,  6300.0 , 8000.0, 10000.0 ,12500.0)

af<-c(0.532,0.506,0.480,
      0.455,0.432,0.409,
      0.387,0.367,0.349,
      0.330,0.315,0.301,
      0.288,0.276,0.267,
      0.259,0.253,0.250,
      0.246,0.244,0.243,
      0.243,0.243,0.242,
      0.242,0.245,0.254,
      0.271,0.301)

Lu<-c(
      -31.6,-27.2,-23,
      -19.1,-15.9,-13,
      -10.3,-8.1,-6.2,
      -4.5,-3.1,-2.0,
      -1.1,-0.4,0,
      0.3,0.5,0,
      -2.7,-4.1,-1.0,
      1.7,2.5,1.2,
      -2.1,-7.1,-11.2,
      -10.7,-3.1
      )

Tf<-c(78.5,68.7,59.5,
      51.1,44.0,37.5,
      31.5,26.5,22.1,
      17.9,14.4,11.4,
      8.6,6.2,4.4,
      3.0,2.2,2.4,
      3.5,1.7,-1.3,
      -4.2,-6,-5.4,
      -1.5,6,12.6,
      13.9,12.3)

Ln<-phon

Af<-4.47*10^(-3)*(10^(0.025*Ln)-1.15)+(0.4*10^((Tf+Lu)/10-9))^af

Lp<-(10/af*log10(Af))-Lu+94
names(Lp)<-freq

return((Lp))

}

eqloudness(40)

plot(eqloudness(40), ylim =c(20,100), type = "l")
for(l in seq(20,80,5) )
    lines(eqloudness(l), lty=2)




    freq<-c(20, 25, 31.5, 40, 50, 63, 80, 100, 125, 160, 200, 250, 315, 400, 500, 630, 800, 1000, 1250, 1600, 2000, 2500, 3150, 4000, 5000, 6300, 8000, 10000, 12500, 16000, 20000)

    freq8<-c(31.5, 63, 125, 250, 500, 1000, 2000, 4000, 8000, 16000, 20000)

#' Frequency nominal central bands
#'
#'
#' @examples
#'
#' freqseq()
#' freqseq(1000,12500)
#' freqseq(bands = "octave")
#'
#' @export

freqseq <- function(from = 20, to = 20000, bands = "thirds")
{
    if(bands == "thirds"){
        ifrom = match(from, freq)
        ito = match(to,freq)
        if(is.na(ifrom))
            stop("'from' is not an available central band")
        if(is.na(ito))
            stop("'to' is not an available central band")

        return(freq[ifrom:ito])

    } else if (bands == "octave") {
        from = 31.5
        ifrom = match(from, freq8)
        ito = match(to,freq8)
        if(is.na(ifrom))
            stop("'from' is not an available octave central band")
        if(is.na(ito))
            stop("'to' is not an available octave central band")

        return(freq8[ifrom:ito])


    } else
        stop("available 'thirds' or 'octave' bands")
}

freqseq()

freqseq(1000,12500)

freqseq(bands = "octave")


Dweights <- read.csv("data-raw/frequency_weighting_61672.csv")

is.data.frame(Dweights)


curvaA<-c(-44.7,-39.4,-34.6,-30.2,-26.2,-22.5,-19.1,-16.1,-13.4,-10.9,-8.6,-6.6,-4.8,-3.2,-1.9,-0.8,0,0.6,1,1.2,1.3,1.2,1,0.5,-0.1,-1.1,-2.5)

cbind(freqseq(25,10000),curvaA)

#' Frequency weighting
#'
#'
#'
#' @export

freqweighting <- function(f, w="A")
{
    i <- match(f, Dweights$Frequency)

    if(w=="A") {
        out <- Dweights$A[i]
        names(out) <- f
        return(out)
    } else if (w == "C")
    {
        return(Dweights$C[i])
    }
}

freqweighting(freqseq(), "A")
freqweighting(freqseq(), "C")


curvaA8<-c(-26.2,-16.1,-8.6,-3.2,0,1.2,1,-1.1)
