#' Day-evening-night noise level computation.
#'
#' The function '\code{lden}' computes the day-evening-night noise level
#' defined by the European Directive END 2002/49/EC.
#' The function takes as input the Lday, Levening and Lnight noise
#' levels and the number of hours of each time period.
#' 
#' The default values for the time periods duration consider the
#' Italian national legislations, with  h_day = 12 hours,
#' h_evening = 2, h_night = 10,
#' corrisponding to daytime 06-20, evening 20-22 and nighttime 22-06.
#' The values can be changed according to each specific national legislation.
#' 
#' 
#' @param Lday Lday is the A-weighted long-term average sound level
#' as defined in ISO 1996-2: 1987,
#' determined over all the day periods of a year
#' @param Levening Levening
#' @param Lnight Lnight
#' @param h_day  Day time period hours
#' @param h_evening Evening hours
#' @param h_night Night hours
#'
#' @references
#' 
#' @examples
#'
#' ## Lden computation using the default Italian time periods
#' lden(63,60,60)
#'
#' ## Lden computation using values for a different legislation
#' lden(63,60,60, h_day=12, h_evening=4, h_night=8)
#'
#' @export

lden <- function(Lday,Levening,Lnight, h_day=14, h_evening=2, h_night =8 ){
    if( (h_day+h_evening+h_night) != 24)
        stop("Error, total hours exceeded 24")
    else{
        Lden <- 10*log10(1/24*(h_day*10^(Lday/10)+h_evening*10^((Levening+5)/10)+h_night*10^((Lnight+10)/10)))
        return(Lden)
        }
}

###########################################################

par_road<-list(
    intercept= -106.97,
    slope= 2.22,
    var_within= 1150.71, 
    var_between = 150.54, 
    slope_var = 0.0023,
    intercept_var = 15.3,
    cov_slop_intercept = -0.146,
    tot_sd = sqrt(1150.71+150.54)
)

par_wtn<-list(
    intercept= -226.88,
    slope= 5.48,
    slope_sd = 0.45,
    slope_var = 0.21,
    intercept_sd = 19.21,
    intercept_var = 369.09,
    cov_slop_intercept = -8.64,
    tot_sd = 47.2
)

par_aircraft<-list(
    intercept= -91.42,
    slope= 2.17,
    var_within= 1187.11, 
    var_between = 77.64, 
    slope_var = 0.0017,
    intercept_var = 10.88,
    cov_slop_intercept = -0.105,
    tot_sd = sqrt(1187.11+77.64)
)

par_railway<-list(
    intercept= -110.09,
    slope= 2.10,
    var_within= 1078.73, 
    var_between = 53.86,   
    slope_var = 0.0071,
    intercept_var = 40.04,
    cov_slop_intercept = -0.478,
    tot_sd = sqrt(1078.73+53.86)
)

par <- list("road"=par_road, "aircraft"= par_aircraft,
           "railway"=par_railway, "wtn" = par_wtn)

sources <- c("road", "aircraft", "railway", "wtn")

###########################################################

miedema_model <- function(Lden, param, C=72, Lden_max = 80)
{
    ## Miedema-Oudshoorn statistical model
    intercept = param$intercept
    slope = param$slope
    tot_sd = param$tot_sd

    if(Lden>Lden_max){
        annoyed<- NA
    } else{
        argument <- (C-intercept-slope*Lden)/tot_sd
        annoyed <-  100*(1-pnorm(argument))
    }
    return(annoyed)
}

## miedema_model(c(60,90),par_road, C=72)

choose_param <- function(type="road")
{
    if(type == "road")
        param = par$road
    else if(type == "aircraft")
        param = par$aircraft
    else if(type == "railway")
        param = par$railway
    else if(type == "wtn")
        param = par$wtn
    else
        stop('Available types: road, aircraft, railway, wtn')
return(param)
}

choose_cutoff <- function(annoyance, C=72 )
{
    if( C == 72 )
    {
        if(annoyance == "HA")
            C = 72
        else if(annoyance == "A")
            C = 50
        else if(annoyance == "LA")
            C = 28
        else
            stop('Available values : HA, A, LA')
    }
    return (C)
}

###########################################################

#' Miedema-Oudshoorn exposure-response curves.
#'
#' The function '\code{pannoyed}' computes the percentage of people
#' highly annoyed, annoyed or little annoyed for different sources,
#' using the Miedema-Oudshoorn exposure-response curves.
#' For transportation noise are used the data from
#' the cumulative analysis by Miedema-Oudshoorn 2001,
#' reported in references.
#' For Wind Turbine Noise the data from Janssen 2011 are used.
#' 
#' '\code{pannoyed}' takes as input both a single Lden value or a
#' vector of Lden values, besides the parameters to choose the
#' source type and annoyance category.
#' The function '\code{\link{lannoyed}}' computes
#' the inverse function of 'pannoyed',
#' obtaining the Lden value corresponding to a given percentage
#' of annoyed people '\code{p}'.
#' To compute the confidence and tolerance interval see the function
#' 'pannoyed_conf'.
#'
#' @param Lden day-evening-night noise level.
#' Can be both a single value or a vector of values.
#' @param type available source type:
#' 'road', 'aircraft', 'railway', 'wtn'.
#' @param annoyance  annoyance category. Use the default 'HA' for the percentage of highly annoyed people, 'A' for annoyed and 'LA' for little annoyed.
#' @param Lden_max maximum admitted Lden value for the specific source.
#' If exceeded the function returns \code{NA}.
#' @param C cut off value for HA, A, LA definitions. Modify only to use different thresholds in the statistical model. If modified the '\code{annoyance}' parameter is ignored.
#'
#' @seealso \code{\link{lannoyed}}, \code{\link{pannoyed_conf}}
#' 
#' @references
#'
#' @examples
#'
#' ## percentage of highly annoyed HA at Lden = 50dB(A)
#' pannoyed(50, type = "road", annoyance = "HA")
#' pannoyed(50, type = "aircraft", annoyance = "HA")
#' pannoyed(50, type = "railway", annoyance = "HA")
#' pannoyed(50, type = "wtn", annoyance = "HA")
#'
#' ## road exposure-response curves for HA, A, LA
#' Lden <- seq(30,70)
#' pannoyed(Lden, type = "road", annoyance = "HA")
#' pHA_road <- pannoyed(Lden, type = "road", annoyance = "HA")
#' pA_road <- pannoyed(Lden, type = "road", annoyance = "A")
#' pLA_road <- pannoyed(Lden, type = "road", annoyance = "LA")
#'
#' plot(pLA_road~Lden, type = "l", col="red",
#' xlab="Lden [dB(A)]", ylab="%p annoyed")
#' lines(pA_road~Lden, col="blue")
#' lines(pHA_road~Lden, col="green")
#'
#' ## comparison of HA exposure-response for different sources
#' Lden <- seq(30,70)
#' p_road <- pannoyed(Lden, type = "road", annoyance = "HA")
#' p_aircraft <- pannoyed(Lden, type = "aircraft", annoyance = "HA")
#' p_railway <- pannoyed(Lden, type = "railway", annoyance = "HA")
#' p_wtn <- pannoyed(Lden, type = "wtn", annoyance = "HA")
#'
#' plot(p_road~Lden, type ="l", col="red",
#' xlab="Lden [dB(A)]", ylab="%HA")                  
#' lines(p_aircraft~Lden, col = "blue")
#' lines(p_railway~Lden, col = "yellow")
#' lines(p_wtn~Lden, col= "green")
#'
#' @export

pannoyed <- function(Lden, type= "road", annoyance = "HA", Lden_max= 80,  C = 72 ){
    
    param <- choose_param( type = type )
    C <- choose_cutoff( annoyance = annoyance, C=C )

    ## apply Miedema-Oudshoorn statistical model

    p <- sapply(Lden, miedema_model, param = param,
                C=C, Lden_max=Lden_max)

    return(p)
}

## percentage of highly annoyed HA at Lden = 50dB(A)
pannoyed(50, type = "road", annoyance = "HA")
pannoyed(50, type = "aircraft", annoyance = "HA")
pannoyed(50, type = "railway", annoyance = "HA")
pannoyed(50, type = "wtn", annoyance = "HA")

## road exposure-response curves for HA, A, LA
Lden <- seq(30,70)
pannoyed(Lden, type = "road", annoyance = "HA")

pHA_road <- pannoyed(Lden, type = "road", annoyance = "HA")
pA_road <- pannoyed(Lden, type = "road", annoyance = "A")
pLA_road <- pannoyed(Lden, type = "road", annoyance = "LA")

plot(pLA_road~Lden, type = "l", col="red", xlab="Lden [dB(A)]", ylab="%p annoyed")
lines(pA_road~Lden, col="blue")
lines(pHA_road~Lden, col="green")

# comparison of HA exposure-response for different sources
Lden <- seq(30,70)
p_road <- pannoyed(Lden, type = "road", annoyance = "HA")
p_aircraft <- pannoyed(Lden, type = "aircraft", annoyance = "HA")
p_railway <- pannoyed(Lden, type = "railway", annoyance = "HA")
p_wtn <- pannoyed(Lden, type = "wtn", annoyance = "HA")

plot(p_road~Lden, type ="l", col="red", xlab="Lden [dB(A)]", ylab="%HA")
lines(p_aircraft~Lden, col = "blue")
lines(p_railway~Lden, col = "yellow")
lines(p_wtn~Lden, col= "green")
###########################################################

## Inverse of the Miedema-Oudshoorn curves.

## function applied only to a single value, not exported

lannoyed_single <- function(p, param, C, Lden_max){
    f_root<-function(x,p){
        annoyed <-miedema_model(Lden=x, param=param, C=C)
        out <- annoyed-p
        return(out)
    }

    root<- uniroot(f_root, interval=c(1, Lden_max), p=p)
    Lden <- root$root

    return(round(Lden,3))
}

#' Inverse of the Miedema-Oudshoorn curves.
#'
#' The function '\code{lannoyed}' computes the Lden noise level
#' corresponding to a given percentage of annoyed people '\code{p}'.
#' It is the inverse function of '\code{\link{pannoyed}}', and
#' can be computed for annoyed 'A', highly annoyed 'HA' or little
#' annoyed 'LA' percentages of people exposed to noise and
#' for different noise sources.
#'  
#' @param p percentage of people annoyed A, highly HA or little LA annoyed.
#' @inheritParams pannoyed
#'
#' @seealso \code{\link{pannoyed}}, \code{\link{pannoyed_conf}}
#'
#' @references
#'
#' @examples
#'
#' ## noise levels at equal percentage %HA=20
#' lannoyed(p = 20, type = "wtn", annoyance = "HA")
#' lannoyed(p = 20, type = "road", annoyance = "HA")
#' lannoyed(p = 20, type = "aircraft", annoyance = "HA")
#' lannoyed(p = 20, type = "railway", annoyance = "HA")
#'
#' ## plotting inverse Miedema-Oudshoorn curves
#' p_HA <- seq(3,30)
#' Lden <- lannoyed(p_HA, type = "wtn", annoyance = "HA")
#'
#' plot(Lden~p_HA)
#' 
#' @export

lannoyed <- function(p, type= "road", annoyance = "HA", Lden_max= 80,  C = 72 ){

    param <- choose_param(type = type)
    C <- choose_cutoff(annoyance = annoyance, C=C)

    Lden <- sapply(p, lannoyed_single, param=param, C=C, Lden_max=Lden_max)
    return(Lden)
}

## noise levels at equal percentage %HA=20 
lannoyed(p = 20, type = "wtn", annoyance = "HA")
lannoyed(p = 20, type = "road", annoyance = "HA")
lannoyed(p = 20, type = "aircraft", annoyance = "HA")
lannoyed(p = 20, type = "railway", annoyance = "HA")

## plotting inverse Miedema-Oudshoorn curves
p_HA <- seq(3,30)
Lden <- lannoyed(p_HA, type = "wtn", annoyance = "HA")

plot(Lden~p_HA)

###########################################################

pannoyed_conf_single_value<- function(Lden, level=0.95, tol=FALSE, Lden_max= 80,  C = 72, param = NULL ){

    intercept <- param$intercept
    slope <- param$slope
    intercept_var <-param$intercept_var
    slope_var <-param$slope_var
    cov_slop_intercept <- param$cov_slop_intercept
    tot_sd <- param$tot_sd
    
    x <- c(1,Lden)    
    b <- c(intercept,slope)

    Sb <- matrix(c(intercept_var,cov_slop_intercept,
                   cov_slop_intercept,slope_var), nrow=2,ncol=2)

    level_less <- (1+level)/2
    k <- qnorm(level_less)

    
    if(tol==FALSE){
        ## confidence interval
        Clu_pos <- x%*%b + k*sqrt(x %*% Sb %*% x)
        Clu_neg <- x%*%b - k*sqrt(x %*% Sb %*% x)

        CL_pos<- 100 * (1- (pnorm( (C-Clu_pos)/tot_sd)))
        CL_neg<- 100 * (1- (pnorm( (C-Clu_neg)/tot_sd)))
    
    } else {
        ## tolerance interval
        sd_city <-sqrt(param$var_between)

        Clu_pos <- x%*%b + k*sqrt(x %*% Sb %*% x+sd_city^2)
        Clu_neg <- x%*%b - k*sqrt(x %*% Sb %*% x+sd_city^2)

        CL_pos<- 100 * (1- (pnorm( (C-Clu_pos)/tot_sd)))
        CL_neg<- 100 * (1- (pnorm( (C-Clu_neg)/tot_sd)))
    }

    if(Lden>Lden_max){
        return(c(conf_low=NA,
                 conf_up=NA
                 ))
    } else {
    return(c(conf_low=CL_neg,
                 conf_up=CL_pos
                 ))

    }
        
}

#' Confidence and tolerance intervals for Miedema-Oudshoorn curves.
#'
#' '\code{pannoyed_conf}' computes the confidence and tolerance
#' intervals for the Miedema-Oudshoorn exposure-response curve
#' calculated with '\code{\link{pannoyed}}'.
#' It returns a dataframe containing the confidence/tolerance intervals
#' and the exposure curve.
#'
#' @inheritParams pannoyed
#' @param level confidence level, default: 0.95
#' @param tol if TRUE returns both the confidence and the tolerance
#' intervals. For wind turbine noise tolerance intervals are not
#' available at the moment.
#'
#' 

pannoyed_conf<- function(Lden, level=0.95, type= "road",  annoyance = "HA", tol=FALSE, Lden_max= 80,  C = 72 ){

    if(type == "wtn" && tol==TRUE)
        stop('Tolerance interval not available for WTN')
    
    param <- choose_param(type = type)
    C <- choose_cutoff(annoyance = annoyance, C=C)

    
    Lconf <- lapply(Lden, pannoyed_conf_single_value,
                        param=param, C=C)

    p <- sapply(Lden, pannoyed,
               type=type, annoyance=annoyance,
               Lden_max=Lden_max, C=C)

    Dconf <- data.frame(Lden, p, do.call(rbind, Lconf))
    colnames(Dconf) <- c("Lden", paste("p",annoyance, sep ="")        , "conf_low", "conf_up")
        
    if(tol==TRUE)
    {
        Ltol <- lapply(Lden, pannoyed_conf_single_value,
                        param=param, C=C, tol=TRUE)
        Dtol <- do.call(rbind, Ltol)
        colnames(Dtol) <- c("tol_low", "tol_up")
        
        Dconf <- cbind(Dconf, Dtol)
                
    }

        
        return(Dconf)
}

## confidence and tolerance interval
pannoyed_conf(60, annoyance = "A", type = "road")
pannoyed_conf(60, annoyance = "A", type = "road", tol = TRUE)


Lden <- seq(30,70)
df <- pannoyed_conf(Lden, annoyance = "HA", type = "road", tol=TRUE)
head(df)

## plotting Miedema-Oudshoorn curves with confidence and tolerance intervals
plot(df$pHA~df$Lden, type = "l", xlab = "Lden [dB(A)]", ylab = "%HA")
lines(df$conf_low~df$Lden, lty = 2)
lines(df$conf_up~df$Lden, lty = 2)
lines(df$tol_low~df$Lden, lty = 3)
lines(df$tol_up~df$Lden, lty = 3)

## plotting with using ggplot
library(ggplot2)
library(reshape2)
## transforming data in long format
df_long <- melt(df, id.vars = "Lden")
head(df_long)

p <- ggplot(df_long, aes(x=Lden, y =value, linetype =variable))
p+geom_line()+ylab("%HA")+xlab("Lden [dB(A)]")

###########################################################

#' Conversion curves for equally annoyed people.
#' 
#' The conversion curves '\code{lconvert}' are linear curves
#' which convert a noise level of a reference source type
#' to the noise level of a different source, having the same percentage
#' of annoyed, highly annoyed or little annoyed people.
#'
#' The function \code{lconvert} directly convert a noise level
#' from a reference source to the noise level of a different source.
#' The default reference source is 'road' traffic noise, but the
#' function can be used to convert levels between any couple of the
#' available sources. The conversion is possible for each annoyance
#' category, 'HA', 'A' or 'LA'.
#' 
#'
#' The conversion curves are presented in the paper of
#' Fredianelli, Carpita, Licitra in references,
#' and are based on the studies of Miedema-Oudshoorn and
#' Janssen et al.
#'
#' @param Lref noise level of the reference source to be converted.
#' @param ref type of reference source, 'road' is the default value.
#' @param type type of the converted noise level, available source type:
#' 'road', 'aircraft', 'railway', 'wtn'.
#' @inheritParams pannoyed
#'
#' @seealso \code{\link{lannoyed}}, \code{\link{pannoyed}} 
#' 
#' @references
#'
#' 
#' @examples
#'
#' ## conversion at equal %HA from road noise level to wtn
#' lconvert(Lref = 70, ref ="road", type = "wtn",
#'                 annoyance = "HA")
#'
#' ## conversion back from wtn to road, in output also the coefficients 
#' lconvert(Lref = 48.91, ref ="wtn", type = "road",
#'                  annoyance = "HA", coeff = TRUE)
#' 
#' # plotting conversion curve road-wtn
#' Lroad <- seq(40,80,1)
#' Lwtn <- lconvert(Lroad, type = "wtn")
#' 
#' plot(Lwtn~Lroad)
#' 
#' ## conversion curves for different sources, with road as reference
#' Lroad <- seq(40,80,1)
#' Lwtn <- lconvert(Lroad, type = "wtn")
#' Lrailway <- lconvert(Lroad, type = "railway")
#' Laircraft <- lconvert(Lroad, type = "aircraft")
#' 
#' plot(Lrailway~Lroad, ylim=c(30,80),
#'      ylab = "Lconv [dB(A)]", xlab= "Lref [dB(A)]")
#' lines(Lwtn~Lroad, lty = 2)
#' lines(Laircraft~Lroad, lty =3)
#'
#' @export

lconvert<- function( Lref, ref = "road", type,
                            annoyance = "HA", C = 72, coeff=FALSE)
{
    param <- choose_param(type = type);
    param_ref <- choose_param(type = ref)

    C <- choose_cutoff(annoyance=annoyance, C=C)
    
    aR <- param_ref$intercept
    bR <- param_ref$slope
    sR <- param_ref$tot_sd

    aW <- param$intercept
    bW <- param$slope
    sW <- param$tot_sd

    ## conversion coefficients
    A <- ( (sW*aR-sR*aW) + (sR-sW)*C )/ (bW*sR)
    B <- (bR*sW)/(bW*sR)
    Lconv <- A+B*Lref

    if(coeff == FALSE)
        return(Lconv)
    else
        return(list(
            Lconv=Lconv,
            A=A,
            B=B)
    )
}

## conversion at equal %HA from road noise level to wtn
lconvert(Lref = 70, ref ="road", type = "wtn",
                 annoyance = "HA")

## conversion back from wtn to road, in output also the coefficients 
lconvert(Lref = 48.91, ref ="wtn", type = "road",
                 annoyance = "HA", coeff = TRUE)

# plotting conversion curve road-wtn
Lroad <- seq(40,80,1)
Lwtn <- lconvert(Lroad, type = "wtn")

plot(Lwtn~Lroad)

## conversion curves for different sources, with road as reference
Lroad <- seq(40,80,1)
Lwtn <- lconvert(Lroad, type = "wtn")
Lrailway <- lconvert(Lroad, type = "railway")
Laircraft <- lconvert(Lroad, type = "aircraft")

plot(Lrailway~Lroad, ylim=c(30,80),
     ylab = "Lconv [dB(A)]", xlab= "Lref [dB(A)]")
lines(Lwtn~Lroad, lty = 2)
lines(Laircraft~Lroad, lty =3)


