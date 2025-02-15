#' Automated Estimation of the Linear Mixed Model with Splines#'
#'
#' The nlive.splines function allows to fit a Linear Mixed Models with the function of time approximated
#' with natural cubic splines or B-splines in the context of longitudinal Gaussian outcomes.
#' This function was designed to be intuitive enough to the less sophisticated users, while using
#' the existing hlme() function from the lcmm R package as well as the existing ns() and bs() functions from the
#' splines  R package.
#' #'
#' CAUTIONS REGARDING THE USE OF THE FUNCTION
#'
#' traj.marg: if "TRUE", this argument automatically plots the estimated marginal trajectories of the longitudinal outcome
#' for the most common profile of covariates, if any (i.e., ref "1" for binary variables and mean values for continuous variables).
#' Thus, users must ensure that continuous variables are centered on the mean.
#'
#'
#' @param dataset data frame containing the variables ID, outcome, time, var.all, and all other var. arguments.
#' @param ID name of the variable representing the grouping structure specified with " (e.g., "ID" representing the unique identifier of participants).
#' @param time name of the variable representing the timescale specified with " (e.g., "time"), which can be negative or positive.
#' @param formula two-sided linear formula object for the fixed-effects in the linear mixed model. The response outcome is on the left of ~ and the covariates are separated by + on the right of ~. By default, an intercept is included. If no intercept, -1 should be the first term included on the right of ~.
#' @param random optional one-sided formula for the random-effects in the linear mixed model. Covariates with a random-effect are separated by +. Default is an intercept. If no intercept, -1 should be the first term included.
#' @param splines "ns" for Natural Cubic Splines. "bs" for Cubic B-splines.
#' @param knots inner knots that define the spline. Typical values are the mean or median for one knot, quantiles for more knots. Default is two equally spaced inner knots (i.e., 33th, 66th percentiles).
#' @param Boundary.knots boundary points at which to impose the natural boundary conditions and anchor the spline basis. Default to the range of the data.
#' @param traj.marg optional logical indicating if the marginal estimated trajectory should be plotted for the most common profile of covariates, if any. Default to FALSE.
#'
#' #' @return An object of from the existing _lcmm_ R package containing the results of the fit of the data
#' by a linear mixed model and with the function of time approximated using natural cubic splines.
#'
#' @author Maude Wagner, Ana W. Capuano, Cecile Proust-lima
#'
#' \email{maude_wagner@@rush.edu}
#'
#' @references
#' Proust-Lima, C., Philipps, V., & Liquet, B. (2017). Estimation of Extended Mixed Models Using Latent Classes and Latent Processes: The R Package lcmm. Journal of Statistical Software, 78(2), 1–56. https://doi.org/10.18637/jss.v078.i02
#' Perperoglou, A., Sauerbrei, W., Abrahamowicz, M. et al. A review of spline function procedures in R. BMC Med Res Methodol 19, 46 (2019). https://doi.org/10.1186/s12874-019-0666-3
#'
#' @examples
#'
#' #### Fitting a linear mixed model with ns splines and 2 inner knots (33th, 66th percentiles)
#' \dontrun{
#' head(dataCog)
#' lmm.ns.fit = nlive.splines(dataCog, ID="ID", outcome="cognition", time="time", splines.type="ns")
#' }
#' #### plot(lmm.ns.fit): diagnostic plots to assess the goodness-of-fit of lmm.ns.fit
#'
#'
#' @import lcmm
#' @import splines
#' @import knitr
#' @import dplyr
#' @import viridis
#' @import sqldf
#' @importFrom lcmm hlme
#' @importFrom splines ns
#' @importFrom splines bs
#' @importFrom dplyr %>%
#' @importFrom dplyr group_by
#' @importFrom dplyr row_number
#' @importFrom dplyr sample_n
#' @importFrom dplyr filter
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 geom_line
#' @importFrom ggplot2 scale_x_continuous
#' @importFrom ggplot2 ggtitle
#' @importFrom ggplot2 labs
#' @importFrom ggplot2 guides
#' @importFrom graphics points
#' @importFrom graphics legend
#' @importFrom graphics plot
#' @importFrom graphics lines
#' @importFrom graphics par
#' @importFrom Rmpfr pmax
#' @importFrom stats quantile
#' @importFrom stats as.formula
#' @importFrom stats na.omit
#' @importFrom sqldf sqldf
#' @importFrom fastDummies dummy_cols
#'
#' @export
########


nlive.splines <- function(dataset, ID, time, formula, random,
                          splines = NULL,
                          knots = NULL,
                          Boundary.knots = NULL,
                          traj.marg = FALSE){

  if (missing(dataset))
    stop("The argument dataset should be specified and defined as a data.frame")
  if (nrow(dataset) == 0)
    stop("Data should not be empty")
  if (missing(ID))
    stop("The argument ID must be specified in any model")
  if (!is.numeric(dataset[[ID]]))
    stop("The argument ID must be numeric")
  if (missing(time))
    stop("The argument time must be specified in any model")
  if (!is.numeric(dataset[[time]]))
    stop("The argument time must be numeric")
  if (missing(formula))
    stop("The argument formula must be specified in any model")
  if (missing(random))
    stop("The argument random must be specified in any model")

  ## dataset
  dataset$ID       = na.omit(dataset[,ID])
  dataset$time     = na.omit(dataset[,time])

  ######################################
  ######      SPLINES BASIS       ######
  ######################################

  # Defaults #
  if (is.null(splines)        == TRUE){splines = "ns"}
  if (is.null(knots)          == TRUE){knots   = c(0.33,0.66)}
  if (is.null(Boundary.knots) == TRUE){Boundary.knots=c(0.05,0.95)}

  # knots + boundaries #
  knots.arg = c(quantile(dataset$time, probs = knots, na.rm=T))
  Boundary.knots.arg = c(quantile(dataset$time, probs = Boundary.knots, na.rm=T))

  # Definition of the splines basis #
  if (splines == "ns"){

  # ns()
  NS = ns(dataset$time, knots = knots.arg, Boundary.knots = Boundary.knots.arg)
  dataset[,((ncol(dataset)+1):(ncol(dataset)+ncol(NS)))]  = NS
  colnames(dataset) = c(colnames(dataset)[1:(ncol(dataset)-ncol(NS))],paste("NS", 1:ncol(NS), sep = ""))
  pos.ns = ncol(dataset) - dim(NS)[2] + 1

  # Formula + Random
  formula.splines = as.formula(gsub(time, paste("(",paste0("NS",1:ncol(NS), collapse = "+"),paste(")")), formula))
  random.splines  = as.formula(gsub(time, paste("(",paste0("NS",1:ncol(NS), collapse = "+"),paste(")")), random))

  }

  else if (splines == "bs") {

  # bs()
  BS = bs(dataset$time, knots = knots.arg)
  dataset[,((ncol(dataset)+1):(ncol(dataset)+ncol(BS)))]  = BS
  colnames(dataset) = c(colnames(dataset)[1:(ncol(dataset)-ncol(BS))],paste("BS", 1:ncol(BS), sep = ""))
  pos.bs = ncol(dataset) - dim(BS)[2] + 1

  # Formula + Random
  formula.splines = as.formula(gsub(time, paste("(",paste0("BS",1:ncol(BS), collapse = "+"),paste(")")), formula))
  random.splines  = as.formula(gsub(time, paste("(",paste0("BS",1:ncol(BS), collapse = "+"),paste(")")), random))

  }

  #print(formula)
  #print(formula.splines)
  # Fitting the models
  ptm           = proc.time()
  model.splines = hlme(formula.splines, random = random.splines, subject = ID, data = dataset)
  cost          = proc.time() - ptm



  #############################################
  ######        SPLINES BASIS PLOT       ######
  #############################################

    datnew = data.frame(time = seq(min(dataset$time), max(dataset$time), length = 100))
    par(mfrow=c(1,1))
    if (splines == "ns"){

    NS = ns(datnew$time, knots = knots.arg, Boundary.knots = Boundary.knots.arg)
    datnew[,((ncol(datnew)+1):(ncol(datnew)+ncol(NS)))] = NS
    colnames(datnew) = c(colnames(datnew)[1:(ncol(datnew)-ncol(NS))],paste("NS", 1:ncol(NS), sep = ""))
    #
    pos.ns = ncol(datnew) - dim(NS)[2] + 1
    min = min(pmin(datnew[,c(pos.ns:ncol(datnew))]))
    max = max(pmax(datnew[,c(pos.ns:ncol(datnew))]))
    #
    plot(datnew$NS1~datnew$time, type='l', lwd=3, col=viridis(dim(NS)[2])[1], las=1, ylim = c(min, max), main = "Natural cubic spline basis", xlab= time, ylab="Value")
    for (j in 1:ncol(NS)-1) lines(datnew[,pos.ns+j] ~ datnew$time, lwd=3, col=viridis(dim(NS)[2])[1+j])

    } else if (splines == "bs"){

      BS = bs(datnew$time, knots = knots.arg)
      datnew[,((ncol(datnew)+1):(ncol(datnew)+ncol(BS)))]  = BS
      colnames(datnew) = c(colnames(datnew)[1:(ncol(datnew)-ncol(BS))],paste("BS", 1:ncol(BS), sep = ""))
      #
      pos.bs = ncol(datnew) - dim(BS)[2] + 1
      min = min(pmin(datnew[,c(pos.bs:ncol(datnew))]))
      max = max(pmax(datnew[,c(pos.bs:ncol(datnew))]))
      #
      plot(datnew$BS1~datnew$time, type='l', lwd=3, col=viridis(dim(BS)[2])[1], las=1, ylim = c(min, max), main = "B-spline basis", xlab=time, ylab="Value")
      for (j in 1:ncol(BS)-1) lines(datnew[,pos.bs+j] ~ datnew$time, lwd=3, col=viridis(dim(BS)[2])[1+j])

    }


  #############################################################
  ######   PLOT Individual observations vs. predictions  ######
  #############################################################
  tab_pred = cbind(as.data.frame(model.splines$pred), dataset$time)
  #
    par(mfrow=c(2,3),mar = c(5.1, 4.1, 4.1, 2.1), xpd=F)
  for (i in 1:6){
    tempo = subset(tab_pred, ID %in% unique(ID)[i])
    plot(tempo$`dataset$time`, tempo$pred_ss,
         ylim = c(min(tempo$obs,tempo$pred_ss),max(tempo$obs,tempo$pred_ss)),
         xlab = time,
         ylab = "outcome",
         main = unique(tempo$ID),
         type = "l",
         col = viridis(7)[i], lwd = 4, las = 1)
    points(tempo$`dataset$time`, tempo$obs, las = 1, cex = 1.5)
    i = i + 1
    legend("bottomleft", legend=c("individual obs","individual pred"), cex = 0.9,
           pch = c(1,NA), lty=c(0,1), lwd = c(1,2), col = c(1,viridis(7)[i]),
           bty = "y")
  }

  for (i in 1:6){
    tempo = subset(tab_pred, ID %in% unique(ID)[i])
    plot(tempo$`dataset$time`, tempo$pred_ss,
         ylim = c(min(tempo$obs,tempo$pred_ss, tempo$pred_m),max(tempo$obs,tempo$pred_ss, tempo$pred_m)),
         xlab = time,
         ylab = "outcome",
         main = unique(tempo$ID),
         type = "l",
         col = viridis(7)[i], lwd = 4, las = 1)
    points(tempo$`dataset$time`, tempo$obs, las = 1, cex = 1.5)
    points(tempo$`dataset$time`, tempo$pred_m, lty = 2, type = "l", las = 1, cex = 1.5)
    legend("bottomleft", legend=c("individual obs","individual pred","marginal pred"), cex = 0.9,
           pch = c(1,NA,NA), lty=c(0,1,2), lwd = c(1,2,1), col = c(1,viridis(7)[i],1),
           bty = "y")
  }


  # output
  print(model.splines)
  cat("----------------------------------------------------\n The program took", round(cost[3],2), "seconds \n")
  return(model.splines)
}

