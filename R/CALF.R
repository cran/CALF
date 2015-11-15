#'@title calf
#'@description Coarse approximation linear function
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable; all other columns contain relevant markers
#'@param nMarkers Maximum number of markers to include in creation of sum
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@examples
#'calf(data = CaseControl, nMarkers = 6)
#'@export
calf <- function(data, nMarkers){
  # setting up some initial values -----------------------------------#
  nVars <- ncol(data) - 1
  dNeg  <- data[ ,2:ncol(data)]
  dNeg  <- dNeg * - 1
  data  <- data.frame(data, dNeg)
  ctrl  <- data[data$case == 0, 2:ncol(data)]
  case  <- data[data$case == 1, 2:ncol(data)]
  indexNegPos <- rep(0, (nVars*2))
  # end of setting up some initial values ----------------------------#

  # initial loop to establish first optimal marker -------------------#
  allP <- numeric()
  for (i in 1:(nVars*2)){
    caseVar    <- case[ ,i]
    ctrlVar    <- ctrl[ ,i]
    p       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
    allP[i] <- p
  }
  # end of initial loop ----------------------------------------------#

  keepMarkers  <- names(case)[which.min(allP)]
  bestP        <- min(allP, na.rm = TRUE)
  keepIndex    <- which.min(allP)
  # second loop to add another marker --------------------------------#

  allP <- numeric()
  casePrev <- case[ ,keepIndex]
  ctrlPrev <- ctrl[ ,keepIndex]
  for (i in 1:(nVars*2)){
    if (i != keepIndex){
      caseVar <- casePrev + case[ ,i]
      ctrlVar <- ctrlPrev + ctrl[ ,i]
      p       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
    } else {
      p <- 1
    }
    allP[i] <- p
  }
  # end of second loop ----------------------------------------------#

  keepMarkers  <- append(keepMarkers, names(case)[which.min(allP)])
  bestP        <- append(bestP, min(allP, na.rm = TRUE))
  keepIndex    <- append(keepIndex, which.min(allP))

  # check if the latest p is lower than the previous p               #
  continue <- bestP[length(bestP)] < bestP[length(bestP)-1]

  # loop for third through nMarker ----------------------------------#

  while (continue == TRUE){
    allP     <- numeric()
    casePrev <- rowSums(case[ ,keepIndex], na.rm = TRUE)
    ctrlPrev <- rowSums(ctrl[ ,keepIndex], na.rm = TRUE)
    for (i in 1:(nVars*2)){
      if (!(i %in% keepIndex)){
        caseVar <- casePrev + case[ ,i]
        ctrlVar <- ctrlPrev + ctrl[ ,i]
        p       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
      } else {
        p <- 1
      }
      allP[i] <- p
    }
    keepMarkers  <- append(keepMarkers, names(case)[which.min(allP)])
    bestP        <- append(bestP, min(allP, na.rm = TRUE))
    keepIndex    <- append(keepIndex, which.min(allP))
    continue     <- bestP[length(bestP)] < bestP[length(bestP)-1]
    # stop the search when it hits the max number of markers
    if (length(keepMarkers) == nMarkers) continue <- FALSE
  }

  indexNegPos[keepIndex] <- ifelse(keepIndex >= nVars, -1, 1)
  finalIndex   <- ifelse(keepIndex <= nVars, keepIndex, keepIndex - nVars)
  finalMarkers <- data.frame(names(case)[finalIndex], indexNegPos[keepIndex])
  names(finalMarkers) <- c("Marker","Weight")
  print.data.frame(finalMarkers, row.names = FALSE)
}


