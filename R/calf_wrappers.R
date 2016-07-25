#'@title calf
#'@description Coarse approximation linear function
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "real"). All other columns contain relevant markers.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "real" for target vector with real numbers.
#'@param margin Real number from 0 to 1. Indicates the amount a potential marker must improve the target criterion (Pearson correlation or p-value) in order to add the marker.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The AUC value for the classification
#'@return rocPlot. A plot object from ggplot2 for the receiver operating curve.
#'@examples
#'calf(data = CaseControl, nMarkers = 6, targetVector = "binary")
#'@export
calf <- function(data,
                 nMarkers,
                 targetVector,
                 margin){
  calf_internal(data,
                nMarkers,
                proportion = NULL,
                randomize  = FALSE,
                targetVector = targetVector,
                times      = 1,
                margin = NULL)
}


#'@title calf_randomize
#'@description Coarse approximation linear function, randomized
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "real"). All other columns contain relevant markers.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param randomize Logical. Indicate TRUE to randomize the case/control status (or real number vector) for each individual. Used to compare results from true data with results from randomized data.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "real" for target vector with real numbers.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param margin Real number from 0 to 1. Indicates the amount a potential marker must improve the target criterion (Pearson correlation or p-value) in order to add the marker.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The AUC value for the classification
#'@return aucHist A histogram of the AUCs across replications.
#'@examples
#'calf_randomize(data = CaseControl, nMarkers = 6, targetVector = "binary", times = 5)
#'@export
calf_randomize <- function(data,
                           nMarkers,
                           randomize  = TRUE,
                           targetVector,
                           times      = 1,
                           margin     = NULL){
  auc        <- numeric()
  finalBest  <- numeric()
  allMarkers <- character()
  count      <- 1
  AUC = NULL
  repeat {
    out <- calf_internal(data,
                         nMarkers,
                         proportion = NULL,
                         randomize  = randomize,
                         targetVector = targetVector,
                         times,
                         margin = margin)
    auc[count] <- out$auc
    selection  <- out$selection
    markers    <- as.character(out$selection[,1])
    finalBest  <- append(finalBest, out$finalBest)
    allMarkers <- as.character((append(allMarkers, markers)))
    if (count == times) break
    count      <- count + 1
  }

  if (times > 1) {
    summaryMarkers <- as.data.frame(table(allMarkers))
    colnames(summaryMarkers) <- c("Marker", "Frequency")
    summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]
    if (targetVector == "binary"){
    auc            <- as.data.frame(auc)
    colnames(auc)  <- "AUC"
    aucHist <- ggplot(auc, aes(AUC)) +
      geom_histogram() +
      ylab("Count") +
      xlab("AUC") +
      scale_x_continuous() +
      theme_bw()
    } else aucHist = NULL
  } else {
    summaryMarkers = NULL
    aucHist        = NULL
  }
  if (times == 1 & targetVector == "binary") {
    rocPlot <- out$rocPlot
  } else {
    rocPlot <- NULL
  }

  est       <- list(selection  = selection,
                    multiple   = summaryMarkers,
                    auc        = auc,
                    randomize  = randomize,
                    targetVec  = targetVector,
                    aucHist    = aucHist,
                    times      = times,
                    finalBest  = finalBest,
                    rocPlot    = rocPlot)
  class(est) <- "calf_randomize"
  return(est)
}


#'@title calf_subset
#'@description Coarse approximation linear function, randomized
#'@param data Matrix or data frame. First column must contain case/control dummy coded variable (if targetVector = "binary"). Otherwise, first column must contain real number vector corresponding to selection variable (if targetVector = "real"). All other columns contain relevant markers.
#'@param nMarkers Maximum number of markers to include in creation of sum.
#'@param proportion Numeric. A value (where 0 < proportion <= 1) indicating the proportion of cases and controls to use in analysis (if targetVector = "binary"). If targetVector = "real", this is just a proportion of the full sample. Used to evaluate robustness of solution. Defaults to 0.8.
#'@param targetVector Indicate "binary" for target vector with two options (e.g., case/control). Indicate "real" for target vector with real numbers.
#'@param times Numeric. Indicates the number of replications to run with randomization.
#'@param margin Real number from 0 to 1. Indicates the amount a potential marker must improve the target criterion (Pearson correlation or p-value) in order to add the marker.
#'@return A data frame containing the chosen markers and their assigned weight (-1 or 1)
#'@return The AUC value for the classification. If multiple replications are requested, this will be a data.frame containing all AUCs across replications.
#'@return aucHist A histogram of the AUCs across replications.
#'@examples
#'calf_subset(data = CaseControl, nMarkers = 6, targetVector = "binary", times = 5)
#'@export

calf_subset <- function(data,
                        nMarkers,
                        proportion = .8,
                        targetVector,
                        times      = 1,
                        margin = NULL){
  auc        <- numeric()
  allMarkers <- character()
  finalBest  <- numeric()
  count      <- 1
  AUC = NULL
  repeat {
    out <- calf_internal(data,
                         nMarkers,
                         proportion = proportion,
                         randomize  = FALSE,
                         targetVector = targetVector,
                         times,
                         margin = margin)
    auc[count] <- out$auc
    selection  <- out$selection
    finalBest  <- append(finalBest, out$finalBest)
    markers    <- as.character(out$selection[,1])
    allMarkers <- as.character((append(allMarkers, markers)))
    if (count == times) break
    count      <- count + 1
  }

  if (times > 1){
    summaryMarkers <- as.data.frame(table(allMarkers))
    colnames(summaryMarkers) <- c("Marker", "Frequency")
    summaryMarkers <- summaryMarkers[order(-summaryMarkers$Frequency),]
    if (targetVector == "binary"){
    auc            <- as.data.frame(auc)
    colnames(auc)  <- "AUC"
    aucHist <- ggplot(auc, aes(AUC)) +
      geom_histogram() +
      ylab("Count") +
      xlab("AUC") +
      scale_x_continuous() +
      theme_bw()
    } else aucHist = NULL
  } else {
    summaryMarkers = NULL
    aucHist        = NULL
  }
  if (times == 1 & targetVector == "binary") {
    rocPlot <- out$rocPlot
  } else {
    rocPlot <- NULL
  }

  est       <- list(selection  = selection,
                    multiple   = summaryMarkers,
                    auc        = auc,
                    proportion = proportion,
                    targetVec  = targetVector,
                    aucHist    = aucHist,
                    times      = times,
                    finalBest  = finalBest,
                    rocPlot    = rocPlot)
  class(est) <- "calf_subset"
  return(est)
}



