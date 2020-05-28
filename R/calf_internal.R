calf_internal <- function(data,
                          nMarkers,
                          randomize  = FALSE,
                          proportion = NULL,
                          times,
                          targetVector = "binary",
                          optimize = "pval",
                          verbose = FALSE){
  
  
  # getting rid of global variable warning -------------------------- #
  x = NULL
  y = NULL
  refx = NULL
  refy = NULL
  
  if (targetVector == "nonbinary")
      optimize <- NULL
  
  # setting up some initial values -----------------------------------#
  if (any(apply(data, 2, is.numeric) == FALSE)) {
    stop("CALF ERROR: Data are not numeric. Please check that data were read in correctly.")
  }
  
  
  nVars <- ncol(data) - 1
  dNeg  <- data[ ,2:ncol(data)]
  dNeg  <- dNeg * - 1
  data  <- data.frame(data, dNeg, check.names = FALSE)
  
  if (nMarkers > nVars){
    stop(paste0("CALF ERROR: Requested number of markers is larger than the number of markers in data set. ",
                "Please revise this value or make sure your data were read in properly."))
  }
  
  if (randomize == TRUE) data[ ,1] <- sample(data[ ,1])
  
  if (!is.null(proportion)) {
    
    if (targetVector == "binary") {
      
      ctrlRows  <- which(data[ ,1] == 0)
      caseRows  <- which(data[ ,1] == 1)
      
      # calculate number of case and control to keep
      
      if(length(proportion) == 2) {
        
        nCtrlKeep <- round(length(ctrlRows)*proportion[1], digits = 0)
        nCaseKeep <- round(length(caseRows)*proportion[2], digits = 0)
        
      } else {
        
        nCtrlKeep <- round(length(ctrlRows)*proportion, digits = 0)
        nCaseKeep <- round(length(caseRows)*proportion, digits = 0)
        
      }
      
      # sample randomly rows of case and control to keep, record rows to keep
      keepRows  <- c(sample(ctrlRows)[1:nCtrlKeep], sample(caseRows)[1:nCaseKeep])
      
      # subset original data to keep these rows
      data      <- data[keepRows, ]
      
    } else {
      
      nDataKeep <- round(nrow(data)*proportion, digits = 0)
      keepRows  <- sample(1:nrow(data))[1:nDataKeep]
      data      <- data[keepRows, ]
      
    }
    
  }
  
  real  <- data[ ,1]
  realMarkers <- data[ , 2:ncol(data)]
  ctrl  <- data[data[ ,1] == 0, 2:ncol(data)]
  case  <- data[data[ ,1] == 1, 2:ncol(data)]
  indexNegPos <- rep(0, (nVars*2))
  # end of setting up some initial values ----------------------------#
  
  # initial loop to establish first optimal marker -------------------#
  allCrit <- numeric()
  for (i in 1:(nVars*2)){
    if (targetVector == "binary"){
      caseVar    <- case[ ,i]
      ctrlVar    <- ctrl[ ,i]
      if (optimize == "pval"){
        crit       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
      } else if (optimize == "auc"){
        crit <- compute.auc(caseVar, ctrlVar)
        crit <- 1/crit
      }
    } else {
      realVar <- realMarkers[ ,i]
      crit    <- suppressWarnings(cor(real, realVar, use = "complete.obs"))
      crit    <- 1/crit
    }
    allCrit[i] <- crit
  }
  allCrit[allCrit < 0] <- NA
  
  # end of initial loop ----------------------------------------------#
  
  keepMarkers  <- names(realMarkers)[which.min(allCrit)]
  bestCrit     <- min(allCrit, na.rm = TRUE)
  keepIndex    <- which.min(allCrit)
  
  if (verbose == TRUE) {
    
    if(targetVector == "binary") {
      
      if (optimize == "pval"){
        cat("Selected:", keepMarkers,
            paste0("p value = ", round(bestCrit, digits = 15), "\n"))
      } else if (optimize == "auc"){
        cat("Selected:", keepMarkers,
            paste0("AUC = ", round((1/bestCrit), digits = 15), "\n"))
      }
      
    } else if (targetVector == "nonbinary") {
        cat("Selected:", keepMarkers,
            paste0("Correlation = ", round((1/bestCrit), digits = 15), "\n"))
    }
  }
  
  # second loop to add another marker --------------------------------#
  if (nMarkers != 1){
    allCrit  <- numeric()
    realPrev <- realMarkers[ ,keepIndex]
    casePrev <- case[ ,keepIndex]
    ctrlPrev <- ctrl[ ,keepIndex]
    for (i in 1:(nVars*2)){
      if (i != keepIndex){
        caseVar <- casePrev + case[ ,i]
        ctrlVar <- ctrlPrev + ctrl[ ,i]
        realVar <- realPrev + realMarkers[ ,i]
        if (targetVector == "binary"){
          if (optimize == "pval"){
            crit       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
          } else if (optimize == "auc"){
            crit <- compute.auc(caseVar, ctrlVar)
            crit <- 1/crit
          }
        } else {
          crit <- suppressWarnings(cor(real, realVar, use = "complete.obs"))
          crit <- 1/crit
        }
      } else {
        crit <- NA
      }
      allCrit[i] <- crit
    }
    # end of second loop ----------------------------------------------#
    
    allCrit[allCrit < 0] <- NA
    
    # check if the latest p is lower than the previous p               #
    continue <- ifelse(bestCrit[length(bestCrit)] > min(allCrit, na.rm = TRUE), TRUE, FALSE)
    
    
    if (continue == TRUE){
      keepMarkers  <- append(keepMarkers, names(realMarkers)[which.min(allCrit)])
      bestCrit     <- append(bestCrit, min(allCrit, na.rm = TRUE))
      keepIndex    <- append(keepIndex, which.min(allCrit))
      
      if (length(keepMarkers) == nMarkers) continue <- FALSE
    }
    
    if (verbose == TRUE) {
      
      if(targetVector == "binary") {

        if (optimize == "pval"){
          cat("Selected:", keepMarkers[length(keepMarkers)],
              paste0("p value = ", round(bestCrit[length(bestCrit)], digits = 15), "\n"))
        } else if (optimize == "auc"){
          cat("Selected:", keepMarkers[length(keepMarkers)],
              paste0("AUC = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
        }
        
      } else if (targetVector == "nonbinary") {
        cat("Selected:", keepMarkers[length(keepMarkers)],
            paste0("Correlation = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
      }
    
    }
    
    # loop for third through nMarker ----------------------------------#
    while (continue == TRUE){
      allCrit  <- numeric()
      casePrev <- rowSums(case[ ,keepIndex], na.rm = TRUE)
      ctrlPrev <- rowSums(ctrl[ ,keepIndex], na.rm = TRUE)
      realPrev <- rowSums(realMarkers[ ,keepIndex], na.rm = TRUE)
      for (i in 1:(nVars*2)){
        if (!(i %in% keepIndex)){
          caseVar <- casePrev + case[ ,i]
          ctrlVar <- ctrlPrev + ctrl[ ,i]
          realVar <- realPrev + realMarkers[ ,i]
          if (targetVector == "binary"){
            if (optimize == "pval"){
              crit       <- t.test(caseVar, ctrlVar, var.equal = FALSE)$p.value
            } else if (optimize == "auc"){
              crit <- compute.auc(caseVar, ctrlVar)
              crit <- 1/crit
            }
          } else {
            crit <- suppressWarnings(cor(real, realVar, use = "complete.obs"))
            crit <- 1/crit
          }
        } else {
          crit <- NA
        }
        allCrit[i] <- crit
      }
      allCrit[allCrit < 0] <- NA
      
      continue <- ifelse(bestCrit[length(bestCrit)] > min(allCrit, na.rm = TRUE),
                         TRUE, FALSE)

      if (continue == TRUE){
        keepMarkers  <- append(keepMarkers, names(realMarkers)[which.min(allCrit)])
        bestCrit     <- append(bestCrit, min(allCrit, na.rm = TRUE))
        keepIndex    <- append(keepIndex, which.min(allCrit))
        continue     <- bestCrit[length(bestCrit)] < bestCrit[length(bestCrit)-1]
        if (verbose == TRUE) {
          
          if(targetVector == "binary") {
            if (optimize == "pval"){
              cat("Selected:", keepMarkers[length(keepMarkers)],
                  paste0("p value = ", round(bestCrit[length(bestCrit)], digits = 15), "\n"))
            } else if (optimize == "auc"){
              cat("Selected:", keepMarkers[length(keepMarkers)],
                  paste0("AUC = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
            }
          } else if (targetVector == "nonbinary") {
            cat("Selected:", keepMarkers[length(keepMarkers)],
                paste0("Correlation = ", round((1/bestCrit[length(bestCrit)]), digits = 15), "\n"))
          }
          
        }
      }
      
      if (length(keepMarkers) == nMarkers) continue <- FALSE
    }
  }
  
  if (verbose == TRUE) cat("\n")
  
  indexNegPos[keepIndex] <- ifelse(keepIndex > nVars, -1, 1)
  finalIndex   <- ifelse(keepIndex <= nVars, keepIndex, keepIndex - nVars)
  finalMarkers <- data.frame(names(case)[finalIndex], indexNegPos[keepIndex], check.names = FALSE)
  names(finalMarkers) <- c("Marker","Weight")
  
  if (targetVector == "nonbinary" || optimize == "auc") {
    finalBestCrit <- 1 / bestCrit[length(bestCrit)]
  } else {
    finalBestCrit <- bestCrit[length(bestCrit)]
  }
  ## AUC -------------------------------------------------------------#
  # create function value for each individual
  if (targetVector == "binary"){
    if (nMarkers != 1 & length(keepIndex) != 1){
      funcValue   <- c(rowSums(case[,c(keepIndex)]), rowSums(ctrl[,c(keepIndex)]))
    } else {
      funcValue   <- c(case[,c(keepIndex)], ctrl[,c(keepIndex)])
    }
    funcValue <- round(funcValue, digits = 8)
    # rank individual function values
    ranks       <- rank(funcValue, ties.method = "average")
    seqCaseCtrl <- c(rep(1, nrow(case)), rep(0, nrow(ctrl)))
    
    # set up plot -----------------------------------------------------#
    all <- data.frame(funcValue,
                      seqCaseCtrl,
                      ranks)
    all <- all[order(all$ranks),]
    all$refx <- seq(0,1,1/(nrow(all)-1))
    all$refy <- seq(0,1,1/(nrow(all)-1))
    initVal  <- all$seqCaseCtrl[1]
    moveRight <- ifelse(initVal == 0, nrow(case), nrow(ctrl))
    moveUp    <- ifelse(initVal == 0, nrow(ctrl), nrow(case))
    # moveLeft
    for (i in 2:nrow(all)){
      all$x[1] <- 0
      all$y[1] <- 0
      if (all$seqCaseCtrl[i] == initVal){
        all$x[i] = all$x[i-1]
        all$y[i] = all$y[i-1] + 1/(moveUp-1)
      } else {
        all$x[i] = all$x[i-1] + 1/(moveRight)
        all$y[i] = all$y[i-1]
      }
    }
    
    # if the plot prints upside-down, switch values for
    # x and y
    n <- round(length(all$refy)/2, digits = 0)
    if (all$refy[n] > all$y[n]){
      all$a <- all$x
      all$b <- all$y
      all$x <- all$b
      all$y <- all$a
    }
    
    rocPlot <- ggplot(all, aes(x = x, y = y)) +
      geom_line(size = 1) +
      geom_line(aes(x = refx, y = refy, colour = "red"), size = 1.5) +
      scale_x_continuous(limits = c(0,1)) +
      theme_bw() +
      theme(legend.position = "none") +
      ylab("True Positive Rate (Sensitivity)") +
      xlab("False Positive Rate (1 - Specificity)")
    # set up plot -----------------------------------------------------#
    
    # compute arguments for AUC
    caseFunc  <- sum(ranks[1:nrow(case)]) - nrow(case)*(nrow(case)+1)/2
    ctrlFunc  <- sum(ranks[(nrow(case)+1):length(ranks)]) - nrow(ctrl)*(nrow(ctrl)+1)/2
    # compute AUC
    auc       <- round(max(ctrlFunc, caseFunc)/(caseFunc + ctrlFunc), digits = 4)
  } else {
    auc     <- NULL
    rocPlot <- NULL
  }
  est       <- list(selection  = finalMarkers,
                    auc        = auc,
                    randomize  = randomize,
                    proportion = proportion,
                    targetVec  = targetVector,
                    rocPlot    = rocPlot,
                    finalBest  = finalBestCrit,
                    optimize   = optimize)
  class(est) <- "calf"
  return(est)
}

compute.auc <- function(caseVar, ctrlVar){
  funcValue <- c(caseVar, ctrlVar)
  funcValue <- round(funcValue, digits = 8)
  ranks     <- rank(funcValue, ties.method = "average")
  caseFunc  <- sum(ranks[1:length(caseVar)]) - length(caseVar)*(length(caseVar)+1)/2
  ctrlFunc  <- sum(ranks[(length(caseVar)+1):length(ranks)]) - length(ctrlVar)*(length(ctrlVar)+1)/2
  auc       <- round(max(ctrlFunc, caseFunc)/(caseFunc + ctrlFunc), digits = 4)
  return(auc)
}