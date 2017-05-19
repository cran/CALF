#'@method print calf
#'@export
print.calf <- function(x, ...){
  if (x$randomize == TRUE) cat("Randomized Output:", "\n", "\n")
  if (!is.null(x$proportion)) cat("Proportion of Data:", x$proportion, "\n", "\n")
  print.data.frame(x$selection, row.names = FALSE, check.names = FALSE)
  if (x$targetVec == "binary") {
    cat("\nAUC:", x$auc)
   if (x$optimize == "pval") cat("\nFinal p-value:", x$finalBest)
  } else {
    cat("\nFinal Correlation:", x$finalBest)
  }
}

#'@method print calf_randomize
#'@export
print.calf_randomize <- function(x, ...){
  if (x$times == 1) {
    cat("Randomized Output Across", x$times, "Replication:", "\n", "\n")
    print.data.frame(x$selection, row.names = FALSE, check.names = FALSE)
    if (x$targetVec == "binary") {
      cat("\nAUC:", x$auc)
      if (x$optimize == "pval") cat("\nFinal p-value:", x$finalBest)
    } else {
      cat("\nFinal Correlation:", x$finalBest)
    }
  } else {
    cat("Randomized Output Across", x$times, "Replications:", "\n", "\n")
    print.data.frame(x$multiple, row.names = FALSE, check.names = FALSE)
    if (x$targetVec == "binary"){
      cat("\n", "\n")
      print.data.frame(x$auc, row.names = F)
    } else {
      cat("\nFinal Correlations:", "\n",
          round(x$finalBest, digits = 4))
    }
  }
}


#'@method print calf_subset
#'@export
print.calf_subset <- function(x, ...){
  if (x$times == 1) {
    cat("Proportion =", x$proportion, "Output Across", x$times, "Replication:", "\n", "\n")
    print.data.frame(x$selection, row.names = FALSE, check.names = FALSE)
    if (x$targetVec == "binary") {
      cat("\nAUC:", x$auc)
      cat("\nFinal p-value:", x$finalBest)
    } else {
      cat("\nFinal Correlation:", x$finalBest)
    }
  } else {
    cat("Proportion =", x$proportion, "Output Across", x$times, "Replications:", "\n", "\n")
    print.data.frame(x$multiple, row.names = FALSE)
    if (x$targetVec == "binary"){
      cat("\n", "\n")
      print.data.frame(x$auc, row.names = F, check.names = FALSE)
    }
  }
}

