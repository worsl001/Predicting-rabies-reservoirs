MatthewsCorrelationCoefficient <- function(confusionMat) {
  # Compute Matthews Correlation Coefficient
  #
  # Args:
  #   confMat: confusion matrix
  #
  # Returns:
  #   MCC 
  #
  TN <- confusionMat[1,1]
  FN <- confusionMat[1,2]
  FP <- confusionMat[2,1]
  TP <- confusionMat[2,2]
  
  num <- as.numeric((TP * TN) - (FP * FN))
  den1 <- as.numeric((TP + FP) * (TP + FN))
  den2 <- as.numeric((TN + FP) * (TN + FN))
  
  MCC <- num/sqrt(den1*den2)
  
}

#---

Confusion2Performance <- function(confusionMat,verbose=FALSE) {
  # Extract model performance from confusionMatrix structure returned by caret
  # package function confusionMatrix(). 
  #
  # MCC, specificity, sensitivity, accuracy and confidence levels for accuracy
  # are returned in a list. 
  #
  # Args:
  #   confusionMat: confusion matrix returned by confusionMatrix() function from
  # caret package
  #
  # Returns:
  #   A list with elements 'prediction' and 'accuracy'
  # 
  # 
  #
  if(missing(confusionMat)) {
    stop("A confusioMatrix class object must be given as parameter for this function.")
  }
  
  if (verbose) {
    print(confusionMat$table)
  }
  
  MCC <- MatthewsCorrelationCoefficient(confusionMat$table)
  SEN <- confusionMat$byClass[1]*100
  SPE <- confusionMat$byClass[2]*100
  ACC <- confusionMat$overall[1]*100
  ACC.low <- confusionMat$overall[3]*100
  ACC.up <- confusionMat$overall[4]*100

  performance <- cbind(ACC,SPE,SEN,MCC)
  colnames(performance) <- c("ACC","SPE","SEN","MCC")
  row.names(performance) <- NULL
  
  accuracy.ci <- cbind (ACC.low, ACC.up)
  colnames(accuracy.ci) <- c("AccuracyLower","AccuracyUpper")
  row.names(accuracy.ci) <- NULL

  
  if (verbose) {
    print(performance)
    print(accuracy.ci)
  }
  
  output <- list(performance,accuracy.ci)
  names(output)<-c("performance","accuracy")

  return(output)
}

#---

EstimatePerformanceCV <- function(oof,repeats) {
  # Computes the average confusion matrix and the standard deviation from the output
  # of a (repeated) cross-validation
  #
  # Args:
  #   oof: predictions for the out-of-fold data (.$pred element from the object 
  # returned by train() function when option savePredictions=TRUE)
  #   repeats: number of repeats of the cross-validation
  #
  # Returns:
  #   A list with the average confusion matrix and the standard deviation when repeated
  # cross-validation is performed.
  # 
  
  #
  if(missing(oof)) {
    stop("The predictions for the out-of-fold data must be given as parameter for this function.")
  }
  
  if(missing(repeats)) {
    stop("It is necessary to inform the number of repetitions for the cross-validation algorithm")
  }
  
  folds <- unique(oof$Resample)
  all.confusion <- NULL
  all.performance <- NULL
  
  #fix repeats names to index oof$Resample vector
  if (repeats >= 10) {
    allrepeats <- 1:repeats
    allrepeats[1:9] <- paste("0",allrepeats[1:9],sep="")
  }
  
  for(ii in 1:repeats) {
    
    resamp <- oof[is.element(oof$Resample,folds[grep(paste("Rep",allrepeats[ii],sep=""),folds)]),]
    oof.conf <- confusionMatrix(resamp$pred,resamp$obs,positive="Positive")
    all.confusion <- c(all.confusion,list(oof.conf$table))
    oof.performance <- Confusion2Performance(oof.conf)
    #print(oof.conf)
    #print(oof.performance)
    all.performance <- c(all.performance,list(oof.performance$performance))

  }
  names(all.confusion) <- paste("Rep",c(1:repeats),sep="")
  
  TP <- NULL
  FP <- NULL
  TN <- NULL
  FN <- NULL
  
  ACC <- NULL
  SPE <- NULL
  SEN <- NULL
  MCC <- NULL
  
  for (ii in 1:length(all.confusion)) {
    TP <- c(TP, all.confusion[[ii]]["Positive","Positive"])
    FN <- c(FN, all.confusion[[ii]]["Negative","Positive"])
    TN <- c(TN, all.confusion[[ii]]["Negative","Negative"])
    FP <- c(FP, all.confusion[[ii]]["Positive","Negative"])
    
    ACC <- c(ACC,all.performance[[ii]][,"ACC"])
    SEN <- c(SEN,all.performance[[ii]][,"SEN"])
    SPE <- c(SPE,all.performance[[ii]][,"SPE"])
    MCC <- c(MCC,all.performance[[ii]][,"MCC"])
  }
  
  # compute mean and deviation for confusion matrices
  TP.mean <- mean(TP)
  TP.sd <- sd(TP)
  if (is.na(TP.sd)) {TP.sd <- 0}
  
  FP.mean <- mean(FP)
  FP.sd <- sd(FP)
  if (is.na(FP.sd)) {FP.sd <- 0}
  
  TN.mean <- mean(TN)
  TN.sd <- sd(TN)
  if (is.na(TN.sd)) {TN.sd <- 0}
  
  FN.mean <- mean(FN)
  FN.sd <- sd(FN)
  if (is.na(FN.sd)) {FN.sd <- 0}

  mean.confusion <- matrix(c(TN.mean,FN.mean,FP.mean,TP.mean),byrow=TRUE,nrow=2,ncol=2,dimnames=list(c("Pred Negative","Pred Positive"),c("Ref Negative","Ref Positive")))
  sd.confusion <- matrix(c(TN.sd,FN.sd,FP.sd,TP.sd),byrow=TRUE,nrow=2,ncol=2,dimnames=list(c("Pred Negative","Pred Positive"),c("Ref Negative","Ref Positive"))) 
  
  # compute mean and deviation for performance metrics
  ACC.mean <- mean(ACC)
  ACC.sd <- sd(ACC)
  SPE.mean <- mean(SPE)
  SPE.sd <- sd(SPE)
  SEN.mean <- mean(SEN)
  SEN.sd <- sd(SEN)
  MCC.mean <- mean(MCC)
  MCC.sd <- sd(MCC)
  
  mean.metrics <- c(ACC.mean, SPE.mean, SEN.mean, MCC.mean)
  names(mean.metrics) <- c("ACC","SPE","SEN","MCC")
  sd.metrics <- c(ACC.sd, SPE.sd, SEN.sd, MCC.sd)
  names(sd.metrics) <- c("ACC","SPE","SEN","MCC")
  
  performance.cv <- list(mean.confusion,sd.confusion,mean.metrics,sd.metrics)
  names(performance.cv) <- c("mean.confusion","deviation.confusion","mean.metrics","deviation.metrics")
  return(performance.cv)
}

#---


EstimateErrorRateCV <- function(oof,repeats) {
  # Computes the average error rate for each class based on (repeated) cross-validation
  #
  # Args:
  #   oof: predictions for the out-of-fold data (.$pred element from the object 
  # returned by train() function when option savePredictions=TRUE)
  #   repeats: number of repeats of the cross-validation
  #
  # Returns:
  #   A list with the average confusion matrix and the standard deviation when repeated
  # cross-validation is performed.
  # 

  if(missing(oof)) {
    stop("The predictions for the out-of-fold data must be given as parameter for this function.")
  }
  
  if(missing(repeats)) {
    stop("It is necessary to inform the number of repetitions for the cross-validation algorithm")
  }
  
  folds <- unique(oof$Resample)
  all.confusion <- NULL
  
  #fix repeats names to index oof$Resample vector
  if (repeats >= 10) {
    allrepeats <- 1:repeats
    allrepeats[1:9] <- paste("0",allrepeats[1:9],sep="")
  }
  
  for(ii in 1:repeats) {
    resamp <- oof[is.element(oof$Resample,folds[grep(paste("Rep",allrepeats[ii],sep=""),folds)]),]
    oof.conf <- confusionMatrix(resamp$pred,resamp$obs,positive="Positive")
    all.confusion <- c(all.confusion,list(oof.conf$table))
    #oof.performance <- Confusion2Performance(oof.conf)
    #print(oof.conf)
    #print(oof.performance)
    #all.performance <- c(all.performance,list(oof.performance$performance))
    
  }
  names(all.confusion) <- paste("Rep",c(1:repeats),sep="")
  
  pos.error <- NULL
  neg.error <- NULL
  
  for (ii in 1:length(all.confusion)) {
    totalExamples <- colSums(all.confusion[[ii]])
    names(totalExamples) <- colnames(all.confusion[[ii]])
    FN <- all.confusion[[ii]]["Negative","Positive"]
    FP <- all.confusion[[ii]]["Positive","Negative"]
    pos.error <- c(pos.error,(FN/totalExamples["Positive"])*100)
    neg.error <- c(neg.error,(FP/totalExamples["Negative"])*100)
  }
  
  
  error.cv <- list(mean(pos.error),sd(pos.error),mean(neg.error),sd(neg.error))
  names(error.cv) <- c("mean.positive","deviation.positive","mean.negative","deviation.negative")
  return(error.cv)
}

#---
