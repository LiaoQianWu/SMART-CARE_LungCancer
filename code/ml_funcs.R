#### Try to implement well-developed package, e.g., 'mlr3' or 'caret', to train ML models
library('randomForest')
library('missForest')
library(xgboost)
library('caret')
library('glmnet')
library('SummarizedExperiment')
library('tidyverse')

# Set plot theme
th <- theme_bw(base_size = 15) +
  theme(axis.title = element_text(face = 'bold'),
        axis.text = element_text(face = 'bold'),
        axis.ticks = element_line(linewidth = 0.8),
        legend.text = element_text(size = 15))


imputeByMF <- function(se) {
  #' Impute missing values using missForest
  #' 
  #' Parameter
  #' se: A SummarizedExperiment object containing a data matrix
  #' 
  #' Return
  #' impuSE: A SummarizedExperiment object containing the imputed data matrix
  
  dat <- t(assay(se))
  impuDat <- missForest(dat, maxiter = 10, verbose = T)$ximp %>%
    t()
  impuSE <- se
  assay(impuSE) <- impuDat
  
  return(impuSE)
}

rmCorrFeats <- function(data, cutoff = 0.8, distance_method = "pearson",
                        cluster_method = "ward.D2") {
  #' Perform agglomerative hierarchical clustering to cluster features into correlated
  #' feature groups and remove most of correlated features in a group from data,
  #' i.e., only one representative of a group is kept
  #' 
  #' Parameter
  #' data: A Feature x Sample data matrix
  #' cutoff: A numeric value specifying where the tree is cut to determine the
  #' number of clusters. The cutoff is applied to the height (distance of a merge)
  #' of the tree and corresponds to feature correlations and should be greater than 0.
  #' distance_method: A character specifying the method for computing a distance
  #' matrix
  #' cluster_method: A character specifying the linkage method for combining clusters
  #' 
  #' Return
  #' A list containing the following components:
  #' reducedData: A matrix where only individual representatives of feature clusters
  #' are kept, i.e., the other correlated features are removed
  #' corrFeats: A list storing clusters of correlated features for latter retrieval
  
  # Check argument
  if (!is(cutoff, 'numeric') | cutoff < 0) {
    stop("Cutoff should be class 'numeric' and greater than 0.")
  }
  
  # Calculate distance matrix
  if (distance_method == 'pearson') {
    distMat <- stats::as.dist(1 - stats::cor(t(data), method = 'pearson',
                                             use = 'pairwise.complete.obs'))
  } else if (distance_method == 'euclidean') {
    distMat <- stats::dist(data, method = 'euclidean')
  } else if (distance_method == 'binary') { #for sparse matrix
    distMat <- stats::dist(data, method = 'binary')
  }
  
  # Perform hierarchical clustering
  hierClust <- stats::hclust(distMat, method = cluster_method)
  clusters <- stats::cutree(hierClust, h = 1 - cutoff)
  reduData <- data[!duplicated(clusters),]
  
  # Record removed correlated feature group for latter retrieval
  corrFeatList <- lapply(rownames(reduData), function(feat) {
    corrFeats <- names(clusters[clusters == clusters[feat]])
    corrFeats[corrFeats != feat]
  })
  names(corrFeatList) <- rownames(reduData)
  
  return(list(reducedData = reduData, corrFeats = corrFeatList))
}

subsetTrainData <- function(x, y, split_method = 'random split', trainSet_ratio = 0.8) {
  #' Subset labeled data into training set in which varied classes are balanced,
  #' and unselected data can be used as test set
  #' 
  #' Parameters
  #' x: A sample x feature matrix containing explanatory variables
  #' y: A response vector containing sample labels that correspond to samples in x
  #' split_method: A character specifying the splitting method to use, which should
  #' be one of 'random split' (default) or 'bootstrap'
  #' trainSet_ratio: A numeric value specifying the ratio of the size of the training
  #' set to that of the total dataset
  #' 
  #' Return
  #' trainIdx: A vector of indices defining the training set and the rest of the
  #' data is used as the test set
  
  # Sanity check
  # stopifnot(nrow(x) == length(y))
  if (!(nrow(x) == length(y))) {
    stop("Row length of x (sample) should be same as length of y.")
  }
  
  # Subset data as training set and balance samples with different labels
  smpIdx <- seq_along(y)
  trainIdx <- lapply(unique(y), function(label) {
    ySub <- smpIdx[y == label]
    if (split_method == 'random split') {
      sample(ySub, size = round(length(ySub)*trainSet_ratio))
    } else if (split_method == 'bootstrap') {
      # Bootstrap data and use OOB samples as test set
      sample(ySub, size = round(length(ySub)*trainSet_ratio), replace = T)
    }
  }) %>% do.call(c, .) %>% sample()
  
  return(trainIdx)
}


runRF <- function(x, y, targetClass, iter = 1, split_method = 'random split',
                  trainSet_ratio = 0.8, ntree = 10000, plot_ROC = F, save_model = F) {
  # This function is currently for binary classification problem. For multi-class
  # classification, data split and bootstrap in random forest has to be reviewed
  # and refined!!!!
  
  #' Perform random forest for binary classification problem. Trained RF model is
  #' evaluated using AUC (Area Under ROC Curve) and feature importance is computed
  #' and retrieved
  #' 
  #' Parameters
  #' x: A sample x feature matrix containing explanatory variables
  #' y: A response vector containing sample labels (classes) that correspond to the samples in x
  #' targetClass: A character or a numeric value indicating the target class of which
  #' the samples are encoded with 1 and the other samples are encoded with 0
  #' iter: A numeric value specifying the number of time to run random forest. Default is 1
  #' split_method: A character specifying the splitting method to use, which should
  #' be one of 'random split' (default) or 'bootstrap'
  #' trainSet_ratio: A numeric value specifying the ratio of the size of the training
  #' set to that of the total dataset. Default is 0.8
  #' ntree: A numeric value specifying the number of trees to grow, which should not
  #' be set to too low to ensure that every input row gets predicted at least a few times.
  #' Default is 10000
  #' plot_ROC: A logical variable indicating whether ROC curve is plotted. Default is FALSE
  #' save_model: A logical variable indicating whether trained RF is saved. Default is FALSE
  #' 
  #' Return
  #' A list containing the following components:
  #' mtry: A vector of tuned best mtry used to train RF models
  #' auc_roc: A vector of computed AUC-ROC scores from trained RF models
  #' MDA, MDG: Matrices containing computed feature importance. Rows and columns
  #' present features and independent trained RF models. MDA and MDG stand for Mean
  #' Decrease Accuracy and Gini
  #' params: A list containing arguments of parameters 'split_method', 'trainSet_ratio',
  #' and 'ntree'
  #' rfRes: A list of objects of class randomForest
  #' yPredList: A list of class predictions of test set samples
  #' yTruthList: A list of ground truth labels of test set samples
  
  # Sanity check
  if (!(nrow(x) == length(y))) {
    stop("Row length of x (sample) should be same as length of y.")
  }
  if (!targetClass %in% unique(y)) {
    stop("Argument for 'targetClass' should be one of classes in y")
  }
  if (iter < 1) {
    stop("Argument for 'iter' should be integer greater than or equal to 1.")
  }
  if (trainSet_ratio < 0 | trainSet_ratio > 1) {
    stop("Argument for 'trainSet_ratio' should be numeric between 0 and 1.")
  }
  
  # Create containers to save results
  mtryList <- c()
  aucList <- c()
  matMDA <- matrix(data = NA, nrow = ncol(x), ncol = iter,
                   dimnames = list(colnames(x), as.character(seq_len(iter))))
  matMDG <- matrix(data = NA, nrow = ncol(x), ncol = iter,
                   dimnames = list(colnames(x), as.character(seq_len(iter))))
  paramList <- list(split_method = split_method, trainSet_ratio = trainSet_ratio,
                    ntree = ntree)
  rfResList <- as.list(rep(NA, iter))
  names(rfResList) <- as.character(seq_len(iter))
  yPredList <- as.list(rep(NA, iter))
  names(yPredList) <- as.character(seq_len(iter))
  yTruthList <- as.list(rep(NA, iter))
  names(yTruthList) <- as.character(seq_len(iter))
  
  # Encode target sample label with 1 and the rest with 0
  y <- ifelse(test = y == targetClass, yes = 1, no = 0)
  
  for (i in seq_len(iter)) {
    # Print progress
    if (i %% 10 == 0) {
      print(paste0('RF', i, ' is building...'))
    }
    
    # Separate data into training and test sets and turn response variable into
    # factor to conduct random forests for classification
    trainIdx <- subsetTrainData(x, y, split_method = split_method,
                                trainSet_ratio = trainSet_ratio)
    x_train <- x[trainIdx,, drop = F]
    y_train <- y[trainIdx] %>% as.factor()
    x_test <- x[-trainIdx,, drop = F]
    y_test <- y[-trainIdx] %>% as.factor()
    
    # Tune mtry
    tuneRes <- try(randomForest::tuneRF(x = x_train, y = y_train, ntreeTry = 1000, stepFactor = 1.5,
                                        improve = 0.01, trace = T, plot = F, dobest = F),
                   silent = T)
    if (!is(tuneRes, 'try-error')) {
      # Retrieve and save mtry with least OOB error
      mtry <- tuneRes[order(tuneRes[, 2]),, drop = F][1, 1]
      mtryList <- c(mtryList, mtry)
      # Run random forest
      rfRes <- randomForest::randomForest(x = x_train, y = y_train, mtry = mtry, ntree = ntree,
                                          xtest = x_test, ytest = y_test,
                                          importance = T, proximity = T, keep.forest = T)
    } else {
      mtryList <- c(mtryList, NA)
      rfRes <- randomForest::randomForest(x = x_train, y = y_train, ntree = ntree,
                                          xtest = x_test, ytest = y_test,
                                          importance = T, proximity = T, keep.forest = T)
    }
    
    # Evaluate model by AUC of ROC curve
    y_pred <- predict(rfRes, newdata = x_test, type = 'prob')
    if (plot_ROC) {
      par(pty = 's')
      rocRes <- try(pROC::roc(response = y_test, predictor = y_pred[, '1'], plot = T,
                              legacy.axes = T, print.auc = T, print.auc.x = 0.4,
                              xlab = 'False positive rate', ylab = 'True positive rate',
                              main = 'ROC Curve for Random Forest',
                              cex.lab = 1.2, cex.main = 1.1, col = '#377eb8', lwd = 4,
                              direction = '<'),
                    silent = T)
      par(pty = 'm')
    } else {
      rocRes <- try(pROC::roc(response = y_test, predictor = y_pred[, '1'], plot = F,
                              direction = '<'),
                    silent = T)
    }
    # Save AUC
    if (!is(rocRes, 'try-error')) {
      aucList <- c(aucList, rocRes$auc)
    } else {
      aucList <- c(aucList, NA)
    }
    
    # Retrieve and save feature importance values from RF object
    matMDA[, i] <- rfRes$importance[, 'MeanDecreaseAccuracy']
    matMDG[, i] <- rfRes$importance[, 'MeanDecreaseGini']
    
    # Save trained RF model
    if (save_model) {
      rfResList[[i]] <- rfRes
    }
    
    # Save predictions and ground truths of test set samples
    y_pred <- as.numeric(y_pred[, '1'] > 0.5) %>%
      as.factor()
    # Make labels interpretable
    names(y_pred) <- rownames(x_test)
    names(y_test) <- rownames(x_test)
    yPredList[[i]] <- y_pred
    yTruthList[[i]] <- y_test
  }
  
  return(list(mtry = mtryList, auc_roc = aucList, MDA = matMDA, MDG = matMDG,
              params = paramList, rfRes = rfResList, y_pred = yPredList, y_truth = yTruthList))
}

iterRF <- function(se, doImputation = F, doInitFeatSelection = F, sigFeatList = NULL,
                   doFeatClustering = F, iter = 100) {
  #' Iteratively train random forest models for stability selection, where input
  #' data goes through data imputation (optional), initial feature selection, feature
  #' grouping, and model training. Input data is imputed using missForest. Initial
  #' feature selection is done with provided significant feature list or by significance
  #' tests using limma, and feature grouping is done with hierarchical clustering
  #' to remove highly correlated features. This function is currently lack of freedom
  #' and very specific to cancer recurrence prediction
  #' 
  #' Parameters
  #' se: A SummarizedExperiment object containing a data matrix and patient recurrence information
  #' doImputation: A logical variable indicating whether to do imputation using
  #' missForest. Default is FALSE
  #' doInitFeatSelection: A logical variable indicating whether to do initial feature
  #' selection. Default is FALSE
  #' sigFeatList: A vector of characters indicating the features for initial feature
  #' selection to obtain a better chance to capture signals in data, for example,
  #' differential analysis results. If NULL, limma will be implemented on the input
  #' data to identify significant features. Default is NULL
  #' doFeatClustering: A logical variable indicating whether to do group highly
  #' correlated features using hierarchical clustering. Default is FALSE
  #' iter: A numeric value specifying the number of models to train. Default is 100
  #' 
  #' Return
  #' A list containing the information of trained models, feature clustering, and
  #' imputed data if parameter 'doImputation' is set to TRUE. They are the returns
  #' of user-defined functions, "runRF", "rmCorrFeats", and "imputeByMF"
  
  # Do imputation using missForest
  if (doImputation) {
    print('---- DATA IS BEING IMPUTED ----')
    se <- imputeByMF(se)
    impuSE <- se
  } else {
    impuSE <- NULL
  }
  
  # Prepare input data for RF
  # Retrieve data matrix from SE object
  datMat <- t(assay(se)) %>%
    tibble::as_tibble(rownames = 'Sample')
  # Retrieve patient recurrence annotations
  recurAnno <- tibble::as_tibble(colData(se), rownames = 'Sample') %>%
    dplyr::select(Sample, Recurrence)
  # Combine all needed information
  inputDat <- dplyr::left_join(recurAnno, datMat, by = 'Sample') %>%
    tibble::column_to_rownames('Sample')
  
  # Do initial feature selection
  if (doInitFeatSelection) {
    if (!is.null(sigFeatList)) {
      inputDat <- dplyr::select(inputDat, c(Recurrence, sigFeatList))
    } else {
      # Identify recurrence-related statistically significant features
      soaRes <- doSOA(se, meta_var = 'Recurrence', use_limma = T)
      inputDat <- dplyr::select(inputDat, c(Recurrence, soaRes$featSigAssoRes$Var1))
    }
    print('---- INITIAL FEATURE SELECTION IS DONE ----')
  }
  
  # Cluster highly correlated features and keep only one representative of each cluster
  if (doFeatClustering) {
    featClusters <- rmCorrFeats(t(inputDat[, -1]), cutoff = 0.9)
    inputDat <- dplyr::select(inputDat, c(Recurrence, rownames(featClusters$reducedData)))
    print('---- FEATURE CLUSTERING IS DONE ----')
  } else {
    featClusters <- NULL
  }
  
  # Run RF
  # Set random seed for reproducible outcomes
  set.seed(42)
  x <- as.matrix(inputDat[, -1])
  y <- inputDat[, 1]
  rfRes <- runRF(x, y, targetClass = 'Yes', iter = iter, split_method = 'random split',
                 trainSet_ratio = 0.8, ntree = 10000, plot_ROC = F, save_model = T)
  
  return(list(rfRes = rfRes, featClusters = featClusters, impuSE = impuSE))
}

doSysTrainRF <- function(rfRes, se, max_numTopFeats = 50) {
  #' Systematically train random forest models with different numbers of top important
  #' features summarized by bootstrapped model training. For instance, models are
  #' trained on top 2 important features 100 times, and mean AUC scores and 95% CI will be reported
  #' 
  #' Parameters
  #' rfRes: An output of function 'iterRF'
  #' se: A SummarizedExperiment object used as the input to function 'iterRF'
  #' max_numTopFeats: A numeric value specifying the maximum number of top important
  #' features allowed, where different combinations of features, from top 2 to top
  #' 'max_numTopFeats', will be used to train models. Default is 50
  #' 
  #' Return
  #' modelPerfTab: A table showing mean AUC and 95% CI of trained models
  
  # Prepare top important feature list
  # Process feature importance scores
  # Rank feature importance scores for all iterations
  rankedFeatImpoTab <- apply(rfRes$rfRes$MDG, 2, function(featImpo) {
    orderIdx <- order(featImpo, decreasing = T)
    featImpoRanks <- seq_along(orderIdx)
    names(featImpoRanks) <- orderIdx
    featImpo[as.numeric(names(featImpoRanks))] <- featImpoRanks
    return(featImpo)
  })
  # Compute median rank of each feature across iterations
  medianFeatImpo <- apply(rankedFeatImpoTab, 1, median)
  # Order features by median ranks
  rankedImpoFeats <- data.frame(Feature = names(medianFeatImpo),
                                ImpoMedian = medianFeatImpo) %>%
    dplyr::arrange(ImpoMedian) %>%
    dplyr::pull(Feature)
  
  # Prepare input data
  x <- t(assay(se))
  y <- colData(se)$Recurrence
  
  # Systematically train models with different combinations of top important features
  numTopImpoFeats <- 2:max_numTopFeats
  # Create containers to save results
  meanAUC <- c()
  lowerAUC <- c()
  upperAUC <- c()
  for (num in numTopImpoFeats) {
    print(paste0('---- TOP ', num, ' FEATURES USED ----'))
    # Subset data
    topImpoFeats <- rankedImpoFeats[1:num]
    # Train random forest model
    xSub <- x[, topImpoFeats, drop = F]
    rfResSystem <- runRF(xSub, y, targetClass = 'Yes', iter = 100, split_method = 'random split',
                         trainSet_ratio = 0.8, ntree = 10000, plot_ROC = F, save_model = F)
    # Save mean AUC and 95% CI
    auc_roc <- rfResSystem$auc_roc
    meanAUC <- c(meanAUC, round(mean(auc_roc), 3))
    ci <- calcCI(auc_roc, bootstrap = T)
    lowerAUC <- c(lowerAUC, ci[1])
    upperAUC <- c(upperAUC, ci[2])
  }
  # Summarize results
  modelPerfTab <- data.frame(FeatComb = paste0('Top', numTopImpoFeats),
                             MeanAUC = meanAUC, CI95 = paste0('[', round(lowerAUC, 2),
                                                              ',', round(upperAUC, 2), ']'))
  
  #### Return result as list for now to fit function 'vizSysTrainModelPerf'
  #### To save scores of test samples, return votes (rfRes$test$votes[, '1']) in
  #### function 'runRF'. Sample order of vote vector will be same as that of rfRes$rfRes$y_pred of certain model
  return(list(summPerformanceTab = modelPerfTab))
}

vizTopImpoFeatsRF <- function(rfRes, featAnno = NULL, fullData = NULL, trainData_smpType = 'Normal',
                              num_p1TopFeats = 15, num_p2TopFeats = NULL) {
  #' Summarize results of trained RF models and visualize top important features.
  #' This function can be further improved and extended
  #' 
  #' Parameters
  #' rfRes: An output of function 'iterRF'
  #' featAnno: A data frame of two columns where the first column stores feature
  #' names of the training data and the second column contains annotations of the
  #' features. Default is NULL
  #' fullData: An SE object of the full version of RF training data for making feature
  #' abundance boxplots (p2), containing a data matrix with all samples, e.g., Normal
  #' and Tumor tissue samples, and sample metadata. If NULL, p2 will not be made.
  #' Default is NULL
  #' trainData_smpType: A character specifying the sample type of RF training data,
  #' which should be 'Normal' (default) or 'Tumor'
  #' num_p1TopFeats: A numeric value specifying the number of top important features
  #' to show in a feature rank plot. Default is 15
  #' num_p2TopFeats: A numeric value specifying the number of top important features
  #' to display for feature abundance boxplots. If NULL, p2 will not be made. Default is NULL
  #' 
  #' Return
  #' A list containing summarized full feature importance table and two ggplot objects,
  #' feature rank and abundance plots
  
  if (!is.null(featAnno)) {
    if (ncol(featAnno) != 2) {
      stop("Column number of annotation data should be 2.")
    }
  }
  if (!is.null(fullData)) {
    if (!trainData_smpType %in% c('Normal', 'Tumor')) {
      stop("Argument for 'trainData_smpType' should be 'Normal' or 'Tumor'.")
    }
  }
  
  # Process feature importance scores
  # Rank feature importance scores for all iterations
  rankedFeatImpoTab <- apply(rfRes$rfRes$MDG, 2, function(featImpo) {
    orderIdx <- order(featImpo, decreasing = T)
    featImpoRanks <- seq_along(orderIdx)
    names(featImpoRanks) <- orderIdx
    featImpo[as.numeric(names(featImpoRanks))] <- featImpoRanks
    return(featImpo)
  })
  # Compute median, SD, and CI of each feature across iterations
  medianFeatImpo <- apply(rankedFeatImpoTab, 1, median)
  sdFeatImpo <- apply(rankedFeatImpoTab, 1, sd)
  lowerCIFeatImpo <- apply(rankedFeatImpoTab, 1, function(featImpoRanks) {
    calcCI(featImpoRanks, bootstrap = T)[1]
  })
  upperCIFeatImpo <- apply(rankedFeatImpoTab, 1, function(featImpoRanks) {
    calcCI(featImpoRanks, bootstrap = T)[2]
  })
  # Prepare feature annotation information
  if (!is.null(featAnno)) {
    colnames(featAnno) <- c('Feature', 'Annotation')
  }
  # Prepare cluster information
  featClustList <- rfRes$featClusters$corrFeats
  for (feat in names(featClustList)) {
    featClustList[[feat]] <- paste(featClustList[[feat]], collapse = '/')
  }
  featClustTab <- data.frame(Feature = names(featClustList), Cluster = unlist(featClustList)) %>%
    dplyr::mutate(Cluster = case_when(Cluster %in% '' ~ NA,
                                      !Cluster %in% '' ~ Cluster))
  # Collect all needed information into a table
  medianRankedFeatImpoTab <- data.frame(Feature = names(medianFeatImpo),
                                        RankMedian = medianFeatImpo,
                                        # RankSD = round(sdFeatImpo, 2),
                                        ImpoLowerCI = round(lowerCIFeatImpo, 1),
                                        ImpoUpperCI = round(upperCIFeatImpo, 1),
                                        CI95 = paste0('[', round(lowerCIFeatImpo, 1),
                                                      ',', round(upperCIFeatImpo, 1), ']')) %>%
    dplyr::arrange(RankMedian) %>%
    dplyr::left_join(featClustTab, by = 'Feature') %>%
    dplyr::relocate(Cluster, .after = Feature) %>%
    as.data.frame()
  medianRankedFeatImpoTab4Viz <- dplyr::mutate(medianRankedFeatImpoTab,
                                               Feature = stringr::str_remove(Feature, ';.+'))
  if (!is.null(featAnno)) {
    medianRankedFeatImpoTab <- dplyr::left_join(medianRankedFeatImpoTab, featAnno, by = 'Feature') %>%
      dplyr::relocate(Annotation, .after = Feature)
    medianRankedFeatImpoTab4Viz <- dplyr::mutate(medianRankedFeatImpoTab,
                                                 Feature = stringr::str_remove(Feature, ';.+'),
                                                 Annotation = stringr::str_remove(Annotation, ';.+'))
  }
  
  # Visualize top important features
  topImpoFeatTab <- medianRankedFeatImpoTab4Viz[seq_len(num_p1TopFeats),]
  if (is.null(featAnno)) {
    p1 <- ggplot(topImpoFeatTab, aes(x=RankMedian, y=factor(Feature, levels = rev(Feature))))
  } else {
    p1 <- ggplot(topImpoFeatTab, aes(x=RankMedian, y=factor(Annotation, levels = rev(Annotation))))
  }
  p1 <- p1 + geom_errorbar(aes(xmin=ImpoLowerCI, xmax=ImpoUpperCI)) +
    geom_point(size = 6) +
    labs(x = 'Median Rank of Feature Importance', y = 'Feature') +
    th
  
  # Visualize abundances of top important features
  if (!is.null(fullData) & !is.null(num_p2TopFeats)) {
    # Prepare data matrix and metadata including both Tumor and Normal samples
    datMat <- assay(fullData)
    rownames(datMat) <- stringr::str_remove(rownames(datMat), ';.+')
    smpAnnoTab <- tibble::as_tibble(colData(fullData), rownames = 'Sample')
    if (trainData_smpType == 'Normal') {
      gpLevel = c('Yes_Normal', 'No_Normal', 'Yes_Tumor', 'No_Tumor')
      condiCol = c('firebrick', 'grey50') #c(Normal, Tumor)
      comparisons = list(c('Yes_Normal', 'No_Normal'), c('Yes_Tumor', 'No_Tumor'))
    } else if (trainData_smpType == 'Tumor') {
      gpLevel = c('Yes_Tumor', 'No_Tumor', 'Yes_Normal', 'No_Normal')
      condiCol = c('grey50', 'firebrick') #c(Normal, Tumor)
      comparisons = list(c('Yes_Tumor', 'No_Tumor'), c('Yes_Normal', 'No_Normal'))
    }
    # Extract top important features and prepare needed information
    topImpoFeatTab <- medianRankedFeatImpoTab4Viz[1:num_p2TopFeats,]
    if (is.null(featAnno)) {
      topImpoFeats <- dplyr::mutate(topImpoFeatTab, Feat4Viz = Feature) %>%
        dplyr::select(Feature, Feat4Viz)
    } else {
      topImpoFeats <- dplr::select(topImpoFeatTab, Feature, Annotation) %>%
        dplyr::mutate(Feat4Viz = paste0(Feature, ' (', Annotation, ')'),
                      Feat4Viz = factor(Feat4Viz, levels = unique(Feat4Viz))) %>%
        dplyr::select(-Annotation)
    }
    topImpoFeatDat <- tibble::as_tibble(datMat[topImpoFeats$Feature,], rownames = 'Feature') %>%
      tidyr::pivot_longer(cols = -'Feature', names_to = 'Sample', values_to = 'Abundance') %>%
      dplyr::left_join(topImpoFeats, by = 'Feature') %>%
      dplyr::left_join(smpAnnoTab, by = 'Sample') %>%
      dplyr::mutate(Recur_Condi = paste0(Recurrence, '_', Condition),
                    Recur_Condi = factor(Recur_Condi, levels = gpLevel))
    p2 <- ggplot(topImpoFeatDat, aes(x=Recur_Condi, y=Abundance, col=Condition, fill=Recurrence)) +
      geom_boxplot(alpha = 1, outlier.shape = NA, linewidth = 1) +
      geom_jitter(position = position_jitter(0.2), size = 2, show.legend = F) +
      ggpubr::stat_compare_means(method = 't.test', paired = F, method.args = list(var.equal = T),
                                 comparisons = comparisons, label = 'p.signif', tip.length = 0.015,
                                 bracket.size = 0.7, size = 4) +
      labs(x = 'Recurrence') +
      scale_color_manual(values = condiCol) +
      scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
      facet_wrap(vars(Feat4Viz), scales = 'free') +
      th +
      theme(strip.text = element_text(size = 13, face = 'bold'),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  } else {
    p2 <- NULL
  }
  
  # Remove columns 'ImpoLowerCI' and 'ImpoUpperCI' for avoiding redundant information
  medianRankedFeatImpoTab <- dplyr::select(medianRankedFeatImpoTab, -c(ImpoLowerCI, ImpoUpperCI))
  
  return(list(fullImpoFeatTab = medianRankedFeatImpoTab, rank = p1, abun = p2))
}




runXGBoost <- function(x, y, targetClass, iter = 1, booster = 'gbtree', nrounds = 1000,
                       maxDepth = 6, eta = 0.3, trainSet_ratio = 0.8, plot_ROC = F,
                       save_model = F) {
  #' Perform XGBoost for binary classification problem. Trained model is evaluated
  #' using AUC (Area Under ROC Curve) and feature importance is retrieved
  #' 
  #' Parameters
  #' x: A sample x feature matrix containing explanatory variables
  #' y: A response vector containing sample labels (classes) that correspond to the samples in x
  #' targetClass: A character or a numeric value indicating the target class of which
  #' the samples are encoded with 1 and the other samples are encoded with 0
  #' iter: A numeric value specifying the number of models to train. Default is 1
  #' booster, nrounds, maxDepth, eta: Parameters as described in function 'xgb.train'
  #' trainSet_ratio: A numeric value specifying the ratio of the size of the training
  #' set to that of the total dataset. The rest of the observations will be evenly
  #' separated into the validation and test sets. Default is 0.8
  #' plot_ROC: A logical variable indicating whether ROC curve is plotted. Default is FALSE
  #' save_model: A logical variable indicating whether trained model is saved. Default is FALSE
  #' 
  #' Return
  #' A list containing the following components:
  #' auc_roc: A vector of computed AUC-ROC scores from trained models
  #' featImpoTab: A table containing summarized feature importance for tree boosted
  #' models and all feature importance for linear boosted models
  #' best_nrounds: A vector of best nrounds used to train models
  #' params: A list containing arguments of parameters 'booster', 'nrounds', 'maxDepth',
  #' 'eta', and 'trainSet_ratio'
  #' modelRes: A list of trained models
  #' yPredList: A list of class predictions of test set samples
  #' yTruthList: A list of ground truth labels of test set samples
  
  # Sanity check
  if (!(nrow(x) == length(y))) {
    stop("Row length of x (sample) should be same as length of y.")
  }
  if (!targetClass %in% unique(y)) {
    stop("Argument for 'targetClass' should be one of classes in y")
  }
  if (iter < 1) {
    stop("Argument for 'iter' should be integer greater than or equal to 1.")
  }
  if (trainSet_ratio < 0 | trainSet_ratio > 1) {
    stop("Argument for 'trainSet_ratio' should be numeric between 0 and 1.")
  }
  
  # Create containers to save results
  bestnroundList <- c()
  aucList <- c()
  if (booster == 'gbtree') {
    featImpoTab <- data.frame(Feature = NA, Gain = NA)
    paramList <- list(booster = booster, nrounds = nrounds, maxDepth = maxDepth,
                      eta = eta, trainSet_ratio = trainSet_ratio)
  } else if (booster == 'gblinear') {
    featImpoTab <- matrix(data = NA, nrow = ncol(x), ncol = iter,
                          dimnames = list(colnames(x), as.character(seq_len(iter))))
    paramList <- list(booster = booster, nrounds = nrounds, maxDepth = NA, eta = NA,
                      trainSet_ratio = trainSet_ratio)
  }
  xgbResList <- as.list(rep(NA, iter))
  names(xgbResList) <- as.character(seq_len(iter))
  yPredList <- as.list(rep(NA, iter))
  names(yPredList) <- as.character(seq_len(iter))
  yTruthList <- as.list(rep(NA, iter))
  names(yTruthList) <- as.character(seq_len(iter))
  
  # Encode target sample label with 1 and the rest with 0
  y <- ifelse(test = y == targetClass, yes = 1, no = 0)
  
  # Prepare parameters for training XGBoost
  if (booster == 'gbtree') {
    parameters <- list(booster = 'gbtree',
                       objective = 'binary:logistic',
                       eval_metric = 'logloss',
                       eval_metric = 'auc',
                       max_depth = maxDepth,
                       eta = eta)
  } else if (booster == 'gblinear') {
    parameters <- list(booster = 'gblinear',
                       objective = 'binary:logistic',
                       eval_metric = 'logloss',
                       eval_metric = 'rmse')
  }
  
  for (i in seq_len(iter)) {
    # Print progress
    if (i %% 10 == 0) {
      print(paste0('XGBoost', i, ' is building...'))
    }
    
    # Separate data into training, validation, and test sets
    trainIdx <- caret::createDataPartition(y, times = 1, p = trainSet_ratio, list = T)
    x_train <- x[trainIdx[[1]],, drop = F]
    y_train <- y[trainIdx[[1]]]
    restIdx <- seq_along(y)[!seq_along(y) %in% trainIdx[[1]]]
    ySub <- y[restIdx]
    validIdx <- caret::createDataPartition(ySub, times = 1, p = 0.5, list = T)
    x_valid <- x[restIdx[validIdx[[1]]],, drop = F]
    y_valid <- y[restIdx[validIdx[[1]]]]
    x_test <- x[restIdx[-validIdx[[1]]],, drop = F]
    y_test <- y[restIdx[-validIdx[[1]]]]
    # Collect split data into xgb.DMatrix objects
    trainDat <- xgb.DMatrix(data = x_train, label = y_train)
    validDat <- xgb.DMatrix(data = x_valid, label = y_valid)
    testDat <- xgb.DMatrix(data = x_test, label = y_test)
    
    # Train XGBoost
    watchlist <- list(train = trainDat, valid = validDat)
    xgboost <- xgb.train(data = trainDat, params = parameters, nrounds = nrounds,
                         watchlist = watchlist, verbose = 1, early_stopping_rounds = 4)
    # Save best nrounds that avoids overfitting
    bestnroundList <- c(bestnroundList, xgboost$best_iteration)
    
    # Evaluate model by AUC of ROC curve
    y_pred <- predict(xgboost, testDat)
    y_truth <- getinfo(testDat, 'label')
    if (plot_ROC) {
      par(pty = 's')
      rocRes <- try(pROC::roc(response = y_truth, predictor = y_pred, plot = T,
                              legacy.axes = T, print.auc = T, print.auc.x = 0.4,
                              xlab = 'False positive rate', ylab = 'True positive rate',
                              main = 'ROC Curve for XGBoost',
                              cex.lab = 1.2, cex.main = 1.1, col = '#377eb8', lwd = 4,
                              direction = '<'),
                    silent = T)
      par(pty = 'm')
    } else {
      rocRes <- try(pROC::roc(response = y_truth, predictor = y_pred, plot = F,
                              direction = '<'),
                    silent = T)
    }
    # Save AUC
    if (!is(rocRes, 'try-error')) {
      aucList <- c(aucList, rocRes$auc)
    } else {
      aucList <- c(aucList, NA)
    }
    
    # Retrieve and save feature importance values from trained model
    if (booster == 'gbtree') {
      featImpo <- xgb.importance(model = xgboost) %>%
        dplyr::select(Feature, Gain)
      featImpoTab <- dplyr::bind_rows(featImpoTab, featImpo)
    } else if (booster == 'gblinear') {
      featImpo <- xgb.importance(model = xgboost) %>%
        dplyr::select(Feature, Weight) %>%
        tibble::column_to_rownames('Feature')
      # Turn feature coefficients into ranks
      featImpo$Weight <- seq_len(nrow(featImpoTab))
      featImpoTab[, i] <- featImpo[rownames(featImpoTab),]
    }
    
    # Save trained model
    if (save_model) {
      xgbResList[[i]] <- xgboost
    }
    
    # Save predictions and ground truths of test set samples
    y_pred <- as.numeric(y_pred > 0.5) %>%
      as.factor()
    y_test <- as.factor(y_test)
    # Make labels interpretable
    names(y_pred) <- rownames(x_test)
    names(y_test) <- rownames(x_test)
    yPredList[[i]] <- y_pred
    yTruthList[[i]] <- y_test
  }
  
  # Summarize feature importance table
  if (booster == 'gbtree') {
    featImpoTab <- featImpoTab[-1,] %>%
      dplyr::group_by(Feature) %>%
      dplyr::summarise(OccurRate = length(Feature)/iter, Occurrence = length(Feature),
                       GainMean = mean(Gain), GainSD = sd(Gain)) %>%
      dplyr::arrange(dplyr::desc(Occurrence))
  } else if (booster == 'gblinear') {
    # Compute median, SD, and 95% CI of ranks of each feature across iterations
    featImpoTab <- lapply(seq_len(nrow(featImpoTab)), function (rowIdx) {
      featRanks <- featImpoTab[rowIdx,,drop = T]
      medianFeatRank <- median(featRanks)
      sdFeatRank <- sd(featRanks)
      ciFeatRank <- calcCI(featRanks, bootstrap = T)
      
      data.frame(Feature = rownames(featImpoTab)[rowIdx],
                 RankMedian = medianFeatRank,
                 # RankSD = sdFeatRank,
                 RankLowerCI = ciFeatRank[1],
                 RankUpperCI = ciFeatRank[2],
                 row.names = c())
    }) %>% dplyr::bind_rows() %>%
      dplyr::arrange(RankMedian)
  }
  
  return(list(auc_roc = aucList, featImpoTab = featImpoTab, best_nrounds = bestnroundList,
              params = paramList, xgbRes = xgbResList, y_pred = yPredList, y_truth = yTruthList))
}

iterXGBoost <- function(se, doImputation = F, doInitFeatSelection = F, sigFeatList = NULL,
                        doFeatClustering = F, booster = 'gbtree', iter = 100, nrounds = 1000) {
  #' Iteratively train XGBoost models for stability selection, where input data
  #' goes through data imputation (optional), initial feature selection, feature
  #' grouping, and model training. Input data is imputed using missForest. Initial
  #' feature selection is done with provided significant feature list or by significance
  #' tests using limma, and feature grouping is done with hierarchical clustering
  #' to remove highly correlated features. This function is currently lack of freedom
  #' and very specific to cancer recurrence prediction
  #' 
  #' Parameter
  #' se: A SummarizedExperiment object containing a data matrix and patient recurrence information
  #' doImputation: A logical variable indicating whether to do imputation using
  #' missForest. Default is FALSE
  #' doInitFeatSelection: A logical variable indicating whether to do initial feature
  #' selection. Default is FALSE
  #' sigFeatList: A vector of characters indicating the features for initial feature
  #' selection to obtain a better chance to capture signals in data, for example,
  #' differential analysis results. If NULL, limma will be implemented on the input
  #' data to identify significant features. Default is NULL
  #' doFeatClustering: A logical variable indicating whether to do group highly
  #' correlated features using hierarchical clustering. Default is FALSE
  #' booster: A character specifying the booster to use, which should be one of
  #' 'gbtree' (default) or 'gblinear'
  #' iter: A numeric value specifying the number of models to train. Default is 100
  #' nrounds: Parameter as described in function 'xgb.train'
  #' 
  #' Return
  #' A list containing the information of trained models, feature clustering, and
  #' imputed data if parameter 'doImputation' is set to TRUE. They are the returns
  #' of user-defined functions, "runXGBoost", "rmCorrFeats", and "imputeByMF"
  
  # Do imputation using missForest
  if (doImputation) {
    print('---- DATA IS BEING IMPUTED ----')
    se <- imputeByMF(se)
    impuSE <- se
  } else {
    impuSE <- NULL
  }
  
  # Prepare input data
  x <- t(assay(se))
  y <- colData(se)$Recurrence
  
  # Do initial feature selection
  if (doInitFeatSelection) {
    if (!is.null(sigFeatList)) {
      x <- x[, colnames(x) %in% sigFeatList]
    } else {
      # Identify recurrence-related statistically significant features
      soaRes <- doSOA(se, meta_var = 'Recurrence', use_limma = T)
      x <- x[, colnames(x) %in% soaRes$featSigAssoRes$Var1]
    }
    print('---- INITIAL FEATURE SELECTION IS DONE ----')
  }
  
  # Cluster highly correlated features and keep only one representative of each cluster
  if (doFeatClustering) {
    featClusters <- rmCorrFeats(t(x), cutoff = 0.9)
    x <- t(featClusters$reducedData)
    print('---- FEATURE CLUSTERING IS DONE ----')
  } else {
    featClusters <- NULL
  }
  
  # Train XGBoost with tree or linear booster
  set.seed(42)
  if (booster == 'gbtree') {
    xgbRes <- runXGBoost(x, y, targetClass = 'Yes', iter = iter, booster = 'gbtree',
                         nrounds = nrounds, maxDepth = 6, eta = 0.3, trainSet_ratio = 0.7,
                         plot_ROC = F, save_model = T)
  } else if (booster == 'gblinear') {
    xgbRes <- runXGBoost(x, y, targetClass = 'Yes', iter = iter, booster = 'gblinear',
                         nrounds = nrounds, trainSet_ratio = 0.7, plot_ROC = F, save_model = T)
  }
  
  return(list(xgbRes = xgbRes, featClusters = featClusters, impuSE = impuSE))
}




runLogisR <- function(x, y, targetClass, regularized_method = 'lasso', cvFold = 10,
                      cvMeasure = 'auc', used_lambda = 'lambda.1se', iter = 1,
                      trainSet_ratio = 0.8, split_method = 'random split', plot_ROC = F,
                      save_model = F) {
  #' Use regularized logistic regression model for binary classification problem
  #' and feature selection. Trained model is evaluated by AUC
  #' 
  #' Parameters
  #' x: A sample x feature matrix containing explanatory variables
  #' y: A response vector containing sample labels that correspond to the samples in x
  #' targetClass: A character or a numeric value indicating the target class of which
  #' the samples are encoded with 1 and the other samples are encoded with 0
  #' regularized_method: A character specifying the regularized method to use, which
  #' should be 'lasso' (default) or 'ridge'
  #' cvFold: A numeric value specifying the number of folds for cross-validation
  #' for optimizing lambda. Default is 10, and it can be as large as the sample size
  #' (leave-one-out cross-validation), but it is not recommended for large datasets
  #' cvMeasure: A character specifying the evaluation metric for cross-validation
  #' used_lambda: A character specifying the tuned lambda to use, which should be
  #' one of 'lambda.1se' (default) or 'lambda.min'
  #' iter: A numeric value specifying the number of models to train. Default is 1
  #' split_method: A character specifying the splitting method to use, which should
  #' be one of 'random split' (default) or 'bootstrap'
  #' trainSet_ratio: A numeric value specifying the ratio of the size of the training
  #' set to that of the total dataset
  #' plot_ROC: A logical variable indicating whether ROC curve is plotted
  #' save_model: A logical variable indicating whether trained model is saved. Default is FALSE
  #' 
  #' Return
  #' A list containing the following components:
  #' coefficient: A matrix containing optimized coefficients. Rows and columns present
  #' features and independent trained models
  #' usedLambda: A vector of used lambda for selected trained models
  #' auc_roc: A vector of computed AUC from trained models
  #' nNonZero: A vector of counts of nonzero coefficients, i.e., selected features
  #' params: A list containing arguments of parameters 'regularized_method', 'cvFold',
  #' 'used_lambda', 'split_method', and 'trainSet_ratio'
  #' lrRes: A list of trained models
  #' y_pred: A list of class predictions of test set samples
  #' y_truth: A list of ground truth labels of test set samples
  #' y_pred_link: A list of linear regression scores of test set samples, which
  #' is the logit transformation of the output of logistic regression, log(p/(1-p))
  
  # Sanity check
  if (!(nrow(x) == length(y))) {
    stop("Row length of x (sample) should be same as length of y.")
  }
  if (!(targetClass %in% unique(y))) {
    stop("Argument for 'targetClass' should be one of classes in y")
  }
  if (!(regularized_method %in% c('lasso', 'ridge'))) {
    stop("Argument for 'regularized_method' should be 'lasso' or 'ridge'.")
  }
  if (!(cvMeasure %in% c('deviance', 'class', 'auc', 'mse'))) {
    stop("Argument for 'cvMeasure' should be one of 'deviance', 'class', 'auc', or 'mse'.")
  }
  if (!(used_lambda %in% c('lambda.1se', 'lambda.min'))) {
    stop("Argument for 'used_lambda' should be 'lambda.1se' or 'lambda.min'.")
  }
  if (iter < 1) {
    stop("Argument for 'iter' should be integer greater than or equal to 1.")
  }
  if (trainSet_ratio < 0 | trainSet_ratio > 1) {
    stop("Argument for 'trainSet_ratio' should be numeric between 0 and 1.")
  }
  
  # Create containers to save results
  matCoeffi  <- matrix(data = NA, nrow = ncol(x), ncol = iter,
                       dimnames = list(colnames(x), as.character(seq_len(iter))))
  lambdaList <- c()
  aucList <- c()
  numNonZeroList <- c()
  paramList <- list(regularized_method = regularized_method, cvFold = cvFold, cvMeasure = cvMeasure,
                    used_lambda = used_lambda, split_method = split_method, trainSet_ratio = trainSet_ratio)
  lrResList <- as.list(rep(NA, iter))
  names(lrResList) <- as.character(seq_len(iter))
  yPredList <- as.list(rep(NA, iter))
  names(yPredList) <- as.character(seq_len(iter))
  yTruthList <- as.list(rep(NA, iter))
  names(yTruthList) <- as.character(seq_len(iter))
  yPredLinkList <- as.list(rep(NA, iter))
  names(yPredLinkList) <- as.character(seq_len(iter))
  
  # Encode target sample label with 1 and the other with 0
  # Will be coerced into factor and target class is the last level in alphabetical order
  y <- ifelse(test = y == targetClass, yes = 1, no = 0)
  
  # Define penalty term for model according to regularized method specified
  if (regularized_method == 'lasso') {
    alpha <- 1
  } else if (regularized_method == 'ridge') {
    alpha <- 0
  }
  
  for (i in seq_len(iter)) {
    # Print progress
    if (i %% 10 == 0) {
      print(paste0('LR model', i, ' is running...'))
    }
    
    # Separate data into training and test sets
    trainIdx <- subsetTrainData(x, y, split_method = split_method,
                                trainSet_ratio = trainSet_ratio)
    x_train <- x[trainIdx,, drop = F]
    y_train <- y[trainIdx]
    x_test <- x[-trainIdx,, drop = F]
    y_test <- y[-trainIdx]
    
    # Train logistic regression model
    lrRes <- glmnet::cv.glmnet(x = x_train, y = y_train, family = 'binomial', type.measure = cvMeasure,
                               nfolds = cvFold, alpha = alpha, standardize = F, intercept = T, nlambda = 100)
    # Save selected lambda of trained model
    lambdaList <- c(lambdaList, round(lrRes[[used_lambda]], 4))
    # Retrieve and save coefficients of fitted model (remove first row intercept)
    matCoeffi[, i] <- as.matrix(coef(lrRes, s = lrRes[[used_lambda]]))[-1, 1]
    
    # Save number of non-zero coefficients of fitted model
    coeffi <- as.matrix(coef(lrRes, s = lrRes[[used_lambda]]))[-1, 1]
    numNonZeroList <- c(numNonZeroList, sum(coeffi != 0))
    
    # Evaluate model by AUC of ROC curve
    y_pred <- predict(lrRes, s = lrRes[[used_lambda]], newx = x_test, type = 'response')
    if (plot_ROC) {
      par(pty = 's')
      rocRes <- try(pROC::roc(response = y_test, predictor = y_pred[, 1], plot = T,
                              legacy.axes = T, print.auc = T, print.auc.x = 0.4,
                              xlab = 'False positive rate', ylab = 'True positive rate',
                              main = 'ROC Curve for Lasso Logistic Regression',
                              cex.lab = 1.2, cex.main = 1.1, col = '#377eb8', lwd = 4,
                              direction = '<'),
                    silent = T)
      par(pty = 'm')
    } else {
      rocRes <- try(pROC::roc(response = y_test, predictor = y_pred[, 1], plot = F,
                              direction = '<'),
                    silent = T)
    }
    # Save AUC
    if (!is(rocRes, 'try-error')) {
      aucList <- c(aucList, rocRes$auc)
    } else {
      aucList <- c(aucList, NA)
    }
    
    # Save trained model
    if (save_model) {
      lrResList[[i]] <- lrRes
    }
    
    # Save predictions and ground truths of test set samples
    y_pred <- as.numeric(y_pred > 0.5) %>%
      as.factor()
    y_test <- as.factor(y_test)
    # Make labels interpretable
    names(y_pred) <- rownames(x_test)
    names(y_test) <- rownames(x_test)
    yPredList[[i]] <- y_pred
    yTruthList[[i]] <- y_test
    
    # Save linear regression scores of test set samples
    y_predLink <- predict(lrRes, s = lrRes[[used_lambda]], newx = x_test, type = 'link')[, 1]
    names(y_predLink) <- rownames(x_test)
    yPredLinkList[[i]] <- y_predLink
  }
  
  return(list(coefficient = matCoeffi, usedLambda = lambdaList, auc_roc = aucList,
              nNonZero = numNonZeroList, params = paramList, lrRes = lrResList,
              y_pred = yPredList, y_truth = yTruthList, y_pred_link = yPredLinkList))
}

iterLogisR <- function(se, doImputation = F, doInitFeatSelection = F, sigFeatList = NULL,
                       doFeatClustering = F, iter = 100, cvFold = 10, cvMeasure = 'auc') {
  #' Iteratively train lasso logistic regression models for stability selection,
  #' where input data goes through data imputation (optional), initial feature selection,
  #' feature grouping, and model training. Input data is imputed using missForest.
  #' Initial feature selection is done with provided significant feature list or
  #' by significance tests using limma, and feature grouping is done with hierarchical
  #' clustering to remove highly correlated features. This function is currently
  #' lack of freedom and very specific to cancer recurrence prediction
  #' 
  #' Parameter
  #' se: A SummarizedExperiment object containing a data matrix and patient recurrence information
  #' doImputation: A logical variable indicating whether to do imputation using
  #' missForest. Default is FALSE
  #' doInitFeatSelection: A logical variable indicating whether to do initial feature
  #' selection. Default is FALSE
  #' sigFeatList: A vector of characters indicating the features for initial feature
  #' selection to obtain a better chance to capture signals in data, for example,
  #' differential analysis results. If NULL, limma will be implemented on the input
  #' data to identify significant features. Default is NULL
  #' doFeatClustering: A logical variable indicating whether to do group highly
  #' correlated features using hierarchical clustering. Default is FALSE
  #' iter: A numeric value specifying the number of models to train. Default is 100
  #' cvFold: A numeric value specifying the number of folds for cross-validation
  #' for optimizing lambda. Default is 10, and it can be as large as the sample size
  #' (leave-one-out cross-validation), but it is not recommended for large datasets
  #' cvMeasure: A character specifying the evaluation metric for cross-validation
  #' 
  #' Return
  #' A list containing the information of trained models, feature clustering, and
  #' imputed data if parameter 'doImputation' is set to TRUE. They are the returns
  #' of user-defined functions, "runLogisR", "rmCorrFeats", and "imputeByMF"
  
  # Do imputation using missForest
  if (doImputation) {
    print('---- DATA IS BEING IMPUTED ----')
    se <- imputeByMF(se)
    impuSE <- se
  } else {
    impuSE <- NULL
  }
  
  # Prepare input data
  x <- t(assay(se))
  y <- colData(se)$Recurrence
  
  # Do initial feature selection
  if (doInitFeatSelection) {
    if (!is.null(sigFeatList)) {
      x <- x[, colnames(x) %in% sigFeatList]
    } else {
      # Identify recurrence-related statistically significant features
      soaRes <- doSOA(se, meta_var = 'Recurrence', use_limma = T)
      x <- x[, colnames(x) %in% soaRes$featSigAssoRes$Var1]
    }
    print('---- INITIAL FEATURE SELECTION IS DONE ----')
  }
  
  # Cluster highly correlated features and keep only one representative of each cluster
  if (doFeatClustering) {
    featClusters <- rmCorrFeats(t(x), cutoff = 0.9)
    x <- t(featClusters$reducedData)
    print('---- FEATURE CLUSTERING IS DONE ----')
  } else {
    featClusters <- NULL
  }
  
  # Train lasso logistic regression model
  set.seed(42)
  lrRes <- runLogisR(x, y, targetClass = 'Yes', iter = iter, regularized_method = 'lasso',
                     cvFold = cvFold, cvMeasure = cvMeasure, used_lambda = 'lambda.min',
                     trainSet_ratio = 0.8, split_method = 'random split',
                     plot_ROC = F, save_model = T)
  
  return(list(lrRes = lrRes, featClusters = featClusters, impuSE = impuSE))
}

doSysTrainLogisR <- function(lrRes, se, max_numTopFeats = 50) {
  #' Systematically train lasso logistic regression models with different numbers
  #' of top important features summarized by bootstrapped model training. For instance,
  #' models are trained on top 2 important features 100 times, and mean AUC scores
  #' and 95% CI will be reported
  #' 
  #' Parameters
  #' lrRes: An output of function 'iterLogisR'
  #' se: A SummarizedExperiment object used as the input to function 'iterLogisR'
  #' max_numTopFeats: A numeric value specifying the maximum number of top important
  #' features allowed, where different combinations of features, from top 2 to top
  #' 'max_numTopFeats', will be used to train models. Default is 50
  #' 
  #' Return
  #' A list containing summarized table showing mean AUC and 95% CI of trained models
  #' and all training results for latter whatever use
  
  # Prepare top important feature list
  # Pinpoint useless models that will be removed
  uselessModels <- which(lrRes$lrRes$nNonZero == 0)
  # Prepare coefficient table
  coefTab <- lrRes$lrRes$coefficient
  if (length(uselessModels) != 0) {
    message(paste(length(uselessModels), 'models with no non-zero coefficient, i.e., useless models.'))
    coefTab <- coefTab[, -uselessModels]
  }
  # Summarize feature selection results into pick rates and rank features by pick rates
  rankedImpoFeats <- as.data.frame(coefTab) %>%
    dplyr::mutate(across(everything(), ~ case_when(.x != 0 ~ 1,
                                                   .x == 0 ~ 0))) %>%
    tibble::rownames_to_column('Feature') %>%
    tidyr::pivot_longer(cols = -'Feature', names_to = 'Model', values_to = 'Pick') %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(PickRate = sum(Pick)/ncol(coefTab)) %>%
    dplyr::filter(PickRate != 0) %>%
    dplyr::arrange(desc(PickRate)) %>%
    dplyr::pull(Feature)
  
  # Prepare input data
  x <- t(assay(se))
  y <- colData(se)$Recurrence
  
  # Systematically train models with different combinations of top important features
  numTopImpoFeats <- 2:max_numTopFeats
  # Create containers to save results
  meanAUC <- c()
  lowerAUC <- c()
  upperAUC <- c()
  trainRes <- list()
  for (num in numTopImpoFeats) {
    print(paste0('---- TOP ', num, ' FEATURES USED ----'))
    # Subset data
    topImpoFeats <- rankedImpoFeats[1:num]
    # Train lasso logistic regression model
    xSub <- x[, topImpoFeats, drop = F]
    lrResSystem <- runLogisR(xSub, y, targetClass = 'Yes', iter = 100, regularized_method = 'lasso',
                             cvFold = lrRes$lrRes$params$cvFold, cvMeasure = lrRes$lrRes$params$cvMeasure,
                             used_lambda = 'lambda.min', trainSet_ratio = 0.8, split_method = 'random split',
                             plot_ROC = F, save_model = F)
    # Save mean AUC and 95% CI
    auc_roc <- lrResSystem$auc_roc
    meanAUC <- c(meanAUC, round(mean(auc_roc), 3))
    ci <- calcCI(auc_roc, bootstrap = T)
    lowerAUC <- c(lowerAUC, ci[1])
    upperAUC <- c(upperAUC, ci[2])
    # Save training results
    trainRes[[paste0('Top', num)]] <- lrResSystem
  }
  # Summarize results
  modelPerfTab <- data.frame(FeatComb = paste0('Top', numTopImpoFeats),
                             MeanAUC = meanAUC, CI95 = paste0('[', round(lowerAUC, 2),
                                                              ',', round(upperAUC, 2), ']'))
  
  return(list(summPerformanceTab = modelPerfTab, allTrainRes = trainRes))
}

vizTopImpoFeatsLR <- function(lrRes, featAnno = NULL, fullData = NULL, trainData_smpType = 'Normal',
                              num_p1TopFeats = 15, num_p2TopFeats = NULL) {
  #' Summarize results of trained logistic regression models and visualize top important
  #' features. This function should be further improved and extended
  #' 
  #' Parameters
  #' lrRes: An output of function 'iterLogisR'
  #' featAnno: A data frame of two columns where the first column stores feature
  #' names of the training data and the second column contains annotations of the
  #' features. Default is NULL
  #' fullData: An SE object of the full version of RF training data for making feature
  #' abundance boxplots (p2), containing a data matrix with all samples, e.g., Normal
  #' and Tumor tissue samples, and sample metadata. If NULL, p2 will not be made.
  #' Default is NULL
  #' trainData_smpType: A character specifying the sample type of RF training data,
  #' which should be 'Normal' (default) or 'Tumor'
  #' num_p1TopFeats: A numeric value specifying the number of top important features
  #' to show in a feature pick rate plot. Default is 15
  #' num_p2TopFeats: A numeric value specifying the number of top important features
  #' to display for feature abundance boxplots. If NULL, p2 will not be made. Default is NULL
  #' 
  #' Return
  #' A list containing summarized full feature importance table and two ggplot objects,
  #' feature pick rate and abundance plots
  
  if (!is.null(featAnno)) {
    if (ncol(featAnno) != 2) {
      stop("Column number of annotation data should be 2.")
    }
  }
  if (!is.null(fullData)) {
    if (!trainData_smpType %in% c('Normal', 'Tumor')) {
      stop("Argument for 'trainData_smpType' should be 'Normal' or 'Tumor'.")
    }
  }
  
  # Summarize feature selection results
  # Pinpoint useless models that will be removed
  uselessModels <- which(lrRes$lrRes$nNonZero == 0)
  # Prepare coefficient table
  coefTab <- lrRes$lrRes$coefficient
  if (length(uselessModels) != 0) {
    coefTab <- coefTab[, -uselessModels]
  }
  # Prepare feature annotation information
  if (!is.null(featAnno)) {
    colnames(featAnno) <- c('Feature', 'Annotation')
  }
  # Prepare cluster information
  featClustList <- lrRes$featClusters$corrFeats
  for (feat in names(featClustList)) {
    featClustList[[feat]] <- paste(featClustList[[feat]], collapse = '/')
  }
  featClustTab <- data.frame(Feature = names(featClustList), Cluster = unlist(featClustList)) %>%
    dplyr::mutate(Cluster = case_when(Cluster %in% '' ~ NA,
                                      !Cluster %in% '' ~ Cluster))
  # Summarize feature selection results into pick rates, combined needed information,
  # and convert it to long table
  impoFeatTab <- as.data.frame(coefTab) %>%
    dplyr::mutate(across(everything(), ~ case_when(.x != 0 ~ 1,
                                                   .x == 0 ~ 0))) %>%
    tibble::rownames_to_column('Feature') %>%
    tidyr::pivot_longer(cols = -'Feature', names_to = 'Model', values_to = 'Pick') %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(PickRate = sum(Pick)/ncol(coefTab)) %>%
    dplyr::filter(PickRate != 0) %>%
    dplyr::arrange(desc(PickRate)) %>%
    dplyr::left_join(featClustTab, by = 'Feature') %>%
    dplyr::relocate(Cluster, .after = Feature) %>%
    as.data.frame()
  impoFeatTab4Viz <- dplyr::mutate(impoFeatTab, Feature = stringr::str_remove(Feature, ';.+'))
  if (!is.null(featAnno)) {
    impoFeatTab <- dplyr::left_join(impoFeatTab, featAnno, by = 'Feature') %>%
      dplyr::relocate(Annotation, .after = Feature)
    impoFeatTab4Viz <- dplyr::mutate(impoFeatTab, Feature = stringr::str_remove(Feature, ';.+'),
                                     Annotation = stringr::str_remove(Annotation, ';.+'))
  }
  
  # Visualize top important features
  topImpoFeatTab <- impoFeatTab4Viz[seq_len(num_p1TopFeats),]
  if (is.null(featAnno)) {
    p1 <- ggplot(topImpoFeatTab, aes(x=PickRate, y=factor(Feature, levels = rev(Feature))))
  } else {
    p1 <- ggplot(topImpoFeatTab, aes(x=PickRate, y=factor(Annotation, levels = rev(Annotation))))
  }
  p1 <- p1 + geom_bar(stat = 'identity', position = position_dodge(), alpha = 0.9) +
    labs(x = 'Feature Pick Rate (~ 1000 models)', y = 'Feature') +
    scale_x_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    th
  
  # Visualize abundances of top important features
  if (!is.null(fullData) & !is.null(num_p2TopFeats)) {
    # Prepare data matrix and metadata including both Tumor and Normal samples
    datMat <- assay(fullData)
    rownames(datMat) <- stringr::str_remove(rownames(datMat), ';.+')
    smpAnnoTab <- tibble::as_tibble(colData(fullData), rownames = 'Sample')
    if (trainData_smpType == 'Normal') {
      gpLevel = c('Yes_Normal', 'No_Normal', 'Yes_Tumor', 'No_Tumor')
      condiCol = c('firebrick', 'grey50') #c(Normal, Tumor)
      comparisons = list(c('Yes_Normal', 'No_Normal'), c('Yes_Tumor', 'No_Tumor'))
    } else if (trainData_smpType == 'Tumor') {
      gpLevel = c('Yes_Tumor', 'No_Tumor', 'Yes_Normal', 'No_Normal')
      condiCol = c('grey50', 'firebrick') #c(Normal, Tumor)
      comparisons = list(c('Yes_Tumor', 'No_Tumor'), c('Yes_Normal', 'No_Normal'))
    }
    # Extract top important features and prepare needed information
    topImpoFeatTab <- impoFeatTab4Viz[1:num_p2TopFeats,]
    if (is.null(featAnno)) {
      topImpoFeats <- dplyr::mutate(topImpoFeatTab, Feat4Viz = Feature) %>%
        dplyr::select(Feature, Feat4Viz)
    } else {
      topImpoFeats <- dplyr::select(topImpoFeatTab, Feature, Annotation) %>%
        dplyr::mutate(Feat4Viz = paste0(Feature, ' (', Annotation, ')'),
                      Feat4Viz = factor(Feat4Viz, levels = unique(Feat4Viz))) %>%
        dplyr::select(-Annotation)
    }
    topImpoFeatDat <- tibble::as_tibble(datMat[topImpoFeats$Feature,], rownames = 'Feature') %>%
      tidyr::pivot_longer(cols = -'Feature', names_to = 'Sample', values_to = 'Abundance') %>%
      dplyr::left_join(topImpoFeats, by = 'Feature') %>%
      dplyr::left_join(smpAnnoTab, by = 'Sample') %>%
      dplyr::mutate(Recur_Condi = paste0(Recurrence, '_', Condition),
                    Recur_Condi = factor(Recur_Condi, levels = gpLevel))
    p2 <- ggplot(topImpoFeatDat, aes(x=Recur_Condi, y=Abundance, col=Condition, fill=Recurrence)) +
      geom_boxplot(alpha = 1, outlier.shape = NA, linewidth = 1) +
      geom_jitter(position = position_jitter(0.2), size = 2, show.legend = F) +
      ggpubr::stat_compare_means(method = 't.test', paired = F, method.args = list(var.equal = T),
                                 comparisons = comparisons, label = 'p.signif', tip.length = 0.015,
                                 bracket.size = 0.7, size = 4) +
      labs(x = 'Recurrence') +
      scale_color_manual(values = condiCol) +
      scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
      facet_wrap(vars(Feat4Viz), scales = 'free') +
      th +
      theme(strip.text = element_text(size = 13, face = 'bold'),
            axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  } else {
    p2 <- NULL
  }
  
  return(list(fullImpoFeatTab = impoFeatTab, pick = p1, abun = p2))
}




vizSysTrainModelPerf <- function(sysTrainRes, numTopFeats = 30) {
  #' Visualize performance of models systematically trained with different numbers
  #' of top important features summarized by bootstrapped model training
  #' 
  #' sysTrainRes: An output of function 'doSysTrainLogisR'
  #' numTopFeats: A numeric value meaning that trained models until with the certain
  #' number of top important features will be shown. Default is 30
  
  summPerfTab <- sysTrainRes$summPerformanceTab
  sysModelTab4Viz <- dplyr::filter(summPerfTab, FeatComb %in% paste0('Top', 2:numTopFeats)) %>%
    dplyr::mutate(lower = stringr::str_extract(CI95, '0.\\d*'),
                  lower = as.numeric(lower),
                  upper = stringr::str_extract(CI95, ',0.\\d*|,1'),
                  upper = stringr::str_remove(upper, '^,'),
                  upper = as.numeric(upper),
                  FeatComb = stringr::str_remove(FeatComb, '^Top'),
                  FeatComb = factor(FeatComb, levels = FeatComb))
  
  ggplot(sysModelTab4Viz, aes(x=FeatComb, y=MeanAUC)) +
    geom_bar(stat = 'identity', position = position_dodge(), alpha = 0.9) +
    geom_errorbar(aes(ymin=lower, ymax=upper), width = 0.5) +
    labs(x = 'Number of top important features', y = 'Mean AUC') +
    th +
    theme(axis.title = element_text(size = 24), axis.text = element_text(size = 22),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}

vizCoreFeats <- function(lrRes, sysTrainRes, se, heatmap_numCoreFeats, scatter_twoCoreFeats,
                         colFeatAnno = NULL, rmUnwantedSmp = NULL, ...) {
  #' Make heatmap and scatterplot to have further investigation of core features
  #' that give rise to models with best performance. Heatmap displays abundances
  #' of core features, where samples are ordered by average linear regression scores
  #' and features are ordered by importance. Scatterplot visualizes predictive power
  #' of selected two core features
  #' 
  #' lrRes: An output of function 'iterLogisR'
  #' sysTrainRes: An output of function 'doSysTrainLogisR'
  #' se: A SummarizedExperiment object used as the input to function 'iterLogisR'
  #' heatmap_numCoreFeats: A numeric value specifying the number of top important
  #' features to display in a heatmap
  #' scatter_twoCoreFeats: A length-2 vector of characters indicating two core features
  #' to use to make a scatterplot. Note that the argument of this parameter depends
  #' on parameter 'colFeatAnno'. If colFeatAnno is NULL, feature names should be
  #' chosen from the original feature space of an input SE object (se). If colFeatAnno
  #' is specified, feature names should be chosen from the annotations of the original
  #' feature space, i.e., the column selected in rowData(se)
  #' colFeatAnno: A character specifying the column name of feature annotations
  #' stored in rowData(se). Default is NULL
  #' rmUnwantedSmp: A character or a vector of characters indicating the sample(s)
  #' to remove from visualization. Default is NULL
  #' ...: Further arguments to be passed to 'pheatmap', e.g., 'fontsize'
  
  # Prepare top important feature list
  # Pinpoint useless models that will be removed
  uselessModels <- which(lrRes$lrRes$nNonZero == 0)
  # Prepare coefficient table
  coefTab <- lrRes$lrRes$coefficient
  if (length(uselessModels) != 0) {
    message(paste(length(uselessModels), 'models with no non-zero coefficient, i.e., useless models.'))
    coefTab <- coefTab[, -uselessModels]
  }
  # Summarize feature selection results into pick rates and rank features by pick rates
  rankedImpoFeats <- as.data.frame(coefTab) %>%
    dplyr::mutate(across(everything(), ~ case_when(.x != 0 ~ 1,
                                                   .x == 0 ~ 0))) %>%
    tibble::rownames_to_column('Feature') %>%
    tidyr::pivot_longer(cols = -'Feature', names_to = 'Model', values_to = 'Pick') %>%
    dplyr::group_by(Feature) %>%
    dplyr::summarise(PickRate = sum(Pick)/ncol(coefTab)) %>%
    dplyr::filter(PickRate != 0) %>%
    dplyr::arrange(desc(PickRate)) %>%
    dplyr::pull(Feature)
  
  # Prepare average linear regression scores of samples
  # Retrieve sample scores of best models
  lrResBest <- sysTrainRes$allTrainRes[[paste0('Top', heatmap_numCoreFeats)]]
  yPredLinkList <- lrResBest$y_pred_link
  # Compute average sample scores
  yPredLinkTab <- data.frame(Sample = names(yPredLinkList[[1]]),
                             Score = yPredLinkList[[1]],
                             row.names = NULL)
  for (i in 2:length(yPredLinkList)) {
    yPredLink <- yPredLinkList[[i]]
    tmp_yPredLinkTab <- data.frame(Sample = names(yPredLinkList[[i]]),
                                   Score = yPredLinkList[[i]],
                                   row.names = NULL)
    yPredLinkTab <- dplyr::bind_rows(yPredLinkTab, tmp_yPredLinkTab)
  }
  yPredLinkTab <- dplyr::group_by(yPredLinkTab, Sample) %>%
    dplyr::summarise(Score = mean(Score)) %>%
    dplyr::arrange(dplyr::desc(Score))
  
  # Make heatmap showing abundances of core features that give rise to best models
  # Prepare data matrix where features are ordered by importance and samples are
  # ordered by average sample scores
  # Remove unwanted samples from visualization if provided
  if (!is.null(rmUnwantedSmp)) {
    yPredLinkTab <- dplyr::filter(yPredLinkTab, !Sample %in% rmUnwantedSmp)
  }
  seSub <- se[rankedImpoFeats[1:heatmap_numCoreFeats], yPredLinkTab$Sample]
  datMatSub <- assay(seSub)
  # Map features to annotations if provided
  if (!is.null(colFeatAnno)) {
    rownames(datMatSub) <- rowData(seSub)[[colFeatAnno]] %>%
      stringr::str_remove(';.+$')
  }
  # Prepare sample annotations (Score and Recurrence) and colors for Recurrence
  smpAnnoTab <- tibble::as_tibble(colData(seSub), rownames = 'Sample') %>%
    dplyr::select(Sample, Recurrence) %>%
    dplyr::left_join(yPredLinkTab, by = 'Sample') %>%
    tibble::column_to_rownames('Sample')
  smpAnnoCols <- list(Recurrence = c(Yes = '#F8766D', No = '#00BFC4'))
  
  heatmap <- pheatmap::pheatmap(datMatSub, annotation_col = smpAnnoTab, scale = 'row',
                                color = colorRampPalette(c('navy', 'white', 'red'))(100),
                                annotation_colors = smpAnnoCols, cluster_rows = F,
                                cluster_cols = F, silent = T, ...) %>%
    ggplotify::as.ggplot()
  
  
  # Make scatterplot visualizing predictive power of two selected core features
  # Prepare long table with two core features selected and needed information for visualization
  recurInfo <- tibble::as_tibble(colData(seSub), rownames = 'Sample') %>%
    dplyr::select(Sample, Recurrence)
  twoFeatTab4Viz <- summExp2df(seSub, row_id = 'Feature', col_id = 'Sample')
  if (!is.null(colFeatAnno)) {
    twoFeatTab4Viz <- dplyr::filter(twoFeatTab4Viz, .data[[colFeatAnno]] %in% scatter_twoCoreFeats) %>%
      dplyr::select(all_of(colFeatAnno), Sample, Value) %>%
      tidyr::pivot_wider(names_from = colFeatAnno, values_from = 'Value')
  } else {
    twoFeatTab4Viz <- dplyr::filter(twoFeatTab4Viz, Feature %in% scatter_twoCoreFeats) %>%
      dplyr::select(Feature, Sample, Value) %>%
      tidyr::pivot_wider(names_from = 'Feature', values_from = 'Value')
  }
  twoFeatTab4Viz <- dplyr::left_join(twoFeatTab4Viz, recurInfo, by = 'Sample') %>%
    dplyr::mutate(Recurrence = factor(Recurrence, levels = c('Yes', 'No')))
  
  scatter <- ggplot(twoFeatTab4Viz, aes(x=.data[[scatter_twoCoreFeats[1]]],
                                        y=.data[[scatter_twoCoreFeats[2]]],
                                        col=Recurrence)) +
    geom_point(size = 6) +
    th +
    theme(axis.title = element_text(size = 24), axis.text = element_text(size = 22),
          legend.title = element_text(size = 24), legend.text = element_text(size = 22))
  
  return(list(heatmap = heatmap, scatterplot = scatter))
}




calcCI <- function(stats, level = 0.95, bootstrap = F) { #iter = 1000
  #' Compute upper and lower bounds of confidence interval. This function provides
  #' two ways of calculating CI: Bootstrapping or CI formula if data is normally
  #' distributed
  #' 
  #' Parameter
  #' stats: A vector of numeric values presenting data or data statistics
  #' level: A numeric value specify the one commonly used confidence level, which
  #' should be one of 0.95 (default), 0.98, 0.99, 0.90, or 0.80
  #' bootstrap: A logical variable indicating whether to use bootstrapping to calculate
  #' CI. Default is FALSE
  #' iter: A numeric value specifying the number of bootstrap iterations
  #' 
  #' Return
  #' A vector containing upper and lower bounds of CI
  
  # Sanity check
  if (!(level %in% c(0.95, 0.98, 0.99, 0.9, 0.8))) {
    stop("Argument for 'level' should be one of 0.95, 0.98, 0.99, 0.90, or 0.80.")
  }
  
  # Define critical value for calculating CI
  if (!bootstrap) {
    criticalVal <- list('0.95' = 1.96, '0.98' = 2.33, '0.99' = 2.58,
                        '0.9' = 1.65, '0.8' = 1.28)
    z <- criticalVal[[as.character(level)]]
  }
  
  # Bootstrap sample and compute mean of each replication
  if (bootstrap) {
    # meanList <- sapply(seq_len(iter), function(i) {
    #   mean(sample(stats, size = length(stats), replace = T))
    # })
    # Locate lower bound by percentile if data is not normally distributed
    p <- (1-level) / 2
    lower <- as.numeric(quantile(stats, p))
    # Locate upper bound
    p <- level + (1-level)/2
    upper <- as.numeric(quantile(stats, p))
  } else {
    # Follow CI formula if data is normally distributed
    lower <- mean(stats) - z * sd(stats)/sqrt(length(stats))
    upper <- mean(stats) + z * sd(stats)/sqrt(length(stats))
  }
  
  return(c(lower = lower, upper = upper))
}

summarizePredPower <- function(pred, truth, auc_roc) {
  #' Summarize performance of trained models, computing average and confidence intervals
  #' of scores obtained from all models
  #' 
  #' pred: A list of class predictions of samples
  #' truth: A list of ground truth labels of samples
  #' auc_roc: A vector of AUC-ROC scores of trained models
  
  # Sanity check
  if (length(pred) != length(truth)) {
    stop("pred and truth should be same length.")
  }
  
  # prepare containers for storing scores
  accu <- c()
  sens <- c()
  spec <- c()
  prec <- c()
  reca <- c()
  f1 <- c()
  # Compute scores from trained models
  for (i in seq_len(length(pred))) {
    # Check if ground truth labels have two levels
    if (length(levels(truth[[i]])) == 2) {
      scores <- caret::confusionMatrix(data = pred[[i]], reference = truth[[i]], positive = '1')
      accu <- c(accu, scores$overall[['Accuracy']])
      sens <- c(sens, scores$byClass[['Sensitivity']])
      spec <- c(spec, scores$byClass[['Specificity']])
      prec <- c(prec, scores$byClass[['Precision']])
      reca <- c(reca, scores$byClass[['Recall']])
      f1 <- c(f1, scores$byClass[['F1']])
    }
  }
  auc <- auc_roc
  # Convert NA score to 0, which means model got bad performance (e.g., model does
  # not have any positive prediction, than precision will be NA)
  prec[is.na(prec)] <- 0
  f1[is.na(f1)] <- 0
  if (any(is.na(auc))) {
    auc <- auc[-which(is.na(auc))]
  }
  
  # Report summarized scores
  # ci <- calcCI(accu, bootstrap = T)
  # cat('Mean Accuracy: ', round(mean(accu), 3), '\n',
  #     '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  # cat('\n')
  # cat('\n')
  # ci <- calcCI(sens, bootstrap = T)
  # cat('Mean Sensitivity: ', round(mean(sens), 3), '\n',
  #     '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  # cat('\n')
  # cat('\n')
  # ci <- calcCI(spec, bootstrap = T)
  # cat('Mean Specificity: ', round(mean(spec), 3), '\n',
  #     '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  # cat('\n')
  # cat('\n')
  # ci <- calcCI(prec, bootstrap = T)
  # cat('Mean Precision: ', round(mean(prec), 3), '\n',
  #     '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  # cat('\n')
  # cat('\n')
  # ci <- calcCI(reca, bootstrap = T)
  # cat('Mean Recall: ', round(mean(reca), 3), '\n',
  #     '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  # cat('\n')
  # cat('\n')
  # ci <- calcCI(f1, bootstrap = T)
  # cat('Mean F1 score: ', round(mean(f1), 3), '\n',
  #     '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  # cat('\n')
  # cat('\n')
  ci <- calcCI(auc, bootstrap = T)
  cat('Mean AUC-ROC score: ', round(mean(auc), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  
  # Show distributions of scores
  # par(mfrow = c(3, 3))
  # hist(accu, xlab = 'Accuracy', main = 'Histogram of Accuracy')
  # hist(sens, xlab = 'Sensitivity', main = 'Histogram of Sensitivity')
  # hist(spec, xlab = 'Specificity', main = 'Histogram of Specificity')
  # hist(prec, xlab = 'Precision', main = 'Histogram of Precision')
  # hist(reca, xlab = 'Recall', main = 'Histogram of Recall')
  # hist(f1, xlab = 'F1 scores', main = 'Histogram of F1 scores')
  hist(auc_roc, xlab = 'AUC-ROC', main = 'Histogram of AUC-ROC')
}
