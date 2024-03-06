#### Try to implement well-developed package, e.g., 'mlr3' or 'caret', to train ML models

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
                  trainSet_ratio = 0.8, ntree = 10000, plot_ROC = F, save_RF = F) {
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
  #' save_RF: A logical variable indicating whether trained RF is saved. Default is FALSE
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
      mtry <- tuneRes[order(tuneRes[, 2]),][1, 1]
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
    if (save_RF) {
      rfResList[[i]] <- rfRes
    }
    
    # Save predictions and ground truths of test set samples
    y_pred <- as.numeric(y_pred > 0.5) %>%
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


runXGBoost <- function(x, y, targetClass, iter = 1, booster = 'gbtree', nrounds = 100,
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
                       eval_metric = 'error',
                       eval_metric = 'logloss',
                       max_depth = maxDepth,
                       eta = eta)
  } else if (booster == 'gblinear') {
    parameters <- list(booster = 'gblinear',
                       objective = 'binary:logistic',
                       eval_metric = 'error',
                       eval_metric = 'logloss')
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
      featImpoTab[, i] <- featImpo[colnames(x),]
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
      dplyr::summarise(Occurrence = length(Feature), ImpoMean = mean(Gain), ImpoSD = sd(Gain)) %>%
      dplyr::arrange(dplyr::desc(Occurrence))
  }
  
  return(list(auc_roc = aucList, featImpoTab = featImpoTab, best_nrounds = bestnroundList,
              params = paramList, xgbRes = xgbResList, y_pred = yPredList, y_truth = yTruthList))
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
  #' iter: A numeric value specifying the number of time to run random forest
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
  #' auc: A vector of computed AUC from trained models
  #' params: A list containing arguments of parameters 'regularized_method', 'cvFold',
  #' 'used_lambda', 'split_method', and 'trainSet_ratio'
  #' yPredList: A list of class predictions of test set samples
  #' yTruthList: A list of ground truth labels of test set samples
  
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
  }
  
  return(list(coefficient = matCoeffi, usedLambda = lambdaList, auc_roc = aucList,
              nNonZero = numNonZeroList, params = paramList, lrRes = lrResList,
              y_pred = yPredList, y_truth = yTruthList))
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
  ci <- calcCI(accu, bootstrap = T)
  cat('Averaged Accuracy: ', round(median(accu), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  cat('\n')
  cat('\n')
  ci <- calcCI(sens, bootstrap = T)
  cat('Averaged Sensitivity: ', round(median(sens), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  cat('\n')
  cat('\n')
  ci <- calcCI(spec, bootstrap = T)
  cat('Averaged Specificity: ', round(median(spec), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  cat('\n')
  cat('\n')
  ci <- calcCI(prec, bootstrap = T)
  cat('Averaged Precision: ', round(median(prec), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  cat('\n')
  cat('\n')
  ci <- calcCI(reca, bootstrap = T)
  cat('Averaged Recall: ', round(median(reca), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  cat('\n')
  cat('\n')
  ci <- calcCI(f1, bootstrap = T)
  cat('Averaged F1 score: ', round(median(f1), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
  cat('\n')
  cat('\n')
  ci <- calcCI(auc, bootstrap = T)
  cat('Averaged AUC-ROC score: ', round(median(auc), 3), '\n',
      '95% CI: [', ci[1], ', ', ci[2], ']', sep = '')
}
