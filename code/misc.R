library(pcaMethods)
library(limma)
library(proDA)
library(sva)
library(SummarizedExperiment)
library(tidyverse)

df2SummExp <- function(df, row_id, col_id, values,
                       row_anno = NULL, col_anno = NULL) {
  #' Convert tidy table to SummarizedExperiment object
  #'
  #' Parameters
  #' df: Input tidy table
  #' row_id: A character variable specifying the column in the tidy table that
  #' will be used as row identifiers in the SE object
  #' col_id: A character variable specifying the column in the tidy table that
  #' will be used as column identifiers in the SE object
  #' values: A character or a vector of characters, specifying the columns that
  #' will be used as measurements (assays)
  #' row_anno: A character or a vector of characters, specifying the columns that
  #' will serve as annotations for row identifiers
  #' col_anno: A character or a vector of characters, specifying the columns that
  #' will serve as annotations for column identifiers
  #'
  #' Return
  #' summExp: A SummarizedExperiment object
  
  # Sanity check
  if (!is(row_id, 'character')) {
    stop("'row_id' has to be character")
  }
  if (!is(col_id, 'character')) {
    stop("'col_id' has to be character")
  }
  
  # Prepare value matrix, converting long data to wide data
  matList <- lapply(values, function(x) {
    dplyr::select(df, all_of(row_id), all_of(col_id), all_of(x)) %>%
      tidyr::pivot_wider(names_from = all_of(col_id), values_from = all_of(x)) %>%
      tibble::column_to_rownames(row_id)
  })
  names(matList) <- values
  
  # Prepare row annotation
  if (!is.null(row_anno)) {
    rowAnno <- df[, c(row_id, row_anno), drop = F]
    rowAnno <- rowAnno[!duplicated(rowAnno[[row_id]]),, drop = F] %>%
      tibble::column_to_rownames(row_id) 
  } else {
    rowAnno <- NULL
    message('Message: There is no row annotation.')
  }
  # Prepare column annotation
  if (!is.null(col_anno)) {
    colAnno <- df[, c(col_id, col_anno), drop = F]
    colAnno <- colAnno[!duplicated(colAnno[[col_id]]),, drop = F] %>%
      tibble::column_to_rownames(col_id)
  } else {
    colAnno <- NULL
    message('Message: There is no column annotation.')
  }
  # Assemble all information into SummarizedExperiment object
  if (!is.null(colAnno)) {
    summExp <- SummarizedExperiment(assays = matList,
                                    rowData = rowAnno,
                                    colData = colAnno)
  } else {
    summExp <- SummarizedExperiment(assays = matList,
                                    rowData = rowAnno)
  }
  return(summExp)
}

summExp2df <- function(summ_exp, assay = NULL,
                       row_id = 'Row_ID', col_id = 'Col_ID') {
  #' Convert SummarizedExperiment object to tidy table
  #' 
  #' Parameters
  #' summ_exp: Input SummarizedExperiment object
  #' assay: A character or an index specifying the assay in the SummarizedExperiment
  #' object that will be used
  #' row_id: A character specifying the name for the row identifiers of the assay
  #' col_id: A character specifying the name for the column identifiers of the assay
  #' 
  #' Return
  #' df: A data frame
  
  # Prepare long data converted from value matrix
  if (is.null(assay)) {
    df <- SummarizedExperiment::assay(summ_exp)
  } else {
    df <- SummarizedExperiment::assays(summ_exp)[[assay]]
  }
  df <- tibble::as_tibble(df, rownames = row_id) %>%
    tidyr::pivot_longer(cols = -all_of(row_id),
                        names_to = col_id,
                        values_to = 'Value')
  # Append row annotations
  rowAnno <- rowData(summ_exp) %>%
    # as.data.frame() %>%
    # tibble::rownames_to_column(row_id) %>%
    tibble::as_tibble(rownames = row_id)
  df <- dplyr::left_join(df, rowAnno, by = row_id)
  # Append column annotations
  colAnno <- colData(summ_exp) %>%
    tibble::as_tibble(rownames = col_id)
  df <- dplyr::left_join(df, colAnno, by = col_id)
  return(df)
}


doPCA <- function(data, pca_method = 'pca', num_PCs = 20) {
  #' Perform conventional PCA or its extensions, probabilistic PCA and Bayesian
  #' PCA, on MS-based data such as proteomics, metabolomics, etc. If missing values
  #' exist, those features with any missingness are removed for running conventional
  #' PCA, leading to loss of information. In this case, PPCA and BPCA can be considered.
  #' 
  #' Parameter
  #' data: A vsn (variance stabilization) normalized Feature x Sample abundance
  #' data matrix since parameter 'scale' is not used for all sorts of PCA methods
  #' pca_method: A character specifying the PCA method to use, which should be one
  #' of 'pca' (default), 'ppca', or 'bpca'
  #' num_PCs: A numeric value specifying the number of PCs to estimate. The preciseness
  #' of the estimation of missing values depends on the number of PCs
  #' 
  #' Return
  #' pcaRes: A list or object containing the components described in 'prcomp' or
  #' 'pcaMethods::pcaRes' documentation
  
  if (pca_method == 'pca') {
    # Run conventional PCA
    # center = T => centering is done by subtracting corresponding column (feature)
    # means of data, which shifts features to be zero centered
    # scale. = T => scaling is done by dividing standard deviation of (centered)
    # columns, which scales features to have unit variance
    # Spot complete features to remove features with any missing values
    completeFeats <- complete.cases(data)
    if (!all(completeFeats)) {
      message("Message: Missing values exist and are removed for running conventional PCA.")
    }
    pcaRes <- prcomp(t(data[completeFeats,]), center = T, scale. = F)
  } else if (pca_method == 'ppca') {
    # Run probabilistic PCA
    # Convert missing values from NaN to NA for performing PCA in 'pcaMethods' package
    data[is.nan(data)] <- NA
    # Remove features fully missing
    rmFeats <- which(apply(data, 1, function(featVec) {all(is.na(featVec))}))
    if (length(rmFeats) != 0) {
      data <- data[-rmFeats,]
    }
    pcaRes <- pcaMethods::pca(t(data), method = 'ppca', nPcs = num_PCs, seed = 123,
                              center = T, scale = 'none')
  } else if (pca_method == 'bpca') {
    # Run Bayesian PCA
    data[is.nan(data)] <- NA
    rmFeats <- which(apply(data, 1, function(featVec) {all(is.na(featVec))}))
    if (length(rmFeats) != 0) {
      data <- data[-rmFeats,]
    }
    pcaRes <- pcaMethods::pca(t(data), method = 'bpca', nPcs = num_PCs,
                              center = T, scale = 'none')
  } else {
    stop("Argument for 'pca_method' should be one of 'pca', 'ppca', or 'bpca'.")
  }
  return(pcaRes)
}


testAssociation <- function(data, metadata, cmn_col, use_limma = F, use_proDA = F,
                            targetVar = NULL, unwantVar = NULL, cor_method = 'pearson') {
  #' Iterate association tests between dependent and independent variables whose
  #' information is stored in data and metadata tables, respectively. Data table
  #' should contain numerical values of samples (e.g., protein abundances or PCs)
  #' and metadata table contains sample's metadata (Patient genders or ages).
  #' Depending on types of independent variables, varied tests are performed. That
  #' is, if independent variable is numerical, correlation is computed; if independent
  #' variable is categorical, t-test or ANOVA is used. Three types of t-test are
  #' provided: Ordinary Student's t-test by lm (default), moderated t-test by limma,
  #' and t-test by proDA.
  #'
  #' Parameters
  #' data: A normalized Sample x Feature numerical data frame
  #' metadata: A data frame containing sample metadata in the column
  #' cmn_col: A character specifying the common column (e.g., sample names) between
  #' data and metadata tables for combining them to match the information
  #' use_limma: A logical variable indicating whether to use limma to conduct moderated
  #' t-test. Default is FALSE
  #' use_proDA: A logical variable indicating whether to use proDA to account for
  #' missing values for obtaining more reliable t-test results. Default is FALSE
  #' targetVar: A character or a vector of characters specifying the metadata variables
  #' whose associations will be tested. If NULL (default), all variables in metadata
  #' will be tested
  #' unwantVar: A character or a vector of characters specifying the metadata variables
  #' that will be accounted for by linear models (regressing out). Default is NULL.
  #' Note that this is usable only for limma and proDA for now
  #' cor_method: A character specifying the correlation method to use, which should
  #' be one of 'pearson' (default), 'spearman', or 'kendall'
  #' 
  #' Return
  #' resTab: A data frame collecting all association test results
  
  # Check arguments
  if (!(is(data, 'data.frame') & is(metadata, 'data.frame'))) {
    stop("Both input tables should be class 'data.frame', 'tbl', or 'tbl_df'.")
  }
  if (!is.character(cmn_col) | length(cmn_col) != 1) {
    stop("Argument for 'cmn_col' should be class 'character' and length-1.")
  }
  if (!(cor_method %in% c('pearson', 'spearman', 'kendall'))) {
    stop("Argument for 'cor_method' should be one of 'pearson', 'spearman', or 'kendall'.")
  }
  if (!(is.logical(use_limma) & is.logical(use_proDA))) {
    stop("Argument for 'use_limma' and 'use_proDA' should be class 'logical'.")
  }
  if (!is.null(targetVar) | !is.null(unwantVar)) {
    if (!all(c(targetVar, unwantVar) %in% colnames(metadata))) {
      stop("Arguments for 'targetVar' and 'unwantVar' should be included in metadata.")
    }
  }
  # Check completeness of metadata and data variable type (for myself aware of everything)
  completeMetadata <- complete.cases(metadata)
  if (!all(completeMetadata)) {
    message("Message: Not all sample's metadata is complete.")
  }
  varTypeData <- apply(data[,-which(colnames(data) == cmn_col)], 2, is.numeric)
  if (!all(varTypeData)) {
    stop("Not all data variables are numeric.")
  }
  
  # Keep all metadata for conducting regressing out
  if (!is.null(unwantVar)) {
    allMetadata <- tibble::column_to_rownames(metadata, cmn_col)
    allMetadata <- allMetadata[data[[cmn_col]],]
  }
  # Combine data and metadata by a common column to ensure matched information
  # Extract metadata variables of interest if 'targetVar' is specified
  if (!is.null(targetVar)) {
    metadata <- dplyr::select(metadata, all_of(cmn_col), all_of(targetVar))
  }
  combinedTab <- dplyr::left_join(data, metadata, by = cmn_col) %>%
    tibble::column_to_rownames(cmn_col)
  # Conduct all possible association tests
  numVar1 <- ncol(data) - 1
  numVar2 <- ncol(metadata) - 1
  resTab <- lapply(seq_len(numVar1), function(i) {
    lapply(seq_len(numVar2), function(j) {
      var1 <- combinedTab[[i]]
      var2 <- combinedTab[[numVar1+j]]
      # Remove samples with missing metadata
      missVal <- which(is.na(var2))
      if (length(missVal) != 0) {
        var1 <- var1[-missVal]
        var2 <- var2[-missVal]
      }
      # Perform corresponding tests depending on types of independent variables
      singleRes <- try(
        if (is.numeric(var2)) {
          # Perform correlation test if independent variable is numerical
          test <- 'Correlation'
          corrRes <- cor.test(var1, var2, method = cor_method,
                              use = "pairwise.complete.obs")
          pValue <- corrRes$p.value
          statistic <- corrRes$estimate
        } else if ((is.character(var2) | is.factor(var2)) &
                   length(unique(var2)) == 2) {
          # Perform Student's t-test if independent variable is categorical and has 2 levels
          if (!(use_limma | use_proDA)) {
            test <- 'T-test'
            tStatRes <- summary(lm(var1 ~ var2))$coefficients
            pValue <- tStatRes[2, 4]
            statistic <- tStatRes[2, 3]
          }
        } else if ((is.character(var2) | is.factor(var2)) &
                   length(unique(var2)) > 2) {
          # Perform ANOVA if independent variable is categorical and has more than 2 levels
          test <- 'ANOVA'
          fStatRes <- summary(aov(var1 ~ var2))
          pValue <- fStatRes[[1]][['Pr(>F)']][1]
          statistic <- fStatRes[[1]][['F value']][1]
        },
        silent = T
      )
      # Create data frame for summarizing results
      if (!is(singleRes, 'try-error') & !is.null(singleRes)) {
        # Prevent error in creating data frame due to differing number of rows,
        # i.e., p-value and statistic are NULL and fail to be retrieved (e.g., 20
        # samples and 20 groups in ANOVA)
        if (is.null(statistic)) {
          pValue <- NA
          statistic <- NA
        }
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = pValue, Stat = statistic, Test = test,
                   stringsAsFactors = F, row.names = c())
      } else if (is(singleRes, 'try-error')) {
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = NA, Stat = NA, Test = test,
                   stringsAsFactors = F, row.names = c())
      }
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  # Use limma to conduct moderated t-test
  if (use_limma) {
    # Prepare data matrix for fitting linear model
    datMat <- combinedTab[, seq_len(numVar1), drop = F] %>%
      t()
    resTab <- lapply(seq_len(numVar2), function(i) {
      metadatVar <- colnames(combinedTab)[numVar1+i]
      metadatVarVec <- combinedTab[[metadatVar]]
      if ((is.character(metadatVarVec) | is.factor(metadatVarVec)) &
          length(unique(metadatVarVec)) == 2) {
        if (is.null(unwantVar)) {
          design <- model.matrix(~ metadatVarVec)
        } else {
          formu <- paste0('~', paste0(c(metadatVar, unwantVar), collapse = '+')) %>%
            as.formula()
          design <- model.matrix(formu, data = allMetadata)
        }
        fit <- limma::lmFit(datMat, design = design)
        fit <- limma::eBayes(fit)
        tStatRes <- limma::topTable(fit, coef = 2, number = Inf)
        data.frame(Var1 = rownames(tStatRes),
                   Var2 = metadatVar,
                   pVal = tStatRes$P.Value, Stat = tStatRes$t,
                   Test = 'T-test (limma)')
      }
    }) %>% dplyr::bind_rows() %>%
      dplyr::bind_rows(resTab)
  }
  
  # Use proDA to account for missing values for more reliable t-test
  if (use_proDA) {
    # Prepare data matrix for fitting linear probabilistic dropout model
    datMat <- combinedTab[, seq_len(numVar1), drop = F] %>%
      t()
    resTab <- lapply(seq_len(numVar2), function(i) {
      metadatVar <- colnames(combinedTab)[numVar1+i]
      metadatVarVec <- combinedTab[[metadatVar]]
      if ((is.character(metadatVarVec) | is.factor(metadatVarVec)) &
          length(unique(metadatVarVec)) == 2) {
        if (is.null(unwantVar)) {
          design <- model.matrix(~ metadatVarVec)
        } else {
          formu <- paste0('~', paste0(c(metadatVar, unwantVar), collapse = '+')) %>%
            as.formula()
          design <- model.matrix(formu, data = allMetadata)
        }
        fit <- proDA(datMat, design = design)
        diffTerm <- proDA::result_names(fit)[2]
        tStatRes <- proDA::test_diff(fit, contrast = diffTerm)
        data.frame(Var1 = tStatRes$name,
                   Var2 = colnames(combinedTab)[numVar1+i],
                   pVal = tStatRes$pval, Stat = tStatRes$t_statistic,
                   Test = 'T-test (proDA)', stringsAsFactors = F)
      }
    }) %>% dplyr::bind_rows() %>%
      dplyr::bind_rows(resTab)
  }
  
  # Compute adjusted p-value
  resTab <- dplyr::mutate(resTab, pValAdj = p.adjust(pVal, method = 'BH'),
                          pVal = as.numeric(formatC(pVal, format = 'g', digits = 3)),
                          Stat = as.numeric(formatC(Stat, format = 'g', digits = 3)),
                          pValAdj = as.numeric(formatC(pValAdj, format = 'g', digits = 3))) %>%
    dplyr::relocate(pValAdj, .after = pVal) %>%
    dplyr::arrange(pVal)
  
  return(resTab)
}


doSOA <- function(summ_exp, meta_var, pca_method = 'pca', num_PCs = 20, alpha = 0.05,
                  num_PCfeats = NULL, do_onlyPCA = F, use_limma = F, use_proDA = F,
                  unwantVar = NULL, level_metaVar = NULL, ...) {
  # To-do: Expand function to e.g., MultiAssayExperiment object
  
  #' Perform single-omics analysis on preprocessed data stored in SE container to
  #' have initial look at data properties and data power in terms of certain scientific
  #' questions, e.g., identifying biomarkers for predicting patient cancer recurrences.
  #' This function mainly includes PCA (multivariate) and univariate association
  #' test between two variables, which attempts to capture significant PCs and features.
  #' 
  #' Parameters
  #' summ_exp: A SummarizedExperiment object containing normalized Feature x Sample
  #' data and metadata
  #' meta_var: A character or a vector of characters indicating the variable(S)
  #' in the sample metadata from an SE object for association tests. Categorical
  #' variable (factor) is used as grouping factor to separate samples into two or
  #' several groups for conducting t-test or ANOVA; numerical variable (covariate)
  #' is used to compute correlation coefficient
  #' pca_method: A character specifying the PCA method to use, which should be one
  #' of 'pca' (default), 'ppca', or 'bpca'
  #' num_PCs: A numeric value specifying the number of PCs to estimate. The preciseness
  #' of the estimation of missing values depends on the number of PCs
  #' alpha: A numeric value specifying the cutoff for statistic significance. Default
  #' is 0.05
  #' num_PCfeats: A numeric value specifying the number of potential PCs' top features
  #' with highest absolute loadings to extract. If NULL (default), the extraction
  #' of lists of features is skipped
  #' do_onlyPCA: A logical variable indicating whether only PCA and its association
  #' tests between PCs and metadata variables are performed. If FALSE (default),
  #' univariate association tests between features (e.g., protein abundances) and
  #' metadata variables are also conducted
  #' use_limma: A logical variable indicating whether to use limma to conduct moderated
  #' t-test for univariate association tests. Default is FALSE
  #' use_proDA: A logical variable indicating whether to use proDA to account for
  #' missing values for obtaining more reliable t-test results for univariate association
  #' tests. Default is FALSE
  #' unwantVar: A character or a vector of characters specifying the metadata variables
  #' that will be accounted for by linear models (regressing out). Default is NULL.
  #' Note that this is usable only for limma and proDA for now
  #' level_metaVar: A vector of characters specifying the level of the variable
  #' 'meta_var' for determining the direction of a comparison. For example, cancer
  #' recurrence got two levels: Yes and No, so 'level_metaVar' will accept c('Yes', 'No')
  #' when 'meta_var' is 'Recurrence', meaning that t-test compares recurrence vs
  #' non-recurrence patient groups (R. - NonR.), and vice versa. Default is NULL
  #' ...: Further arguments to be passed to 'testAssociation', e.g., correlation
  #' method 'cor_method'
  #' 
  #' Return
  #' A list containing the following components:
  #' data: A data matrix from the SE object
  #' dataCorrect: A data matrix that is adjusted for unwanted variance according
  #' to parameter 'unwantVar'
  #' smpMetadata: A table storing sample metadata, e.g., patient genders, ages, etc.
  #' pcaRes: A complete PCA result obtained using 'prcomp' or 'pcaMethods::pca'
  #' pcSigAssoRes: A table showing significant association test results between PCs
  #' and metadata variables (defined by parameter 'meta_var'). Percentages after
  #' PC are variance explained
  #' pcTopFeatTab: A list containing lists of top features (determined by parameter
  #' 'num_PCfeats') of potential PCs
  #' pcTab: A PC table containing sample representations and metadata. Percentages
  #' after PC are variance explained
  #' featSigAssoRes: A table showing significant association test results between
  #' features and metadata variables (defined by parameter 'meta_var')
  #' featAssoRes: A table showing all association test results between features and
  #' metadata variables (defined by parameter 'meta_var')
  
  # Check arguments
  if (!is(summ_exp, 'SummarizedExperiment')) {
    stop("This function takes only 'SummarizedExperiment' object for now.")
  }
  if (!is(meta_var, 'character') | !all(meta_var %in% colnames(colData(summ_exp)))) {
    stop("Argument for 'meta_var' should be class 'character' and included in sample metadata.")
  }
  if (!is.null(unwantVar)) {
    if (!is(unwantVar, 'character') | !all(unwantVar %in% colnames(colData(summ_exp)))) {
      stop("Argument for 'unwantVar' should be class 'character' and included in sample metadata.")
    }
  }
  if (!(pca_method %in% c('pca', 'ppca', 'bpca'))) {
    stop("Argument for 'pca_method' should be one of 'pca', 'ppca', or 'bpca'.")
  }
  if (!is(do_onlyPCA, 'logical')) {
    stop("Argument for 'do_onlyPCA' should be class 'logical'.")
  }
  
  # Extract data matrix from SE object
  datMat <- SummarizedExperiment::assay(summ_exp) %>%
    as.matrix()
  # Extract sample metadata from SE object and ensure 'Sample' column exists and
  # contains sample names of data matrix for further association tests
  smpMetadat <- colData(summ_exp)
  if (!'Sample' %in% colnames(smpMetadat)) {
    smpMetadat <- tibble::as_tibble(smpMetadat, rownames = 'Sample')
  } else {
    smpMetadat <- tibble::as_tibble(smpMetadat, rownames = 'tmp_Sample') %>%
      dplyr::select(-Sample) %>%
      dplyr::rename(Sample = tmp_Sample)
  }
  # Level metadata variable of interest if parameter 'level_metaVar' is specified
  # for determining direction of comparison
  if (!is.null(level_metaVar)) {
    smpMetadat[[meta_var]] <- factor(smpMetadat[[meta_var]], levels = rev(level_metaVar))
  }
  # Correct data for unwanted variance
  if (!is.null(unwantVar)) {
    # Define design matrix to keep particular effects, e.g., treatments, and adjust
    # unwanted effects
    formu <- paste0('~', paste0(c(meta_var, unwantVar), collapse = '+')) %>%
      as.formula()
    design <- model.matrix(formu, data = smpMetadat)
    # Determine number of columns that is of interest in design matrix
    wantedColCount <- 0
    for (i in seq_len(length(meta_var))) {
      smpMetaVar <- smpMetadat[[meta_var[i]]]
      clas <- class(smpMetaVar)
      if (!is.numeric(clas)) {
        wantedColCount <- wantedColCount + length(unique(smpMetaVar)) - 1
      } else {
        wantedColCount <- wantedColCount + 1
      }
    }
    wantedDesign <- design[, seq_len(wantedColCount+1)] #+1 due to intercept term
    unwantedDesign <- design[, -seq_len(wantedColCount+1)]
    # Do correction
    datMatCorrect <- limma::removeBatchEffect(datMat, design = wantedDesign,
                                              covariates = unwantedDesign)
  } else {
    datMatCorrect <- NULL
  }
  
  # Perform PCA
  if (is.null(unwantVar)) {
    pcaRes <- doPCA(datMat, pca_method = pca_method, num_PCs = num_PCs)
  } else {
    pcaRes <- doPCA(datMatCorrect, pca_method = pca_method, num_PCs = num_PCs)
  }
  # Conduct association test between PCs and factors
  if (pca_method == 'pca') {
    pcTab <- pcaRes$x %>%
      tibble::as_tibble(rownames = 'Sample')
    varExplained <- pcaRes$sdev^2 / sum(pcaRes$sdev^2)
  } else if (pca_method == 'ppca' | pca_method == 'bpca') {
    pcTab <- pcaRes@scores %>%
      tibble::as_tibble(rownames = 'Sample')
    varExplained <- pcaRes@sDev^2 / sum(pcaRes@sDev^2)
  }
  newPCs <- paste0('PC', seq_along(varExplained), ' (', round(varExplained*100, 1), '%)')
  colnames(pcTab) <- c('Sample', newPCs)
  metaTab <- dplyr::select(smpMetadat, c('Sample', all_of(meta_var)))
  pcSigAssoRes <- testAssociation(pcTab, metaTab, cmn_col = 'Sample', use_limma = F,
                                  use_proDA = F, targetVar = NULL, unwantVar = NULL,
                                  cor_method = 'pearson') %>%
    dplyr::filter(pVal <= alpha)
  # Create complete PC table containing sample representations and metadata
  pcTab <- dplyr::left_join(pcTab, smpMetadat, by = 'Sample')
  # Extract top features with highest absolute loadings of significant PCs
  if (!is.null(num_PCfeats)) {
    sigPCs <- stringr::str_remove(pcSigAssoRes[['Var1']], ' \\(.+\\)')
    pcTopFeatTab <- lapply(sigPCs, function(pc) {
      if (pca_method == 'pca') {
        featLoad <- pcaRes$rotation[, pc]
      } else if (pca_method == 'ppca' | pca_method == 'bpca') {
        featLoad <- pcaRes@loadings[, pc]
      }
      topFeatLoadIdx <- order(abs(featLoad), decreasing = T)[1:num_PCfeats]
      topFeatLoad <- featLoad[topFeatLoadIdx]
      topFeatIDs <- names(topFeatLoad)
      data.frame(Feature = topFeatIDs,
                 Loading = as.numeric(formatC(topFeatLoad, format = 'g', digits = 3)),
                 row.names = c())
    })
    names(pcTopFeatTab) <- sigPCs
  } else {
    pcTopFeatTab <- NULL
  }
  
  # Conduct univariate association test between features and factors
  if (do_onlyPCA == F) {
    featTab <- as_tibble(t(datMat), rownames = 'Sample')
    featAssoRes <- testAssociation(featTab, smpMetadat, cmn_col = 'Sample',
                                   use_limma = use_limma, use_proDA = use_proDA,
                                   targetVar = meta_var, unwantVar = unwantVar, ...)
    # Extract features that can significantly separate groups of samples
    featSigAssoRes <- dplyr::filter(featAssoRes, pVal <= alpha)
  } else {
    featAssoRes <- NULL
    featSigAssoRes <- NULL
  }
  
  # Convert metadata variable of interest (factor) back to character in output tables
  if (!is.null(level_metaVar)) {
    smpMetadat[[meta_var]] <- as.character(smpMetadat[[meta_var]])
    pcTab[[meta_var]] <- as.character(pcTab[[meta_var]])
  }
  
  return(list(data = datMat, dataCorrect = datMatCorrect, smpMetadata = smpMetadat,
              pcaRes = pcaRes, pcSigAssoRes = pcSigAssoRes, pcTopFeatTab = pcTopFeatTab,
              pcTab = pcTab, featSigAssoRes = featSigAssoRes, featAssoRes = featAssoRes))
}


doSVA <- function(summExp, wantedVar, unwantedVar = NULL, numSV_method = 'be',
                  asso_metaVar = NULL) {
  #' Do surrogate variable analysis to identify latent factors that explain unknown
  #' variance and correct data using identified SVs
  #' 
  #' Parameters
  #' summExp: An input SummarizedExperiment object containing the data matrix to
  #' do SVA and sample metadata to do association tests (optional)
  #' wantedVar: A character or a vector of characters indicating the variable(S)
  #' in the sample metadata from an SE object, the sources of variation to keep
  #' unwantedVar: A character or a vector of characters indicating the variable(S)
  #' in the sample metadata from an SE object, the sources of variation to adjust.
  #' Default is NULL, which adjusts all identified sources of variation, except
  #' variance of interest specified by the parameter 'wantedVar'
  #' numSV_method: A character specifying the approach to estimate the number of
  #' SVs, which should be 'be' (default) or 'leek'
  #' asso_metaVar: A character or a vector of characters indicating the variable(S)
  #' in the sample metadata from an SE object for association tests between SVs
  #' and sample metadata variables. Default is NULL
  #' 
  #' Return
  #' A list containing the following components:
  #' summExp_SVs: A SummarizedExperiment object identical to the input SE object,
  #' besides identified SVs also included in the sample metadata slot
  #' sigAssoTab: A table showing significant associations between SVs and sample
  #' metadata variables
  #' correctDatMat: A data matrix accounted for SVs
  
  # Do sanity check
  if (!is(summExp, 'SummarizedExperiment')) {
    stop("This function takes only 'SummarizedExperiment' object for now.")
  }
  if (!is(wantedVar, 'character') | !all(wantedVar %in% colnames(colData(summExp)))) {
    stop("Argument for 'wantedVar' should be class 'character' and included in sample metadata.")
  }
  if (!is.null(unwantedVar)) {
    if (!is(unwantedVar, 'character') | !all(unwantedVar %in% colnames(colData(summExp)))) {
      stop("Argument for 'unwantedVar' should be class 'character' and included in sample metadata.")
    }
  }
  if (!is.null(asso_metaVar)) {
    if (!is(asso_metaVar, 'character') | !all(asso_metaVar %in% colnames(colData(summExp)))) {
      stop("Argument for 'asso_metaVar' should be class 'character' and included in sample metadata.")
    }
  }
  
  # Prepare data matrix and sample metadata table for creating model matrices
  datMat <- as.matrix(assay(summExp))
  # Keep only complete features
  keptFeats <- apply(datMat, 1, function(featVec) {
    !any(is.na(featVec))
  })
  datMat <- datMat[keptFeats,]
  metadatTab <- colData(summExp) %>%
    tibble::as_tibble(rownames = 'Sample')
  # Create full model matrix including both adjustment variables and variables of interest
  formu <- paste0('~', paste0(wantedVar, collapse = '+')) %>%
    as.formula()
  modelMat <- model.matrix(formu, data = metadatTab)
  # Create null model matrix including only adjustment variables
  if (is.null(unwantedVar)) {
    modelMat0 <- model.matrix(~ 1, data = metadatTab)
  } else {
    formu <- paste0('~', paste0(unwantedVar, collapse = '+')) %>%
      as.formula()
    modelMat0 <- model.matrix(formu, data = metadatTab)
  }
  
  # Perform SVA
  # Identify number of latent factors that need to be estimated
  numSV = sva::num.sv(datMat, modelMat, method = numSV_method)
  if (numSV == 0) {
    stop("No SV identified...")
  }
  spsUtil::quiet({
    svObj = sva(datMat, modelMat, modelMat0, n.sv = numSV)
  })
  
  # Prepare surrogate variable table
  svTab <- svObj$sv
  rownames(svTab) <- metadatTab$Sample
  colnames(svTab) <- paste0('SV', seq_len(ncol(svTab)))
  svTab <- tibble::as_tibble(svTab, rownames = 'Sample')
  # Do association tests between SVs and patient metadata
  if (!is.null(asso_metaVar)) {
    # Prepare patient metadata table
    metaTab <- dplyr::select(metadatTab, Sample, all_of(asso_metaVar))
    # Do association tests
    sigAssoTab <- testAsso(svTab, metaTab, cmn_col = 'Sample') %>%
      dplyr::filter(pVal <= 0.05)
  } else {
    sigAssoTab <- NULL
  }
  
  # Include surrogate variables into SE object
  longTab <- summExp2df(summExp, row_id = 'Feature', col_id = 'Sample')
  longTab <- dplyr::left_join(longTab, svTab, by = 'Sample')
  summExp <- df2SummExp(longTab, row_id = 'Feature', col_id = 'Sample', values = 'Value',
                        row_anno = colnames(rowData(summExp)),
                        col_anno = c(colnames(colData(summExp)), paste0('SV', seq_len(ncol(svTab)-1))))
  assay(summExp) <- as.matrix(assay(summExp))
  
  # Correct data using identified SVs
  datMatCorrect <- doBatchCorrection(summExp, wantedVar = wantedVar,
                                     unwantedVar = paste0('SV', seq_len(ncol(svTab)-1)))
  
  return(list(summExp_SVs = summExp, sigAssoTab = sigAssoTab, correctDatMat = datMatCorrect))
}


doBatchCorrection <- function(summExp, wantedVar, unwantedVar) {
  #' Regress out unwanted effects from data using function limma::removeBatchEffect
  #' 
  #' Parameters
  #' summExp: An input SummarizedExperiment object containing the data matrix to correct
  #' wantedVar: A character or a vector of characters indicating the variable(S)
  #' in the sample metadata from an SE object, the sources of variation to keep
  #' unwantedVar: A character or a vector of characters indicating the variable(S)
  #' in the sample metadata from an SE object, the sources of variation to adjust
  #' 
  #' Return
  #' A data matrix corrected for unwanted variance
  
  # Do sanity check
  if (!is(summExp, 'SummarizedExperiment')) {
    stop("This function takes only 'SummarizedExperiment' object for now.")
  }
  if (!is(wantedVar, 'character') | !all(wantedVar %in% colnames(colData(summExp)))) {
    stop("Argument for 'wantedVar' should be class 'character' and included in sample metadata.")
  }
  if (!is(unwantedVar, 'character') | !all(unwantedVar %in% colnames(colData(summExp)))) {
    stop("Argument for 'unwantedVar' should be class 'character' and included in sample metadata.")
  }
  
  # Prepare data matrix
  datMat <- as.matrix(assay(summExp))
  # Define design matrix including wanted effects to preserve and unwanted effects to adjust
  metadatTab <- colData(summExp) %>%
    tibble::as_tibble(rownames = 'Sample')
  formu <- paste0('~', paste0(c(wantedVar, unwantedVar), collapse = '+')) %>%
    as.formula()
  design <- model.matrix(formu, data = metadatTab)
  # Determine number of columns that is of interest in design matrix
  wantedColCount <- 0
  for (i in seq_len(length(wantedVar))) {
    smpMetaVar <- metadatTab[[wantedVar[i]]]
    clas <- class(smpMetaVar)
    if (!is.numeric(clas)) {
      wantedColCount <- wantedColCount + length(unique(smpMetaVar)) - 1
    } else {
      wantedColCount <- wantedColCount + 1
    }
  }
  wantedDesign <- design[, seq_len(wantedColCount+1)] #+1 due to intercept term
  unwantedDesign <- design[, -seq_len(wantedColCount+1)]
  # Do correction
  datMatCorrect <- limma::removeBatchEffect(datMat, design = wantedDesign,
                                            covariates = unwantedDesign)
  
  return(datMatCorrect)
}


testAsso <- function(tabX, tabY, cmn_col, cor_method = 'pearson', var_equal = T) {
  #' Systematically perform association tests between variables from two input tables
  #' 
  #' Parameters
  #' cor_method: A character specifying the correlation method to use, which should
  #' be one of 'pearson' (default), 'spearman', or 'kendall'
  #' var_equal: A logical variable indicating whether to treat the variances of
  #' two groups as being equal. If TRUE (default), Student's t-test is performed,
  #' otherwise Welch's t-test is used
  
  # Combine two input tables by a common column to ensure matched information
  combinedTab <- dplyr::left_join(tabX, tabY, by = cmn_col) %>%
    dplyr::select(-all_of(cmn_col))
  # Conduct all possible association tests
  numVar1 <- ncol(tabX) - 1
  numVar2 <- ncol(tabY) - 1
  resTab <- lapply(seq_len(numVar1), function(i) {
    lapply(seq_len(numVar2), function(j) {
      var1 <- combinedTab[[i]]
      var2 <- combinedTab[[numVar1+j]]
      # Stop when any variable has only one level
      if (length(unique(var1)) == 1 | length(unique(var2)) == 1) {
        stop("Please remove categorical variables with only one level.")
      }
      # Perform corresponding tests depending on variable types
      singleRes <- try(
        if (is.numeric(var1) & is.numeric(var2)) {
          # Perform correlation test if both variables are numerical
          # use = 'complete.obs' (casewise) => remove all cases with missing data
          # and calculate correlation on remaining data
          # use = 'pairwise.complete.obs' => remove cases with missing data for each
          # correlation calculation
          test <- 'Correlation'
          corrRes <- cor.test(var1, var2, method = cor_method,
                              use = "pairwise.complete.obs")
          pValue <- corrRes$p.value
          statistic <- corrRes$estimate
        } else if (!is.numeric(var1) & !is.numeric(var2)) {
          # Perform chi-square test if both variables are categorical
          #### Use Fisher's exact test for small sample size? ####
          test <- 'Chi-square'
          chisqRes <- chisq.test(var1, var2)
          pValue <- chisqRes$p.value
          statistic <- chisqRes$statistic
        } else if ((is.numeric(var1) & !is.numeric(var2)) &
                   length(unique(var2)) == 2) {
          # Perform t-test if one variable is numerical and one is categorical and
          # has 2 levels
          test <- 'T-test'
          tStatRes <- t.test(var1 ~ var2, paired = F, var.equal = var_equal)
          pValue <- tStatRes$p.value
          statistic <- tStatRes$statistic
        } else if ((!is.numeric(var1) & is.numeric(var2)) &
                   length(unique(var1)) == 2) {
          # Perform t-test if one variable is numerical and one is categorical and
          # has 2 levels
          test <- 'T-test'
          tStatRes <- t.test(var2 ~ var1, paired = F, var.equal = var_equal)
          pValue <- tStatRes$p.value
          statistic <- tStatRes$statistic
        } else if ((is.numeric(var1) & !is.numeric(var2)) &
                   length(unique(var2)) > 2) {
          # Perform ANOVA if one variable is numerical and one is categorical and
          # has more than 2 levels
          test <- 'ANOVA'
          anovaRes <- summary(aov(var1 ~ var2))
          pValue <- anovaRes[[1]][['Pr(>F)']][1]
          statistic <- anovaRes[[1]][['F value']][1]
        } else if ((!is.numeric(var1) & is.numeric(var2)) &
                   length(unique(var1)) > 2) {
          # Perform ANOVA if one variable is numerical and one is categorical and
          # has more than 2 levels
          test <- 'ANOVA'
          anovaRes <- summary(aov(var2 ~ var1))
          pValue <- anovaRes[[1]][['Pr(>F)']][1]
          statistic <- anovaRes[[1]][['F value']][1]
        },
        silent = T
      )
      # Create data frame for summarizing results
      if (!is(singleRes, 'try-error')) {
        # Prevent error in creating data frame due to differing number of rows,
        # i.e., p-value and statistic are NULL and fail to be retrieved (e.g., 20
        # samples and 20 groups in ANOVA)
        if (is.null(statistic)) {
          pValue <- NA
          statistic <- NA
        }
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = pValue, Stat = statistic, Test = test,
                   stringsAsFactors = F, row.names = c())
      } else {
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = NA, Stat = NA, Test = test,
                   stringsAsFactors = F, row.names = c())
      }
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows() %>%
    dplyr::arrange(pVal)
  
  # Compute adjusted p-value
  resTab <- dplyr::mutate(resTab, pValAdj = p.adjust(pVal, method = 'BH')) %>%
    dplyr::relocate(pValAdj, .after = pVal)
  return(resTab)
}


prepRankedFeatList <- function(statRes, featAnno = NULL, to_gene = F) {
  #' This function is to conduct more reliable enrichment analysis, which tidy up
  #' significance testing results. Remove extra inferred proteins and annotated
  #' genes of same features (remove whatever after ';') and check whether there
  #' is any protein with duplicated annotated gene if 'to_gene' is set to TRUE.
  #' If there is, keep one with highest magnitude of statistic. In the end, feature
  #' list is ranked.
  #' 
  #' Parameters
  #' statRes: A table containing significance testing results, where the row names
  #' are features and the first column contains statistics with whatever column name.
  #' Note that it does not matter if excess inferred features are removed
  #' featAnno: A data frame holding original features (i.e., proteins) as the row
  #' names and corresponding feature annotations (i.e., gene) in the first column
  #' with whatever column name. Note that it does not matter if excess inferred
  #' features are removed, and the feature space must not be the same as that of 'statRes'
  #' to_gene: A logical value indicating whether the row names of 'statRes' are mapped
  #' from proteins to genes
  #' 
  #' Return
  #' rankedList: A named vector in which features are ranked
  
  #### Can provide parameter specifying increasing or decreasing for ordering provided
  #### stats (e.g., t-values, p-values, or -log10(p))
  
  #### Can combine statRes and featAnno tables into a table, so later manipulation
  #### does not have to be using index, which is very prone to error sometimes.
  #### Also, avoid working on row names, because R may modify them without telling
  #### you (e.g., row names containing '-') or some tidyverse function removes row names
  
  # Keep first inferred protein if there are many
  rownames(statRes) <- stringr::str_remove(rownames(statRes), ';.+')
  # Sanity check
  if (any(duplicated(rownames(statRes)))) {
    stop("There are duplicated features in 'statRes' when whatever after ';' is removed.")
  }
  
  # Change protein names to gene names. Those proteins with duplicated annotated
  # genes or no annotated gene will be removed
  if (to_gene) {
    # Keep first annotated gene if there are many
    rownames(featAnno) <- stringr::str_remove(rownames(featAnno), ';.+')
    featAnno[[1]] <- stringr::str_remove(featAnno[[1]], ';.+')
    # Remove proteins with duplicated annotated genes (keep one with highest magnitude
    # of statistic)
    
    #### Can order gene list based on stats and remove latter duplicated ones, so
    #### nuch code can be skipped
    
    genes <- featAnno[[1]]
    dupGenes <- unique(genes[duplicated(genes)])
    rmFeatIdxList <- numeric(0)
    for (g in dupGenes) {
      rmFeats <- rownames(featAnno)[genes %in% g]
      rmFeatIdx <- which(rownames(statRes) %in% rmFeats)
      rmFeatIdx <- rmFeatIdx[-which.max(abs(statRes[rmFeatIdx,]))]
      rmFeatIdxList <- c(rmFeatIdxList, rmFeatIdx)
    }
    # Remove proteins with no annotated gene
    rmFeats <- rownames(featAnno)[is.na(genes)]
    rmFeatIdx <- which(rownames(statRes) %in% rmFeats)
    rmFeatIdxList <- unique(c(rmFeatIdxList, rmFeatIdx))
    # Filter out unwanted proteins
    if (length(rmFeatIdxList) != 0) {
      statRes <- statRes[-rmFeatIdxList,, drop = F]
    }
    # Change protein names to gene names
    pro2gene <- featAnno[[1]]
    names(pro2gene) <- rownames(featAnno)
    rownames(statRes) <- pro2gene[rownames(statRes)]
  }
  
  # Rank feature
  rankedList <- statRes[order(statRes[[1]], decreasing = T),]
  names(rankedList) <- rownames(statRes)[order(statRes[[1]], decreasing = T)]
  
  return(rankedList)
}


doEA <- function(rankedList, categoryDB, subcategoryDB = NULL, numSigGS = 10, rmPrefixGS = NULL) {
  #' Do human gene set enrichment analysis. Visualization in this function is specific
  #' to recurrence prediction task for now, assuming positive ES score corresponds
  #' to recurrence patient group
  #' 
  #' Parameters
  #' rankedList: A named vector in which features are ranked (descending), which
  #' can be the return of the function 'prepRankedFeatList'
  #' categoryDB: A character specifying the MSigDB collection, e.g., 'H'
  #' subcategoryDB: A character specifying the subcategory in the specified MSigDB
  #' collection by the parameter 'categoryDB', e.g., 'GO:BP' in C5. Default is NULL
  #' numSigGS: A numeric value specifying the number of overrepresented terms shown
  #' in the output plot. Default is 10
  #' rmPrefixGS: A regular expression specifying the prefix(es) of terms to remove
  #' for better visualization, e.g., '^HALLMARK_' for H or '^GOCC_|^GOMF_' for C5.
  #' Default is NULL
  #' 
  #' Return
  #' A list containing the following components:
  #' gseaRes: A gseaResult object returned by 'clusterProfiler::GSEA'
  #' gseaRes4Plot: A table containing overrepresented terms, normalized enrichment
  #' scores, adjusted p-values, etc. for visualization
  #' plotTopSigGS: A ggplot object visualizing overrepresented terms
  
  # Retrieve annotated gene sets
  msigTab <- msigdbr(species = 'Homo sapiens', category = categoryDB, subcategory = subcategoryDB) %>%
    dplyr::select(gs_name, human_gene_symbol)
  # Run GSEA
  gseaRes <- clusterProfiler::GSEA(geneList = rankedList, TERM2GENE = msigTab,
                                   minGSSize = 10, maxGSSize = 500,
                                   pvalueCutoff = 0.05, pAdjustMethod = 'BH',
                                   by = 'fgsea', eps = 0)
  
  # Plot enrichment analysis results
  if (nrow(gseaRes@result) != 0) {
    gseaRes4Plot <- dplyr::select(gseaRes@result, Description, setSize, NES, pvalue, p.adjust, leading_edge) %>%
      dplyr::arrange(p.adjust) %>%
      dplyr::slice_head(n = numSigGS) %>%
      dplyr::mutate(leading_edge = stringr::str_extract(leading_edge, '\\d\\d'),
                    leading_edge = as.numeric(leading_edge),
                    Recurrence = ifelse(test = NES > 0, yes = 'Up', no = 'Down'),
                    Recurrence = factor(Recurrence, levels = c('Up', 'Down'))) %>%
      dplyr::rename(gene_ratio = leading_edge)
    #
    if (!is.null(rmPrefixGS)) {
      gseaRes4Plot <- dplyr::mutate(gseaRes4Plot, Description = stringr::str_remove_all(Description, rmPrefixGS))
    }
    midGrad <- mean(gseaRes4Plot$p.adjust)
    
    plotTopSigGS <- ggplot(gseaRes4Plot, aes(x=NES, y=reorder(Description, NES), fill=p.adjust, col=Recurrence)) +
      geom_bar(stat = 'identity', width = 0.7, size = 2) +
      scale_fill_gradient2(low='#5B1C00', high='#FFFCA5', mid = '#E86900',
                           midpoint = midGrad, name = 'Adjusted pVal') +
      scale_color_manual(values = c(Up = 'red', Down = 'darkblue')) +
      labs(x = 'Normalized Enrichment Score', y = 'Pathway') +
      theme_minimal()
  } else {
    stop("No gene sets overrepresented...")
  }
  
  return(list(gseaRes = gseaRes, gseaRes4Plot = gseaRes4Plot, plotTopSigGS = plotTopSigGS))
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








####################################################################################

subData <- function(summ_exp, split_by, factors, suffixes, smp_anno) {
  #' IF SUFFIXES IS SPECIFIED AND THEN?
  #' Parameter
  #' split_by: A character specifying the column of sample metadata
  #' factors: A vector of characters indicating factors
  #' suffixes: A vector of characters indicating suffixes to be removed, corresponding to factors
  #' smp_anno: A character or a vector of characters specifying the column of sample
  #' metadata
  
  # Subset samples
  smpAnno <- as.data.frame(colData(summ_exp))
  se1 <- summ_exp[, smpAnno[[split_by]] == factors[1]]
  se2 <- summ_exp[, smpAnno[[split_by]] == factors[2]]
  # Modify column names for matching two matrices
  colnames(se1) <- stringr::str_remove(colnames(se1), suffixes[1])
  colnames(se2) <- stringr::str_remove(colnames(se2), suffixes[2])
  # Retrieve data matrix
  Mat1 <- SummarizedExperiment::assay(se1)
  Mat2 <- SummarizedExperiment::assay(se2)
  # Subtract tumor matrix by normal matrix
  if (identical(colnames(Mat1), colnames(Mat2)) &
      identical(rownames(Mat1), rownames(Mat2))) {
    diffMat <- Mat1 - Mat2
  } else {
    stop("Align samples or features first.")
  }
  # Retrieve samples metadata for creating SE object
  smpAnno <- dplyr::select(as.data.frame(colData(se1)), all_of(smp_anno))
  diffSE <- SummarizedExperiment(assays = diffMat, colData = smpAnno)
  
  return(list(SE1 = se1, SE2 = se2, diffSE = diffSE))
}


normTest <- function(data) {
  #' Perform Shapiro-Wilk's normality test on row of input data matrix
  #' 
  #' Parameter
  #' data: A Feature x Sample data matrix
  #' 
  #' Return
  #' normRes: A data frame storing the normality of features
  
  # Check argument
  if (!is(data, 'matrix')) {
    stop("Data should be class 'matrix'. Please change it with 'as.matrix()'.")
  }
  
  normRes <- sapply(seq_len(nrow(data)), function(i) {
    featVals <- data[i,]
    normality <- try(shapiro.test(featVals), silent = T)
    if (!is(normality, 'try-error')) {
      normality$p.value > 0.05
    } else {
      NA
      }}) %>% data.frame(Feature = rownames(data), Normality = .)
  return(normRes)
}