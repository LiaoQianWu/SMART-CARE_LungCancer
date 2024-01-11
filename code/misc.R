library(pcaMethods)
library(limma)
library(proDA)
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
    pcaRes <- pcaMethods::pca(t(data), method = 'ppca', nPcs = num_PCs, seed = 123,
                              center = T, scale = 'none')
  } else if (pca_method == 'bpca') {
    # Run Bayesian PCA
    data[is.nan(data)] <- NA
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
                  unwantVar = NULL, ...) {
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
  #' ...: Further arguments to be passed to 'testAssociation', e.g., correlation
  #' method 'cor_method'
  #' 
  #' Return
  #' A list containing the following components:
  #' data: A data matrix from the SE object
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
  
  # Perform PCA
  pcaRes <- doPCA(datMat, pca_method = pca_method, num_PCs = num_PCs)
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
  
  return(list(data = datMat, smpMetadata = smpMetadat, pcaRes = pcaRes,
              pcSigAssoRes = pcSigAssoRes, pcTopFeatTab = pcTopFeatTab, pcTab = pcTab,
              featSigAssoRes = featSigAssoRes, featAssoRes = featAssoRes))
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
                  trainSet_ratio = 0.8, ntree = 10000, plot_ROC = F) {
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
  #' iter: A numeric value specifying the number of time to run random forest
  #' split_method: A character specifying the splitting method to use, which should
  #' be one of 'random split' (default) or 'bootstrap'
  #' trainSet_ratio: A numeric value specifying the ratio of the size of the training
  #' set to that of the total dataset
  #' ntree: A numeric value specifying the number of trees to grow, which should not
  #' be set to too low to ensure that every input row gets predicted at least a few times
  #' plot_ROC: A logical variable indicating whether ROC curve is plotted
  #' 
  #' Return
  #' A list containing the following components:
  #' mtry: A vector of tuned best mtry used to train RF models
  #' auc: A vector of computed AUC from trained RF models
  #' MDA, MDG: Matrices containing computed feature importance. Rows and columns
  #' present features and independent trained RF models. MDA and MDG stand for Mean
  #' Decrease Accuracy and Gini
  #' params: A list containing arguments of parameters 'split_method', 'trainSet_ratio',
  #' and 'ntree'
  
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
                                          # xtest = x_test, ytest = y_test,
                                          importance = T, proximity = T, keep.forest = T)
    } else {
      mtryList <- c(mtryList, NA)
      rfRes <- randomForest::randomForest(x = x_train, y = y_train, ntree = ntree,
                                          # xtest = x_test, ytest = y_test,
                                          importance = T, proximity = T, keep.forest = T)
    }
      
      # Evaluate model by AUC of ROC curve
      y_pred <- predict(rfRes, newdata = x_test, type = 'prob')
      if (plot_ROC) {
        par(pty = 's')
        rocRes <- pROC::roc(response = y_test, predictor = y_pred[, '1'], plot = T,
                            legacy.axes = T, print.auc = T, print.auc.x = 0.4,
                            xlab = 'False positive rate', ylab = 'True positive rate',
                            main = 'ROC Curve for Random Forest',
                            cex.lab = 1.2, cex.main = 1.1, col = '#377eb8', lwd = 4,
                            direction = '<')
        par(pty = 'm')
      } else {
        rocRes <- pROC::roc(response = y_test, predictor = y_pred[, '1'], plot = F,
                            direction = '<')
      }
      # Save AUC
      aucList <- c(aucList, rocRes$auc)
      
      # Retrieve and save feature importance values from RF object
      matMDA[, i] <- rfRes$importance[, 'MeanDecreaseAccuracy']
      matMDG[, i] <- rfRes$importance[, 'MeanDecreaseGini']
  }
  
  return(list(mtry = mtryList, auc = aucList, MDA = matMDA, MDG = matMDG, params = paramList))
}


runLogisR <- function(x, y, targetClass, regularized_method = 'lasso', cvFold = 10,
                      cvMeasure = 'auc', used_lambda = 'lambda.1se', iter = 20,
                      trainSet_ratio = 0.8, split_method = 'random split', plot_ROC = F) {
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
  #' 
  #' Return
  #' A list containing the following components:
  #' coefficient: A matrix containing optimized coefficients. Rows and columns present
  #' features and independent trained models
  #' usedLambda: A vector of used lambda for selected trained models
  #' auc: A vector of computed AUC from trained models
  #' params: A list containing arguments of parameters 'regularized_method', 'cvFold',
  #' 'used_lambda', 'split_method', and 'trainSet_ratio'
  
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
  paramList <- list(regularized_method = regularized_method, cvFold = cvFold, cvMeasure = cvMeasure,
                    used_lambda = used_lambda, split_method = split_method, trainSet_ratio = trainSet_ratio)
  
  ####
  numNonZeroList <- c()
  
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
                               nfolds = cvFold, alpha = alpha, standardize = F, intercept = T)
    # Save selected lambda of trained model
    lambdaList <- c(lambdaList, round(lrRes[[used_lambda]], 4))
    # Retrieve and save coefficients of fitted model
    matCoeffi[, i] <- as.matrix(coef(lrRes, s = lrRes[[used_lambda]]))[-1, 1]
    
    ####
    coeffi <- as.matrix(coef(lrRes, s = lrRes[[used_lambda]]))[-1, 1]
    numNonZeroList <- c(numNonZeroList, sum(coeffi != 0))
    
    # Evaluate model by AUC of ROC curve
    y_pred <- predict(lrRes, s = lrRes[[used_lambda]], newx = x_test, type = 'response')
    if (plot_ROC) {
      par(pty = 's')
      rocRes <- pROC::roc(response = y_test, predictor = y_pred[, 1], plot = T,
                          legacy.axes = T, print.auc = T, print.auc.x = 0.4,
                          xlab = 'False positive rate', ylab = 'True positive rate',
                          main = 'ROC Curve for Random Forest',
                          cex.lab = 1.2, cex.main = 1.1, col = '#377eb8', lwd = 4,
                          direction = '<')
      par(pty = 'm')
    } else {
      rocRes <- pROC::roc(response = y_test, predictor = y_pred[, 1], plot = F,
                          direction = '<')
    }
    # Save AUC
    aucList <- c(aucList, rocRes$auc)
  }
  
  return(list(coefficient = matCoeffi, usedLambda = lambdaList, auc = aucList,
              nNonZero = numNonZeroList, params = paramList))
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