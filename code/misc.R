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
    rowAnno <- df[,c(row_id, row_anno)]
    rowAnno <- rowAnno[!duplicated(rowAnno[[row_id]]),] %>%
      tibble::column_to_rownames(row_id) 
  } else {
    rowAnno <- NULL
    # print('There is no row annotation.')
  }
  # Prepare column annotation
  if (!is.null(col_anno)) {
    colAnno <- df[,c(col_id, col_anno)]
    colAnno <- colAnno[!duplicated(colAnno[[col_id]]),] %>%
      tibble::column_to_rownames(col_id)
  } else {
    colAnno <- NULL
    # print('There is no column annotation.')
  }
  # Assemble all information into SummarizedExperiment object
  summExp <- SummarizedExperiment(assays = matList,
                                  rowData = rowAnno,
                                  colData = colAnno)
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


PCA <- function(data) {
  #' Perform PCA (using prcomp) on input data matrix
  #' 
  #' Parameter
  #' data: A Feature x Sample data matrix
  #' 
  #' Return
  #' pcaRes: A list containing the components described in prcomp documentation
  
  # Check completeness of data
  naFeats <- complete.cases(data)
  # Run PCA if there is no missing value
  if (all(naFeats)) {
    # center = T => centering is done by subtracting corresponding column (feature)
    # means of data, which shifts features to be zero centered
    # scale. = T => scaling is done by dividing standard deviation of (centered)
    # columns, which scales features to have unit variance
    pcaRes <- prcomp(t(data), center = T, scale. = F)
  } else {
    # Remove features with missing values
    dataSub <- data[naFeats,]
    pcaRes <- prcomp(t(dataSub), center = T, scale. = F)
  }
  return(pcaRes)
}


tTestIter <- function(data, metadata, cmn_col) {
  #' Iterate t-test on features or representations grouped by varied categorical
  #' variables. For instance, t-test is systematically performed on tissue samples
  #' grouped by sample conditions (i.e, tumor and normal) across PCs
  #'
  #' Parameters
  #' data: A data frame containing numeric values of samples in the column
  #' metadata: A data frame containing sample metadata in the column
  #' cmn_col: A character specifying the common column (e.g., sample names) between
  #' data and metadata tables for combining information from them
  #' 
  #' Return
  #' resTab: A data frame collecting all t-test results
  #' 
  #' Try to simply silently return NA for errors due to limited number (one in
  #' each group) or zero variation of observations when function is called multiple
  #' time as part of automated procedure
  #' my.t.test.p.value <- function(...) {
  #'   obj <- try(t.test(...), silent=TRUE)
  #'   if (is(obj, 'try-error')) return(NA) else return(obj$p.value)
  #' }
  
  # Check arguments
  if (!(is(data, 'data.frame') & is(metadata, 'data.frame'))) {
    stop("Both tables should be class 'data.frame', 'tbl', or 'tbl_df'.")
  }
  # Check completeness of metadata tables
  naVals <- complete.cases(metadata)
  if (!all(naVals)) {
    stop("Metadata table is not complete. Please remove samples with missing values.")
  }
  
  # Combine data and metadata tables to ensure matched information
  combinedTab <- dplyr::left_join(data, metadata, by = cmn_col) %>%
    dplyr::select(-all_of(cmn_col))
  # Iterate all possible association tests
  numVar1 <- ncol(data) - 1
  numVar2 <- ncol(metadata) - 1
  # CHECK combn FUNCTION THAT CAN GENERATE ALL POSSIBLE COMBINATIONS
  resTab <- lapply(seq(numVar1), function(i) {
    lapply(seq(numVar2), function(j) {
      var1 <- combinedTab[[i]]
      var2 <- combinedTab[[numVar1+j]]
      
      # Perform t-test and extract p-value and group means
      # var.equal = F (default) => Welch's t-test
      # var.equal = T => Student's t-test
      tStatRes <- try(t.test(var1 ~ factor(var2), paired = F, var.equal = T),
                      silent = T)
      if (!is(tStatRes, 'try-error')) {
        pVal <- tStatRes$p.value
        estimates <- tStatRes$estimate
        tStats <- tStatRes$statistic
        
        # Create data frame for summarizing results 
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = pVal,
                   tStats = tStats,
                   Group1 = stringr::str_remove(names(estimates)[1],
                                                'mean in group '),
                   Group2 = stringr::str_remove(names(estimates)[2],
                                                'mean in group '),
                   Mu1 = estimates[1],
                   Mu2 = estimates[2],
                   stringsAsFactors = F)
      } else {
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = NA, tStats = NA, Group1 = NA, Group2 = NA,
                   Mu1 = NA, Mu2 = NA, stringsAsFactors = F)
        }
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows() %>%
    dplyr::arrange(pVal)
  
  # Compute adjusted p-value
  resTab <- dplyr::mutate(resTab, pValAdj = p.adjust(pVal, method = 'BH'))
  rownames(resTab) <- c()
  return(resTab)
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


SOA <- function(summ_exp, fac, fac4uni = NULL, num_feats = NULL) {
  # To-do: Expand function to e.g., MultiAssayExperiment object and association
  # test for more than two groups e.g., ANOVA
  
  #' Perform single-omics analysis on preprocessed data stored in SE container to
  #' identify potential biomarkers for certain scientific questions, e.g., predicting
  #' patient cancer recurrences. This function mainly includes univariate association
  #' test and PCA, which attempts to capture significant features and PCs.
  #' 
  #' Parameters
  #' summ_exp: A SummarizedExperiment object containing data and metadata
  #' fac: A character or a vector of characters indicating the sample metadata from
  #' a SE object, separating samples into the groups for conducting association
  #' test (t-test)
  #' fac4uni: Similar to 'fac'. If NULL (default), all single-omics analysis uses
  #' 'fac'. Specify 'fac4uni' only if different factors are used for univariate
  #' analysis of each feature.
  #' num_feats: A numeric value specifying the number of potential PCs' top features
  #' with highest absolute loadings to extract. If NULL (default), the extraction
  #' of lists of features is skipped
  #' 
  #' Return
  #' A list containing the following components:
  #' data: A data matrix from the SE object
  #' smpMetadata: A table storing sample metadata, e.g., cancer recurrence or not
  #' normRes: A table indicating the normality of features
  #' tFeatSigRes: A table showing significant association test results between
  #' features and our questions (defined by parameter 'fac')
  #' tFeatRes: A table showing all association test results between features and
  #' our questions (defined by parameter 'fac')
  #' pcaRes: A complete PCA result obtained using 'prcomp'
  #' tPCASigRes: A table showing significant association test results between PCs
  #' and our questions (defined by parameter 'fac')
  #' pcTopFeatTab: A list containing lists of top features (determined by parameter
  #' 'num_feats') of potential PCs
  #' pcTab: A PC table containing sample representations and metadata
  
  # Check arguments
  if (!is(summ_exp, 'SummarizedExperiment')) {
    stop("This function takes only 'SummarizedExperiment' object for now.")
  }
  if (!is(fac, 'character')) {
    stop("Argument for 'fac' should be class 'character'.")
  }
  if (is.null(fac4uni)) {
    fac4uni <- fac
  } else {
    if(!is(fac4uni, 'character')) {
      stop("Argument for 'fac4uni' should be class 'character'.")
    }
  }
  
  # Retrieve data matrix from SE object
  datMat <- SummarizedExperiment::assay(summ_exp) %>%
    as.matrix()
  # Retrieve sample metadata from SE object and ensure 'Sample' column exists and
  # contains sample names of data matrix for further operations (tTestIter)
  smpMetadat <- colData(summ_exp)
  if (!'Sample' %in% colnames(smpMetadat)) {
    smpMetadat <- tibble::as_tibble(smpMetadat, rownames = 'Sample')
  } else {
    smpMetadat <- tibble::as_tibble(smpMetadat, rownames = 'tmp_Sample') %>%
      dplyr::select(-Sample) %>%
      dplyr::rename(Sample = tmp_Sample)
  }
  
  # Perform PCA
  pcaRes <- PCA(datMat)
  # Conduct association test between PCs and factors
  pcTab <- pcaRes$x %>%
    tibble::as_tibble(rownames = 'Sample')
  metaTab <- dplyr::select(smpMetadat, c('Sample', all_of(fac)))
  tPCASigRes <- tTestIter(pcTab, metaTab, cmn_col = 'Sample') %>%
    dplyr::filter(pVal <= 0.05) %>%
    dplyr::select(-c(Group1, Group2, Mu1, Mu2))
  # Create complete PC table containing sample representations and metadata
  pcTab <- dplyr::left_join(pcTab, smpMetadat, by = 'Sample')
  # Extract top features with highest absolute loadings of significant PCs
  if (!is.null(num_feats)) {
    sigPCs <- tPCASigRes[['Var1']]
    pcTopFeatTab <- lapply(sigPCs, function(pc) {
      featLoad <- pcaRes$rotation[, pc]
      topFeatLoad <- sort(abs(featLoad), decreasing = T)[1:num_feats]
      topFeatIDs <- names(topFeatLoad)
      data.frame(Feature = topFeatIDs,
                 Loading = topFeatLoad,
                 row.names = c())
    })
    names(pcTopFeatTab) <- sigPCs
  } else {
    pcTopFeatTab <- NULL
  }
  
  # Check normality of features for parametric association test
  normRes <- normTest(datMat)
  # Conduct association test between features and factors
  featTab <- as_tibble(t(datMat), rownames = 'Sample')
  metaTab <- dplyr::select(smpMetadat, c('Sample', all_of(fac4uni)))
  tFeatRes <- tTestIter(featTab, metaTab, cmn_col = 'Sample')
  # Extract features that can significantly separate groups of samples
  tFeatSigRes <- dplyr::filter(tFeatRes, pVal <= 0.05) %>%
    dplyr::select(-c(Group1, Group2, Mu1, Mu2))
  
  return(list(data = datMat, smpMetadata = smpMetadat, normRes = normRes,
              tFeatSigRes = tFeatSigRes, tFeatRes = tFeatRes, pcaRes = pcaRes,
              tPCASigRes = tPCASigRes, pcTopFeatTab = pcTopFeatTab, pcTab = pcTab
  ))
}


effSize <- function(mofaObj, gp = NULL, fac, abs_es = T) {
  #' Compute effect size of factor learned from MOFA model, which can potentially
  #' predict or separate patients that (will) suffer cancer recurrences or not
  #' 
  #' Parameters
  #' mofaObj: A trained MOFA model containing learned factors and metadata
  #' gp: A character specifying the group name of interest to retrieve the factor
  #' table from the trained MOFA model. If there is only one group defined, this
  #' parameter can be ignored
  #' fac: An index or a character (e.g., 1 or 'Factor1') specifying the factor of
  #' interest to retrieve the factor vector from the factor table
  #' abs_es: A logical value indicating
  #' 
  #' Return
  #' es: A numerical value indicating the effect size
  
  if (length(MOFA2::get_factors(mofaObj)) > 1) {
    if (is.null(gp)) {
      print('The MOFA model has more than 1 group. 
            Please specify a group of interest.')
    }
    # Retrieve factor vector
    facTab <- tibble::as_tibble(MOFA2::get_factors(mofaObj)[[gp]][, fac],
                                rownames = 'sample')
    # Retrieve sample recurrence annotations
    metaTab <- tibble::as_tibble(samples_metadata(mofaObj)) %>%
      dplyr::filter(grepl(gp, group)) %>%
      dplyr::select(sample, Recurrence)
  } else if (length(MOFA2::get_factors(mofaObj)) == 1) {
    # Retrieve factor vector
    facTab <- tibble::as_tibble(MOFA2::get_factors(mofaObj)[[1]][, fac],
                                rownames = 'sample')
    # Retrieve sample recurrence annotations
    metaTab <- tibble::as_tibble(samples_metadata(mofaObj)) %>%
      dplyr::select(sample, Recurrence)
  }
  # Compute mean and standard deviation
  statTab <- dplyr::left_join(facTab, metaTab, by = 'sample') %>%
    dplyr::group_by(Recurrence) %>%
    dplyr::summarise(Mean = mean(value), Std = sd(value), Median = median(value))
  if (abs_es) {
    # Compute absolute effect size
    es <- abs(statTab$Median[1] - statTab$Median[2])
  } else {
    # Compute Cohen's D
    es <- abs(statTab$Mean[1] - statTab$Mean[2]) /
      sqrt((statTab$Std[1]^2 + statTab$Std[2]^2) / 2)
  }
  return(es)
}


modiFeatNames <- function(mofaObj, view, feat_anno = NULL, view2 = NULL) {
  #' Modify feature names from MOFA model by removing suffixes for cleaner plots.
  #' If feat_anno is not null, features are converted to corresponding annotations,
  #' e.g., proteins to genes
  #' 
  #' Parameters
  #' mofaObj: A trained MOFA model
  #' view: A character or a vector of characters indicating the view in which
  #' feature names are retrieved
  #' feat_anno: A data frame with old features as the row names containing
  #' corresponding gene annotations
  #' view2: A character or a vector of characters indicating the view in which
  #' feature names are converted to corresponding annotations. This parameter has
  #' to be specified if feat_anno is specified
  #' 
  #' Return
  #' mofaObj: A MOFA model with modified feature names
  
  if (length(view) > 1) {
    for (v in view) {
      # Retrieve weight matrix
      W <- MOFA2::get_weights(mofaObj)[[v]]
      # Remove suffixes from feature names
      rownames(W) <- stringr::str_remove(rownames(W), '_.*')
      # Convert feature names to corresponding annotations
      if (!is.null(feat_anno)) {
        if (v %in% view2) {
          geneAnno <- feat_anno[[1]]
          names(geneAnno) <- rownames(feat_anno)
          rownames(W) <- geneAnno[rownames(W)]
          }
        }
      # Replace old weight matrix with new one with modified feature names
      mofaObj@expectations$W[[v]] <- W
      }
    } else {
      W <- MOFA2::get_weights(mofaObj)[[view]]
      rownames(W) <- stringr::str_remove(rownames(W), '_.*')
      if (!is.null(feat_anno)) {
        if (view == view2) {
          geneAnno <- feat_anno[[1]]
          names(geneAnno) <- rownames(feat_anno)
          rownames(W) <- geneAnno[rownames(W)]
        }
      }
      mofaObj@expectations$W[[view]] <- W
    }
  return(mofaObj)
}


data_heatmap <- function(mofa_obj, factor, view, group = 1, features,
                         feat_anno = NULL, smp_anno = NULL, smp_cluster = NULL,
                         rm_smp_suffix = NULL, ...) {
  #' Plot a heatmap to display those features with high weights (molecular
  #' signatures) in the input data
  #' 
  #' Parameters
  #' mofa_obj: An trained MOFA object
  #' factor: An integer specifying the factor of interest
  #' view: A character specifying the view for retrieving the corresponding data
  #' group: A character specifying the group for retrieving the corresponding data
  #' features: An integer specifying the number of features with highest weights
  #' (absolute values) to extract
  #' feat_anno: A data frame with proteins as the row names containing corresponding
  #' gene annotations
  #' smp_anno: A data frame with samples as the row names containing sample
  #' annotations, e.g, patient conditions
  #' smp_cluster: A character specifying the sample annotations used to cluster
  #' samples in heatmap for easier heatmap interpretation
  #' rm_smp_suffix: A character specifying the sample suffix that will be removed
  #' 
  #' Return
  #' P: A heatmap
  
  # Retrieve weight matrix from MOFA object
  W <- MOFA2::get_weights(mofa_obj)[[view]]
  # Remove suffixes from feature names for matching gene annotations
  rownames(W) <- stringr::str_remove(rownames(W), paste0('_', view))
  # Map feature names from proteins to corresponding genes
  if (!is.null(feat_anno)) {
    geneAnno <- feat_anno[[1]]
    names(geneAnno) <- rownames(feat_anno)
    rownames(W) <- geneAnno[rownames(W)]
  }
  # Extract features with highest weights
  weightVec <- W[, paste0('Factor', factor)]
  featOrder <- order(abs(weightVec), decreasing = T)[seq_len(features)]
  topFeat <- names(weightVec[featOrder])
  # Remove NA features, or all NA features will be extracted later
  if (any(is.na(topFeat))) {
    featOrder <- featOrder[-which(is.na(topFeat))]
  }
  
  # Retrieve input data from MOFA object
  D <- MOFA2::get_data(mofa_obj)[[view]][[group]]
  # Remove suffixes from feature names for cleaner plots
  rownames(D) <- stringr::str_remove(rownames(D), paste0('_', view))
  # Map features names from proteins to corresponding genes
  if (!is.null(feat_anno)) {
    rownames(D) <- geneAnno[rownames(D)]
  }
  # Prepare data for plotting heatmap
  if (!is.null(smp_anno)) {
    # Tidy up sample annotations
    if (!identical(sort(rownames(smp_anno)), sort(colnames(D)))) {
      smp_anno <- tibble::rownames_to_column(smp_anno, 'sample') %>%
        dplyr::filter(sample %in% colnames(D)) %>%
        tibble::column_to_rownames('sample')
    }
    # Order samples by specific annotation for easier heatmap interpretation
    if (!is.null(smp_cluster)) {
      smpAnno <- smp_anno[[smp_cluster]]
      smpOrder <- c()
      for (clust in unique(smpAnno)) {
        smpOrder <- append(smpOrder, rownames(smp_anno)[which(smpAnno == clust)])
      }
      # Subset input data
      DSub <- D[featOrder, smpOrder] #rownames(D) %in% topFeat
    } else {
      # Subset input data
      DSub <- D[featOrder,]
    }
  } else {
    # Subset input data
    DSub <- D[featOrder,]
  }
  # Remove suffix of samples
  if(!is.null(rm_smp_suffix)) {
    colnames(DSub) <- stringr::str_remove(colnames(DSub), rm_smp_suffix)
  }
  
  # Make heatmap
  P <- pheatmap::pheatmap(DSub, annotation_col = smp_anno, ...)
  return(P)
}


unique_features <- function(weight_mat, factor, feat_anno, g_col) {
  # ADD SANITY CHECK FOR FEATURES FROM WEIGHT MATRIX AND ANNOTATION TABLE
  
  #' Remove extra inferred proteins from same features and keep single feature
  #' with highest weight if duplicated
  #' 
  #' Parameters
  #' weight_mat: A weight matrix obtained from the trained MOFA model
  #' factor: An integer specifying the factor of interest
  #' feat_anno: A feature metadata containing gene annotations of proteins
  #' g_col: A character specifying the column of gene annotations
  #' 
  #' Return
  #' W: A SummarizedExperiment object containing the modified weight matrix and
  #' feature metadata
  
  # Specify factor of interest
  fac <- paste0('Factor', factor)
  
  # Remove extra inferred proteins from same features (Pick the first one)
  modiFeat <- stringr::str_remove(rownames(feat_anno), ';.*')
  names(modiFeat) <- rownames(feat_anno)
  # Pinpoint duplicated protein features and keep the one with highest weight
  dupFeatIdx <- which(duplicated(modiFeat))
  
  # CHECK THIS BLOCK
  if (length(dupFeatIdx) != 0) {
    dropFeatList <- c()
    for (i in dupFeatIdx) {
      feat <- modiFeat[i]
      featIdx <- which(modiFeat == feat)
      dropFeatList <- append(dropFeatList,
                             featIdx[-which.max(weight_mat[featIdx, fac])])
    }
    # Filter out unwanted protein features and corresponding gene annotations
    # Features
    modiFeat <- modiFeat[-dropFeatList]
    # Weight matrix
    weightMatModi <- weight_mat[-dropFeatList,]
    rownames(weightMatModi) <- modiFeat
    # Feature annotations
    if (ncol(feat_anno) == 1) {
      modiGeneAnno <- stringr::str_remove(feat_anno[-dropFeatList,], ';.*')
      rowAnnoModi <- data.frame(modiGeneAnno, row.names = modiFeat)
      colnames(rowAnnoModi) <- g_col
    } else {
      rowAnnoModi <- feat_anno[-dropFeatList,]
      modiGeneAnno <- rowAnnoModi[, g_col]
      modiGeneAnno <- stringr::str_remove(modiGeneAnno, ';.*')
      rowAnnoModi[, g_col] <- modiGeneAnno
    }
  } else {
    # Weight matrix
    weightMatModi <- weight_mat
    rownames(weightMatModi) <- modiFeat[rownames(weightMatModi)]
    # Feature annotations
    if (ncol(feat_anno) == 1) {
      modiGeneAnno <- stringr::str_remove(feat_anno[, g_col], ';.*')
      rowAnnoModi <- data.frame(modiGeneAnno, row.names = modiFeat)
      colnames(rowAnnoModi) <- g_col
    } else {
      rowAnnoModi <- feat_anno
      modiGeneAnno <- rowAnnoModi[, g_col]
      modiGeneAnno <- stringr::str_remove(modiGeneAnno, ';.*')
      rowAnnoModi[, g_col] <- modiGeneAnno
    }
  }
  
  # Create SummarizedExperiment object to include gene annotations for proteins
  W <- SummarizedExperiment(weightMatModi, rowData = rowAnnoModi)
  return(W)
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
  
  return(list(reducedData = reduData,
              corrFeats = corrFeatList))
}


corrIter <- function(tabX, tabY, cmn_col) {
  #' 
  
  # Combine Table X and Table Y to ensure matched information
  combinedTab <- dplyr::left_join(tabX, tabY, by = cmn_col) %>%
    dplyr::select(-all_of(cmn_col))
  
  # Iterate all possible association tests
  numVar1 <- ncol(tabX) - 1
  numVar2 <- ncol(tabY) - 1
  resTab <- lapply(seq(numVar1), function(i) {
    lapply(seq(numVar2), function(j) {
      var1 <- combinedTab[[i]]
      var2 <- combinedTab[[numVar1+j]]
      
      # Count number of non-NA value pairs, which means none of pair values from
      # pair features is NA
      nonNACount <- 0
      for (p in seq_len(length(var1))) {
        if (!is.na(var1[p]) & !is.na(var2[p])) {
          nonNACount <- nonNACount + 1
        }
      }
      
      # Compute correlations and extract p-value
      if (nonNACount >= 3) {
        # use = 'complete.obs' (casewise) => remove all cases with missing data
        # and calculate correlation on remaining data
        # use = 'pairwise.complete.obs' => remove cases with missing data for each
        # correlation calculation
        corrRes <- stats::cor.test(var1, var2, method = 'pearson',
                                   use = "pairwise.complete.obs")
        corr <- corrRes$estimate
        pVal <- corrRes$p.value
        
        # Create data frame for summarizing results 
        data.frame(Var1 = colnames(combinedTab)[i],
                   Var2 = colnames(combinedTab)[numVar1+j],
                   pVal = pVal,
                   Corr = corr,
                   stringsAsFactors = F)
      }
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  # Compute adjusted p-value
  resTab <- dplyr::mutate(resTab, pValAdj = p.adjust(pVal, method = 'BH'))
  rownames(resTab) <- c()
  return(resTab)
}


# To-do: Extend the function to different classes of variables
associationTest <- function(tabX, tabY, cmn_col) {
  #' Conduct association test between columns (variables) from two tables, e.g.,
  #' Table X contains dependent variables and Table Y contains independent variables.
  #' One-way ANOVA can be used to assess whether an independent variable has a
  #' real impact on dependent variables.
  #' 
  #' Parameters:
  #' 
  
  # Combine Table X and Table Y to ensure matched information
  combinedTab <- dplyr::left_join(tabX, tabY, by = cmn_col) %>%
    dplyr::select(-all_of(cmn_col))
  
  # Iterate all possible association tests
  numVar1 <- ncol(tabX) - 1
  numVar2 <- ncol(tabY) - 1
  resTab <- lapply(seq(numVar1), function(i) {
    lapply(seq(numVar2), function(j) {
      var1 <- combinedTab[[i]]
      var2 <- combinedTab[[numVar1+j]]
      
      # Perform ANOVA and extract p-value
      pVal <- summary(aov(var1 ~ factor(var2)))[[1]][['Pr(>F)']][1]
      
      # Create data frame for summarizing results 
      data.frame(Var1 = colnames(combinedTab)[i],
                 Var2 = colnames(combinedTab)[numVar1+j],
                 pVal = pVal,
                 stringsAsFactors = F)
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows() %>%
    dplyr::arrange(pVal)
  
  return(resTab)
}