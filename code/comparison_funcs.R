doFeatAssoTest <- function(data, meta_var, subset = NULL, prot_gene_tbl = NULL,
                           method = 'proDA', unwantVar = NULL) {
  #' Conduct univariate association test (t-test for now) between each feature and
  #' sample metadata variable of interest
  #' 
  #' Parameters
  #' data: A SummarizedExperiment object that contains the data matrix and metadata
  #' meta_var: A character specifying the sample metadata variable that has two
  #' levels for conducting t-test
  #' subset: A vector of length 2 containing characters specifying the sample metadata
  #' variable and one of its level to split the dataset. Default is NULL, meaning
  #' the whole dataset is used
  #' prot_gene_tbl: A table containing proteins and their encoding genes for providing
  #' gene names in a resulting table. Default is NULL. The column for proteins should
  #' be named 'Proteins' and for genes should be named 'Genes'. Furthermore, only
  #' the first protein and gene should be kept if there are more than one
  #' method: A character specifying the method used to conduct t-test, which should
  #' be one of 'proDA' (default) and 'limma'
  #' unwantVar: A character or a vector of characters specifying the metadata variables
  #' that will be accounted for by linear models (regressing out). Default is NULL
  #' 
  #' Return
  #' assoTestRes: A table of feature significance test and statistical results
  
  # Check arguments
  # meta_var should be length of 1 and in sample metadata
  # subset should be length of 2 and in sample metadata
  # method should be one of 'proDA' and limma''
  
  # Subset samples if specified
  if (!is.null(subset)) {
    subsetIdx <- which(colData(data)[[subset[1]]] == subset[2])
    data <- data[, subsetIdx]
  }
  datMat <- SummarizedExperiment::assay(data)
  smpAnno <- tibble::as_tibble(colData(data), rownames = 'Sample')
  
  # Perform analysis
  # Prepare model design matrix
  if (is.null(unwantVar)) {
    design <- model.matrix(~ smpAnno[[meta_var]])
  } else {
    formu <- paste0('~', paste0(c(meta_var, unwantVar), collapse = '+')) %>%
      as.formula()
    design <- model.matrix(formu, data = smpAnno)
  }
  if (method == 'proDA') {
    # Fit linear probabilistic dropout model to normalized data
    fit <- proDA(datMat, design = design)
    diffTerm <- proDA::result_names(fit)[2]
    assoTestRes <- proDA::test_diff(fit, contrast = diffTerm, sort_by = 'pval') %>%
      dplyr::select(name, pval, adj_pval, diff) %>%
      dplyr::rename(Feature = name, pVal = pval, pValAdj = adj_pval, logFC = diff)
  } else if (method == 'limma') {
    fit <- limma::lmFit(datMat, design = design)
    fit <- limma::eBayes(fit)
    assoTestRes <- limma::topTable(fit, coef = 2, number = Inf) %>%
      tibble::as_tibble(rownames = 'Feature') %>%
      dplyr::select(Feature, P.Value, adj.P.Val, logFC) %>%
      dplyr::rename(pVal = P.Value, pValAdj = adj.P.Val)
  }
  assoTestRes <- dplyr::mutate(assoTestRes,
                               pVal = as.numeric(formatC(pVal, format = 'g', digits = 3)),
                               pValAdj = as.numeric(formatC(pValAdj, format = 'g', digits = 3)))
  # Include feature 
  if (!is.null(prot_gene_tbl)) {
    assoTestRes <- dplyr::mutate(assoTestRes, Gene = plyr::mapvalues(Feature,
                                                                     from = prot_gene_tbl$Proteins,
                                                                     to = prot_gene_tbl$Genes,
                                                                     warn_missing = F)) %>%
      dplyr::relocate(Gene, .after = Feature)
  }
  
  return(assoTestRes)
}


summarizeFeatAssoRes <- function(data1, data2, name1, name2, meta_var, subset = NULL,
                                 prot_gene_tbl = NULL, method = 'proDA',
                                 unwantVar1 = NULL, unwantVar2 = NULL,
                                 alpha = 0.05, num_topSigFeats = 6, plot_pVal = T,
                                 plot_topSigFeats = T, plot_logFC = T) {
  #' Conduct univariate analysis (t-test for now) between features and sample metadata
  #' variable of interest, and combine and tidy up significance test and statistical
  #' results (p-values, top significant features, and log fold changes) into tables for making plots
  #' 
  #' Parameters
  #' data1, data2: SummarizedExperiment objects that contain the data matrices and metadata
  #' name1, name2: A character indicating the names of datasets (data1 and data2)
  #' in the plot legend
  #' meta_var: A character specifying the sample metadata variable that has two
  #' levels for grouping samples to conduct t-test
  #' subset: A vector of length 2 containing characters specifying the sample metadata
  #' variable and one of its level to split the dataset. Default is NULL, meaning
  #' the whole dataset is used
  #' prot_gene_tbl: A table containing proteins and their encoding genes for providing
  #' and labeling gene names in plots. Default is NULL. The column for proteins
  #' should be named 'Proteins' and for genes should be named 'Genes'. Furthermore,
  #' only the first protein and gene should be kept if there are more than one
  #' method: A character specifying the method used to conduct t-test, which should
  #' be one of 'proDA' (default) and 'limma'
  #' unwantVar1, unwantVar2: A character or a vector of characters specifying the
  #' metadata variables (data1 and data2) that will be accounted for by linear models
  #' (regressing out). Default is NULL
  #' alpha: A numeric value specifying the significance level for summarizing p-value
  #' information. Default is 0.05
  #' num_topSigFeats: A numeric value specifying the number of top significant features
  #' to visualize. Default is 6
  #' plot_pVal: A logical variable indicating whether to generate the summarized
  #' p-value table for plotting. Default is TRUE
  #' plot_topSigFeats: A logical variable indicating whether to generate the summarized
  #' table of top significant features based on 'data1' for plotting. Default is TRUE
  #' plot_logFC: A logical variable indicating whether to generate the summarized
  #' log(FC) table for plotting. Default is TRUE
  #' 
  #' Return
  #' A list containing the following components:
  #' tbl_featAsso1: A table containing association test results of 'data1'
  #' tbl_featAsso2: A table containing association test results of 'data2'
  #' tbl_pVal: A table containing summarized significance test results, labeling
  #' information, etc. for making plots
  #' tbl_topSigFeats: A table containing abundances of top significant features,
  #' grouping information, etc. for making plots
  #' tbl_logFC: A table containing summarized log(FC) results, coloring information,
  #' etc. for making plots
  
  # Extract common features between two datasets
  # Select first protein as representative if there are more than one
  rownames(data1) <- stringr::str_remove(rownames(data1), ';.*')
  rownames(data2) <- stringr::str_remove(rownames(data2), ';.*')
  commonFeats <- intersect(rownames(data1), rownames(data2))
  
  # Compute associations between features and sample metadata variable of interest
  # and keep only results of common features
  assoTestRes1 <- doFeatAssoTest(data1, meta_var = meta_var, subset = subset,
                                 prot_gene_tbl = prot_gene_tbl, method = method,
                                 unwantVar = unwantVar1)
  assoTestRes_dat1 <- dplyr::filter(assoTestRes1, Feature %in% commonFeats)
  assoTestRes2 <- doFeatAssoTest(data2, meta_var = meta_var, subset = subset,
                                 prot_gene_tbl = prot_gene_tbl, method = method,
                                 unwantVar = unwantVar2)
  assoTestRes_dat2 <- dplyr::filter(assoTestRes2, Feature %in% commonFeats)
  
  # Summarize resulting p-Values
  if (plot_pVal) {
    # Prepare feature significance information
    tmp_assoTestRes_dat1 <- dplyr::select(assoTestRes_dat1, Feature, pVal) %>%
      dplyr::mutate(Sig_dat1 = dplyr::case_when(pVal <= alpha ~ 'Yes',
                                                pVal > alpha ~ 'No')) %>%
      dplyr::rename(pVal_dat1 = pVal)
    tmp_assoTestRes_dat2 <- dplyr::select(assoTestRes_dat2, Feature, pVal) %>%
      dplyr::mutate(Sig_dat2 = dplyr::case_when(pVal <= alpha ~ 'Yes',
                                                pVal > alpha ~ 'No')) %>%
      dplyr::rename(pVal_dat2 = pVal)
    # Combine association test results and all needed information into a table
    combinedAssoTestRes_pVal <- dplyr::left_join(tmp_assoTestRes_dat1, tmp_assoTestRes_dat2,
                                                 by = 'Feature') %>%
      # Include coloring information
      dplyr::mutate(Significance = dplyr::case_when(
        Sig_dat1 == 'Yes' & Sig_dat2 == 'Yes' ~ 'Sig. in Both',
        Sig_dat1 == 'Yes' & Sig_dat2 == 'No' ~ paste('Sig. in', name1),
        Sig_dat1 == 'No' & Sig_dat2 == 'Yes' ~ paste('Sig. in', name2),
        Sig_dat1 == 'No' & Sig_dat2 == 'No' ~ 'Not Sig.'
      ),
      Significance = factor(Significance, levels = c('Sig. in Both',
                                                     paste('Sig. in', name1),
                                                     paste('Sig. in', name2),
                                                     'Not Sig.')))
    # Include labeling information
    if (is.null(prot_gene_tbl)) {
      combinedAssoTestRes_pVal <- dplyr::mutate(combinedAssoTestRes_pVal,
                                                Label = dplyr::case_when(
                                                  Significance == 'Sig. in Both' ~ Feature
                                                ))
    } else {
      combinedAssoTestRes_pVal <- dplyr::mutate(combinedAssoTestRes_pVal,
                                                Gene = plyr::mapvalues(
                                                  Feature,
                                                  from = prot_gene_tbl$Proteins,
                                                  to = prot_gene_tbl$Genes,
                                                  warn_missing = F
                                                ),
                                                Label = dplyr::case_when(
                                                  Significance == 'Sig. in Both' ~ Gene
                                                ))
    }
    # Reorder points in graph
    significance <- combinedAssoTestRes_pVal$Significance
    combinedAssoTestRes_pVal <- combinedAssoTestRes_pVal[c(which(significance == 'Not Sig.'),
                                                           which(significance == 'Sig. in Both'),
                                                           which(significance == paste('Sig. in', name1)),
                                                           which(significance == paste('Sig. in', name2))),]
  } else {
    combinedAssoTestRes_pVal <- NULL
  }
  
  # Summarize top significant features
  # Check if there is any common significant feature
  sigFeats1 <- dplyr::filter(assoTestRes_dat1, pVal <= alpha) %>%
    dplyr::arrange() %>%
    dplyr::pull(Feature)
  sigFeats2 <- dplyr::filter(assoTestRes_dat2, pVal <= alpha) %>%
    dplyr::arrange() %>%
    dplyr::pull(Feature)
  sigFeats <- intersect(sigFeats1, sigFeats2)
  if (plot_topSigFeats & (length(sigFeats) != 0)) {
    # Retrieve top common significant features based on 'data1'
    topSigFeats <- intersect(sigFeats1, sigFeats2)
    if (length(topSigFeats) > num_topSigFeats) {
      topSigFeats <- topSigFeats[1:num_topSigFeats]
    }
    # Extract abundances of top common significant features and combine all needed
    # information into a table
    # Subset samples if specified
    if (!is.null(subset)) {
      data1 <- data1[, which(colData(data1)[[subset[1]]] == subset[2])]
      data2 <- data2[, which(colData(data2)[[subset[1]]] == subset[2])]
    }
    # Prepare sample grouping information
    smpGpInfo1 <- tibble::as_tibble(colData(data1), rownames = 'Sample') %>%
      dplyr::select(Sample, all_of(meta_var))
    smpGpInfo2 <- tibble::as_tibble(colData(data2), rownames = 'Sample') %>%
      dplyr::select(Sample, all_of(meta_var))
    # Prepare abundance information for top common significant features
    topSigFeatTab1 <- t(SummarizedExperiment::assay(data1)) %>%
      tibble::as_tibble(rownames = 'Sample') %>%
      dplyr::select(c(Sample, topSigFeats)) %>%
      tidyr::pivot_longer(cols = -'Sample', names_to = 'Feature', values_to = 'Abundance') %>%
      dplyr::left_join(smpGpInfo1, by = 'Sample') %>%
      dplyr::mutate(DataName = name1)
    topSigFeatTab2 <- t(SummarizedExperiment::assay(data2)) %>%
      tibble::as_tibble(rownames = 'Sample') %>%
      dplyr::select(c(Sample, topSigFeats)) %>%
      tidyr::pivot_longer(cols = -'Sample', names_to = 'Feature', values_to = 'Abundance') %>%
      dplyr::left_join(smpGpInfo2, by = 'Sample') %>%
      dplyr::mutate(DataName = name2)
    # Combine all information
    combinedTopSigFeatTab <- dplyr::bind_rows(topSigFeatTab1, topSigFeatTab2) %>%
      # Include new grouping information
      dplyr::mutate(newGroup = paste0(.data[[meta_var]], '_', DataName))
    # Create new feature names if feature annotations are provided
    if (!is.null(prot_gene_tbl)) {
      combinedTopSigFeatTab <- dplyr::mutate(combinedTopSigFeatTab,
                                             Gene = plyr::mapvalues(Feature,
                                                                    from = prot_gene_tbl$Proteins,
                                                                    to = prot_gene_tbl$Genes,
                                                                    warn_missing = F),
                                             newFeature = paste0(Feature, ' (', Gene, ')'),
                                             newFeature = factor(newFeature,
                                                                 levels = unique(newFeature)))
    }
  } else {
    combinedTopSigFeatTab <- NULL
  }
  
  # Summarize resulting log(FC)
  if (plot_logFC) {
    # Prepare feature log(FC) information
    tmp_assoTestRes_dat1 <- dplyr::select(assoTestRes_dat1, Feature, logFC) %>%
      dplyr::mutate(Dir_dat1 = dplyr::case_when(logFC > 0 ~ 'Up',
                                                logFC < 0 ~ 'Down')) %>%
      dplyr::rename(logFC_dat1 = logFC)
    tmp_assoTestRes_dat2 <- dplyr::select(assoTestRes_dat2, Feature, logFC) %>%
      dplyr::mutate(Dir_dat2 = dplyr::case_when(logFC > 0 ~ 'Up',
                                                logFC < 0 ~ 'Down')) %>%
      dplyr::rename(logFC_dat2 = logFC)
    # Combine association test results and all needed information into a table
    combinedAssoTestRes_logFC <- dplyr::left_join(tmp_assoTestRes_dat1, tmp_assoTestRes_dat2,
                                                  by = 'Feature') %>%
      # Include coloring information
      dplyr::mutate(Direction = dplyr::case_when(
        Dir_dat1 == 'Up' & Dir_dat2 == 'Up' ~ 'Up in Both',
        Dir_dat1 == 'Down' & Dir_dat2 == 'Down' ~ 'Down in Both',
        Dir_dat1 == 'Up' & Dir_dat2 == 'Down' ~ 'Inconsistent',
        Dir_dat1 == 'Down' & Dir_dat2 == 'Up' ~ 'Inconsistent'
      ),
      Direction = factor(Direction, levels = c('Up in Both', 'Down in Both', 'Inconsistent')))
    # Include gene annotations
    if (!is.null(prot_gene_tbl)) {
      combinedAssoTestRes_logFC <- dplyr::mutate(combinedAssoTestRes_logFC,
                                                 Gene = plyr::mapvalues(Feature,
                                                                        from = prot_gene_tbl$Proteins,
                                                                        to = prot_gene_tbl$Genes,
                                                                        warn_missing = F))
    }
    # Reorder points in graph
    direction <- combinedAssoTestRes_logFC$Direction
    combinedAssoTestRes_logFC <- combinedAssoTestRes_logFC[c(which(direction == 'Inconsistent'),
                                                             which(direction == 'Down in Both'),
                                                             which(direction == 'Up in Both')),]
  } else {
    combinedAssoTestRes_logFC <- NULL
  }
  
  return(list(tbl_featAsso1 = assoTestRes1, tbl_featAsso2 = assoTestRes2,
              tbl_pVal = combinedAssoTestRes_pVal, tbl_topSigFeats = combinedTopSigFeatTab,
              tbl_logFC = combinedAssoTestRes_logFC))
}


calcFeatCorr <- function(data1, data2, show_nFeats = F) {
  #' Systematically compute feature correlations between two PROTEOMICS (for now) datasets
  #' 
  #' Parameters
  #' data1, data2: SummarizedExperiment objects that contain data matrices or data matrices
  #' 
  #' Return
  #' corrRes: A table of correlations of overlapped features between two datasets
  
  if (all(is(data1, 'SummarizedExperiment'), is(data2, 'SummarizedExperiment'))) {
    # Extract data matrix
    dat1 <- SummarizedExperiment::assay(data1)
    dat2 <- SummarizedExperiment::assay(data2)
  } else if (all(is.matrix(data1), is.matrix(data2))) {
    dat1 <- data1
    dat2 <- data2
  } else {
    stop("Please input either SE objects containing data matrices or data matrices")
  }
  
  # Make two matrices share same sample and feature spaces and in same row and column order
  # Keep only common samples
  commonSmps <- intersect(colnames(dat1), colnames(dat2))
  dat1 <- dat1[, commonSmps]
  dat2 <- dat2[, commonSmps]
  # Keep only common features
  # Select first protein as representative if there are more than one
  rownames(dat1) <- stringr::str_remove(rownames(dat1), ';.*')
  rownames(dat2) <- stringr::str_remove(rownames(dat2), ';.*')
  commonFeats <- intersect(rownames(dat1), rownames(dat2))
  dat1 <- dat1[commonFeats,]
  dat2 <- dat2[commonFeats,]
  if (show_nFeats) {
    cat('Features in data1: ', nrow(data1), '\n', 'Features in data2: ', nrow(data2), '\n',
        'Features in common: ', length(commonFeats), sep = '')
  }
  
  # Compute feature correlations between two proteomics datasets
  # Do sanity check
  if (identical(rownames(dat1), rownames(dat2)) &
      identical(colnames(dat1), colnames(dat2))) {
    corrs <- sapply(seq_len(nrow(dat1)), function(i) {
      corr <- cor.test(dat1[i,], dat2[i,], method = 'pearson', use = "pairwise.complete.obs")
      corr$estimate
    })
    corrRes <- data.frame(Feature = rownames(dat1), Correlation = corrs) %>%
      dplyr::arrange(Correlation)
  }
  
  return(corrRes)
}
