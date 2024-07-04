library(vsn)
library(limma)
library(visdat)
library(ggplotify)
library(pheatmap)

doPreprocessing <- function(data, feat, smp, val, featAnno = NULL, smpAnno = NULL,
                            do_featFilt = F, cutoff = 0.67, viz_miss = F, bins = 50,
                            save_path = NULL) {
  #' Preprocess tidy long data containing information about features and samples,
  #' i.e., identifiers, abundances, and metadata. Preprocessing in this function
  #' mainly consists of two parts: Feature filtering and data normalization. Preprocessed
  #' data and its metadata is stored in SummarizedExperiment object for further use.
  #' 
  #' Parameters
  #' data: A tidy long table that must have three columns containing feature and
  #' sample identifiers and feature abundances. Feature and sample metadata is optional.
  #' feat, smp, val: A character specifying the column in the tidy long table, which
  #' will be used as row and column identifiers and matrix values in an SE object.
  #' featAnno, smpAnno: A character or a vector of characters specifying the columns
  #' in the tidy long table, which will be used as row and column annotations in
  #' an SE object. Defaults are NULL.
  #' do_featFilt: A logical indicating whether feature filtering is conducted. Default
  #' is FALSE.
  #' cutoff: A numeric specifying the cutoff for feature filtering. Default is 0.67,
  #' which indicates features quantified in less than 67% of samples are removed.
  #' viz_miss: A logical indicating whether data missingness is visualized. Default
  #' is FALSE.
  #' bins: A numeric specifying the number of bins in both vertical and horizontal
  #' directions for the hexagonal heatmap. Default is 50. 
  #' save_path: A character specifying the path to save preprocessed SE objects.
  #' '.rds' in a file name should be omitted. Default is NULL.
  #' 
  #' Return
  #' A list containing the following components:
  #' ori.data: An SE object containing the original data matrix and the feature
  #' and sample metadata,
  #' ori.data.dist: A ggplot object for visualizing the distribution of the original data,
  #' ori.smp.mean.var: A ggplot object for visualizing the sample mean-variance
  #' relationship of the original data,
  #' ori.feat.mean.var: A ggplot object for visualizing the feature mean-variance
  #' relationship of the original data,
  #' ori.data.miss: A ggplot object for visualizing the missingness of the original data,
  #' and so forth. 'filt', 'vsn', and 'medi' represent filtered and vsn and median
  #' normalized data.
  
  # Do sanity check
  if (!is(feat, 'character')) {
    stop("Argument for 'feat' should be class 'character'.")
  }
  if (!is(smp, 'character')) {
    stop("Argument for 'smp' should be class 'character'.")
  }
  
  # Set plot theme
  th <- theme_bw(base_size = 15) +
    theme(axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'),
          axis.ticks = element_line(linewidth = 0.8),
          legend.text = element_text(size = 15))
  th4MissViz <- theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0,
                                                 size = 11, face = 'bold'),
                      axis.text.y = element_blank(),
                      axis.title.y = element_text(size = 13, face = 'bold'),
                      legend.text = element_text(size = 12, face = 'bold'))
  
  # Convert original long data to SE object
  seOri <- df2SummExp(data, row_id = feat, col_id = smp, values = val,
                      row_anno = featAnno, col_anno = smpAnno)
  # Assign seOri to se, so that data for normalization can be up-to-date
  se <- seOri
  # Visualize missingness of original data
  if (viz_miss) {
    exprMat <- SummarizedExperiment::assay(seOri)
    oriMiss <- visdat::vis_miss(exprMat, cluster = T, sort_miss = T,
                                warn_large_data = F) +
      labs(y = 'Features') +
      th4MissViz
  } else {
    oriMiss <- NULL
  }
  # Visualize distribution of original data
  oriDist <- ggplot(data, aes(x=.data[[smp]], y=.data[[val]])) +
    geom_boxplot() +
    scale_y_log10() +
    labs(x = 'Sample', y = 'Abundance') +
    th + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # Visualize sample mean-variance relationship of original data
  # Compute sample means and standard deviations
  smpStats <- dplyr::group_by(data, .data[[smp]]) %>%
    dplyr::summarise(Mean = mean(.data[[val]], na.rm = T), SD = sd(.data[[val]], na.rm = T))
  oriSmpMeanVar <- ggplot(smpStats, aes(x=Mean, y=SD)) +
    geom_point(size = 4) +
    ggpubr::stat_cor(aes(label=after_stat(r.label)), method = 'pearson', size = 6) +
    th
  # Visualize feature mean-variance relationship of original data
  exprMat <- SummarizedExperiment::assay(seOri) %>%
    as.matrix()
  oriFeatMeanVar <- vsn::meanSdPlot(exprMat, ranks = T, bins = bins, plot = F)$gg +
    labs(x = 'Rank of mean', y = 'SD') +
    th
  
  # Do feature filtering
  if (do_featFilt) {
    # Remove features quantified in less than certain proportion of samples
    rmFeats <- dplyr::group_by(data, .data[[feat]]) %>%
      # Compute proportion of observed data points of each group (feature)
      dplyr::summarise(frac_nonNA = round(sum(!is.na(.data[[val]])) / length(.data[[val]]), 2)) %>%
      dplyr::ungroup() %>%
      # Keep features observed in less than certain proportion of samples to remove
      dplyr::filter(frac_nonNA < cutoff) %>%
      dplyr::pull(.data[[feat]])
    dataFilt <- dplyr::filter(data, !.data[[feat]] %in% rmFeats)
    
    # Convert filtered long data to SE object
    seFilt <- suppressMessages(
      df2SummExp(dataFilt, row_id = feat, col_id = smp, values = val,
                 row_anno = featAnno, col_anno = smpAnno)
      )
    # Assign seFilt to se, so that data normalization is on filtered data
    se <- seFilt
    if (viz_miss) {
      # Visualize missingness of filtered data
      exprMat <- SummarizedExperiment::assay(seFilt)
      filtMiss <- visdat::vis_miss(exprMat, cluster = T, sort_miss = T,
                                   warn_large_data = F) +
        labs(y = 'Features') +
        th4MissViz
    } else {
      filtMiss <- NULL
    }
    # Visualize distribution of filtered data
    filtDist <- ggplot(dataFilt, aes(x=.data[[smp]], y=.data[[val]])) +
      geom_boxplot() +
      scale_y_log10() +
      labs(x = 'Sample', y = 'Abundance') +
      th + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    seFilt <- NULL
    filtDist <- NULL
    filtMiss <- NULL
  }
  
  # Do data normalization
  # Perform VSN
  exprMat <- as.matrix(SummarizedExperiment::assay(se))
  fit <- suppressMessages(vsnMatrix(exprMat))
  seVsn <- se
  SummarizedExperiment::assay(seVsn) <- predict(fit, exprMat)
  
  # Convert SE object containing vsn normalized data to long data for plotting
  dataVsn <- summExp2df(seVsn, assay = val, row_id = feat, col_id = smp)
  # Visualize distribution of vsn normalized data
  vsnDist <- ggplot(dataVsn, aes(x=.data[[smp]], y=Value)) +
    geom_boxplot() +
    labs(x = 'Sample', y = 'Vsn normalized abundance') +
    th + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # Visualize sample mean-variance relationship of vsn normalized data
  # Compute sample means and standard deviations
  smpStats <- dplyr::group_by(dataVsn, .data[[smp]]) %>%
    dplyr::summarise(Mean = mean(Value, na.rm = T), SD = sd(Value, na.rm = T))
  vsnSmpMeanVar <- ggplot(smpStats, aes(x=Mean, y=SD)) +
    geom_point(size = 4) +
    ggpubr::stat_cor(aes(label=after_stat(r.label)), method = 'pearson', size = 6) +
    th
  # Visualize feature mean-variance relationship of vsn normalized data
  exprMat <- SummarizedExperiment::assay(seVsn)
  vsnFeatMeanVar <- vsn::meanSdPlot(exprMat, ranks = T, bins = bins, plot = F)$gg +
    labs(x = 'Rank of mean', y = 'SD') +
    th
  # Save vsn normalized data
  if (!is.null(save_path)) {
    saveRDS(seVsn, paste0(save_path, 'Vsn.rds'))
  }
  
  # Perform median normalization, for computing log2(FC) later
  exprMat <- as.matrix(SummarizedExperiment::assay(se)) %>%
    limma::normalizeBetweenArrays(exprMat, method = 'scale') %>%
    log2()
  seMedi <- se
  SummarizedExperiment::assay(seMedi) <- exprMat
  
  # Convert SE object containing median normalized data to long data for plotting
  dataMedi <- summExp2df(seMedi, assay = val, row_id = feat, col_id = smp)
  # Visualize distribution of median normalized data
  mediDist <- ggplot(dataMedi, aes(x=.data[[smp]], y=Value)) +
    geom_boxplot() +
    labs(x = 'Sample', y = 'Median normalized abundance') +
    th + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # Save median normalized data
  if (!is.null(save_path)) {
    saveRDS(seMedi, paste0(save_path, 'Medi.rds'))
  }
  
  return(list(ori.data = seOri, ori.data.dist = oriDist, ori.smp.mean.var = oriSmpMeanVar,
              ori.feat.mean.var = oriFeatMeanVar, ori.data.miss = oriMiss,
              filt.data = seFilt, filt.data.dist = filtDist, filt.data.miss = filtMiss,
              vsn.data = seVsn, vsn.data.dist = vsnDist, vsn.smp.mean.var = vsnSmpMeanVar,
              vsn.feat.mean.var = vsnFeatMeanVar, medi.data = seMedi, medi.data.dist = mediDist))
}


doSingleOmicsAnalysis <- function(summ_exp, soa_metaVar, soa_unwantVar = NULL,
                                  num_sigFeats = 6, pca_method = 'pca', num_PCs = 20,
                                  num_PCfeats = NULL, use_limma = F, use_proDA = F,
                                  alpha = 0.05, pca_metaVar = 'all', feat_anno = NULL,
                                  feat_conv = F, plot_title = NULL, ...) {
  #' Do single-omics analysis of preprocessed data stored in SE container. Analyses
  #' in this function mainly consist of two parts: Main SOA and association tests
  #' between PCs and sample metadata variables. Main SOA focuses only on question
  #' of most interest and includes all kinds of analyses (univariate and multivariate)
  #' in function 'doSOA'; associations of PCs is meant to overview sources of variance
  #' in data.
  #'
  #' Parameters
  #' summ_exp: A SummarizedExperiment object containing normalized Feature x Sample
  #' data and metadata
  #' soa_metaVar: A character indicating the variable in the sample metadata from
  #' an SE object for association tests in the main SOA. Categorical variable (factor)
  #' is used as grouping factor to separate samples into two or several groups for
  #' conducting t-test or ANOVA; numerical variable (covariate) is used to compute
  #' correlation coefficient. Note that only one sample metadata variable can be
  #' specified at a time
  #' soa_unwantVar: A character or a vector of characters specifying the metadata
  #' variables that will be accounted for by linear models (regressing out). Default
  #' is NULL. Note that this is usable only for limma and proDA for now
  #' num_sigFeats: A numeric value specifying the number of top significant features
  #' to visualize through boxplots. Default is 6 and this parameter can be set to
  #' NULL if distribution of top significant features is not of interest
  #' pca_method: A character specifying the PCA method to use, which should be one
  #' of 'pca' (default), 'ppca', or 'bpca'
  #' num_PCs: A numeric value specifying the number of PCs to estimate. The preciseness
  #' of the estimation of missing values depends on the number of PCs
  #' num_PCfeats: A numeric value specifying the number of potential PCs' top features
  #' with highest absolute loadings to extract. If NULL (default), the extraction
  #' of lists of features is skipped
  #' use_limma: A logical variable indicating whether to use limma to conduct moderated
  #' t-test for univariate association tests. Default is FALSE
  #' use_proDA: A logical variable indicating whether to use proDA to account for
  #' missing values for obtaining more reliable t-test results for univariate association
  #' tests. Default is FALSE
  #' alpha: A numeric value specifying the cutoff for statistic significance. Default
  #' is 0.05
  #' pca_metaVar: A character or a vector of characters indicating the variable(s)
  #' in the sample metadata from an SE object for association tests between PCs
  #' and sample metadata variables to overview sources of variance in data. Default
  #' is 'all', which uses all variables in sample metadata
  #' feat_anno: A data frame holding original features (e.g., protein UniProt IDs)
  #' as the row names and corresponding feature annotations (e.g., protein name,
  #' gene, etc.) in columns with whatever column names for providing extra feature
  #' information in analysis result tables or converting original feature names
  #' to the feature annotations of the first column in exploratory plots. Default
  #' is NULL
  #' feat_conv: A logical variable indicating whether to convert original feature
  #' names to the feature annotations stored in the first column of a feature annotation
  #' table provided by parameter 'feat_anno' in exploratory plots. Default is FALSE
  #' plot_title: A character specifying the title of plots. Default is NULL
  #' ...: Further arguments to be passed to 'pheatmap', e.g., 'show_rownames'
  #' 
  #' Return
  #' A list containing the following components:
  #' sig.feat.tab: A table showing significant features identified from the main SOA
  #' sig.feat.heat: A heatmap visualizing molecular signatures of significant features
  #' from the main SOA
  #' top.sig.feat.dist: Boxplots visualizing comparisons of top significant feature
  #' abundances between sample groups
  #' pval.hist: A histogram of p-values from the main SOA univariate analysis, showing
  #' the power of the data
  #' sig.pc.tab: A table showing significant PCs identified from the main SOA
  #' sig.pc.dist: Boxplots visualizing comparisons of significant PCs from the main
  #' SOA between sample groups
  #' pc.top.feat.heat: A heatmap visualizing molecular signatures of top features
  #' of significant PCs from the main SOA
  #' pc.top.feat.tab: A table showing top features of significant PCs from the main SOA
  #' SOA.res: A list containing all analysis results of the main SOA. More details
  #' can be found in the function 'doSOA'
  #' data.var.source: A table showing significant associations between PCs and sample
  #' metadata variables
  #' Note that all plots are ggplot objects
  
  # Do sanity check
  if (length(soa_metaVar) != 1) {
    stop("Argument for 'soa_metaVar' should be a character of length 1.")
  }
  if (!all(pca_metaVar %in% 'all')) {
    if (!is(pca_metaVar, 'character') | !all(pca_metaVar %in% colnames(colData(summ_exp)))) {
      stop("Argument for 'pca_metaVar' should be class 'character' and included in sample metadata.")
    }
  }
  if (!is.null(feat_anno)) {
    if (!identical(sort(rownames(feat_anno)), sort(rownames(summ_exp)))) {
      stop("Row names of 'feat_anno' should be identical to feature space of 'summ_exp'.")
    }
  }
  if (!is(feat_conv, 'logical')) {
    stop("Argument for 'feat_conv' should be class 'logical'.")
  }
  if (!is.null(plot_title)) {
    if (!is.character(plot_title)) {
      stop("Argument for 'plot_title' should be class 'character'.")
    }
  }
  
  # Set plot theme
  th <- theme_bw(base_size = 15) +
    theme(axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'),
          axis.ticks = element_line(linewidth = 0.8),
          legend.text = element_text(size = 15))
  
  # Do PCA and association tests between learned PCs and sample metadata variables
  # to overview sources of variance in data
  if (all(pca_metaVar %in% 'all')) {
    metaVarAll <- colnames(colData(summ_exp))
    datVarSour <- doSOA(summ_exp, meta_var = metaVarAll, do_onlyPCA = T, alpha = alpha,
                        pca_method = pca_method, num_PCs = num_PCs)$pcSigAssoRes
  } else {
    datVarSour <- doSOA(summ_exp, meta_var = pca_metaVar, do_onlyPCA = T, alpha = alpha,
                        pca_method = pca_method, num_PCs = num_PCs)$pcSigAssoRes
  }
  
  # Perform main SOA
  soaRes <- doSOA(summ_exp, meta_var = soa_metaVar, pca_method = pca_method,
                  num_PCs = num_PCs, num_PCfeats = num_PCfeats, alpha = alpha,
                  use_limma = use_limma, use_proDA = use_proDA, unwantVar = soa_unwantVar)
  # Retrieve needed information (check return of 'doSOA' in misc.R for more details)
  # Data matrix and sample metadata for making plots
  if (is.null(soa_unwantVar)) {
    datMat <- soaRes$data
  } else {
    datMat <- soaRes$dataCorrect
  }
  smpAnnoTab <- soaRes$smpMetadata
  # Univariate analysis results
  featAssoTab <- soaRes$featAssoRes
  featSigAssoTab <- soaRes$featSigAssoRes
  # Append feature annotations to significant association result table if specified
  if (!is.null(feat_anno)) {
    featAnno <- tibble::as_tibble(feat_anno, rownames = 'Var1')
    featSigAssoTab <- dplyr::left_join(featSigAssoTab, featAnno, by = 'Var1')
  }
  # Multivariate analysis (PCA) results
  pcTab <- soaRes$pcTab
  pcSigAssoTab <- soaRes$pcSigAssoRes
  if (!is.null(num_PCfeats)) {
    pcTopFeatTab <- soaRes$pcTopFeatTab
    # Append feature annotations to table of top features with highest absolute loadings
    # of significant PC(s) if specified
    if (!is.null(feat_anno)) {
      featAnno2 <- dplyr::rename(featAnno, Feature = Var1)
      for (i in seq_along(pcTopFeatTab)) {
        pcTopFeatTab[[i]] <- dplyr::left_join(pcTopFeatTab[[i]], featAnno2, by = 'Feature')
      }
    }
  } else {
    pcTopFeatTab <- NULL
  }
  
  # Explore univariate analysis results
  # Plot molecular signatures in input data through heatmap
  # Order significant features according to t-statistics
  featOrder <- dplyr::arrange(featSigAssoTab, dplyr::desc(Stat))$Var1
  # Cluster samples by metadata variable of interest
  smpOrder <- order(smpAnnoTab[[soa_metaVar]])
  # Subset and arrange data
  datMatSub <- datMat[featOrder, smpOrder]
  # Convert original feature names to provided feature annotations if specified
  if (!is.null(feat_anno) & feat_conv) {
    rownames(datMatSub) <- plyr::mapvalues(rownames(datMatSub),
                                           from = featAnno[[1]],
                                           to = featAnno[[2]],
                                           warn_missing = F)
  }
  # Prepare sample annotation table
  smpAnno <- dplyr::select(smpAnnoTab, Sample, all_of(soa_metaVar)) %>%
    tibble::column_to_rownames('Sample')
  # Make heatmap
  # Change parameter 'plot_title' from NULL to NA for pheatmap if NULL is specified
  if (is.null(plot_title)) {
    main <- NA
  } else {
    main <- plot_title
  }
  featSigHeat <- pheatmap(datMatSub, annotation_col = smpAnno, scale = 'row', #row scaling is across columns
                          color = colorRampPalette(c('navy', 'white', 'red'))(100),
                          cluster_rows = F, cluster_cols = F, main = main, silent = T, ...) %>%
    # Convert pheatmap object to ggplot object
    ggplotify::as.ggplot()
  
  # Visualize data of top 6 significant features through boxplots
  # Extract top significant features and needed information
  if (!is.null(num_sigFeats)) {
    if (!is.null(feat_anno) & feat_conv) {
      featConvVar <- colnames(featAnno)[2]
      topSigFeats <- featSigAssoTab[1:num_sigFeats,] %>%
        dplyr::select(Var1, pVal, all_of(featConvVar)) %>%
        dplyr::mutate(newVar1 = paste0(.data[[featConvVar]], '\n(p = ', round(pVal, 4), ')'),
                      newVar1 = factor(newVar1, levels = unique(newVar1)))
    } else {
      topSigFeats <- featSigAssoTab[1:num_sigFeats,] %>%
        dplyr::select(Var1, pVal) %>%
        dplyr::mutate(newVar1 = paste0(Var1, '\n(p = ', round(pVal, 4), ')'),
                      newVar1 = factor(newVar1, levels = unique(newVar1)))
    }
    topSigFeatDat <- tibble::as_tibble(datMat[topSigFeats$Var1,], rownames = 'Var1') %>%
      tidyr::pivot_longer(cols = -'Var1', names_to = 'Sample', values_to = 'Abundance') %>%
      dplyr::left_join(topSigFeats, by = 'Var1') %>%
      dplyr::left_join(smpAnnoTab, by = 'Sample')
    # Make boxplots
    topFeatSigDist <- ggplot(topSigFeatDat, aes(x=.data[[soa_metaVar]], y=Abundance,
                                                col=.data[[soa_metaVar]], fill=.data[[soa_metaVar]])) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.2), size = 2, show.legend = F) +
      labs(title = plot_title) +
      scale_color_manual(values=c('#00BFC4', '#F8766D')) +
      scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
      facet_wrap(vars(newVar1), scales = 'free') +
      th + theme(strip.text = element_text(size = 13, face = 'bold'))
  } else {
    topFeatSigDist <- NULL
  }
  
  # Assess data power by p-value histogram
  pValHist <- ggplot(featAssoTab, aes(x=pVal)) +
    geom_histogram(breaks = seq(0, 1, 0.1), color = 'black', fill = 'grey80') +
    scale_x_continuous(breaks = c(seq(0, 1, 0.2))) +
    labs(x = 'P-value', y = 'Frequency', title = plot_title) +
    th
  
  # Explore multivariate analysis (PCA) results
  # Visualize significant PC(s) by boxplot
  pcSigDist <- lapply(pcSigAssoTab$Var1, function(pc) {
    ggplot(pcTab, aes(x=.data[[soa_metaVar]], y=.data[[pc]],
                      col=.data[[soa_metaVar]], fill=.data[[soa_metaVar]])) +
      geom_boxplot(alpha = 0.7, outlier.shape = NA) +
      geom_jitter(position = position_jitter(0.2), size = 2, show.legend = F) +
      labs(title = plot_title) +
      scale_color_manual(values=c('#00BFC4', '#F8766D')) +
      scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
      ggpubr::stat_compare_means(method = 't.test', paired = F,
                                 method.args = list(var.equal = T),
                                 size = 6, show.legend = F) +
      th
  })
  names(pcSigDist) <- stringr::str_remove(pcSigAssoTab$Var1, ' \\(.+\\)')
  
  # Plot molecular signatures in input data through heatmap
  if (!is.null(num_PCfeats)) {
    pcTopFeatHeat <- lapply(names(pcTopFeatTab), function(pc) {
      # Order top features of significant PC(s) according to loadings
      featOrder <- dplyr::arrange(pcTopFeatTab[[pc]], dplyr::desc(Loading))$Feature
      # Subset and arrange data
      datMatSub <- datMat[featOrder, smpOrder]
      # Convert original feature names to provided feature annotations if specified
      if (!is.null(feat_anno) & feat_conv) {
        rownames(datMatSub) <- plyr::mapvalues(rownames(datMatSub),
                                               from = featAnno[[1]],
                                               to = featAnno[[2]],
                                               warn_missing = F)
      }
      # Make heatmap
      pheatmap(datMatSub, annotation_col = smpAnno, scale = 'row',
               color = colorRampPalette(c('navy', 'white', 'red'))(100),
               cluster_rows = F, cluster_cols = F, main = main, silent = T) %>%
        # Convert pheatmap object to ggplot object
        ggplotify::as.ggplot()
    })
    names(pcTopFeatHeat) <- names(pcTopFeatTab)
  } else {
    pcTopFeatHeat <- NULL
  }
  
  return(list(sig.feat.tab = featSigAssoTab, sig.feat.heat = featSigHeat,
              top.sig.feat.dist = topFeatSigDist, pval.hist = pValHist,
              sig.pc.tab = pcSigAssoTab, sig.pc.dist = pcSigDist,
              pc.top.feat.heat = pcTopFeatHeat, pc.top.feat.tab = pcTopFeatTab,
              SOA.res = soaRes, data.var.source = datVarSour))
}
