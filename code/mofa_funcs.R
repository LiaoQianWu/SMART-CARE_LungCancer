mofa_summExp2df <- function(summ_exp, smp_type, data_acqui, data_modal) {
  #' Convert SummarizedExperiment objects to long data for creating MOFA object
  #' through 'create_mofa_from_df' to include metadata
  #' 
  #' Parameters
  #' summ_exp: An SE object containing a data matrix and metadata
  #' smp_type: A character specifying the sample type, which should be either 'Plasma' or 'Tissue'
  #' data_acqui: A character specifying the data acquisition method, which could
  #' be 'DIA', 'DDA', 'Targeted', 'Untargeted', and so on
  #' data_modal: A character specifying the data modality, which could be 'Proteomics',
  #' 'Metabolomics', and so on
  #' 
  #' Return
  #' longDat: A table (long data) containing at least five columns, 'sample', 'feature',
  #' 'value', 'view', 'group'. Extra columns may exist if sample metadata is provided
  
  # Convert SE objects to long data
  #### SPECIFIC SITUATION FOR UNTARGETED METABOLOMICS AND LIPIDOMICS
  featMetadat <- stringr::str_replace(colnames(rowData(summ_exp)), '/', '.')
  longDat <- summExp2df(summ_exp = summ_exp, row_id = 'feature', col_id = 'sample') %>%
    # Remove feature metadata and sample identifiers
    # dplyr::select(-c(colnames(rowData(summ_exp)), Identifier)) %>%
    dplyr::select(-c(featMetadat, Identifier)) %>% ####
    dplyr::rename(value = Value)
  if (smp_type == 'Plasma') {
    longDat <- dplyr::mutate(longDat,
                             group = dplyr::case_when(Condition == 'Baseline' ~ 'Baseline Group',
                                                      Condition != 'Baseline' ~ 'Follow-up Group'),
                             view = paste(data_acqui, 'Plasma', data_modal)) %>%
      dplyr::select(-Condition)
  } else if (smp_type == 'Tissue') {
    longDat <- dplyr::mutate(longDat,
                             group = 'Baseline Group',
                             view = dplyr::case_when(Condition == 'Tumor' ~ paste(data_acqui, 'Tumor', data_modal),
                                                     Condition == 'Normal' ~ paste(data_acqui, 'Normal', data_modal)),
                             sample = stringr::str_replace(sample, '_TG|_NG', '_P_B')) %>%
      dplyr::select(-Condition)
  }
  
  return(longDat)
}


trainMOFA <- function(data_comb, view_data = T, train_mofa = T, num_factors = 10, save_path) {
  #' Create MOFA object through 'create_mofa_from_df' to include metadata, train
  #' MOFA model, and save trained model
  #' 
  #' Parameters
  #' data_comb: A list of long data
  #' save_path: A character specifying the path to save trained model. The saved
  #' data formats, '.hdf5' and '.rds', in the file names should be omitted.
  
  # Create MOFA object
  # Concatenate long data
  summTab <- dplyr::bind_rows(data_comb)
  mofaObj <- MOFA2::create_mofa_from_df(summTab, extract_metadata = T)
  # Overview data in MOFA object
  if (view_data) {
    print(MOFA2::plot_data_overview(mofaObj))
  }
  
  # Train MOFA model
  if (train_mofa) {
    # Train MOFA model
    data_opts <- MOFA2::get_default_data_options(mofaObj)
    model_opts <- MOFA2::get_default_model_options(mofaObj)
    model_opts$num_factors <- num_factors
    train_opts <- MOFA2::get_default_training_options(mofaObj)
    train_opts$convergence_mode <- 'slow'
    
    mofaObj <- MOFA2::prepare_mofa(mofaObj,
                                   data_options = data_opts,
                                   model_options = model_opts,
                                   training_options = train_opts)
    
    mofaObj <- MOFA2::run_mofa(mofaObj,
                               outfile = paste0(save_path, '.hdf5'),
                               use_basilisk = T)
    # Save trained MOFA model
    saveRDS(mofaObj, paste0(save_path, '.rds'))
  }
  
  return(mofaObj)
}


mofa_vizSigFactor <- function(mofaObj, smpGroup = 'Baseline Group', alpha = 0.05, show = T) {
  #' Pinpoint significant cancer recurrence-related factors in trained MOFA model
  #' and visualize data through boxplots
  
  # Conduct association tests between learnt factors and cancer recurrence
  facTab <- tibble::as_tibble(MOFA2::get_factors(mofaObj)[[smpGroup]], rownames = 'sample')
  metaTab <- tibble::as_tibble(MOFA2::samples_metadata(mofaObj)) %>%
    dplyr::filter(group == smpGroup) %>%
    dplyr::select(sample, Recurrence)
  sigAssoRes <- testAsso(facTab, metaTab, cmn_col = 'sample') %>%
    dplyr::filter(pVal <= alpha)
  if (show) {
    print(sigAssoRes)
  }
  
  # Make boxplots of significant recurrence-related factors
  # Prepare information table and container for making and storing plots
  tab4Plot <- dplyr::left_join(facTab, metaTab, by = 'sample')
  plotList <- list()
  if (show) {
    for (i in seq_len(nrow(sigAssoRes))) {
      fac <- sigAssoRes$Var1[i]
      p <- ggplot(tab4Plot, aes(x=Recurrence, y=.data[[fac]], col=Recurrence, fill=Recurrence)) +
        geom_boxplot(alpha = 0.3, outlier.shape = NA) +
        geom_jitter(position = position_jitter(0.2), show.legend = F) +
        scale_color_manual(values=c('#00BFC4', '#F8766D')) +
        scale_fill_manual(values=c('#00BFC4', '#F8766D')) +
        ggpubr::stat_compare_means(method = 't.test', paired = F,
                                   method.args = list(var.equal = T),
                                   show.legend = F) +
        th
      print(p)
      # Store plot in list
      plotList[[i]] <- p
      names(plotList)[i] <- fac
    }
  }
  
  return(list(sigAssoRes = sigAssoRes, plotList = plotList, tab4Plot = tab4Plot))
}


mofa_rmFeatSuffix <- function(mofaObj, view, feat_anno = NULL) {
  #' Remove suffixes from feature names in trained MOFA model for more concise plots.
  #' Feature names are suffixed with views if duplicated across different views.
  #' If feat_anno is not null, features are converted to corresponding feature
  #' annotations, e.g., proteins to genes
  #' 
  #' Parameters
  #' mofaObj: A trained MOFA model
  #' view: A character or a vector of characters indicating the view from which
  #' feature names are retrieved
  #' feat_anno: A list of data frames named corresponding to specified 'view', where
  #' each table holds original features (e.g., proteins) as the row names and corresponding
  #' feature annotations (e.g., genes) in the first column with whatever column name
  #' 
  #' Return
  #' mofaObj: A MOFA model with modified feature names
  
  # Sanity check
  if (!is.null(feat_anno)) {
    if (!all(names(feat_anno) %in% view)) {
      stop("Feature annotation tables provided in a list should be named corresponding to argument for 'view'.")
    }
  }
  
  for (v in view) {
    # Retrieve weight matrix
    W <- MOFA2::get_weights(mofaObj)[[v]]
    # Remove suffix from feature names, i.e., _view
    rownames(W) <- stringr::str_remove(rownames(W), paste0('_', v))
    # Map feature names to corresponding annotations
    if (!is.null(feat_anno)) {
      featAnnoTab <- feat_anno[[v]]
      featAnno <- featAnnoTab[[1]]
      names(featAnno) <- rownames(featAnnoTab)
      rownames(W) <- featAnno[rownames(W)]
    }
    # Replace original weight matrix with new one with modified feature names
    mofaObj@expectations$W[[v]] <- W
  }
  
  return(mofaObj)
}


mofa_plotDataHeatmap <- function(mofaObj, factor, view, group = 1, num_feats,
                                 feat_anno = NULL, smp_anno = NULL, smp_cluster = NULL,
                                 rm_smpSuffix = NULL, ...) {
  #' Plot a heatmap to display features with highest weights (molecular signature)
  #' in input data. Comparing to MOFA2::plot_data_heatmap, this function helps remove
  #' suffix from feature names, resulting in more concise plot
  #' 
  #' Parameters
  #' mofaObj: An trained MOFA object
  #' factor: An integer specifying the factor of interest
  #' view: A character specifying the view for retrieving the corresponding data
  #' group: A character specifying the group for retrieving the corresponding data.
  #' Default is 1, indicating the first group
  #' num_feats: An integer specifying the number of features with highest absolute
  #' weights to extract
  #' feat_anno: A data frame holding original features (e.g., proteins) as the row
  #' names and corresponding feature annotations (e.g., gene) in the first column
  #' with whatever column name
  #' smp_anno: A data frame holding samples as the row names and sample annotations
  #' in columns, e.g, patient suffering cancer recurrence or not. Note that multiple
  #' sample annotations could be provided.
  #' smp_cluster: A character specifying the sample annotation used to cluster samples
  #' in the heatmap for easier heatmap interpretation
  #' rm_smpSuffix: A character specifying the sample suffix to remove
  #' 
  #' Return
  #' P: A heatmap
  
  # Extract factor features with highest absolute weights
  # Retrieve weight matrix from MOFA object
  W <- MOFA2::get_weights(mofaObj)[[view]]
  weightVec <- W[, paste0('Factor', factor)]
  featIdx <- order(abs(weightVec), decreasing = T)[seq_len(num_feats)]
  
  # Retrieve input data from MOFA object
  D <- MOFA2::get_data(mofaObj)[[view]][[group]]
  # Remove null samples
  nullSmp <- sapply(seq_len(ncol(D)), function(i) {
    if(all(is.na(D[, i]))) {T} else {F}
  })
  D <- D[, !nullSmp]
  # Remove suffix from feature names for more concise plots and matching feature
  # annotations if 'feat_anno' is specified
  rownames(D) <- stringr::str_remove(rownames(D), paste0('_', view))
  # Map feature names to corresponding annotations
  if (!is.null(feat_anno)) {
    featAnno <- feat_anno[[1]]
    names(featAnno) <- rownames(feat_anno)
    rownames(D) <- featAnno[rownames(D)]
  }
  # Prepare data for plotting heatmap
  if (!is.null(smp_anno)) {
    # Keep only sample annotations of interest
    if (!identical(sort(rownames(smp_anno)), sort(colnames(D)))) {
      smp_anno <- smp_anno[colnames(D),, drop = F]
    }
    # Order samples by specific annotation for easier heatmap interpretation
    if (!is.null(smp_cluster)) {
      tmp_smpAnno <- smp_anno[[smp_cluster]]
      smpOrder <- lapply(unique(tmp_smpAnno), function(level) {
        rownames(smp_anno)[which(tmp_smpAnno == level)]
      }) %>% do.call(c, .)
      # Subset input data
      DSub <- D[featIdx, smpOrder, drop = F]
    } else {
      # Subset input data
      DSub <- D[featIdx,, drop = F]
    }
  } else {
    # Subset input data
    DSub <- D[featIdx,, drop = F]
  }
  # Remove suffix of samples
  if(!is.null(rm_smpSuffix)) {
    colnames(DSub) <- stringr::str_remove(colnames(DSub), rm_smpSuffix)
  }
  
  # Make heatmap
  P <- pheatmap::pheatmap(DSub, annotation_col = smp_anno, ...)
  
  return(P)
}


mofa_keepUniFeats <- function(mofaObj, factor, view, feat_anno, to_genes = F) {
  #' Remove extra inferred proteins/genes of same features and check whether there
  #' are duplicated proteins/genes. If yes, keep one with highest weight. This function
  #' is to conduct more reliable enrichment analysis
  #' 
  #' Parameters
  #' mofaObj: An trained MOFA object
  #' factor: An integer specifying the factor of interest
  #' view: A character specifying the view for retrieving the corresponding weight matrix
  #' feat_anno: A data frame holding original features (i.e., proteins) as the row
  #' names and corresponding feature annotations (i.e., gene) in the first column
  #' with whatever column name
  #' to_genes: A logical value indicating whether the row names of the weight matrix
  #' are mapped from proteins to genes
  #' 
  #' Return
  #' seW: A SummarizedExperiment object containing the weight matrix with a modified
  #' feature space and feature metadata
  
  # Retrieve proteomics weight matrix from trained MOFA model
  W <- MOFA2::get_weights(mofaObj)[[view]]
  rownames(W) <- stringr::str_remove(rownames(W), paste0('_', view))
  colnames(feat_anno) <- 'Genes'
  
  # Sanity check
  if(!identical(rownames(W), rownames(feat_anno))) {
    stop("Feature order of specified view shold be same as that of provided feature annotation table.")
  }
  
  # Pinpoint and remove duplicated protein features from weight matrix
  # Keep first inferred protein if there are many
  modiFeat <- stringr::str_remove(rownames(W), ';.*')
  dupFeats <- unique(modiFeat[duplicated(modiFeat)])
  if (length(dupFeats) == 0) {
    # Keep first inferred protein if there are many
    # weight matrix
    rownames(W) <- modiFeat
    # feature annotation table
    rownames(feat_anno) <- stringr::str_remove(rownames(feat_anno), ';.*')
    feat_anno[[1]] <- stringr::str_remove(feat_anno[[1]], ';.*')
  } else {
    dropFeatList <- c()
    for (feat in dupFeats) {
      dropFeatIdx <- which(modiFeat == feat)
      # to keep one with highest weight
      dropFeatIdx <- dropFeatIdx[-which.max(abs(W[dropFeatIdx, paste0('Factor', factor)]))]
      dropFeatList <- c(dropFeatList, dropFeatIdx)
    }
    # Filter out unwanted duplicated protein features and corresponding gene annotations
    # weight matrix
    modiFeat <- modiFeat[-dropFeatList]
    W <- W[-dropFeatList,, drop = F]
    rownames(W) <- modiFeat
    # feature annotation table
    feat_anno <- feat_anno[-dropFeatList,, drop = F]
    rownames(feat_anno) <- stringr::str_remove(rownames(feat_anno), ';.*')
    feat_anno[[1]] <- stringr::str_remove(feat_anno[[1]], ';.*')
  }
  
  # Create SummarizedExperiment object to include all information (weight matrix
  # and feature annotations)
  if (!to_genes) {
    seW <- SummarizedExperiment(W, rowData = feat_anno)
  } else {
    # Pinpoint and remove duplicated annotated genes
    featAnno <- feat_anno[[1]]
    dupFeats <- unique(featAnno[duplicated(featAnno)])
    dropFeatList <- c()
    for (feat in dupFeats) {
      dropFeatIdx <- which(featAnno == feat)
      # to keep one with highest weight
      dropFeatIdx <- dropFeatIdx[-which.max(abs(W[dropFeatIdx, paste0('Factor', factor)]))]
      dropFeatList <- c(dropFeatList, dropFeatIdx)
    }
    # Append all NA features to remove
    naFeatIdx <- which(is.na(featAnno))
    dropFeatList <- unique(c(dropFeatList, naFeatIdx))
    # Filter out unwanted duplicated genes and corresponding proteins
    # feature annotation table (proteins)
    feat_anno <- data.frame(Proteins = rownames(W)[-dropFeatList],
                            row.names = featAnno[-dropFeatList])
    # weight matrix
    W <- W[-dropFeatList,, drop = F]
    featAnno <- featAnno[-dropFeatList]
    rownames(W) <- featAnno
    
    seW <- SummarizedExperiment(W, rowData = feat_anno)
  }
  
  return(seW)
}


mofa_rankFeatList <- function(seW, factor) {
  # Extract weight matrix
  W <- SummarizedExperiment::assay(seW)
  # Rank feature weights
  rankedGeneList <- W[, factor]
  names(rankedGeneList) <- rownames(W)
  rankedGeneList <- sort(rankedGeneList, decreasing = T)
  
  return(rankedGeneList)
}


SYMBOL2ENTREZID <- function(seW, show_quality = F) {
  # This function can be generalized by turning 'keytype' and 'columns' into parameters
  
  #' Map SummarizedExperiment object's row names from gene symbols to Entrez IDs
  #' for conducting GSEA. Input data can be from output of 'mofa_keepUniFeats' function
  #' with 'to_genes' parameter is set to TRUE
  #' 
  #' Parameters
  #' seW: An SE object that contains a weight matrix whose features are gene symbols
  #' show_quality: A logical value specifying whether the mapping quality is shown
  #' 
  #' Return
  #' seW: An SE object containing the weight matrix with Entrez IDs as row names
  
  # Create annotation mapping table
  # AnnotationDbi::columns(hs)
  # AnnotationDbi::keys(hs, keytype = 'SYMBOL')
  nomenTab <- AnnotationDbi::select(hs, keys = rownames(seW),
                                    keytype = 'SYMBOL', columns = 'ENTREZID') #multiVals = 'first'
  
  # Check mapping quality from gene symbols to Entrez IDs
  if (show_quality) {
    # Show number of protein-coding genes fail to map to Entrez IDs
    print(paste(sum(is.na(nomenTab$ENTREZID)), 'out of', nrow(nomenTab),
                'protein-coding genes do not map to Entrez IDs.'))
    # Show duplicated genes (could be duplicated genes in our feature space or a gene
    # with duplicated mappings, e.g., 'HBD' maps to two Entrez IDs)
    dupGenes <- unique(nomenTab$SYMBOL[duplicated(nomenTab$SYMBOL)])
    print(paste(length(dupGenes), 'genes have multiple mapping problems.'))
    print(nomenTab[nomenTab$SYMBOL %in% dupGenes,])
  }
  
  # Filter out duplicated and NA mappings
  nomenTab <- dplyr::filter(nomenTab, !duplicated(SYMBOL), !is.na(ENTREZID))
  
  # Update new feature annotations
  # Keep gene symbols in metadata table, or good practice should specify then as row names
  rowData(seW) <- tibble::as_tibble(rowData(seW), rownames = 'SYMBOL') %>%
    dplyr::left_join(nomenTab, by = 'SYMBOL')
  # Change row names to Entrez IDs and remove rows of NA Entrez IDs for conducting GSEA
  seW <- seW[-which(is.na(rowData(seW)$ENTREZID)),]
  rownames(seW) <- rowData(seW)$ENTREZID
  
  return(seW)
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
