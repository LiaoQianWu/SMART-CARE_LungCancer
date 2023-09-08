library(vsn)
library(limma)
library(visdat)
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

doPreprocessing <- function(data, feat, smp, val, featAnno = NULL, smpAnno = NULL,
                            do_featFilt = F, cutoff = 0.67, viz_miss = F, save_path = NULL) {
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
  #' save_path: A character specifying the path to save preprocessed SE objects.
  #' '.rds' in a file name should be omitted. Default is NULL.
  #' 
  #' Return
  #' A list containing the following components:
  #' ori.data: An SE object containing the original data matrix and the feature
  #' and sample metadata,
  #' ori.data.dist: A ggplot object for visualizing the distribution of the original data,
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
  
  # Convert original long data to SE object
  seOri <- df2SummExp(data, row_id = feat, col_id = smp, values = val,
                      row_anno = featAnno, col_anno = smpAnno)
  # Assign seOri to se, so that data for normalization can be up-to-date
  se <- seOri
  # Visualize missingness of original data
  if (viz_miss) {
    exprMat <- SummarizedExperiment::assay(seOri)
    oriMiss <- visdat::vis_miss(exprMat, cluster = T, sort_miss = T) +
      labs(y = 'Features') +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 11, face = 'bold'),
            axis.text.y = element_blank(),
            axis.title.y = element_text(size = 13, face = 'bold'),
            legend.text = element_text(size = 12, face = 'bold'))
  } else {
    oriMiss <- NULL
  }
  # Visualize distribution of original data
  oriDist <- ggplot(data, aes(x=.data[[smp]], y=.data[[val]])) +
    geom_boxplot() +
    scale_y_log10() +
    labs(x = 'Sample', y = 'Abundance') +
    theme_bw(base_size = 15) +
    theme(axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'),
          axis.ticks = element_line(linewidth = 0.8),
          legend.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  
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
    seFilt <- df2SummExp(dataFilt, row_id = feat, col_id = smp, values = val,
                         row_anno = featAnno, col_anno = smpAnno)
    # Assign seFilt to se, so that data normalization is on filtered data
    se <- seFilt
    if (viz_miss) {
      # Visualize missingness of filtered data
      exprMat <- SummarizedExperiment::assay(seFilt)
      filtMiss <- visdat::vis_miss(exprMat, cluster = T, sort_miss = T) +
        labs(y = 'Features') +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 0, size = 11, face = 'bold'),
              axis.text.y = element_blank(),
              axis.title.y = element_text(size = 13, face = 'bold'),
              legend.text = element_text(size = 12, face = 'bold'))
    } else {
      filtMiss <- NULL
    }
    # Visualize distribution of filtered data
    filtDist <- ggplot(dataFilt, aes(x=.data[[smp]], y=.data[[val]])) +
      geom_boxplot() +
      scale_y_log10() +
      labs(x = 'Sample', y = 'Abundance') +
      theme_bw(base_size = 15) +
      theme(axis.title = element_text(face = 'bold'),
            axis.text = element_text(face = 'bold'),
            axis.ticks = element_line(linewidth = 0.8),
            legend.text = element_text(size = 15),
            axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  } else {
    seFilt <- NULL
    filtDist <- NULL
    filtMiss <- NULL
  }
  
  # Do data normalization
  # Perform VSN
  exprMat <- as.matrix(SummarizedExperiment::assay(se))
  fit <- vsnMatrix(exprMat)
  seVsn <- se
  SummarizedExperiment::assay(seVsn) <- predict(fit, exprMat)
  
  # Convert SE object containing vsn normalized data to long data for plotting
  dataVsn <- summExp2df(seVsn, assay = val, row_id = feat, col_id = smp)
  # Visualize distribution of vsn normalized data
  vsnDist <- ggplot(dataVsn, aes(x=.data[[smp]], y=Value)) +
    geom_boxplot() +
    labs(x = 'Sample', y = 'Vsn normalized abundance') +
    theme_bw(base_size = 15) +
    theme(axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'),
          axis.ticks = element_line(linewidth = 0.8),
          legend.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
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
    theme_bw(base_size = 15) +
    theme(axis.title = element_text(face = 'bold'),
          axis.text = element_text(face = 'bold'),
          axis.ticks = element_line(linewidth = 0.8),
          legend.text = element_text(size = 15),
          axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
  # Save median normalized data
  if (!is.null(save_path)) {
    saveRDS(seMedi, paste0(save_path, 'Medi.rds'))
  }
  
  return(list(ori.data = seOri, ori.data.dist = oriDist, ori.data.miss = oriMiss,
              filt.data = seFilt, filt.data.dist = filtDist, filt.data.miss = filtMiss,
              vsn.data = seVsn, vsn.data.dist = vsnDist, medi.data = seMedi, medi.data.dist = mediDist))
}
