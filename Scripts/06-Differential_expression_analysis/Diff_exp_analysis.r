
################################################################################
##                                                                            
##  06-Diff_exp_analysis.r                                                          
##                                                                            
##  1. Exploratory analysis
##
##  This program takes the tables of absolute counts from each project and
##  performs a Principal Component Analysis (PCA) for each of the stress
##  events considered in that project. From the results of this analysis,
##  it takes the coordinates generated for each sample from the values of the
##  first three principal components and calculates the Euclidean distances
##  between samples of the same condition or INTRA-group (e.g. treated1 -
##  treated2) and the distances between samples of different conditions or
##  INTER-group (e.g. treated1 - control1). Then, it performs a Mann-whitney-
##  wilcoxon test to check if there are differences between the INTRA and
##  INTER-group distances previously calculated. 
##
##  2. Differential expression analysis
##
##  Then, the program performs a differential expression analysis using
##  DESeq2. The absolute counts tables contain a group of control samples and
##  different treatment samples to which they are related. The differential
##  expression analysis is performed considering all possible combinations of
##  control vs treated (c_vs_t1, c_vs_t2, etc), so the program returns a result
##  table for each of them. The results table contains all the information
##  provided by the results() function of DESeq2 together with the
##  log2FoldChange and lfcSE from lfcShrink. In addition to the raw data
##  obtained in the analysis, this script also provides tables with those
##  sequences with an adjusted p-value lower than 0.05.
##                                                                            
##                                                                            
##  Author: Antonio Gonzalez Sanchez                                          
##  Date: 10/01/2023                                                          
##                                                                            
###############################################################################


################################### MODULES ####################################

suppressMessages(library(argparse))
suppressMessages(library(DESeq2))
suppressMessages(library(ff))
suppressMessages(library(htmltools))
suppressMessages(library(plotly))
suppressMessages(library(stringr))
suppressMessages(library(tibble))
suppressMessages(library(tidyverse))


################################## FUNCTIONS ###################################

#' Get the command line arguments
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values

get_arguments <- function() {
  
  # create parser object
  parser <- ArgumentParser(prog = '06-Diff_exp_analysis.r',
                           description = '
    This program takes the tables of absolute counts from each project and
    performs a Principal Component Analysis (PCA) for each of the stress
    events considered in that project. From the results of this analysis,
    it takes the coordinates generated for each sample from the values of
    the first three principal components and calculates the Euclidean
    distances between samples of the same condition or INTRA-group (e.g.
    treated1 - treated2) and the distances between samples of different
    conditions or INTER-group (e.g. treated1 - control1). Then, it performs
    a Mann-Whitney-Wilcoxon test to check if there are differences between
    the INTRA and INTER-group distances previously calculated. Then, the
    program performs a differential expression analysis using DESeq2. The
    absolute counts tables contain a group of control samples and different
    treatment samples to which they are related. The differential expression
    analysis is performed considering all possible combinations of control
    vs treated (c_vs_t1, c_vs_t2, etc), so the program returns a result table
    for each of them. The results table contains all the information provided
    by the results() function of DESeq2 together with the log2FoldChange and
    lfcSE from lfcShrink. In addition to the raw data obtained in the analysis,
    this script also provides tables with those sequences with an adjusted
    p-value lower than 0.05.',
                           formatter_class = 'argparse.RawTextHelpFormatter')
  
  required <- parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--input',
                        type = 'character',
                        help = 'Project directory path.',
                        required = TRUE)
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Differential expression analysis output directory path. If it does not exist, it will be created',
                        required = TRUE)
  required$add_argument('-e', '--exploratory',
                        type = 'character',
                        help = 'Exploratory analysis output directory path. If it does not exist, it will be created',
                        required = TRUE)
  parser$add_argument('-a', '--alpha',
                      default = 0.05,
                      type = 'double',
                      help = 'Alpha significance level. Default is 0.05')
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('input', 'output', 'exploratory', 'alpha')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  # Check if the input directory exists
  if (!dir.exists(args$input)) {
    stop('Error. The input directory does not exist.')
  }
  
  return(args)
}


#' Save a csv file in a DESeqdataset
#' 
#' This function receives the path to a comma-separated CSV file, reads it,
#' extracts information from the column names and creates a DESeqDataSet object.
#'
#' @param path_file Path to a comma-separated .csv file
#' @return  A DESeqDataSet
#' @examples 
#' csv_to_deseq_dataset("/home/minimind/Desktop/Results/arth/PRJNA277424_1.csv")

csv_to_deseq_dataset <- function(path_file) {

  # Get and prepare counts df
  counts_df <- read.csv(path_file, sep = ",", header = TRUE)
  
  # The df is not empty
  if (nrow(counts_df) > 0) {
  
    # Convert df to matrix
    counts_matrix <- data.matrix(counts_df[, -1])
    rownames(counts_matrix) <- counts_df[, 1]
    
    # Get ColData
    sample_v <- c()
    condition_v <- c()
    cols <- colnames(counts_matrix)
    for (col in cols) {
      
      # Create condition names
      col_elements <- strsplit(col, "_")
      condition_elements <- col_elements[[1]][-3]
      condition <- paste(condition_elements, collapse = "_")
      
      # Save control condition in variable
      if (col_elements[[1]][1] == 'control') {
        control = condition
      }
      
      # Save data in vector
      sample_v <- c(sample_v, col)
      condition_v <- c(condition_v, condition)
    }
    
    # Save sample and condition in dataframe
    coldata <- data.frame(sample = sample_v,
                          condition = condition_v,
                          row.names = 'sample')
    coldata$condition <- as.factor(coldata$condition)

    # Create DESeqDataSet
    dds <- suppressMessages(DESeqDataSetFromMatrix(countData = counts_matrix,
                                                   colData = coldata,
                                                   design = ~ condition))
    
    # Set control condition as reference
    dds$condition <- relevel(dds$condition, ref = control)
  
  # The df is empty
  } else {
    dds <- -1
  }
  return(dds)
}

#' Row variance
#' 
#' This function calculates the variance of each row in a matrix or data frame
#' "a" using the R function "apply".
#'
#' @param a A DataFrame or Matrix object
#' @return A vector with the variances of each row of the DataFrame or Matrix.
#' @examples
#' row_var(df)

row_var <- function(a) {
  apply(a, 1, var)
}


#' Euclidean distance
#' 
#' This function calculates the Euclidean distance between two points using
#' the formula: \code{sqrt(sum((a - b) ^ 2))}
#'
#' @param a A numeric vector
#' @param b A numeric vector
#' @return Euclidean distance between a and b
#' @examples
#' euclidean_dist(c(3, 15, 4), c(23,11,8))

euclidean_dist <- function(a, b) {
  return(sqrt(sum((a - b) ^ 2)))
}


#' Exploratory Analysis
#' 
#' This function receives a DESeqDataSet object and performs an individual
#' exploratory analysis for each of the stress events it contains. In this
#' analysis, a variance stabilisation transformation (VST) is performed on the
#' absolute counts data to subsequently perform a principal component analysis
#' (PCA). To discriminate whether the samples are grouped by conditions
#' (control, treated...), the Mann-Whitney-Wilcoxon test (MWW) is performed to
#' check i there are significant differences between intra-group (replicates of
#' the same condition e.g. control1 vs control2) and inter-group (samples that
#' are associated to different conditions e.g. control vs treated) distances.
#' The function returns a vector whose first 6 elements indicate the variance
#' explained by each of the first 6 principal components and the last element
#' represents the p-value obtained in the MWW test.
#'
#' @param dds A DESeqDataSet object
#' @param path_out Path of the output directory to save the generated plots
#'                 (PCA and barplot of the variance explained by each component).
#' @param file_name Name for output files
#' @return  A numeric vector with the variance explained by the first 6
#'          principal components and de p-value obtained in the MWW test.
#' @examples
#' euclidean_dist(dds,/home/minimind/Desktop/Results, PRJNA277424_1)
#'

exploratory_analysis <- function(dds, path_out_dir, file_name) {
  
  # Get sample names
  sample_names <- colnames(dds)
  
  # Obtain the list of conditions.
  condition_samples_list <- list() 
  for (sample in sample_names) {
    sample_elements <- unlist(strsplit(sample, "_"))
    sample_elements <- sample_elements[-3]
    condition_name <- paste(sample_elements, collapse = "_")
    condition_samples_list[[condition_name]] <- c(condition_samples_list[[condition_name]], sample)
  }
  
  # Iterate through the list of conditions
  ea_results <- list()
  sub_file <- 1
  for (condition in names(condition_samples_list)) {
    # Select only treated conditions (e.g. treated_heat_37ºC_6h_seedling_15d)
    if (startsWith(condition, "treated")) {
      
      # Get subfile name
      sub_file_name <- paste0(file_name, "_", sub_file)
      
      # Create output directory
      path_out <- paste0(path_out_dir,"/", sub_file_name)
      dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
      
      # Get control and treated sample names (Create sample column)
      control_samples_names <- condition_samples_list[grep("^control", names(condition_samples_list))][[1]]
      treated_samples_names <- condition_samples_list[[condition]]
      samples_v <- c(control_samples_names, treated_samples_names)
      
      # Get condition column
      control_condition_names <- rep(names(condition_samples_list)[grep("^control", names(condition_samples_list))][1]
                                     , length(control_samples_names))
      treated_condition_names <- rep(condition, length(treated_samples_names))
      condition_v <- c(control_condition_names, treated_condition_names)
      
      # Save sample and condition in dataframe (colData)
      coldata <- data.frame(sample = samples_v,
                            condition = condition_v,
                            row.names = 'sample')
      coldata$condition <- as.factor(coldata$condition)
      
      # Divide the count table to use only the information from the samples of interest.
      col_index <- match(samples_v, colnames(counts(dds)))
      sub_counts_matrix <- counts(dds)[, col_index, drop = FALSE]
      
      # Create a new DESeqDataSet using the samples of interest.
      sub_dds <- suppressMessages(DESeqDataSetFromMatrix(countData = sub_counts_matrix,
                                                         colData = coldata,
                                                         design = ~ condition))
      
      # Get comparison name
      comp <- paste0("condition_",
                     condition,
                     "_vs_",
                     names(condition_samples_list[grep("^control", names(condition_samples_list))]))
      
      ### 1.1 VST NORMALIZATION
      # Absolute counts normalization for meand vs variance plot (blind = FALSE)
      deseqds <- suppressMessages(DESeq2::estimateSizeFactors(sub_dds))
      assay(deseqds, 'counts.norm.VST.false') <- as.data.frame(assay(varianceStabilizingTransformation(deseqds, blind = FALSE)))
      
      ### 1.2 MEAN VS VARIANCE PLOT
      png(file = paste(path_out, '/', sub_file_name, '_meanvsvar.png', sep = ''),
          width     = 3.25,
          height    = 3.25,
          units     = "in",
          res       = 1200,
          pointsize = 4)
      par(mfrow = c(1, 2))
      plot(log10(rowMeans(assay(deseqds, 'counts')) + 1),
           log10(row_var(assay(deseqds, 'counts')) + 1),
           xlab = expression('Log'[10] ~ 'Mean count'),
           ylab = expression('Log'[10] ~ 'Variance'),
           main = 'Counts')
      plot(rowMeans(assay(deseqds, 'counts.norm.VST.false')),
           row_var(assay(deseqds, 'counts.norm.VST.false')),
           xlab = 'Mean count',
           ylab = 'Variance',
           main = 'VST')
      dev.off()
      
      ### 2.1 PRINCIPAL COMPONENT ANALYSIS (PCA)
      # Recalculate vst, this time with blind = TRUE for a fully unsupervised calculation
      assay(deseqds, 'counts.norm.VST.true') <- as.data.frame(assay(varianceStabilizingTransformation(deseqds, blind = TRUE)))
      
      # Perform the PCA
      pca_res <- prcomp(x=t(assay(deseqds,'counts.norm.VST.true')), rank. = 6)
      pca_df <- as.data.frame(pca_res$x)
      
      # Create the conditions column for coloring the plot
      samples <- rownames(pca_df)
      condition_col <- c()
      for (sample in samples){
        sample_v <- strsplit(sample, split = '_')
        sample_v <- sample_v[[1]][-c(3)]
        condition <- paste(sample_v, collapse = '_')
        condition_col <- c(condition_col, condition)
      }
      
      # Bind condition column to PCA counts dataframe (to specify sample color)
      pca_df <- cbind(pca_df, condition_col)
      
      # Create interactive plot
      plot <- plot_ly(pca_df,
                      x = pca_df[, 1],
                      y = pca_df[, 2],
                      z = pca_df[, 3],
                      type = 'scatter3d',
                      mode = 'markers',
                      color = ~pca_df$condition_col,
                      colors = 'Paired') %>%
        layout(scene = list(xaxis = list(title = 'PC1'),
                            yaxis = list(title = 'PC2'),
                            zaxis = list(title = 'PC3')),
               showlegend = TRUE,
               legend = list(font = list(size = 20)))
      
      # Create html file with interactive plot
      htmlwidgets::saveWidget(widget = plot,
                              file = paste(path_out, '/', sub_file_name, '.html', sep=''),
                              selfcontained = FALSE)
      
      ### 2.2 EUCLIDEAN AND INTRA-/INTER-GROUP DISTANCES
      
      # Create Null matrix
      pca_comp <- pca_res$x
      distance_matrix <-  matrix(data = rep(0, nrow(pca_comp) * nrow(pca_comp)),
                                 nrow = nrow(pca_comp),
                                 ncol = nrow(pca_comp)
      )
      rownames(distance_matrix) <- rownames(pca_comp)
      colnames(distance_matrix) <- rownames(pca_comp)
      
      # INTRA- and INTER-group distance vectors
      intra_group <- c()
      inter_group <- c()
      
      # Complete matrix and vectors with the corresponding distances
      for (i in 1:nrow(pca_comp)) {
        for (j in 1:nrow(pca_comp)) {
          # Save distances in matrix
          distance_matrix[i,j] <- euclidean_dist(pca_comp[i,][1:3], pca_comp[j,][1:3]);
          # Prevent the presence of repeated distances in the vectors.
          if (j > i) {
            # INTRA-group distances
            if (condition_col[i] == condition_col[j]) {
              intra_group <- c(intra_group, euclidean_dist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
            }
            # INTER-group distances
            else{
              inter_group <- c(inter_group, euclidean_dist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
            }
          }
        }
      }
      
      ### 2.3 MANN-WHITNEY-WILCOXON TEST
      mww_res <- wilcox.test(inter_group, intra_group, paired = FALSE)
      p_value <- mww_res[3][[1]]
      
      ### 2.4 PROPORTION OF VARIANCE EXPLAINED BY EACH COMPONENT
      # Create plot
      png(file = paste(path_out, '/', sub_file_name, '_variance.png', sep = ''),
          width     = 3.25,
          height    = 3.25,
          units     = "in",
          res       = 1200,
          pointsize = 4)
      data <- head(round(pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2) * 100, 2), n = 6)
      plot <- barplot(data, ylim = c(0, max(data) * 1.2),
                      las = 2,
                      names.arg = colnames(pca_res$x),
                      ylab = '% variance explained')
      text(x = plot, y = data, labels = data, pos = 3)
      dev.off()
      
      # Complete vector until it has 6 elements.
      # If there are no more principal components, add 0
      n_ceros <- 6 - length(data)
      ceros <- rep(0.00, n_ceros)
      data_six <- c(data, ceros)
      
      ### RESULTS
      # Return % variance explained and MWW p-value
      ea_results[[comp]] <- c(sub_file_name, data_six, p_value)
      
      sub_file <- sub_file + 1
    }
  }
  
  ### RESULTS
  # Return subfile name, % variance explained and MWW p-value (List)
  return(ea_results)
  
}


##################################### MAIN #####################################

# Get programm arguments
args <- get_arguments()

# Save the rest of the arguments in variables
path_project <- args$input
path_out <- args$output
path_out_ea <- args$exploratory
alpha_value <- args$alpha

path_project <- "/home/antonio/Escritorio/pruebas/pruebas_diffexpanalysis/Results/04-Projects_divided_by_subprojects/arth/PRJNA277424"
path_out <- "/home/antonio/Escritorio/pruebas/pruebas_diffexpanalysis/Results/06-DiffExpAnalysis_res"
path_out_ea <- "/home/antonio/Escritorio/pruebas/pruebas_diffexpanalysis/Results/05-PCA"
alpha_value <- 0.05

# Summary dataframe
sum <- data.frame()

# Exploratory analysis results table
ea_df <- data.frame()

# Get species and project name
project <- basename(path_project)
species <- basename(dirname(path_project))
  
# Create output paths
path_raw_out <- paste(path_out, '01-DEA_raw', species, project, sep = '/')
path_sig_out <- paste(path_out, '02-DEA_sig', species, project, sep = '/')

# List project files
files_list <- list.files(path = path_project)
cat('Analysing project', project, '(', species, ')\n')

# Iterate project files
for (file in files_list) {
  # Files paths
  path_project_file <- paste(path_project, file, sep = '/')
  
  # New file out (without extension)
  new_file <- gsub('.csv', '', file)
  
  # Create DeseqDataSet
  dds_outer <- csv_to_deseq_dataset(path_project_file)
  
  # If the count table is empty...
  if (!class(dds_outer) == 'DESeqDataSet') {
    next
  }
  
  # Exploratory analysis
  path_out_ea_pro <- paste(path_out_ea, species, project, new_file, sep = "/")
  ea_results <- exploratory_analysis(dds_outer, path_out_ea_pro, new_file)
    
  # Create output directories
  if (!dir.exists(path_raw_out)) {
    dir.create(path_raw_out, recursive = TRUE, showWarnings = FALSE)
    dir.create(path_sig_out, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Differential expression analysis
  dds_outer <- suppressMessages(DESeq(dds_outer))
  
  # Save results
  results_names <- resultsNames(dds_outer)
  for (comp in results_names){
    if (comp != 'Intercept') {
      
      # Get subproject file id (PRJNA277424_3_2)
      comp_subproject_id <- ea_results[[comp]][1]
      
      # Save the exploratory analysis results in data.frame
      ea_df <- rbind(ea_df, c(species, ea_results[[comp]]))
      
      # Get results
      res_deseq <- results(dds_outer, name = comp, alpha = alpha_value)
      resLFC <- suppressMessages(lfcShrink(dds_outer, coef = comp, res = res_deseq))
      
      # Add Shrunken LFC to results
      final_res <- na.omit(cbind(res_deseq,
                                 Shrunkenlog2FoldChange = resLFC$log2FoldChange,
                                 ShrunkenlfcSE = resLFC$lfcSE))
      
      # Extract significant differentially expressed miRNAs
      final_res_sig <- final_res %>%
        data.frame() %>%
        rownames_to_column(var = 'seq') %>%
        as_tibble() %>%
        dplyr::filter(padj <= alpha_value)
      
      ### RESULTS TABLES
      # Save raw final results
      final_res <- rownames_to_column(as.data.frame(final_res), var = 'seq')
      path_file_raw_out <- paste(path_raw_out, paste(comp_subproject_id, '_dea_raw.csv', sep = ''), sep = '/')
      write.csv(final_res, file = path_file_raw_out, quote = FALSE, row.names = FALSE)
      
      # Save significant final results
      path_file_sig_out <- paste(path_sig_out, paste(comp_subproject_id, '_dea_sig.csv', sep = ''), sep = '/')
      write.csv(as.data.frame(final_res_sig), file = path_file_sig_out, quote = FALSE, row.names = FALSE)
      
      ### SUMMARY TABLE ###
      # Get experiment
      file_elements <- strsplit(comp_subproject_id, split = '_')
      exp <- paste(file_elements[[1]][2], file_elements[[1]][3], sep = '.')
      
      # Significant miRNAs (p.adj < 0.05)
      final_res_sig <- as.data.frame(final_res_sig)
      sig <- nrow(final_res_sig)
      
      # Add data to dataframe
      sum <- rbind(sum, c(species, project, exp, comp, sig, nrow(final_res)))
    }
  }
  # end results loop
}
# end files loop
cat(project, '(', species, ') Done!\n')

# Save summary file
if (nrow(sum) > 0) {
  colnames(sum) <- c('species', 'Project', 'Experiment', 'Comparison', 'Padj<0.05', 'Total')
  sum_sorted <- sum[order(sum$Project, sum$Experiment), ]
  write.csv(as.data.frame(sum_sorted),
            file = paste(path_out, '/', project, '_summary.csv', sep = ''),
            quote = FALSE,
            row.names = FALSE)
}

# Save the proportion of variance explained by each component and the p-value obtained in the MWW test in a .csv file
if (nrow(ea_df) > 0) {
  colnames(ea_df) <- c('species', 'Project', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'P-value(MWW)')
  ea_df[is.na(ea_df)] <- 0
  ea_df_sorted <-  ea_df[order(ea_df$Project), ]
  write_csv(ea_df_sorted, paste(path_out_ea, '/', project, '_ea_summary_table.csv', sep = ''))
}
