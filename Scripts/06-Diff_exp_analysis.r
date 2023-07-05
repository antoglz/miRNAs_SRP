
################################################################################
##                                                                            
##  06-Diff_exp_analysis.r                                                          
##                                                                            
##  1. Exploratory analysis  
##
##  This programm takes the absolute count tables (INNER) of each project to
##  perform a Principal Component Analysis. From the results of this analysis,
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
##  If there are differences between INTRA and INTER-group distances, the
##  programm takes the absolute counts tables of the project (OUTER) to
##  perform a Differential Expression Analysis using DESeq2. The absolute
##  counts tables contain a group of control samples and different treatment
##  samples to which they are related. The differential expression analysis
##  is performed considering all possible combinations of control vs treated
##  (c_vs_t1, c_vs_t2, etc), so the program returns a result table for each
##  of them. The results table contains all the information provided by the
##  results() function of DESeq2 together with the log2FoldChange and lfcSE
##  from lfcShrink. In addition to the raw data obtained in the analysis,
##  this script also provides tables with those sequences with an adjusted
##  p-value lower than 0.05.             
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
#' 
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values
#' 

getArguments <- function(){
  
  # create parser object
  parser <- ArgumentParser(prog='06-DifExpAnalysis.r',
                           description= '
                           
   This programm takes the absolute count tables (INNER) of each project to
   perform a Principal Component Analysis. From the results of this analysis,
   it takes the coordinates generated for each sample from the values of the
   first three principal components and calculates the Euclidean distances
   between samples of the same condition or INTRA-group (e.g. treated1 -
   treated2) and the distances between samples of different conditions or
   INTER-group (e.g. treated1 - control1). Then, it performs a Mann-whitney-
   wilcoxon test to check if there are differences between the INTRA and
   INTER-group distances previously calculated. If there are differences
   between INTRA and INTER-group distances, the programm takes the absolute
   counts tables of the project (OUTER) to perform a Differential Expression
   Analysis using DESeq2. The absolute counts tables contain a group of
   control samples and different treatment samples to which they are related.
   The differential expression analysis is performed considering all possible
   combinations of control vs treated (c_vs_t1, c_vs_t2, etc), so the program
   returns a result table for each of them. The results table contains all
   the information provided by the results() function of DESeq2 together with
   the log2FoldChange and lfcSE from lfcShrink. In addition to the raw data
   obtained in the analysis, this script also provides tables with those
   sequences with an adjusted p-value lower than 0.05.',
                           formatter_class= 'argparse.RawTextHelpFormatter')
  
  required = parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-i', '--input',
                      type = 'character',
                      help = 'species directory path.',
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
#' csv2DESeqDataSet("/home/minimind/Desktop/Results/arth/PRJNA277424_1.csv")
#'

csv2DESeqDataSet <- function(path_file) {
  
  # Get and prepare counts df
  counts_df <- read.csv(path_file, sep=",", header = TRUE)
  
  # The df is not empty
  if(nrow(counts_df) > 0){
  
    # Convert df to matrix
    counts_matrix <- data.matrix(counts_df[,-1])
    rownames(counts_matrix) <- counts_df[,1]
    
    # Get ColData
    sample_v <- c()
    condition_v <- c()
    cols <- colnames(counts_matrix)
    for (col in cols){
      
      # Create condition names
      col_elements <- strsplit(col, "_")
      condition_elements <- col_elements[[1]][-3]
      condition <- paste(condition_elements, collapse = "_")
      
      # Save control condition in variable
      if (col_elements[[1]][1] == 'control'){
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
    
    # Create DESeqDataSet (Inner)
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
#' rowVar(df)
#' 

rowVar <- function(a){
  apply(a,1,var)
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
#' euclideanDist(c(3, 15, 4), c(23,11,8))
#' 

euclideanDist <- function(a, b){
  return(sqrt(sum((a - b) ^ 2)))
}


#' Exploratory Analysis
#' 
#' This function receives a DESeqDataSet object and performs an exploratory
#' analysis. In this analysis, a variance stabilisation transformation (VST) is
#' performed on the absolute counts data to subsequently perform a principal
#' component analysis (PCA). To discriminate whether the samples are grouped by
#' conditions (control, treated...), the Mann-Whitney-Wilcoxon test (MWW) is
#' performed to check i there are significant differences between intra-group
#' (replicates of the same condition e.g. control1 vs control2) and inter-group
#' (samples that are associated to different conditions e.g. control vs treated)
#' distances. The function returns a vector whose first 6 elements indicate the
#' variance explained by each of the first 6 principal components and the last
#' element represents the p-value obtained in the MWW test.
#'
#' @param dds A DESeqDataSet object
#' @param path_out Path of the output directory to save the generated plots
#'                 (PCA and barplot of the variance explained by each component).
#' @param file_name Name for output files
#' @return  A numeric vector with the variance explained by the first 6
#'          principal components and de p-value obtained in the MWW test.
#' @examples
#' euclideanDist(dds,/home/minimind/Desktop/Results, PRJNA277424_1)
#'

exploratoryAnalysis <- function(dds, path_out, file_name){
  
  # Create output directory
  dir.create(path_out, recursive = TRUE, showWarnings = FALSE)
  
  ### 1.1 VST NORMALIZATION 
  # Absolute counts normalization
  deseqds <- suppressMessages(DESeq2::estimateSizeFactors(dds))
  assay(deseqds, 'counts.norm.VST') <- as.data.frame(assay(varianceStabilizingTransformation(deseqds, blind=F)))
  
  ### 1.2 MEAN VS VARIANCE PLOT
  png(file = paste(path_out, '/', file_name, '_meanvsvar.png', sep=''),
      width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)
  par(mfrow=c(1,2))
  plot(log10(rowMeans(assay(deseqds,'counts'))+1),
       log10(rowVar(assay(deseqds, 'counts'))+1),
       xlab=expression('Log'[10]~'Mean count'),
       ylab=expression('Log'[10]~'Variance'),
       main='Counts')
  plot(rowMeans(assay(deseqds,'counts.norm.VST')),
       rowVar(assay(deseqds, 'counts.norm.VST')),
       xlab='Mean count',
       ylab='Variance',
       main='VST')
  dev.off()
  
  ### 2.1 PRINCIPAL COMPONENT ANALYSIS (PCA)
  
  # Perform the PCA
  pca_res <- prcomp(x=t(assay(deseqds,'counts.norm.VST')), rank. = 6)
  pca_df <- as.data.frame(pca_res$x)
  
  # Create the conditions column for coloring the plot
  samples = rownames(pca_df)
  condition_col = c()
  for (sample in samples){
    sample_v <-strsplit(sample, split = '_')
    sample_v <- sample_v[[1]][-c(3)]
    condition <- paste(sample_v, collapse = '_')
    condition_col <- c(condition_col, condition)
  }
  
  # Bind condition column to PCA counts dataframe (to specify sample color)
  pca_df <- cbind(pca_df, condition_col)
  
  # Create interactive plot
  plot <- plot_ly(pca_df,
                  x = pca_df[,1],
                  y = pca_df[,2],
                  z = pca_df[,3],
                  type='scatter3d',
                  mode='markers',
                  color=~pca_df$condition_col,
                  colors = 'Paired') %>%
    layout(scene = list(xaxis = list(title = 'PC1'),
                        yaxis = list(title = 'PC2'),
                        zaxis = list(title = 'PC3')),
           showlegend = TRUE,
           legend = list(font = list(size = 20)))
  
  # Create html file with interactive plot
  htmlwidgets::saveWidget(widget = plot,
                          file = paste(path_out, '/', file_name, '.html', sep=''),
                          selfcontained = FALSE)
  
  ### 2.2 EUCLIDEAN AND INTRA-/INTER-GROUP DISTANCES
  
  # Create Null matrix
  pca_comp = pca_res$x
  distance_matrix <-  matrix(data = rep(0,nrow(pca_comp) * nrow(pca_comp)),
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
    for (j in 1:nrow(pca_comp)){
      # Save distances in matrix
      distance_matrix[i,j] <- euclideanDist(pca_comp[i,][1:3], pca_comp[j,][1:3]);
      # Prevent the presence of repeated distances in the vectors.
      if (j > i){
        # INTRA-group distances
        if (condition_col[i] == condition_col[j]) {
          intra_group <- c(intra_group, euclideanDist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
        }
        # INTER-group distances
        else{
          inter_group <- c(inter_group, euclideanDist(pca_comp[i,][1:3], pca_comp[j,][1:3]))
        }}}}
  
  ### 2.3 MANN-WHITNEY-WILCOXON TEST
  mww_res <- wilcox.test(inter_group, intra_group, paired = FALSE )
  p_value <- mww_res[3][[1]]
  
  ### 2.4 PROPORTION OF VARIANCE EXPLAINED BY EACH COMPONENT 
  # Create plot
  png(file = paste(path_out, '/', file_name, '_variance.png', sep=''),
      width     = 3.25,
      height    = 3.25,
      units     = "in",
      res       = 1200,
      pointsize = 4)
  data <- head(round(pca_res$sdev ^ 2 / sum(pca_res$sdev ^ 2) * 100,2), n=6)
  plot <- barplot(data, ylim = c(0, max(data) * 1.2),
                  las=2,
                  names.arg=colnames(pca_res$x),
                  ylab='% variance explained')
  text(x = plot, y = data, labels = data, pos = 3)
  dev.off()
  
  # Complete vector until it has 6 elements. If there are no more principal components, add 0
  n_ceros <- 6 - length(data)
  ceros <- rep(0.00, n_ceros)
  data_six <- c(data, ceros)
  
  ### RESULTS
  # Return % variance explained and MWW p-value
  return(c(data_six, p_value))
  
}


##################################### MAIN #####################################

# Get programm arguments
args <- getArguments()

# Save the rest of the arguments in variables
path_species <- args$input
path_out <- args$output
path_out_ea <- args$exploratory
alpha_value <- args$alpha

# Summary dataframe
sum <- data.frame()

# Exploratory analysis results table
ea_df <- data.frame()

# Paths
path_inner <- paste(path_species, 'Inner_joined_tables', sep='/')
path_outer <- paste(path_species, 'Outer_joined_tables', sep='/')

# List species files
projects_list <- list.files(path = path_outer)
species <- basename(path_species)

# Iterate species projects
for (project in projects_list){
  
  # Create output paths
  path_raw_out <- paste(path_out, '01-DEA_raw', species, project, sep = '/')
  path_sig_out <- paste(path_out, '02-DEA_sig', species, project, sep = '/')
  path_nosig_out <- paste(path_out, '03-DEA_nosig', species, project, sep = '/')
  
  # List project files
  path_inner_project <- paste(path_inner, project, sep = '/')
  path_outer_project <- paste(path_outer, project, sep = '/')
  files_list <- list.files(path = path_outer_project)
  cat('Analysing project', project, '(', species, ')\n')
  
  # Iterate project files
  for (file in files_list){
    # Files paths
    path_inner_file <- paste(path_inner_project, file, sep = '/')
    path_outer_file <- paste(path_outer_project, file, sep = '/')
    
    # New file out (without extension)
    new_file <- gsub('.csv', '', file)
    
    # Create DeseqDataSet (Inner and outer)
    dds_inner <- csv2DESeqDataSet(path_inner_file)
    dds_outer <- csv2DESeqDataSet(path_outer_file)
    
    # If the count table is empty... 
    if (!class(dds_inner) == 'DESeqDataSet' | !class(dds_outer) == 'DESeqDataSet') {
      next
    }
    
    # Exploratory analysis
    path_out_ea_pro <- paste(path_out_ea, species, project, new_file, sep="/")
    ea_results <- exploratoryAnalysis(dds_inner, path_out_ea_pro, new_file)
    ea_df <- rbind(ea_df, c(species, new_file, ea_results))
    mww_p_value <- tail(ea_results, 1)
    
    # If the p-value of the Man-whitney-wilcoxon test is equal or less than 0.05
    if (mww_p_value <= 0.05){
      
      # Create output directories
      if (!dir.exists(path_raw_out)) {
        dir.create(path_raw_out, recursive = TRUE, showWarnings = FALSE)
        dir.create(path_sig_out, recursive = TRUE, showWarnings = FALSE)
        dir.create(path_nosig_out, recursive = TRUE, showWarnings = FALSE)
      }
      
      # Differential expression analysis
      dds_outer <- suppressMessages(DESeq(dds_outer))
      
      # Save results
      results_names <- resultsNames(dds_outer)
      i = 1
      for (comp in results_names){
        if (comp != 'Intercept'){
          
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
          
          # Extracting non-significant miRNAs.
          final_res_nosig <- final_res %>% 
            data.frame() %>%
            rownames_to_column(var = 'seq') %>%
            as_tibble() %>%
            dplyr::filter(padj > alpha_value)
          
          ### RESULTS TABLES
          # Save raw final results
          final_res <- rownames_to_column(as.data.frame(final_res), var = 'seq')
          path_file_raw_out <- paste(path_raw_out, paste(new_file,'_', i,'_dea_raw.csv', sep=''), sep='/')
          write.csv(final_res, file=path_file_raw_out, quote=FALSE, row.names = FALSE)
          
          # Save significant final results
          path_file_sig_out <- paste(path_sig_out,paste(new_file,'_', i,'_dea_sig.csv', sep=''),sep='/')
          write.csv(as.data.frame(final_res_sig), file=path_file_sig_out, quote=FALSE, row.names = FALSE)
          
          # Save significant final results
          path_file_nosig_out <- paste(path_nosig_out,paste(new_file,'_', i,'_dea_nosig.csv', sep=''),sep='/')
          write.csv(as.data.frame(final_res_nosig), file=path_file_nosig_out, quote=FALSE, row.names = FALSE)
          
          ### SUMMARY TABLE ###
          # Get experiment
          file_elements <- strsplit(new_file, split = '_')
          exp <- paste(file_elements[[1]][2], i, sep = '.')
          
          # Significant miRNAs (p.adj < 0.05)
          final_res_sig <- as.data.frame(final_res_sig)
          sig <- nrow(final_res_sig)
          
          # Add data to dataframe
          sum <- rbind(sum, c(species, project, exp, comp, sig, nrow(final_res)))
          
          i = i + 1
        }
      }
      # end results loop
    }
  }
  # end files loop
  cat(project,'(', species, ') Done!\n')
}
# end projects loop


# Save summary file
if (nrow(sum) > 0){
  colnames(sum) <- c('species', 'Project', 'Experiment', 'Comparison', 'Padj<0.05', 'Total')
  write.csv(as.data.frame(sum),
            file=paste(path_out, '/', species,'_summary.csv', sep=''),
            quote=FALSE,
            row.names = FALSE)
}

# Save the proportion of variance explained by each component and the p-value obtained in the MWW test in a .csv file
if (nrow(ea_df) > 0){
  colnames(ea_df) <- c('species', 'Project','PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'P-value(MWW)')
  ea_df[is.na(ea_df)] <- 0
  write_csv(ea_df, paste(path_out_ea,'/', species, '_ea_summary_table.csv', sep=''))
}



