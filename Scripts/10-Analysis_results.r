
################################################################################
##                                                                            
##  10-Analysis_results.r                                                          
##                                                                            
##  This program uses different results obtained by running some of the
##  pipeline programs to generate the "nodes" and "edges" tables needed
##  to build the miRNA stress response network with Cytoscape. In addition,
##  it provides information as a table or plot on the classification of
##  miRNAs according to their range of stress response, the number of
##  differentially expressed sequences obtained and how many of them have
##  been annotated as miRNAs and an overview of the miRNA-mediated stress
##  response in the different plant species analyzed (Heatmap). 
##                                                                            
##  Author: Antonio Gonzalez Sanchez                                          
##  Date: 16/02/2021                                                         
##                                                                            
###############################################################################


################################### MODULES ####################################

suppressMessages(library(argparse))
suppressMessages(library(Ckmeans.1d.dp))
suppressMessages(library(ComplexHeatmap))
suppressMessages(library(dplyr))
suppressMessages(library(ggplot2))
suppressMessages(library(grid))
suppressMessages(library(gridExtra))
suppressMessages(library(ggthemes))

################################## FUNCTIONS ###################################

#' Get the command line arguments
#' 
#' This function parse the command line arguments entered into the program.
#'
#' @return List with the argument values
#' 

getArguments <- function(){
  
  # create parser object
  parser <- ArgumentParser(prog='10-Analysis_results.r ',
                           description= '
  This program uses different results obtained by running some of the
  pipeline programs to generate the "nodes" and "edges" tables needed
  to build the miRNA stress response network with Cytoscape. In addition,
  it provides information as a table or plot on the classification of
  miRNAs according to their range of stress response, the number of
  differentially expressed sequences obtained and how many of them have
  been annotated as miRNAs and an overview of the miRNA-mediated stress
  response in the different plant species analyzed (Heatmap). ',
                           formatter_class= 'argparse.RawTextHelpFormatter')
  
  required = parser$add_argument_group('required arguments')
  
  # specify our desired options 
  # by default ArgumentParser will add an help option 
  required$add_argument('-t', '--preabstable',
                        type = 'character',
                        help = 'Path of the CSV file representing the presence-absence table of miRNAs in each stress event analyzed.',
                        required = TRUE)
  required$add_argument('-i', '--ids',
                        type = 'character',
                        help = 'Path of the CSV file representing the identifier assigned to each stress event.',
                        required = TRUE)
  required$add_argument('-s', '--spnames',
                        type = 'character',
                        help = 'Path of the CSV file representing the name of the species analyzed and the abbreviation used for each of them.',
                        required = TRUE)
  required$add_argument('-d', '--diffexp',
                        type = 'character',
                        help = 'Path of the CSV file representing the number of differentially expressed small RNAs in each stress event analyzed.',
                        required = TRUE)
  required$add_argument('-m', '--mirnasannot',
                        type = 'character',
                        help = 'Path of the CSV file representing the number of differentially expressed small RNAs annotated as miRNAs in each stress event analyzed.',
                        required = TRUE)
  required$add_argument('-o', '--output',
                        type = 'character',
                        help = 'Output directory path',
                        required = TRUE)
  
  # Arguments list
  args <- parser$parse_args()
  
  #  Check for missing arguments
  expected_arguments <- c('preabstable', 'ids', 'spnames', 'diffexp', 'mirnasannot', 'output')
  if (any(sapply(args, is.null))) {
    empty_args <- names(args[sapply(args, is.null)])
    error_message <- paste('\n\tError. Unspecified argument:', empty_args, sep = ' ')
    stop(error_message)
  }
  
  return(args)
}


#' Species and stresses proportion
#' 
#' La función analyzes the proportion of species and stresses in which
#' each miRNA is represented.
#'
#' @param pre_abs_table Presence-absence table of miRNAs in the different stress
#'                      events analyzed.
#' @return  List of tables with the proportion of species and stresses in
#'          which miRNAs are represented.
#' @examples 
#' SpeciesAndStressesProportion(pre_abs_table)
#'

SpeciesAndStressesProportion <- function(pre_abs_table){
  
  #### Get species and stresses vectors
  species <- c()
  stresses <- c()
  # Iterate identifiers
  for (i in 1:ncol(pre_abs_table)){
    
    # Split identifier
    event <- colnames(pre_abs_table)[i]
    event_elements = strsplit(event, ".", fixed =TRUE)
    
    # Save species id
    specie <- gsub('X','',event_elements[[1]][1])
    species <- c(species, specie)
    
    # Save stress id
    stress <- event_elements[[1]][2]
    stresses <- c(stresses, stress)
  }
  
  # Add species and stresses rows to initial dataframe
  pre_abs_table <- rbind(pre_abs_table, species, stresses)
  rownames(pre_abs_table)[nrow(pre_abs_table) - 1] <- "species"
  rownames(pre_abs_table)[nrow(pre_abs_table)] <- "stresses"
  
  ## Species
  specie_v <- unique(as.numeric(pre_abs_table[which(rownames(pre_abs_table) == "species"), ]))
  
  # Create output specie table
  table_mod_sp <- data.frame(matrix(0, nrow = nrow(pre_abs_table) - 2, ncol = length(specie_v)))
  colnames(table_mod_sp) <- specie_v
  rownames(table_mod_sp) <- head(rownames(pre_abs_table), -2)
  
  # Iterate specie vector
  for (i in 1:length(specie_v)){
    
    # Get specie pre_abs_table
    sp <- specie_v[i]
    table_specie <- head(as.data.frame(pre_abs_table[,(pre_abs_table[nrow(pre_abs_table) - 1,]) == sp]), -2)
    colnames(table_specie) <- colnames(pre_abs_table)[which((pre_abs_table[nrow(pre_abs_table) - 1,]) == sp)]
    rownames(table_specie) <- head(rownames(pre_abs_table), -2)
    
    # Iterate table_specie rows
    for (j in 1:nrow(table_specie)) {
      # Get miRNA name
      miRNA <- rownames(table_specie)[j]
      # Check if any of the columns associated with the stress contains a 1
      if (any(table_specie[j, ] == 1)) {
        table_mod_sp[miRNA, sp] <- 1
      } else {
        table_mod_sp[miRNA, sp] <- 0
      }
    }
  }
  
  Number_of_sp <- rowSums(table_mod_sp[])
  total_sp <- sp_col <- rep(ncol(table_mod_sp) ,times=length(Number_of_sp))
  sp_df <- data.frame(row.names = rownames(table_mod_sp),
                      de_in_sp_number = Number_of_sp, Total_sp_number = total_sp,
                      proportion_sp = (Number_of_sp * 100) / total_sp)
  
  ## Stress
  stress_v <- sort(as.numeric(unique(sapply(strsplit(colnames(pre_abs_table), "\\."), function(x) x[2]))))
  
  # Create output stress table
  table_mod_st <- data.frame(matrix(0, nrow = nrow(pre_abs_table) - 2, ncol = length(stress_v)))
  colnames(table_mod_st) <- stress_v
  rownames(table_mod_st) <- head(rownames(pre_abs_table), -2)
  
  # Iterate specie vector
  for (i in stress_v){
    
    # Get specie table
    selected_cols <- colnames(pre_abs_table)[sapply(strsplit(colnames(pre_abs_table), "\\."), function(x) x[2] == as.character(i))]
    table_stress <- as.data.frame(head(pre_abs_table[, selected_cols], -2))
    colnames(table_stress) <- selected_cols
    rownames(table_stress) <- head(rownames(pre_abs_table), -2)
    
    # Iterate table_specie rows
    for (j in 1:nrow(table_stress)) {
      # Get miRNA name
      miRNA <- rownames(table_stress)[j]
      # Check if any of the columns associated with the stress contains a 1
      if (any(table_stress[j, ] == 1)) {
        table_mod_st[miRNA, as.character(i)] <- 1
      } else {
        table_mod_st[miRNA, as.character(i)] <- 0
      }
    }
  }
  
  Number_of_st <- rowSums(table_mod_st[])
  total_st <- st_col <- rep(ncol(table_mod_st) ,times=length(Number_of_st))
  st_df <- data.frame(row.names = rownames(table_mod_st), de_in_st_number = Number_of_st, Total_st_number = total_st)
  st_df <- data.frame(row.names = rownames(table_mod_st),
                      de_in_st_number = Number_of_st, Total_st_number = total_st,
                      proportion_st = (Number_of_st * 100) / total_st)
  
  # Create output_list
  output_list <- list(species = sp_df, stress = st_df)

  return(output_list)
  
}


#' Check for bias due to overrepresentation of species and stresses
#' 
#' This function generates two scatter plots to study the possible bias due
#' to overrepresentation of certain species or stresses in the analysis. These
#' plots represent the miRNAs ordered from most general to least general
#' according to their range of stress response and the proportion of
#' species/stresses in which each is represented. In addition, the function
#' also returns the tables used to generate both plots.
#'
#' @param pre_abs_table Presence-absence table of miRNAs in the different stress
#'                      events analyzed.
#' @param miRNA_class_table Table with the classification of miRNAs according
#'                          to their range of stress response.
#' @param save_table_dir Output directory path for the output tables.
#'                       Optional (Default = "")
#' @param save_plot_dir Output directory path for the output plots. Optional
#'                      (Default = "")
#' @return  List of plots representing the proportion of species and stresses
#'          in which miRNAs are represented and the tables used to generate them.
#' @examples
#' CheckSpeciesAndStressesBias(pre_abs_table, miRNA_class_table)
#' CheckSpeciesAndStressesBias(pre_abs_table, miRNA_class_table,
#'                             save_table_dir= "/home/dir/tables")
#' CheckSpeciesAndStressesBias(pre_abs_table, miRNA_class_table,
#'                             save_plot_dir= "/home/dir/plots")
#' CheckSpeciesAndStressesBias(pre_abs_table, miRNA_class_table,
#'                             save_table_dir= "/home/dir/tables",
#'                             save_plot_dir= "/home/dir/plots")
#'

CheckSpeciesAndStressesBias <- function(pre_abs_table, miRNA_class_table, save_table_dir="", save_plot_dir=""){
  
  # Create Output list
  output_list <- list()
  
  ## 1. Prepare table
  # Proportion of species and stresses in which each miRNA is represented.
  proportion_results <- SpeciesAndStressesProportion(pre_abs_table)
  
  # Join results tables
  table_res_final <- merge(miRNA_class_list$table,  proportion_results$species, by = "row.names")
  rownames(table_res_final ) <- table_res_final$Row.names
  table_res_final$Row.names <- NULL
  table_res_final <- merge(table_res_final,  proportion_results$stress, by = "row.names")
  rownames(table_res_final ) <- table_res_final$Row.names
  table_res_final$Row.names <- NULL
  
  # Add a new column with the generality order
  table_res_final <- table_res_final %>%
    mutate(Order = dense_rank(desc(per)))
  
  # Sort table by generality order
  table_res_final_sort <- table_res_final %>% arrange(Order)
  
  # Create output_table
  output_table <- table_res_final_sort[, c("miRNA_fam", "Order", "clusters",
                                           "de_in_sp_number", "Total_sp_number",
                                           "proportion_sp", "de_in_st_number",
                                           "Total_st_number", "proportion_st")]
  output_list$table <- output_table
  
  ## 2. Create plots
  plots_list <- list()
  # Species plot
  plot_sp <- ggplot(table_res_final_sort, aes(x = Order, y = proportion_sp, color= as.factor(clusters))) +
    geom_point() +
    scale_color_manual(values = c("#CE5A57", "#6C91BF", "#E1B16A"),
                       labels = c("Específicos", "Intermediarios", "Generales")) +
    theme_few()  +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10)) +
    labs(y = "Proporción de especies (%)", x = "Familias de miARNs en orden de generalidad")
  
  # Save plot_sp in list
  plots_list$plot_sp <- plot_sp
  
  # Stresses plot
  plot_st <- ggplot(table_res_final_sort, aes(x = Order, y = proportion_st, color= as.factor(clusters))) +
    geom_point() +
    scale_color_manual(values = c("#CE5A57", "#6C91BF", "#E1B16A"),
                       labels = c("Específicos", "Intermediarios", "Generales")) +
    theme_few()  +
    theme(legend.title = element_blank(),
          legend.text = element_text(size = 10)) +
    labs(y = "Proporción de agentes inductores de estrés (%)", x = "Familias de miARNs en orden de generalidad")
  
  # Save plot_st in list
  plots_list$plot_st <- plot_st
  
  # Save plot_list in output_list
  output_list$plots <- plots_list
  
  ## 3. Save table and plots
  if (nchar(save_table_dir) != 0){
    write.csv2(output_table,paste(save_table_dir, "sp_st_bias_table.csv", sep="/"), quote = FALSE, row.names = FALSE)
  }
  if (nchar(save_plot_dir) != 0){
    ggsave(paste(save_plot_dir, "species_proportion.png", sep="/"), plot_sp, dpi = 400)
    ggsave(paste(save_plot_dir, "stresses_proportion.png", sep="/"), plot_st, dpi = 400)
  }
  
  return(output_list)
  
}


#' Classification of miRNAs according to their range of stress response.
#' 
#' This function ranks the miRNAs according to their range of stress response.
#' This is represented in a bar plot where the length of the bars represents
#' the proportion of stress events in which each miRNA is represented.
#'
#' @param pre_abs_table Presence-absence table of miRNAs in the different stress
#'                      events analyzed.
#' @param save_table_dir Output directory path for the output tables.
#'                       Optional (Default = "")
#' @param save_plot_dir Output directory path for the output plots. Optional
#'                      (Default = "")
#' @return List with the plot representing the distribution of miRNA families
#'         classified according to their range of stress response and the table
#'         used to generate it.
#' @examples
#' miRNAClasification(pre_abs_table)
#' miRNAClasification(pre_abs_table, save_table_dir= "/home/dir/tables")
#' miRNAClasification(pre_abs_table, save_plot_dir= "/home/dir/plots")
#' miRNAClasification(pre_abs_table, save_table_dir= "/home/dir/tables",
#'                    save_plot_dir= "/home/dir/plots")
#'

miRNAClasification <- function(pre_abs_table, save_table_dir="", save_plot_dir="") {
  
  ## 1. Create summary table
  # Results list
  output_list <- list()
  
  # Create output table columns.
  Number_of_events <- rowSums(pre_abs_table)
  num_max_oc <- rep(ncol(pre_abs_table), times=nrow(pre_abs_table))
  per_oc <- (Number_of_events / num_max_oc) * 100
  
  # Classify miRNAs according to their stress response range.
  results <- Ckmeans.1d.dp(Number_of_events, 3)
  clusters <- results$cluster
  
  # Create output table
  output_table <- data.frame("miRNA_fam" = rownames(pre_abs_table),
                             "Occurrences" = Number_of_events,
                             "Num_Max_Occurrences" = num_max_oc,
                             "Percentage_of_occurrence" = per_oc,
                             "clusters" = clusters)
  # Delete rows with 0 occurrences
  output_table <-output_table[output_table$Occurrences > 0,]
  
  # Sort table
  output_table <- output_table[order(output_table$Occurrences, decreasing = TRUE),]
  
  # Add miRNA distribution column.
  per <- sapply(nrow(output_table):1, function(i) i * 100 / nrow(output_table))
  output_table <- cbind(output_table, per)
  
  # Save output_table in output_list
  output_list$table <- output_table
  
  ### 2. Create plot
  # Main plot
  plot_main <- ggplot(output_table, aes(x = per, y = Percentage_of_occurrence, fill = as.factor(clusters))) +
    geom_bar(stat = "identity") +
    coord_flip() +
    scale_fill_manual(values = c("#CE5A57", "#6C91BF", "#E1B16A"),
                      labels = c("Específicos", "Intermediarios", "Generales")) +
    theme_few()  +
    theme(strip.text = element_text(face = "italic"),
          legend.title = element_blank()) +
    guides(fill = guide_legend(reverse = TRUE)) +
    labs(y = "Proporción de eventos de estrés (%)", x = "Distribución de familias de miARNs (%)")
  
  # Secondary plot
  plot_zoom <- ggplot(output_table[output_table$clusters == 3, ], aes(x = reorder(miRNA_fam, Percentage_of_occurrence), Percentage_of_occurrence)) +
    geom_bar(stat = "identity", fill = "#E1B16A") +
    coord_flip() +
    theme_few()  +
    theme(legend.position = "none",
          axis.text.y = element_text(
            face = ifelse(output_table$miRNA_fam == "miRN485", "bold", "plain")
          )) +
    labs(x = "",y = "")
  
  # Join plots
  zoom_x <- c(20, 80)
  zoom_y <- c(20, 80)
  
  combined_plot <- plot_main +
    annotation_custom(ggplotGrob(plot_zoom), xmin = zoom_x[1], xmax = zoom_x[2], ymin = zoom_y[1], ymax = zoom_y[2]) +
    geom_rect(xmin = zoom_x[1] -1, xmax = zoom_x[2] + 1, ymin = zoom_y[1] + 1, ymax = zoom_y[2] + 1, color = "#666666", fill = NA, size = 0.2) +
    geom_rect(xmin = 97, xmax = 101, ymin = 0, ymax = 80, color = "#666666", fill = NA, size = 0.2) +
    geom_segment(aes(x = 97, xend = zoom_x[2] + 1, y = 0, yend = zoom_y[1] + 1), color = "#666666", size = 0.2) + 
    geom_segment(aes(x = 97, xend = zoom_x[2] + 1, y = 80, yend = zoom_y[2] + 1), color = "#666666", size=0.2)
  
  # Save plot in output_list
  output_list$plot <- combined_plot
  
  ## 3. Save table and plots
  if (nchar(save_table_dir) != 0){
    write.csv2(output_table,paste(save_table_dir, "ocurrences_table.csv", sep="/"), quote = FALSE, row.names = FALSE)
  }
  if (nchar(save_plot_dir) != 0){
    ggsave(paste(save_plot_dir, "barplot.png", sep="/"), combined_plot, dpi = 600)
  }
  
  return(output_list)
}


#' Summary of differential expression analysis.
#' 
#' This function graphically represents the number of differentially expressed
#' sequences in each of the stress events analyzed and how many of them have
#' been annotated as miRNAs.
#'
#' @param desRNAs_table Table with the number of sequences differentially
#'                      expressed in each stress event.
#' @param miRNA_annotated_table Table with the number of differentially
#'                              expressed sequences annotated as miRNAs in
#'                              each stress event.
#' @param save_plot_dir Output directory path for the output plots. Optional
#'                      (Default = "")
#' @examples
#' DEASummaryPlot(desRNAs_table, miRNA_annotated_table)
#' DEASummaryPlot(desRNAs_table, miRNA_annotated_table,
#'                save_plot_dir="/home/dir/plots")
#'

DEASummaryPlot <- function(desRNAs_table, miRNA_annotated_table, save_plot_dir=""){
  
  ## 1. Prepare the table
  # Creating a common identifier between the two tables
  sig_srnas_table$Complete_id <- paste(sig_srnas_table$Project, sub("\\.", "_", sig_srnas_table$Experiment), sep = "_")
  
  # Join the tables
  final_table <- merge(sig_srnas_table, subset(miRNAs_ident_table, select = -c(species, Annot_precursor)), by.x = "Complete_id", by.y = "Experiment", all = FALSE)
  
  # Sort the table in descending order based on the "Padj0.05" column
  final_table_sig <- final_table %>% arrange(desc(Padj.0.05))
  
  # Create a new column indicating the order position
  final_table_sig <- final_table_sig %>% mutate(Orden = row_number())

  # Set the resolution and size of the PNG plot
  png(paste(save_plot_dir, "dea_sum.png", sep="/"), width = 10, height = 8, units = "in", res = 700)
  
  # Create plot of differentially expressed sequences
  plot1 <- ggplot(data=final_table_sig, aes(x=Orden, y = Padj.0.05)) +
    geom_bar(stat = "identity", position=position_dodge(), fill="#6C91BF") +
    theme_few()  +
    theme(plot.margin = margin(30,20,10,40)) +
    labs(y = "Secuencias diferencialmente expresadas", x = "Eventos de estrés")
  
  ## 2. Create the plot
  # Create plot of sequences annotated as miRNAs
  plot2 <- ggplot(data=final_table_sig, aes(x=Orden, y = Annot_miRNAs)) +
    geom_bar(stat = "identity", position=position_dodge(), fill="#6C91BF") +
    theme_few()  +
    theme(plot.margin = margin(30,20,10,40)) +
    labs(y = "miARNs anotados", x = "Eventos de estrés") 
  
  # Join the plots
  # Arrange the plots in a grid with 2 rows and 1 column
  combined_plot <- grid.arrange(plot1, plot2, nrow = 2, ncol = 1)
  grid.text("A", x = unit(0.020, "npc"), y = unit(0.97, "npc"), gp = gpar(fontsize = 16, fontface = "bold"))
  grid.text("B", x = unit(0.020, "npc"), y = unit(0.47, "npc"), gp = gpar(fontsize = 16, fontface = "bold"))
  dev.off()
}


#' Create Heatmap of stress-reactive miRNAs
#' 
#' This function creates a heatmap with the differentially expressed miRNAs
#' in each independent stress event. This heatmap is divided according to the
#' species analyzed in order to compare the miRNA expression patterns of each
#' of them.
#'
#' @param pre_abs_table Presence-absence table of miRNAs in the different stress
#'                      events analyzed.
#' @param save_plot_dir Output directory path for the output plots. Optional
#'                      (Default = "")
#' @return Heatmap with de miRNA expression pattern of each species.
#' @examples
#' miRNAsStressEventsHeatmap(pre_abs_table)
#' miRNAsStressEventsHeatmap(pre_abs_table, save_plot_dir="/home/dir/plots")
#'

miRNAsStressEventsHeatmap <- function(pre_abs_table, save_plot_dir="") {
  
  # Prepare presence-absence matrix
  matrix_heatmap <- as.matrix(pre_abs_table, stringsAsFactors =FALSE)
  mat_num <- matrix(as.numeric(matrix_heatmap),
                    ncol = ncol(matrix_heatmap))
  colnames(mat_num) <- colnames(matrix_heatmap)
  rownames(mat_num) <- rownames(matrix_heatmap)
  
  # Obtain species identifiers
  species <- sapply(strsplit(colnames(mat_num), "\\."), `[`, 1)
  
  # Creating a list of submatrices using lapply
  submatrix_list <- lapply(unique(species), function(num) mat_num[, species == num])
  species_v <- unique(as.numeric(gsub("X", "", species)))
  species_names <- unique(info_table$Specie_name[info_table$Specie_id %in% species_v])
  names(submatrix_list) <- species_names
  
  # Create heatmap list
  heatmaps <- lapply(names(submatrix_list), function(x) Heatmap(as.matrix(submatrix_list[[x]]),
                                                                col = c("#FFF2AE", "#CE5A57"),
                                                                row_names_gp = gpar(fontsize = 6),
                                                                border = TRUE,
                                                                show_row_dend = FALSE,
                                                                show_column_dend = FALSE,
                                                                row_title = NULL,
                                                                show_row_names = FALSE,
                                                                show_heatmap_legend = FALSE,
                                                                show_column_names = FALSE,
                                                                column_title_side = "bottom",
                                                                column_title_rot = 90,
                                                                column_title = gt_render(paste0("*", x, "*"))
  ))
  
  # Concatenate heatmaps horizontally
  heatmaps_conc <- heatmaps[[1]]
  for (i in 2:length(heatmaps)) {
    heatmaps_conc <- heatmaps_conc + heatmaps[[i]]
  }
  
  # Save Heatmap in png format
  if (nchar(save_plot_dir) != 0){
    # Save heatmap as png
    png(paste(save_plot_dir, "heatmap.png", sep="/"), width = 10000, height = 7000, res = 700)
    draw(heatmaps_conc)
    dev.off()
  }
  
  return(heatmaps_conc)
}


#' Create Cytoscape network tables
#' 
#' This function generates the tables needed to build the miRNA-mediated stress
#' response network with Cytoscape: nodes and edges. The nodes table represents
#' the nodes of the network (each miRNA) and the classification to which each
#' one belongs (General=3, Intermediate=2 or Especific=1). On the other hand,
#' the edges table indicates which miRNAs are connected and the weight of the
#' connections. The connections establish that miRNAs are jointly represented
#' in at least one stress event, while the weight of connections establishes
#' the degree of co-occurrence.
#'
#' @param pre_abs_table Presence-absence table of miRNAs in the different stress
#'                      events analyzed.
#' @param miRNA_class_table Table with the classification of miRNAs according
#'                          to their range of stress response.
#' @param save_table_dir Output directory path for the output tables.
#' @return List with the output tables
#' @examples
#' CreateCytoscapeTables(pre_abs_table, save_table_dir)
#' CreateCytoscapeTables(pre_abs_table, save_table_dir= "/home/dir/tables")
#'

CreateCytoscapeTables <- function(pre_abs_table, miRNA_class_table, save_table_dir="") {
  
  ## 1. Create node dataframe
  cat("Preparing nodes and edges tables", sep="\n")
  nodes <- data.frame(id = miRNA_class_table$miRNA_fam,
                      score = as.integer(miRNA_class_table$Occurrences),
                      group = miRNA_class_table$clusters)
  
  ## 2. Create edges dataframe
  comb <- combn(rownames(pre_abs_table), 2)
  edges <- data.frame()
  for (i in 1:ncol(comb)) {
    miRNA1 <- comb[1, i]
    miRNA2 <- comb[2, i]
    events <- as.numeric(sum(pre_abs_table[miRNA1, ] == 1 & pre_abs_table[miRNA2, ] == 1))
    edges <- rbind(edges, c(miRNA1, miRNA2, events))
  }
  # Use column names recognizable by Cytoscape
  colnames(edges) <- c("source", "target", "weight")
  edges$weight <- as.numeric(edges$weight)
  cat("Nodes and edges tables prepared", sep="\n")
  
  # Remove interactions that do not exist
  edges <- edges[edges$weight != 0, ]
  
  ## 3. Save tables
  if (nchar(save_table_dir) != 0){
    write.csv2(nodes,paste(save_table_dir, "nodes.csv", sep="/"), quote = FALSE, row.names = FALSE)
    write.csv2(edges,paste(save_table_dir, "edges.csv", sep="/"), quote = FALSE, row.names = FALSE)
  }
  
  return(list(nodes=nodes, edges=edges))

}


##################################### MAIN #####################################

# Get programm arguments
args <- getArguments()

# Save the arguments in variables
pre_abs_table_path <- args$preabstable
ids_path <- args$ids
species_name_path <- args$spnames
sRNAs_sig_path <- args$diffexp
miRNAs_ident_path <- args$mirnasannot
output_directory_path <- args$output


pre_abs_table_path="/home/antonio/Escritorio/Antonio_G/Antonio_TFM/Resultados/09-Presence_absence_table_prueba/presence_absence_table_sig.csv"
ids_path="/home/antonio/Escritorio/Antonio_G/Antonio_TFM/Resultados/09-Presence_absence_table/ids_table.csv"
species_name_path="/home/antonio/Escritorio/Antonio_G/Antonio_TFM/Additional_info/Metaanalysis_miRNA/10-miRNAsAnnotation/species_id.csv"
sRNAs_sig_path="/home/antonio/Escritorio/Antonio_G/Antonio_TFM/Resultados/Tablas/num_sRNAs_difexp.csv"
miRNAs_ident_path="/home/antonio/Escritorio/Antonio_G/Antonio_TFM/Resultados/Tablas/identificacion_miRNAs.csv"
output_directory_path="/home/antonio/Escritorio/10-Analysis_results"
## 1. Create output directories
# Create output directory paths
path_out_tables <- paste(output_directory_path, "01-Tables", sep="/")
path_out_plots <- paste(output_directory_path, "02-Plots", sep="/")

# Create output directories if they do not exist
if (!dir.exists(path_out_tables)) {
  dir.create(path_out_tables, recursive = TRUE)
}
if (!dir.exists(path_out_plots)) {
  dir.create(path_out_plots, recursive = TRUE)
}


## 2. Read files

pre_abs_table <- read.csv(pre_abs_table_path, sep = ",", header = TRUE, row.names = 1)
ids_table <- read.csv(ids_path, sep = ",", header = TRUE)
sp_names <- read.csv(species_name_path, sep = ",", header = FALSE)
sig_srnas_table <- read.csv(sRNAs_sig_path,sep=";")
miRNAs_ident_table <- read.csv(miRNAs_ident_path,sep=",")


## 3. Create table with information about files

# Merge ids_table and sp_names tables
info_table <- merge(ids_table, sp_names, by.x = "Species", by.y = "V1")

# Select and create columns of interest
colnames(info_table)[which(colnames(info_table) == "V3")] <- "Specie_name"
colnames(info_table)[which(colnames(info_table) == "V2")] <- "Specie_3l_id"
info_table$Specie_id <- sapply(strsplit(as.character(info_table$ID), "\\."), function(x) x[1])


## 4. Prepare presence-absence table

# Count the number of miRNA families per stress event
sum_cols <- colSums(pre_abs_table)

# Select those events with more than 5 annotated families
pre_abs_table <- pre_abs_table[, sum_cols >=5]


## 5. Classification of miRNAs acording to their stress response range
cat("Classifying miRNAs according to their range of stress response...", sep="\n")
miRNA_class_list <- miRNAClasification(pre_abs_table,
                                       save_table_dir = path_out_tables,
                                       save_plot_dir = path_out_plots)
cat("Done!", sep="\n")

## 6. Check species and stress bias
cat("Testing biases due to overrepresentation of species and stresses...", sep="\n")
sp_st_bias <- CheckSpeciesAndStressesBias(pre_abs_table,
                                          miRNA_class_list$table,
                                          save_table_dir = path_out_tables,
                                          save_plot_dir = path_out_plots)
cat("Done!", sep="\n")

## 7. Differentialy expressed sRNAs annotated as miRNAs
cat("Representing the number of differentially expressed sRNAs and how many of them are miRNAs...", sep="\n")
dea_sum <- DEASummaryPlot(sig_srnas_table,
                          miRNAs_ident_path,
                          save_plot_dir = path_out_plots)
cat("Done!", sep="\n")

## 7. Heatmap
cat("Creating Heatmap of stress-reactive miRNAs.", sep="\n")
heatmap <- miRNAsStressEventsHeatmap(pre_abs_table,
                                     save_plot_dir = path_out_plots)
cat("Done!", sep="\n")

## 8.Create Cytoscape tables
cat("Creating Cytoscape tables to construct the miRNA-mediated stress response network in plants.", sep="\n")
Cytotables <- CreateCytoscapeTables(pre_abs_table,
                                    miRNA_class_list$table,
                                    save_table_dir = path_out_tables)
cat("Done!", sep="\n")




