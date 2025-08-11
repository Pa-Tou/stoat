# Load necessary libraries
# install.packages("ggplot2")
# install.packages("dplyr")
# install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)

generate_boxplots <- function(phenotype_file, dir_path, output_path) {
  
  # Load phenotype data (FID, IID, PHENO)
  pheno_data <- read.table(phenotype_file, header = TRUE, sep = "\t")
  
  # List all files in the provided directory path (N-table files)
  n_table_files <- list.files(path = dir_path, full.names = TRUE)
  
  for (n_table_file in n_table_files) {
    
    base_name <- tools::file_path_sans_ext(basename(n_table_file))
    raw_header <- readLines(n_table_file, n = 1)
    header_cols <- strsplit(raw_header, "\t")[[1]]
    
    n_table <- read.table(n_table_file, header = FALSE, sep = "\t", skip = 1, stringsAsFactors = FALSE)
    colnames(n_table) <- header_cols
    
    # Rename first column to IID
    colnames(n_table)[1] <- "IID"
    
    merged_data <- n_table %>%
      left_join(pheno_data, by = "IID")
    
    genotype_columns <- setdiff(colnames(merged_data), c("IID", "PHENO"))
    
    # Function to classify genotype based on normalized allele values
    classify_genotype <- function(row) {
      values <- as.numeric(row[genotype_columns])
      if (all(is.na(values))) return(NA)
      top_two <- sort(values, decreasing = TRUE)[1:2]
      top_indices <- order(values, decreasing = TRUE)[1:2]
      
      tol <- 1e-6
      
      if (abs(top_two[1] - 1) < tol) {
        return(paste0(genotype_columns[top_indices[1]], "/", genotype_columns[top_indices[1]]))
      } else if (abs(top_two[1] - 0.5) < tol && abs(top_two[2] - 0.5) < tol) {
        return(paste0(genotype_columns[top_indices[1]], "/", genotype_columns[top_indices[2]]))
      } else {
        return("Other")
      }
    }
    
    merged_data$Genotype <- apply(merged_data, 1, classify_genotype)
    merged_data <- merged_data %>% filter(!is.na(PHENO) & !is.na(Genotype))
    
    genotype_counts <- merged_data %>%
      group_by(Genotype) %>%
      summarise(count = n(), .groups = "drop")
    
    merged_data <- merged_data %>%
      left_join(genotype_counts, by = "Genotype") %>%
      mutate(Genotype = paste0(gsub("/", "\n", Genotype), "\n(", count, ")"))
    
    # Plot
    p <- ggplot(merged_data, aes(x = Genotype, y = PHENO)) +
      geom_violin(fill = "cadetblue3", alpha = 0.3) +
      geom_boxplot(width = 0.2, outlier.size = 2, outlier.colour = "red", alpha = 0.5, fill = "darkcyan") +
      labs(x = "Genotype", y = "Phenotype", title = paste("Boxplot -", base_name)) +
      theme_minimal() +
      theme(
        axis.title.x = element_text(size = 14, color = "cadetblue4"),
        axis.title.y = element_text(size = 14, color = "cadetblue4"),
        plot.title = element_text(size = 16, color = "cadetblue4", face = "bold", hjust = 0.5)
      )

    output_file <- file.path(output_path, paste0(base_name, "_boxplot.jpeg"))
    ggsave(output_file, plot = p, device = "jpeg", width = 8, height = 6, dpi = 300)
    message("Saved plot for ", base_name, " as ", output_file)
  }
}

# Command line arguments
args <- commandArgs(trailingOnly = TRUE)
dir_path <- args[which(args == "-d") + 1]
phenotype_file <- args[which(args == "-p") + 1]
output_path <- args[which(args == "-o") + 1]

generate_boxplots(phenotype_file, dir_path, output_path)
print("Program finished")

# Rscript box_plox.R -d ../output/regression -p ../data/quantitative/phenotype.tsv -o ../output/plots
