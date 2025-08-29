library(ggplot2)
library(readr)  # safer reader

plot_pvalue_hist <- function(file_path, p_threshold = 0.1, bin = 500) {
  # Read the TSV file safely as all characters
  data <- read_tsv(file_path,
                   col_types = cols(.default = "c"))  # read everything as character
  
  # Check if P column exists
  if (!"P" %in% colnames(data)) {
    stop("Column 'P' not found in the data.")
  }
  
  # Convert P column to numeric safely
  data$P <- as.numeric(data$P)
  
  # Remove NAs and filter by threshold
  data_filtered <- data[!is.na(data$P) & data$P <= p_threshold, ]
  
  if (nrow(data_filtered) == 0) {
    stop("No valid P-values found within the specified threshold.")
  }
  
  # Plot histogram
  ggplot(data_filtered, aes(x = P)) +
    geom_histogram(bins = bin, fill = "skyblue", color = "black") +
    theme_minimal() +
    labs(title = paste("Distribution of P-values (0 -", p_threshold, ")"),
         x = "P-value",
         y = "Frequency")
}

# Example usage:
plot_pvalue_hist("output_droso/quantitative_table_vcf.tsv")
