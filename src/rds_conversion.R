# Load required libraries
library(Seurat)
library(SeuratDisk, lib.loc="/home/jovyan/R")
library(tools, lib.loc="/home/jovyan/R")

# Define the function with the argparse argument
rds_conversion <- function(rdsFile) {
  # Parse the command line arguments
  args <- parse_args(rdsFile)
  
  # Print the help message if --help was passed
  if (args$help) {
    cat("Usage: rds_conversion [options] <input.rds>\n")
    cat("Options:\n")
    cat("  -h, --help                Show this help message and exit\n")
    cat("  <input.rds>            Input RDS file\n")
    return
  }
  
  # Check if the input file exists
  if (!file.exists(args$input)) {
    stop("Input file does not exist.")
  }
  
  # Run the rest of the script
  rds_conversion_main(args$input)
}

# Define the main function
rds_conversion_main <- function(input) {
  # Load the Seurat Clustered object
  seurat_clustered_dataset.object <- readRDS(file = input)
  
  # Create the output file names
  output_file_name <- file.path(file.dirname(input), "dropdown_data.txt")
  h5seurat_name <- file.path(file.dirname(input), paste0(basename(input), ".h5Seurat"))
  
  # Perform the RDS -> H5AD conversion
  suppressMessages(suppressWarnings({
    SaveH5Seurat(seurat_clustered_dataset.object, filename = h5seurat_name, overwrite = TRUE)
    Convert(h5seurat_name, dest = "h5ad", overwrite = TRUE)
  }))
  
  # Extract the metadata names
  meta_data_names <- as.list(colnames(seurat_clustered_dataset.object@meta.data))
  
  # Find the indices of the SeuratClustering iteration numbers
  matching_indices <- grep("snn", meta_data_names)
  snn_strings <- meta_data_names[matching_indices]
  
  # Remove the SeuratClustering iteration numbers from the list
  meta_data_names <- meta_data_names[!(meta_data_names %in% c("seurat_clusters", snn_strings))]
  
  # Order the remaining metadata names
  meta_data_names_ordered <- c("seurat_clusters", snn_strings, meta_data_names)
  
  # Write the ordered metadata names to a text file
  writeLines(unlist(meta_data_names_ordered), output_file_name)
  
  # Print a success message
  print(paste0(h5seurat_name, ' and ', output_file_name, ' successfully written to the same folder as this notebook!'))
}

# Make the script executable
if (interactive()) {
  rds_conversion()
} else {
  # Run the script with the input file as a command line argument
  rds_conversion(commandArgs(trailingOnly = TRUE)[1])
}

