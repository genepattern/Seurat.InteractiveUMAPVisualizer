#!/usr/bin/env R
# Rscript /opt/genepatt/rds_conversion.R --input /data/seurat_preprocessed_dataset.clustered.rds
# Load required libraries
library(argparse)
library(Seurat)
library(SeuratDisk)
library(tools)

# Parse command line arguments
parser <- argparse::ArgumentParser(description = "Convert Seurat RDS files to H5AD format")
parser$add_argument("--input", "-i", required = TRUE, help = "Input Seurat RDS file")
args <- parser$parse_args()

rds_conversion <- function(rdsFile) {
    print("converting now ...")
    file_name <- rdsFile
    
    base_name <- basename(file_name)
    # Get name of clustered RDS file without extension, used to name the dropdown_data.txt file
    file_name_without_extension <- file_path_sans_ext(base_name)
    
    # Read in Seurat Clustered object
    seurat_clustered_dataset.object <- readRDS(file = file_name)
    
    # Create file name for h5Seurat object
    h5seurat_name <-paste0(file_name_without_extension, ".h5Seurat")
    
    # Quick fix to prevent missing data; this apparently sometimes happens as a result from changes to metadata columns in the conversion
    # (source: https://github.com/mojaveazure/seurat-disk/issues/23)
    i <- sapply(seurat_clustered_dataset.object@meta.data, is.factor)
    seurat_clustered_dataset.object@meta.data[i] <- lapply(seurat_clustered_dataset.object@meta.data[i], as.character)
    
    # RDS -> H5AD conversion
    suppressMessages(suppressWarnings({SaveH5Seurat(seurat_clustered_dataset.object, filename = h5seurat_name, overwrite = TRUE)
        Convert(h5seurat_name, dest = "h5ad", overwrite = TRUE)
    }))
    # Extract list of metadata names from Seurat Clustered object
    meta_data_names <- as.list(colnames(seurat_clustered_dataset.object@meta.data))

    # Use grep to find indices of strings matching SeuratClustering iterations
    matching_indices <- grep("snn", meta_data_names)
    # Extract those matching strings from the list
    snn_strings <- meta_data_names[matching_indices]
    # Move those strings to beginning of list (by removing then adding at beginning of list)
    meta_data_names <- meta_data_names[!(meta_data_names %in% c("seurat_clusters", snn_strings))]
    meta_data_names_ordered <- c("seurat_clusters", snn_strings, meta_data_names)
    # Specify the name of the text file containing metadata 
    output_file <- paste0(file_name_without_extension, "_dropdown_data.txt")

    # Convert the list to character vector and write to the file
    writeLines(unlist(meta_data_names_ordered), output_file)
    
    # Ideally this message would be displayed via UIOutput, but I had trouble getting that to work, so just printing for now
    print(paste0(h5seurat_name, ' and ', output_file, ' successfully written to the same folder as this notebook!'))
    
    
    #gp_output <- gpUIOutput(id="output_id", title="Output Title")
    return
}


rds_conversion(args$input)