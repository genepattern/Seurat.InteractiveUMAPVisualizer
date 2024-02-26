# Import packages
import anndata as ad
import pandas as pd
import numpy as np
import plotly.express as px
import os
import warnings
import sys
import argparse

# Prevent printing any errors to console
warnings.filterwarnings('ignore')

# Get file names in the current working directory
curr_dir = os.getcwd()
in_directory = os.listdir(curr_dir)
# Filter out only the files (excluding directories)
files = [file for file in in_directory if os.path.isfile(os.path.join(curr_dir, file))]

content_list = []

# Iterate through the list of files to find text file containing metadata columns, used to populate the corresponding dropdown
for file_name in files:
    if file_name.endswith("_dropdown_data.txt"):
        file_path = os.path.join(os.getcwd(), file_name)
        try:
            with open(file_path, 'r') as file:
                # Read in file contents, line by line
                content_list = file.read().splitlines()
        except FileNotFoundError:
            print(f"Missing a file required to run this cell: {file_name}")

color_scale_options = ["categorical", "continuous"]
file_options = [file for file in os.listdir(os.getcwd()) if file.endswith(".h5ad")]

def determine_coloring(umap_df, target_var, color_scale):
    """ Determine categorical vs continuous coloring, in case user's choice needs to be overridden"""
    # if Categorical coloring & large # of unique numerical values (e.g. >= 40), visualizing is not advisable
    if np.unique(umap_df[target_var]).size >= 40 and color_scale == "categorical":
        try:
            umap_df[target_var] = umap_df[target_var].apply(pd.to_numeric)
        except ValueError:  # if values are non-numerical, plot somewhat confusing categorical graph
            raise Exception(f"WARNING: The selected column has a large number of unique non-numeric entries. Plotting will be done with a continuous color scale.")
        else:
            color_scale = "continuous"
            raise Exception(f"WARNING: The selected column has a large number of unique numeric entries. Color scaling was changed from categorical to continuous for better interpretability.")
    return umap_df



def plot_umap(clustered_file, selected_column, target_gene, color_scale):
    """Show interactive UMAP plot, colored to represent either one of Seurat metadata columns or the expression level of a specified gene

    Parameters:
    clustered_file (str): Name of H5AD file containing AnnData object, which has gene expression data and UMAP loadings + clustering data from Seurat.Clustering
    selected_column (str): Selected metadata column (loaded into UI builder cell using "_dropdown_data.txt" file)
    target_gene(str): Gene for which expression data is to be shown
    color_scale(str): Either categorical or continuous

    Returns:
    None, displays interactive UMAP plot

   """
    # Check for file format compatibility
    if clustered_file.endswith("h5ad") == False:
        raise Exception("<p class='alert alert-warning'>ERROR: Selected file must be of type: h5ad</p>")
    
    is_gene = False # boolean to store whether UMAP plot is for visualizing metadata or gene expression level
    target_var = '' # var to store target variable name (gene or col)
    
    # Either a metadata column or target gene is to be selected
    if str(target_gene) == '' and str(selected_column) == '': # no parameter selected
        raise Exception("<p class='alert alert-warning'>ERROR: Select an input parameter (either a metadata column or a gene name)</p>")
    elif str(target_gene) != '' and str(selected_column) != '': # two parameters selected
        raise Exception("<p class='alert alert-warning'>ERROR: Select either a metadata column or a gene name as an input parameter</p>")
    elif str(selected_column) == '': # if metadata column not selected -> gene expression plot
        is_gene = True
        target_var = target_gene
    else: # metadata column selected
        target_var = selected_column
    
    # Identify cluster-colored plots for later step
    cluster_coloring = False
    if "RNA_snn_res" in selected_column or selected_column == "seurat_clusters":
        cluster_coloring = True
    
    # Read in clustered H5AD object
    adata = ad.read_h5ad(clustered_file)

    # Find UMAP in annotation of observations (allowing for some deviation in naming) -- as a precaution in case of change in naming convention
    axis_array_keys = list(adata.obsm)
    umap_key = ''
    for key in axis_array_keys:
        if 'umap' in key.lower():
            umap_key = key
    
    # Access the UMAP coordinates and Seurat clusters from the AnnData object
    umap_coords = adata.obsm[umap_key]
    cluster_assignments = adata.obs["seurat_clusters"]
    
    unique_clusters = np.unique(cluster_assignments)
    cluster_cell_counts = {cluster: np.sum(cluster_assignments == cluster) for cluster in unique_clusters}
    
    # Create a DataFrame with UMAP coordinates and gene expression values
    umap_df = pd.DataFrame(umap_coords, columns=['UMAP1', 'UMAP2'])
    
    # For UMAP plots colored by gene expression:
    if is_gene:
        try: # Access gene expression data, assuming storage in .X attribute
            ## reshaped to length 
            gene_expression = adata[:, target_var].X.reshape(len(adata))
        except KeyError: # if gene cannot be identified in gene expression matrix
            raise KeyError("<p class='alert alert-warning'>ERROR: The gene name was not found.</p>")
        
        # Identify the clusters and the number of cells ine ach cluster
        clusters = [cluster_assignments.iloc[entry_index] for entry_index in range(len(umap_coords))]
        cells_in_cluster = [cluster_cell_counts[cluster] for cluster in clusters]

        # Add gene expression data to dataframe
        umap_df[target_var] = gene_expression
        
        # Override color scheme if necessary
        umap_df = determine_coloring(umap_df, target_var, color_scale)
        
        # Create a scatter plot colored by the gene expression values
        fig = px.scatter(umap_df, x='UMAP1', y='UMAP2', color=target_var,
                        hover_data={'UMAP1': False, 'UMAP2': False,
                                    'UMAP cluster': clusters,
                                     'Number of cells in cluster': cells_in_cluster
                                    })
        fig.update_layout(title=f"UMAP Scatterplot Colored by Expression of {target_var} Gene", title_x = 0.5)
    
    # UMAP plot colored by metadata column
    else:
        # Access specified metadata column
        selected_col = adata.obs[target_var] 
        
        # Loop through to map together cluster assignments and UMAP coordinates through UMAP indices
        selected_col_looped = {}
        cluster_assignments_looped = {}
        umap_indices = []
        for entry_index in range(len(umap_coords)):
            umap_coordinates_entry = tuple(umap_coords[entry_index])
            selected_col_looped[umap_coordinates_entry] = selected_col.iloc[entry_index]
            cluster_assignments_looped[umap_coordinates_entry] = cluster_assignments.iloc[entry_index]
            umap_indices.append(umap_coordinates_entry)
    
        # Add information to dataframe
        umap_df["indices"] = umap_indices
        umap_df["seurat_clusters"] = umap_df["indices"].map(cluster_assignments_looped)
        umap_df["cells_in_cluster"] = umap_df["seurat_clusters"].map(cluster_cell_counts)
        if target_var not in umap_df:
            umap_df[target_var] = umap_df["indices"].map(selected_col_looped)
  
        # Override color scheme if necessary
        umap_df = determine_coloring(umap_df, target_var, color_scale)
    
        # Sort UMAP in order of selected column (in order for legend to reflect order)
        umap_df = umap_df.sort_values(by=[target_var], ignore_index=True)
        
        # Change data types for plotly to plot continuous/categorical color scheme accordingly
        if color_scale == "categorical":
            umap_df[target_var] = umap_df[target_var].apply(str)
        elif color_scale == "continuous":
            umap_df[target_var] = umap_df[target_var].apply(pd.to_numeric)
        
        # If cluster coloring, show # of cells per cluster; if not, then omit from hover data
        if cluster_coloring:
            fig = px.scatter(umap_df, x='UMAP1', y='UMAP2', color=target_var,
                    hover_data={'UMAP1': False, 'UMAP2': False,
                                selected_column:True, "seurat_clusters": True, "cells_in_cluster":True
                                })
        else:
            fig = px.scatter(umap_df, x='UMAP1', y='UMAP2', color=target_var,
                    hover_data={'UMAP1': False, 'UMAP2': False,
                                selected_column:True, "seurat_clusters": True, "cells_in_cluster":False
                                })
        
        fig.update_layout(title=f"UMAP Scatterplot ({target_var} clustering)", title_x = 0.5)

    # Customize to adjust size, adjust margins for axis labels
    fig.update_layout(width=750, height=750, margin=dict(l=65, r=65, t=100, b=65))
    # Customize axes to hide tick marks, maintain scale ratio
    fig.update_xaxes(tickfont=dict(color="rgba(0,0,0,0)"))
    fig.update_yaxes(tickfont=dict(color="rgba(0,0,0,0)"),scaleanchor="x",scaleratio=1)
    # fig.show()
    return fig.to_html(full_html=False, include_plotlyjs='cdn')
    # return fig








