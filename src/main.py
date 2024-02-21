import argparse
import os
import pandas as pd
from plotting import plot_umap


def main():
    parser = argparse.ArgumentParser(description='Convert Seurat RDS files to H5AD format and plot UMAP')
    parser.add_argument('--input', '-i', required=True, help='Input Seurat clustered RDS file')
    parser.add_argument('--target_genes', '-t', required=False, help='	List of genes to visualize. If you write multiple genes they must be separated by a comma and a space (for example: MS4A1, CD8A, CD14).')
    parser.add_argument('--color_choice', '-c', required=False, help="Color Choice, either Continuous or Categorical", default="continuous", choices=['categorical', 'continuous'])
    args = parser.parse_args()

    # Run rds_conversion function
    os.system(f"Rscript /opt/genepatt/rds_conversion.R --input {args.input}")

    

    directory = '/'  # replace with actual directory path
    found_file = None
    for root, dirs, files in os.walk(directory):
        for file in files:
            if '_dropdown_data.txt' in file:
                found_file = os.path.join(root, file)
                break

    clusterh5ad = None
    for root, dirs, files in os.walk(directory):
        for file in files:
            if '.h5ad' in file:
                clusterh5ad = os.path.join(root, file)
                break

        
    if found_file:
        print('Found file:', found_file)
        selected_columns = found_file
    else:
        raise Exception('Could not find file with "_dropdown_data.txt" in its name.')
    
    if clusterh5ad:
        print('Found file:', clusterh5ad)
        clustered_file = clusterh5ad
    else:
        raise Exception('Could not find file with ".h5ad" in its name.')
    

    with open(selected_columns, 'r') as f:
        lines = f.readlines()
    # Convert the list of strings into a list of integers
    selected_columns = [(line.strip()) for line in lines]

    color_scale = args.color_choice

    ## plotting: 
    with open('index.html', 'a') as file:
        if args.target_genes:
            print("Will plot " + args.target_genes)
            for gene in args.target_genes.split(','):
                try:
                    file.write(plot_umap(clustered_file, "", gene, color_scale))
                except Exception as e:
                    print(e)
                    continue

        
        for col in selected_columns:
            try:
                file.write(plot_umap(clustered_file, col, "", color_scale))
            except Exception as e:
                print(e)
                continue

    print('finished')



if __name__ == '__main__':
    main()




## python3 /opt/genepatt/main.py -i /data/seurat_preprocessed_dataset.clustered.rds
## docker build --platform linux/amd64 -t test . && docker run --rm -it -v ~/Downloads:/data test bash
    
