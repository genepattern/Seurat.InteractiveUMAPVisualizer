import argparse
import os
import pandas as pd
from plotting import plot_umap


def main(args):
    

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

    figs = []
    columns_in = []
    genes_in = []
    with open('index.html', 'a') as file:
    
    
        for col in selected_columns:
            try:
                viz = plot_umap(clustered_file, col, "", color_scale)
    #             file.write("hello")
                file.write(f'<div id={col}>')
                file.write(viz)
                file.write(f"</div>")
                columns_in.append(col)
            except Exception as e:
                print(e)
                continue
        for gene in args.target_genes.split(','):
            try:
                viz = plot_umap(clustered_file, "", gene, color_scale)
                file.write(f'<div id={gene}>')
                file.write(viz)
                file.write(f"</div>")
                genes_in.append(gene)
            except Exception as e:
                print(e)
                continue
        file.write("<div> Select Columns  ")
        file.write("<select id='mySelect'> Select Column to Display")
        for col in columns_in:
            file.write(f"<option value='{col}'>{col}</option>")
        file.write("</select>")
        file.write("</div>")
        
        file.write("<div> Select Genes  ")
        file.write("<select id='mySelectGenes'>")
        file.write(f"<option value='-- SELECT A GENE --'>-- SELECT A GENE --</option>")
        for gene in genes_in:
            file.write(f"<option value='{gene}'>{gene}</option>")
        file.write("</select>")
        file.write("</div>")
    #     file.write(f'<script>const divElements_genes = {genes_in}; // array of ids for the divs to toggle</script>')
        file.write(f'<script>const divElements = {columns_in + genes_in}; // array of ids for the divs to toggle</script>')
        file.write(
            """
        <script>
        function display() {
        for (let i = 1; i < divElements.length; i++) {
        try{
            console.log(divElements[i]);
            const divElement = document.getElementById(divElements[i]);
            divElement.style.display = "none";
        }catch(err){
        console.log(err.message);
        }
        }
        }
        display()
        const selectElement = document.getElementById("mySelect");
        // listen for changes to the dropdown menu
        selectElement.addEventListener("change", () => {
        document.getElementById("mySelectGenes").selectedIndex = 0;
        try{
            const selectedOption = selectElement.options[selectElement.selectedIndex].value;
            // loop through all the divs and toggle their visibility
        for (let i = 0; i < divElements.length; i++) {
            const divElement = document.getElementById(divElements[i]);
            if (divElements[i] === selectedOption) {
            divElement.style.display = "block";
            } else {
            divElement.style.display = "none";
            }
        }
        
        }
        catch(err){
            console.log(err.message);
        }
        });
        const selectElement2 = document.getElementById("mySelectGenes");
        // listen for changes to the dropdown menu
        selectElement2.addEventListener("change", () => {
        
        try{
            const selectedOption = selectElement2.options[selectElement2.selectedIndex].value;
            // loop through all the divs and toggle their visibility
        for (let i = 0; i < divElements.length; i++) {
            const divElement = document.getElementById(divElements[i]);
            if (divElements[i] === selectedOption) {
            divElement.style.display = "block";
            } else {
            divElement.style.display = "none";
            }
        }
        }
        catch(err){
            console.log(err.message);
        }
        });
        </script>
        """)

    print('finished')



if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Convert Seurat RDS files to H5AD format and plot UMAP')
    parser.add_argument('--input', '-i', required=True, help='Input Seurat clustered RDS file')
    parser.add_argument('--target_genes', '-t', required=True, help='List of genes to visualize. If you write multiple genes they must be separated by a comma and a space (for example: MS4A1, CD8A, CD14).')
    parser.add_argument('--color_choice', '-c', required=False, help="Color Choice, either Continuous or Categorical", default="continuous", choices=['categorical', 'continuous'])
    args = parser.parse_args()

    main(args)




## python3 /opt/genepatt/main.py -i /data/seurat_preprocessed_dataset.clustered.rds
## docker build --platform linux/amd64 -t test . && docker run --rm -it -v ~/Downloads:/data test bash
    
