FROM satijalab/seurat:5.0.0

USER root
RUN apt-get update && \
    apt-get install -y python3-pip && \
    apt-get install -y  python3-dev && \ 
    apt-get install -y python3-rpy2

# # Install other dependencies
RUN pip install pandas numpy matplotlib anndata plotly



# # Install remotes


# # Install tools


# # Load remotes
RUN R -e "install.packages(c('argparse'), repos='https://cran.rstudio.com')"
# RUN R -e "install.packages(c('ggplot2', 'dplyr', 'tidyr'), repos='https://cran.rstudio.com')"
# # Install SeuratDisk from GitHub
# RUN mkdir -p /home/joyvan/R
# RUN chmod 777 /home/joyvan/R
# RUN sudo R -e "remotes::install_github('mojaveazure/seurat-disk', lib='/home/jovyan/R')"


RUN mkdir -p /opt/genepatt

COPY src/* /opt/genepatt/
RUN chmod 777 /opt/genepatt
# Run the command to start the app
#CMD ["python", "app.py"]
