# Use an official Python runtime as a parent image
FROM continuumio/miniconda3:4.9.2

# Set the working directory in the container to /app
WORKDIR /app

# Create and activate a new conda environment
RUN conda create -y -n leiden python=3.8

# Update pip
RUN pip install --upgrade pip

ENV PATH /opt/conda/envs/leiden/bin:$PATH

SHELL ["conda", "run", "-n", "leiden", "/bin/bash", "-c"]

RUN conda update -n base -c defaults conda
#RUN pip install --upgrade pip

# Install python packages

RUN pip install anndata
RUN pip install numpy
RUN pip install pandas
RUN pip install scipy
RUN pip install scikit-learn
RUN pip install matplotlib
RUN pip install scanpy
RUN pip install squidpy
RUN pip install leidenalg
RUN pip install jupyter

WORKDIR /app

# Make port 8888 available to the world outside this container
EXPOSE 8888

# Run Jupyter Notebook when the container launches
CMD ["/opt/conda/envs/leiden/bin/jupyter", "notebook", "--no-browser", "--ip=0.0.0.0", "--allow-root", "--NotebookApp.token=''", "--NotebookApp.password=''", "--notebook-dir=/app"]
