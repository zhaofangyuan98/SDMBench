# SDMBench
## Introduction
SDMBench is a python package that provides the metrics function to benchmark methods for indentifing spatial domain among Spatial transcriptome datasets.

## Installation
1. Clone the source code.
```
git clone https://github.com/zhaofangyuan98/SDMBench.git
cd SDMBench
```
2. Create a conda environment and activate it.
```
conda env create -n SDMBench --file SDMBench.yml
conda activate SDMBench
```
3. Install SDMBench as a dependency or third-party package with pip:
```
pip install .
```

After the SDMBench package is successfully installed, you can now analysis your data. [***The tutorial***](https://github.com/zhaofangyuan98/SDMBench/tree/main/Tutorial) serves as an illustrative example that demonstrates the process of benchmarking new methods against existing ones. You can run the jupyter notebook of tutorial.ipynb to reproduce it, and new_method.txt provides a format about your results file.
