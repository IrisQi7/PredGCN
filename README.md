# PredGCN

![framework_v2 1](https://github.com/IrisQi7/PredGCN/assets/67795742/90bb006d-07e4-4b41-a694-8437977e8d28)



Pruning-enabled Gene-Cell Net(PredGCN) is an automatic framework for scRNA-seq data annotation. It incorporates a Coupled Gene-Cell Net (CGCN) to enable representation learning and information storage. PredGCN integrates a Gene Splicing Net (GSN) and a Cell Stratification Net (CSN), employing a pruning operation (PrO) to dynamically tackle the complexity of heterogeneous cell identification. Among them, GSN leverages multiple statistical and hypothesis-driven feature extraction methods to selectively assemble genes with specificity for scRNA-seq data while CSN unifies elements based on diverse region demarcation principles, exploiting the representations from GSN and precise identification from different regional homogeneity perspectives. Furthermore, we develop a multi-objective Pareto pruning operation(Pareto PrO) to expand the dynamic capabilities of CGCN, optimizing the sub-network structure for accurate cell type annotation.

# Installation Guide
       
Before starting the installation, we recommend starting by creating a local conda environment and installing the following packages in that environment in order. PredGCN mainly depends on the following Python packages:
      
python:

    python=3.7.12
    pycaret=2.3.10
    scanpy
    rpy2
    r-devtools
    
After installation of the python package, switch to the R environment and execute the following R commands:
        
R:

    library(devtools)
    install_github("suke18/FEAST")

### Usage

After completing the installation process, follow the commands below to use PredGCN:

    cd PredGCN
    mkdir MLmodels
    mkdir tmp
    python run.py

# Webserver

![image](https://github.com/IrisQi7/PredGCN/blob/master/webserver_screenshoot.png)

As PredGCN is developed using Python and R programming languages, the installation process may be complicated and time-consuming. For user convenience, a dedicated PredGCN webserver has been established. It is now publicly accessible via <a href="https://www.aibio-lab.com/PredGCN/index/">link</a>.

# Data Availability

To comprehensively evaluate PredGCN across various contexts, we curated several real scRNA-seq datasets spanning different species (human and mouse), organs (pancreas and PBMC), and sequencing platforms(including SMARTer, 10x Genomics, inDrop, CEL-seq2, and Smart-seq2).The data sets we used can be download in <a href="https://figshare.com/articles/dataset/scCPEP/22333150">datas</a>.

# License
This project is covered under the MIT License.
