# PredGCN

<img width="525" alt="framework_screenshot" src="https://github.com/IrisQi7/PredGCN/assets/67795742/07400e1a-58fa-46e2-b582-0d445e8a4bda">

Pruning-enabled Gene-Cell Net(PredGCN) is an automatic framework for scRNA-seq data annotation. It incorporates a Coupled Gene-Cell Net (CGCN) to enable representation learning and information storage. PredGCN integrates a Gene Splicing Net (GSN) and a Cell Stratification Net (CSN), employing a pruning operation (PrO) to dynamically tackle the complexity of heterogeneous cell identification. Among them, GSN leverages multiple statistical and hypothesis-driven feature extraction methods to selectively assemble genes with specificity for scRNA-seq data while CSN unifies elements based on diverse region demarcation principles, exploiting the representations from GSN and precise identification from different regional homogeneity perspectives. Furthermore, we develop a multi-objective Pareto pruning operation(Pareto PrO) to expand the dynamic capabilities of CGCN, optimizing the sub-network structure for accurate cell type annotation.

# Installation Guide
       
Before running the tool, we recommend creating an environment locally using '''conda''' first, then you need to install the following packages in your environment, PredGCN mainly depends on the following Python packages:
      
python:

    python=3.7.12
    pycaret=2.3.10
    scanpy
    rpy2
    r-devtools
    
After installing the python package, switch to '''R''' environment, and install the following tool:
        
R:

    library(devtools)
    install_github("suke18/FEAST")

# Usage

    cd PredGCN
    mkdir MLmodels
    mkdir tmp
    python run.py

# Webserver

<img width="359" alt="webserver_screenshoot" src="https://github.com/IrisQi7/PredGCN/assets/67795742/6cfab30f-019a-4fea-b857-dbac5f6ffc0a">

As PredGCN is developed using Python and R programming languages, the installation process may be complicated and time-consuming. For user convenience, a dedicated PredGCN webserver has been established. It is now publicly accessible via <a href="https://www.aibio-lab.com/PredGCN/index/">link</a>.

# Data Availability

The data sets we used can be download in <a href="https://figshare.com/articles/dataset/scCPEP/22333150">datas</a>.

# License
This project is covered under the MIT License.
