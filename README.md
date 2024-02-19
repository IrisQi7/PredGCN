# PredGCN

# Installation Guide
       
Before running the tool, we recommend creating an environment locally using ‘’‘conda’‘’ first, then you need to install the following packages in your environment, PredGCN mainly depends on the followinb Python packages:
      
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

# Data Availability

The data sets we used can be download in <a href="https://figshare.com/articles/dataset/scCPEP/22333150">datas</a>.

# License
This project is covered under the MIT License.
