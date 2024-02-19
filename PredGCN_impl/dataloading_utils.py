import pandas as pd
import numpy as np
import anndata
import os, sys
import scCPEP.method_utils as method_utils

def prepare_data(traindata, obs_train, testdata):
    df_train = traindata.T
    df_test = testdata.T

    obs_train.columns = ["cell.type"]
    obs_train.index = df_train.columns

    var_train = pd.DataFrame(data=df_train.index, index=df_train.index)
    var_train.columns = ['gene_symbols']
    var_test = pd.DataFrame(data=df_test.index, index=df_test.index)
    var_test.columns = ['gene_symbols']

    train_adata = anndata.AnnData(X=df_train.T, obs=obs_train, var=var_train)
    train_adata.obs_names_make_unique(join="-")
    train_adata.var_names_make_unique(join="-")
    test_adata = anndata.AnnData(X=df_test.T, var=var_test)
    test_adata.obs_names_make_unique(join="-")
    test_adata.var_names_make_unique(join="-")
    
    return train_adata, test_adata

def load_data(datapath, datalabel=None):
    df = pd.read_csv(datapath, index_col=0)
    df = df.T
    obs = pd.read_csv(datalabel, index_col=0)

    obs.columns = ["cell.type"]
    obs.index = df.columns
    var = pd.DataFrame(data=df.index, index=df.index)
    var.columns = ['gene_symbols']

    adata = anndata.AnnData(X=df.T, obs=obs, var=var)
    adata.obs_names_make_unique(join="-")
    adata.var_names_make_unique(join="-")

    return adata



def preprocessing_data(train_adata, test_adata, lognorm=True):

    train_adata = method_utils.process_adata(train_adata, lognorm)
    test_adata = method_utils.process_adata(test_adata, lognorm)

    common_feas = set(train_adata.var_names.tolist()).intersection(set(test_adata.var_names.tolist()))
    common_feas = list(common_feas)

    train_adata = train_adata[:, common_feas]
    test_adata = test_adata[:, common_feas]

    return train_adata, test_adata

    