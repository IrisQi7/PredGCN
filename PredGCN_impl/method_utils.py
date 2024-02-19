import os, sys
import numpy as np
import pandas as pd
import scipy
import scanpy as sc
import anndata

import rpy2.robjects as robjects
from rpy2.robjects.packages import importr

from sklearn import metrics


def process_adata(adata, min_genes=10, min_cells=10,
                  lognorm=True,
                  celltype_label="cell.type"):
    adata.var_names = [i.upper() for i in list(adata.var_names)]#avod some genes having lower letter

    adata.var_names_make_unique()
    adata.obs_names_make_unique()

    sc.pp.filter_cells(adata, min_genes=min_genes)

    sc.pp.filter_genes(adata, min_cells=min_cells)

    Gene1Pattern="ERCC"
    Gene2Pattern="MT-"
    id_tmp1 = np.asarray([not str(name).startswith(Gene1Pattern) for name in adata.var_names], dtype=bool)
    id_tmp2 = np.asarray([not str(name).startswith(Gene2Pattern) for name in adata.var_names], dtype=bool)
    id_tmp = np.logical_and(id_tmp1, id_tmp2)
    adata._inplace_subset_var(id_tmp)

    if adata.shape[0] < 3 or adata.shape[1] < 3:
        return None

    if lognorm == 1:
        sc.pp.normalize_per_cell(adata)
        sc.pp.log1p(adata)
    
    if celltype_label in adata.obs:
        cells = adata.obs.dropna(subset=[celltype_label]).index.tolist()
        adata = adata[cells]
    return adata


def generate_tmp(train_adata, result_dir):
    tmp_adata = train_adata.copy()
    if scipy.sparse.issparse(tmp_adata.X) or isinstance(tmp_adata.X, pd.DataFrame):
        tmp_data = tmp_adata.X.toarray()
    else:
        tmp_data = tmp_adata.X
    tmp_df = pd.DataFrame(data=tmp_data, index=tmp_adata.obs_names, columns=tmp_adata.var_names).T

    if not os.path.exists(result_dir):
        os.makedirs(result_dir)
        
    tmp_df_path = result_dir + os.sep + "tmp_counts.csv"
    
    tmp_df.to_csv(tmp_df_path)

    cell_annots = tmp_adata.obs["cell.type"].tolist()
    cell_annots_path = result_dir + os.sep + "tmp_cell_annots.txt"
    with open(cell_annots_path, 'w') as f:
        for cell_annot in cell_annots:
            f.write("%s\n" % cell_annot)
            f.flush()
    f.close()
    del tmp_adata
    return tmp_df_path, cell_annots_path

DV_PATH = "Rimpl/doDV.R"
DD_PATH = "Rimpl/doDD.R"
chisq_PATH = "Rimpl/doChisSquared.R"
BI_PATH = "Rimpl/doBI.R"
FTEST_PATH = "Rimpl/doFtest.R"


def feature_selection(method, train_adata, test_adata, result_dir, tmp_df_path, cell_annots_path, celltype_label="cell.type"):

    topN = 50
    pSig = 0.001
    gene_no = 2000

    if "Ftest" == method:
        os.system("Rscript --vanilla " + FTEST_PATH + " " + tmp_df_path + " " +
                  cell_annots_path + " " + str(gene_no))
    if "DV" == method:
        os.system("Rscript --vanilla " + DV_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)
    if "DD" == method:
        os.system("Rscript --vanilla " + DD_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)
    if "chisq" == method:
        os.system("Rscript --vanilla " + chisq_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)
    if "BI" == method:
        os.system("Rscript --vanilla " + BI_PATH + " " + tmp_df_path + " " +
                  cell_annots_path)

    if train_adata is None or test_adata is None:
        return None, None

    sub_feature_file = result_dir + os.sep + method + "_features.txt"
    with open(sub_feature_file) as f:
        features = f.read().splitlines()
    sub_features = [x.upper() for x in features]

    sub_genes = set(sub_features).intersection(set(train_adata.var_names.tolist()))
    train_adata = train_adata[:, list(sub_genes)]

    features = set(train_adata.var_names.tolist()).intersection(set(test_adata.var_names.tolist()))
    features = list(features)
    features.sort()
    print("Number of", method, "features:", len(features))

    with open(result_dir + os.sep + method + "_features.txt", 'w') as f:
        for feature in features:
            f.write("%s\n" % feature)
            f.flush()
    f.close()

    train_adata = train_adata[:, features]
    test_adata = test_adata[:, features]
    return train_adata, test_adata, features


