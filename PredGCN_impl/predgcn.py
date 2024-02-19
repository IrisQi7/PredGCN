import pandas as pd
import PredGCN_impl.dataloading_utils as dataloading_utils
import PredGCN_impl.method_utils as method_utils
import PredGCN_impl.pareto_ensemble_pruning as pareto_ensemble_pruning
import numpy as np
import scipy
import gc
import os
import math
import rpy2.robjects as robjects
import time
import scanpy as sc
import zipfile

from sklearn.preprocessing import OneHotEncoder
from pycaret.classification import *

from sklearn.naive_bayes import MultinomialNB
from sklearn.svm import LinearSVC
from sklearn.calibration import CalibratedClassifierCV

from sklearn.metrics import f1_score
from sklearn.metrics import accuracy_score

import anndata


def mainFunc(args):
    traindata = pd.read_csv(args.TrainDataPath, index_col=0, sep=',')
    train_labels = pd.read_csv(args.TrainLabelsPath, header=0, index_col=None, sep=',')
    testdata = pd.read_csv(args.TestDataPath, index_col=0, sep=',')

    train_adata, test_adata = dataloading_utils.prepare_data(traindata, train_labels, testdata)
    train_adata, test_adata = dataloading_utils.preprocessing_data(train_adata, test_adata, lognorm=args.lognorm)
    tmp_df_path, cell_annots_path = method_utils.generate_tmp(train_adata, str(args.result_dir))

    celltype_cols = "cell.type"
    enc_train = OneHotEncoder(handle_unknown='ignore')
    y_train = enc_train.fit_transform(train_adata.obs[[celltype_cols]]).toarray()

    trueClass = y_train.argmax(1)
    num_celltype = len(np.unique(y_train.argmax(1)))
    model_count = 0

    model_num = 6

    fs_count = 0
    for fs in args.feature_selection:
        fs_count = fs_count + 1
        train_sub, test_sub, features = method_utils.feature_selection(fs, train_adata, test_adata,
                                                                       str(args.result_dir), tmp_df_path,
                                                                       cell_annots_path)
        trainXY = train_sub.to_df()
        trainXY['cell.type'] = y_train.argmax(1)
        model_setup = setup(data=trainXY, target='cell.type', preprocess=False, silent=True,
                            session_id=args.random_seed)

        best_models = compare_models(n_select=model_num, include=[MultinomialNB(alpha=0.01),
                                                                  CalibratedClassifierCV(LinearSVC()),
                                                                  'mlp',
                                                                  'knn',
                                                                  'rbfsvm',
                                                                  'rf'])

        save_path = os.path.join(args.result_dir, 'MLmodels')
        if not os.path.exists(save_path):
            os.makedirs(save_path)
        os.chdir(save_path)
        for i in range(min(model_num, len(best_models))):

            model_count = model_count + 1

            if 1 == model_count:
                test_feas = pd.DataFrame({str(model_count): features})
            else:
                test_feas = pd.concat([test_feas, pd.DataFrame({str(model_count): features})], axis=1)  # column

            model = best_models[i]
            pred = predict_model(model, data=trainXY)
            predClass = pred[['Label']]

            predScore = pred[['Score']].to_numpy()
            if 0 == i and 1 == fs_count:
                predictClass = predClass
            else:
                predictClass = np.append(predictClass, predClass, axis=1)

            save_model(model, 'model_' + str(model_count))
        gc.collect()

    ensembleResults = pareto_ensemble_pruning.startTrain(predictClass, trueClass)

    start = time.time()
    ind = np.where(ensembleResults == 1)[0]  # need to add 1(start from 0)
    for test in range(len(ind)):
        model_ind = ind[test] + 1  # selected models
        # load model
        model = load_model('model_' + str(model_ind))
        fea = test_feas.iloc[:, ind[test]].tolist()
        fea = [elem for elem in fea if elem == elem]
        testX = test_adata[:, fea].to_df()
        pred = predict_model(model, data=testX)
        predLabel = pred[['Label']]
        testScore = pred[['Score']].to_numpy()
        unknown = np.where(testScore < args.rejection)
        predLabel.iloc[unknown[0].tolist()] = num_celltype
        if 0 == test:
            testPred = predLabel
            tScore = pred[['Score']]
        else:
            testPred = np.append(testPred, predLabel, axis=1)
            tScore = np.append(tScore, pred[['Score']], axis=1)

    test_index = test_adata.obs_names
    testP = pd.DataFrame(testPred, index=test_index)
    predResult = pd.DataFrame.mode(testP, axis=1)[0].astype('int64')
    cellTypes_train = enc_train.categories_[0]
    with_unknown = np.append(cellTypes_train, 'unknown')
    predResult_type = with_unknown[predResult.astype(np.int_)]
    end = time.time()

    predLabs = pd.DataFrame(predResult_type, columns=['pred'])
    testingTime = str(end - start)
    testTimeFile = str(args.result_dir) + os.sep + " testTimeFile.txt"
    with open(testTimeFile, 'a') as f:
        f.write("%s\n" % testingTime)
        f.flush()
    f.close()

    test_rows = testdata.index.tolist()
    preLabs_row = pd.concat([pd.DataFrame(test_rows, columns=['cell']), predLabs], axis=1)

    tsScore = pd.DataFrame(tScore)
    predLabs.to_csv(str(args.OutputDir) + os.sep + "predLabels_" + str(args.rejection) + ".csv", index=False)
    ensembleResults.to_csv(str(args.OutputDir) + os.sep + "ensembleResult_" + str(args.rejection) + ".csv", index=False)
    preLabs_row.to_csv(str(args.OutputDir) + os.sep + "predLabels_row_" + str(args.rejection) + ".csv", index=False)
    tsScore.to_csv(str(args.OutputDir) + os.sep + "testScore_" + str(args.rejection) + ".csv", index=False)

    df_testPred = pd.DataFrame(testPred)
    df_testPred.to_csv(str(args.OutputDir) + os.sep + "testPred_" + str(args.rejection) + ".csv", index=False)
    
    df_ct = pd.DataFrame(with_unknown)
    df_ct.to_csv(str(args.OutputDir) + os.sep + "celltypeFile.csv", index=False)

    predictC = pd.DataFrame(predictClass)
    predictC.to_csv(str(args.OutputDir) + os.sep + "predClass_" + str(args.rejection) + ".csv", index=False)
