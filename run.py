import argparse
import random
# from pathlib import Path
import numpy as np
from PredGCN_impl.predgcn import mainFunc


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--random_seed", type=int, default=21)
    parser.add_argument("--feature_selection", default=["DV", "DD", "chisq", "BI", "Ftest"])
    parser.add_argument("--rejection", type=float, default=0)
    parser.add_argument("--confidence_threshold", type=float, default=1)
    parser.add_argument("--lognorm", default=False)
    parser.add_argument("--OutputDir")
    parser.add_argument("--TrainDataPath")
    parser.add_argument("--TrainLabelsPath")
    parser.add_argument("--TestDataPath")
    parser.add_argument("--TestLabelsPath")
    parser.add_argument("--result_dir")
    parser.add_argument('--task_id')

    params = parser.parse_args()

    random.seed(params.random_seed)
    np.random.seed(params.random_seed)
 
    mainFunc(params)


    

     
    