import pandas as pd
import numpy as np
import os
import csv
import glob
import sys
from functools import partial
from joblib import Parallel, delayed

def decimal(x):
        y = '{:.12f}'.format(x)
        return y

def average(path, sim, stat, center, window):
        # average stats over all SNPs (everything except HAF and H12)
        file = f"{path}/stats/center_{center}/window_{window}/norm/{sim}_c{center}_w{window}_norm.{stat}"
        if not os.path.exists(file):
                print(f"no {stat} file at center {center} and window size {window} for {sim}")
        else:
                try:
                        values=pd.read_csv(file, sep=",", quotechar='"', skipinitialspace=True) # NA values read as NaN
                except pd.errors.EmptyDataError:
                        print(f"empty {stat} file at center {center} and window size {window} for {sim}, exiting")
                        values = pd.DataFrame([["NaN","NaN","NaN","NaN","NaN"]], columns=["freq_bins","position","stat","DAF","stat_norm"])
                        sys.exit(1)
                avgStat = values['stat_norm'].mean() # excludes NaN values
                if pd.isna(avgStat):
                        pass
                return avgStat

def writeFV(fV, path, type):
        # first write list of features to first line
        os.makedirs(f"{path}/fvs/", exist_ok=True)
        fV.to_csv(f"{path}/fvs/{type}.fv", mode='w', header=True, index=False, na_rep="nan")

def fvWrap(argsDict, simFile, type):
        sim = os.path.basename(simFile)
        iter = int(simFile.split("_")[-1])
        path = f"{argsDict['outputDir']}/training_data/{type}_data"
        statFiles = glob.glob(f"{argsDict['outputDir']}/training_data/{type}_data/stats/center_*/window_*/norm/{type}_data_{iter}_*_norm.*")
        fV = []
        features = []
        if os.path.exists(f"{argsDict['outputDir']}/training_data/{type}_data/array_{type}.txt") and os.path.getsize(f"{argsDict['outputDir']}/training_data/{type}_data/array_{type}.txt") > 0:
                params = pd.read_csv(f"{argsDict['outputDir']}/training_data/{type}_data/array_{type}.txt", sep="\s+", quotechar='"', skipinitialspace=True, header=0)                
                paramsIter = params.loc[params['#ArrayIndex'] == iter].values.tolist()[0] # get the row with index matching iter
                paramsIter[0] = int(paramsIter[0])
                if len(paramsIter) == 4:
                        featureFile = pd.DataFrame(np.empty((0, 1155+4))) # total number of features: 11*5*21 + 4 for params
                elif len(paramsIter) == 8:
                        featureFile = pd.DataFrame(np.empty((0, 1155+8))) # total number of features: 11*5*21 + 8 for params
                else:
                        print(f"Parameter array file should have 4 or 8 columns (#ArrayIndex, NE, MU, RHO), but {argsDict['outputDir']}/training_data/{type}_data/array_{type}.txt has len(paramsIter)")
                        sys.exit(1)

                fV.extend(paramsIter)
        else:
                print(f"No parameter array at {argsDict['outputDir']}/training_data/{type}_data/array_{type}.txt, creating {type} feature vector without this information.")
        for stat in ["ihs","iSAFE","nsl","DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                for center in range(argsDict['minCenter'], argsDict['maxCenter'] + argsDict['distCenters'], argsDict['distCenters']):
                        for window in [50000, 100000, 200000, 500000, 1000000]:
                                avgStat=average(path, sim, stat, center, window)
                                if not argsDict['keepStats']:
                                        os.remove(f"{path}/stats/center_{center}/window_{window}/norm/{sim}_c{center}_w{window}_norm.{stat}")
                                if not pd.isna(avgStat):
                                        avgStat=decimal(avgStat)
                                        if stat in ["HAF"]:
                                                fV.append(float(avgStat))
                                        else:
                                                fV.append(float(avgStat))
                                else:
                                        fV.append(float("nan"))
        featureFile.loc[iter]=fV
        return featureFile

def parallelFV(argsDict, simFiles, type):
        featureFile = Parallel(n_jobs=argsDict["numJobs"], verbose=100, backend="multiprocessing")(delayed(fvWrap)(argsDict, i, type) for i in simFiles)
        return pd.concat(featureFile)

def main(argsDict, simFiles, type):
        fullFeatureFile = parallelFV(argsDict, simFiles, type)
        features=[]
        if len(fullFeatureFile.columns) == 1155+4:
                params=["iter","NE","THETA","RHO"]
        elif len(fullFeatureFile.columns) == 1155+8:
                params=["iter","NE","SAF","SEL","TIME","EAF","THETA","RHO"]
        elif len(fullFeatureFile.columns == 1155):
                params=[]
        else:
                print(f"wrong number of parameters")
        features.extend(params)
        for stat in ["ihs","iSAFE","nsl","DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                for center in range(500000, 710000, 10000):
                        for window in [50000, 100000, 200000, 500000, 1000000]:
                                features.append(f"{stat}_w{window}_c{center}")
        fullFeatureFile.columns=features
        numNans = fullFeatureFile.isnull().sum(axis=1)
        numDropFV = sum(numNans > 115)
        fullFeatureFile.drop(fullFeatureFile[numNans > 115].index, inplace=True) # dump fvs with more than 10% nans
        print(f"Removing {numDropFV} feature vectors because more than 10% of statistics are nan")
        fullFeatureFile=fullFeatureFile.fillna(0) # replace remaining nans with 0s
        writeFV(fullFeatureFile, argsDict['outputDir'], type)


if __name__=="__main__":
        main()
