#!/bin/python3.6

import argparse
import os
import pandas as pd
import numpy as np
import csv
import glob
import sys
from functools import partial

### note that sim file iterations should be designated "<simfilename>_<iter>.<suffix>"

def decimal(x):
        y = '{:.12f}'.format(x)
        return y

def average(windowDir, stat, center, window):
        # average stats over all SNPs (everything except HAF and H12)
        classWindow = os.path.basename(windowDir)
        normFile = f"{windowDir}/stats/center_{center}/window_{window}/norm/{classWindow}_c{center}_w{window}_norm.{stat}"
        try:
                values = pd.read_csv(normFile, sep=",", quotechar='"', skipinitialspace=True) # NA values read as NaN
        except (FileNotFoundError,pd.errors.EmptyDataError):
                values=pd.DataFrame([["NaN","NaN","NaN","NaN","NaN"]],columns=["freq_bins","position","stat","DAF","stat_norm"])
        avg_stat=values['stat_norm'].mean() # excludes NaN values
        return avg_stat

def writeFV(fV, outDir):
        # first write list of features to first line
        fV.to_csv(f"{outDir}.fv", index=False)

def fvWrap(argsDict, windowDir):
        classWindow = os.path.basename(windowDir)
        #statFiles = glob.glob(f"{windowDir}/stats/center_*/window_*/norm/{classWindow}_*_norm.*")
        
        features=[]
        feature_vec=[]
        for stat in ["ihs","iSAFE","nsl","DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                for center in range(argsDict['minCenter'], argsDict['maxCenter'] + argsDict['distCenters'], argsDict['distCenters']):
                        for window in [50000, 100000, 200000, 500000, 1000000]:
                                avg_stat=average(windowDir, stat, center, window)
                                if not pd.isna(avg_stat):
                                        avg_stat=decimal(avg_stat)
                                        feature_vec.append(float(avg_stat))
                                else:
                                        feature_vec.append(float("nan"))
        return feature_vec                

def main(argsDict, windowDir):
        featureVec = fvWrap(argsDict, windowDir)
        features=[]
        for stat in ["ihs","iSAFE","nsl","DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                for center in range(argsDict['minCenter'], argsDict['maxCenter'] + argsDict['distCenters'], argsDict['distCenters']):
                        for window in [50000, 100000, 200000, 500000, 1000000]:
                                features.append(f"{stat}_w{window}_c{center}")
        featureDF = pd.DataFrame([featureVec], columns=features)
        numNans = featureDF.isnull().sum(axis=1)
        if numNans[0] <= 115: # fewer than 10% of stats are nans
                featureDF = featureDF.fillna(0)
                windowName = os.path.basename(windowDir)
                outDir = f"{argsDict['outputDir']}/classification/fvs/{windowName}"
                writeFV(featureDF, outDir)
        else:
                print(f"Not enough statistics were calculated for {windowDir}, could not create feature vector")
if __name__=="__main__":
        main()
