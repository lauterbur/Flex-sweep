#!/bin/python3

# this script normalizes the data from the sweep simulations
# using itself

#1. read in all derived allele frequencies and <stat> scores into a pandas table
#2. bin sites by 2% allele frequency
#3. calculate mean and SD for each bin
#4. normalize for each bin ( (value-mean)/SD )

import os
import pandas as pd
import numpy as np
import sys

def readWindow(argsDict, stat, simName, simPath, center):
        # make a pandas table with columns of:
        ## position
        ## statistic value
        ## derived allele frequency
        values=pd.DataFrame()
        stats=pd.DataFrame()
        if stat in ["DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                center = int(argsDict['locusLength']/2) # these statistics always have window calculated from the locus center
                file=f"{simPath}/stats/{simName}_c{center}.{stat}"
                if not os.path.exists(file):
                        print(f"Normalizing: Can't read {file}, does not exist")
#                        sys.exit(1)
                else:
                        if stat == "DIND":
                                values=pd.read_csv(file, sep=' ', names=["position","DIND","DAF"], skiprows=1, engine='python') # read file, these have headers to replace
                                values=values[['position','DIND','DAF']]
                                values=values.rename(columns={'DIND':'stat'})
                        if stat == "hDo":
                                values=pd.read_csv(file, sep=' ', names=["position","hDo","DAF"], skiprows=1, engine='python')
                                values=values[['position','hDo','DAF']] 
                                values=values.rename(columns={'hDo':'stat'})
                        if stat == "hDs":
                                values=pd.read_csv(file, sep=' ', names=["position","hDs","DAF"], skiprows=1, engine='python')
                                values=values[['position','hDs','DAF']] 
                                values=values.rename(columns={'hDs':'stat'})
                        if stat == "hf":
                                values=pd.read_csv(file, sep=' ', names=["position","hf","DAF"], skiprows=1, engine='python')
                                values=values[['position','hf','DAF']] 
                                values=values.rename(columns={'hf':'stat'})
                        if stat == "lf":
                                values=pd.read_csv(file, sep=' ', names=["position","lf","DAF"], skiprows=1, engine='python')
                                values=values[['position','lf','DAF']] 
                                values=values.rename(columns={'lf':'stat'})
                        if stat == "S":
                                values=pd.read_csv(file, sep=' ', names=["position","S","DAF"], skiprows=1, engine='python')
                                values=values[['position','S','DAF']] 
                                values=values.rename(columns={'S':'stat'})
                        if stat == "HAF":
                                values=pd.read_csv(file, sep=' ', names=["position","HAF"], skiprows=0, engine='python')
                                values=values[['position','HAF']] 
                                values['DAF']=1 # because there needs to be a value for DAF to get bin
                                values=values.rename(columns={'HAF':'stat'})
                        if stat == "H12":
                                values=pd.read_csv(file, sep='\s+', names=["position","H12"], skiprows=0, engine='python')
                                values=values[['position','H12']] 
                                values['DAF']=1 # because there needs to be a value for DAF to get bin
                                values=values.rename(columns={'H12':'stat'})
        elif stat in ["ihs","iSAFE","nsl"]:
                file=f"{simPath}/stats/center_{center}/{simName}_c{center}.{stat}"
                if not os.path.exists(file):
                        print(f"Normalizing: Can't read {file}, does not exist")
#                        sys.exit(1)
                else:
                        if stat == "ihs":
                                values=pd.read_csv(file, sep='\s', names=["position","index","DAF","iHH_0","iHH_1","iHS","standardized_iHS"], skiprows=1, engine='python') # read ihs file, these have headers to replace
                                values=values[['position','iHS','DAF']]
                                values=values.rename(columns={'iHS':'stat'})
                        if stat == "iSAFE":
                                values=pd.read_csv(file, sep='\t', names=["position","iSAFE","DAF"], skiprows=1, engine='python')
                                values=values[['position','iSAFE','DAF']]        
                                values=values.rename(columns={'iSAFE':'stat'})
                        if stat == "nsl":
                                values=pd.read_csv(file, sep='\s', names=["position","DAF","nSL"], engine='python') # no header
                                values=values[['position','nSL','DAF']] 
                                values=values.rename(columns={'nSL':'stat'})
        return values

def bin(values, stat, center, freq=.02):
        # bin sites by x% allele frequency
        try:
                values['freq_bins'] = pd.cut(x=values['DAF'], bins=np.arange(0,1+freq,freq), include_lowest=True, precision=2) #, labels=np.arange(0,1,freq))
        except:
                print(f"Normalization bin values for {stat} at center {center} do not exist")
                sys.exit(1)
        return values

def cut(values, center, window):
        # cut window stats to only SNPs within the window around center
        lower=center-window/2
        upper=center+window/2
        cut_values=values[(values['position'] >= lower) & (values['position'] <= upper)]
        return cut_values

def normalizeWindow(argsDict, binsDir, values, path, runname, stat, center, window):
        # get expected value (mean) and standard deviation from neutral data
        neutral=pd.read_csv(f"{binsDir}/neutral_data_bins.{stat}",index_col='freq_bins')
       # now normalize using those bins values
        values.index=values['freq_bins']
        values['freq_bins']=values['freq_bins'].astype(str)
        expected_col=values['freq_bins'].map(neutral['expected'])
        stdev_col=values['freq_bins'].map(neutral['stdev'])
        try:
                stdev_col.iloc[pd.Index(stdev_col).get_loc(0)]=1 # when stdev is 0, means no variation, so scaled value -> 0
        except KeyError:
                pass
        if stat in ["HAF","H12"]:
                values['stat_norm']=values['stat'] # don't normalize HAF and H12
        else:
                values['stat_norm']=(values['stat']-expected_col)/stdev_col
        values=values.drop(['freq_bins'],axis=1)
        if stat not in ["HAF","H12"]:
                cut_values=cut(values, center, window)
                if os.path.exists(f"{path}/stats/center_{center}/window_{window}/norm/"):
                        cut_values.to_csv(f"{path}/stats/center_{center}/window_{window}/norm/{runname}_c{center}_w{window}_norm.{stat}", na_rep="NA")
                else:
                        print(f"{path}/stats/center_{center}/window_{window}/norm/ doesn't exist?")
        elif stat in ["HAF","H12"]:
                cut_values=values
                cut_values.to_csv(f"{path}/stats/center_{center}/window_{window}/norm/{runname}_c{center}_w{window}_norm.{stat}", na_rep="NA")
        return cut_values

def main(argsDict, binsDir, simFile):
        simName = os.path.splitext(os.path.basename(simFile))[0]
        simPath = os.path.dirname(simFile)
        for center in range(argsDict['minCenter'], argsDict['maxCenter'] + argsDict['distCenters'], argsDict['distCenters']):
                for stat in ["ihs","iSAFE","nsl", "DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                        values = readWindow(argsDict, stat, simName, simPath, center)
                        if values is not None:
                                binnedValues=bin(values,stat,center)
                                for window in [50000, 100000, 200000, 500000, 1000000]:
                                        normValues=normalizeWindow(argsDict, binsDir, binnedValues, simPath, simName, stat, center, window)
                        else:
                                print(f"no {stat} files for center = {center}")
                        if not argsDict['keepStats']:
                                if stat in ["DIND","hDo","hDs","hf","lf","S","HAF","H12"]:
                                        os.remove(f"{simPath}/stats/{simName}_c{center}.{stat}")
                                elif stat in ["ihs","iSAFE","nsl"]:
                                        os.remove(f"{simPath}/stats/center_{center}/{simName}_c{center}.{stat}")
if __name__=="__main__":
        main()
