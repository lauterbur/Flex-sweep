#!/bin/python3

# this script normalizes the data from the neutral simulations
# using itself

#1. read in all derived allele frequencies and <stat> scores into a pandas table
#2. bin sites by 2% allele frequency
#3. calculate mean and SD for each bin
#4. normalize for each bin ( (value-mean)/SD )

import pandas as pd
import argparse
import glob
import numpy as np
import os

def read(argsDict, stat, fileList):
        # make a pandas table with columns of:
        ## position
        ## statistic value
        ## derived allele frequency
        values = pd.DataFrame()
#        center = int(argsDict["locusLength"]/2)
        for file in fileList:
                df = pd.read_csv(file, sep='\s', names=["position","stat","DAF"], engine='python', skiprows=0)
                values = pd.concat([values,df], ignore_index=True) # concatenate with previous files read
        values = values.loc[values['stat'] == stat]
        if stat == "HAF" or stat == "H12":
                values['DAF']=1 # because there needs to be a value for DAF to get bins
        return values

def bin(values, freq=.02):
        # bin sites by x% allele frequency
        values['freq_bins'] = pd.cut(x=values['DAF'], bins=np.arange(0,1+freq,freq), include_lowest=True, precision=2) #, labels=np.arange(0,1,freq))
        return values

def normalize(values, argsDict, stat):
        # get expected value (mean) and standard deviation
        expected = values.groupby('freq_bins')[['stat']].mean()
        stdev = values.groupby('freq_bins')[['stat']].std()
        bins = pd.concat([expected,stdev], axis=1)
        bins.columns = ['expected','stdev']
             # save these for reference later when normalizing sweep simulations
        bins.to_csv(f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins/neutral_data_bins.{stat}", na_rep="NA")
        # now normalize using those bins values
# this doesn't seem to do anything at this point
#        expected_col = values['freq_bins'].map(bins['expected'])
#        expected_col = pd.to_numeric(expected_col)
#        stdev_col = values['freq_bins'].map(bins['stdev'])
#        stdev_col = pd.to_numeric(stdev_col)
#        try:
#                stdev_col.iloc[pd.Index(stdev_col).get_loc(0)] = 1 # when stdev is 0, means no variation, so scaled value -> 0
#        except KeyError:
#                pass
#        values['stat_norm'] = (values['stat']-expected_col)/stdev_col
#        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/norm", exist_ok=True)
#        return values

def main(argsDict, stat):
        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins", exist_ok=True)
        fileList = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/neutral_data_*_fulllocus.stats")
        values = read(argsDict, stat)
        bin_values = bin(values)
        norm_values = normalize(bin_values, argsDict, stat)

if __name__=="__main__":
        main()
