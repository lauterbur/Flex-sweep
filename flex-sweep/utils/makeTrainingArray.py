#!/bin/python3

import numpy as np
import sys
import random
import os
from scipy import stats

def parseConfig(configFile):
        print(f"Parsing config file at {configFile}") 
        config = open(configFile, "r")
        lines = [line for line in config.readlines() if line.strip()]
        configDict={}
        params = ["NE","MU","RHO","SEL","TIME","SAF","EAF","CHANGESIZE","CHANGETIME"]
        param_dists = ["NE_DIST","MU_DIST","RHO_DIST","SEL_DIST","TIME_DIST","SAF_DIST","EAF_DIST"]
        count=0
        for line in lines:
                count+=1
                if not line.strip().startswith("#"):
                        try:
                                param = line.split('=',1)[0].strip()
                                if param in params or param in param_dists:
                                        configDict[param] = [x.strip() for x in line.split('=',1)[1].strip().split('#',1)[0].split(',')]
                                else:
                                        print(f"{param} is not an appropriate parameter, check spelling.")
                                        sys.exit(1)
                        except(IndexError):
                                print(f"{param} does not have expected format in config file. Are you missing an '='?")
                                sys.exit(1)

        if len(configDict) == 16 and "CHANGESIZE" in configDict and "CHANGETIME" in configDict or len(configDict) == 14 and "CHANGESIZE" not in configDict and "CHANGETIME" not in configDict: # should have all parameters
                if "CHANGESIZE" in configDict and "CHANGETIME" in configDict:
                        realParams = params
                else:
                        realParams = params[:-2]
                for param in realParams: # CHANGESIZE and CHANGETIME don't currently have distributions
                        if param == "CHANGESIZE" or param == "CHANGETIME":
                                for i in configDict[param]:
                                        if float(i) < 0:
                                                print(f"Parameter values < 0 are invalid, check {param}, you provided {i}")
                                                sys.exit(1)
                        else:
                                param_val = configDict[param]
                                dist = configDict[param+'_DIST']
                                if len(configDict[param]) == 1:
                                        if float(configDict[param][0]) < 0:
                                                print(f"Parameter values < 0 are invalid, check {param}")
                                                sys.exit(1)
                                        elif float(configDict[param][0]) > 1:
                                                if param == "MU" or param == "RHO" or param == "SEL" or param == "SAF" or param == "EAF":
                                                        print(f"{param} cannot be greater than 1, you provided {configDict[param][0]}")
                                                        sys.exit(1)
                                        elif float(configDict[param][0]) < 1:
                                                if param == "NE":
                                                        print(f"You can't have {param} less than 1!")
                                                        sys.exit(1)
                                elif len(configDict[param]) > 1:
                                        for i in configDict[param]:
                                                if float(i) < 0:
                                                        print(f"Parameter values < 0 are invalid, check {param}, you provided {i}")
                                                        sys.exit(1)
                                                elif float(i) > 1:
                                                        if param == "MU" or param == "RHO" or param == "SEL" or param == "SAF" or param == "EAF":
                                                                print(f"{param} cannot be greater than 1, you provided {i}")
                                                                sys.exit(1)
                                                elif float(i) < 1:
                                                        if param == "NE":
                                                                print(f"You can't have {param} less than 1!")
                                                                sys.exit(1)
                                if dist == "fixed" or dist == "exponential":
                                        if len(configDict[param]) > 1:
                                                print(f"You have specified a {dist} distribution for {param} but specified more than one value: {param_val}, please modify config file appropriately")
                                                sys.exit(1)
                                elif dist == "uniform" or dist == "normal":
                                        if len(configDict[param]) != 2:
                                                print(f"You have specified a {dist} distribution for {param} but specified more than one value: {param_val}, please modify config file appropriately")
                                                sys.exit(1)
        else:
                print(f"One or more parameters are not defined, please double check {configFile}. You provided: {configDict}")
                sys.exit(1)
        if "CHANGESIZE" in configDict and "CHANGETIME" in configDict:
                if len(configDict["CHANGETIME"]) != len(configDict["CHANGESIZE"]):
                        print(f"CHANGETIME and CHANGESIZE must be the same length, but have lenths {len(configDict['CHANGETIME'])} and {len(configDict['CHANGESIZE'])}")
                        sys.exit(1)
        return configDict

def makeNeutral(numberSims, configDict, outputDir):
        if "CHANGESIZE" in configDict and "CHANGETIME" in configDict:
                print(f"Making array for {numberSims} neutral simulations with the following parameters:\n \
                                \tNe: {configDict['NE_DIST']} with {configDict['NE']}\n \
                                \tper-bp mutation rate: {configDict['MU_DIST']} with {configDict['MU']}\n \
                                \tper-bp recombination rate: {configDict['RHO_DIST']} with {configDict['RHO']}\n \
                                \tpopulation sizes relative to Ne, in order backward in time: {configDict['CHANGESIZE']}\n \
                                \tpopulation size change times, in generations back in time: {configDict['CHANGETIME']}")
        else:
                print(f"Making array for {numberSims} neutral simulations with the following parameters:\n \
                                \tNe: {configDict['NE_DIST']} with {configDict['NE']}\n \
                                \tper-bp mutation rate: {configDict['MU_DIST']} with {configDict['MU']}\n \
                                \tper-bp recombination rate: {configDict['RHO_DIST']} with {configDict['RHO']}\n \
                                \tno demographic changes")
        paramDict = {}
        for param in ["NE", "MU", "RHO"]:
                dist = configDict[f"{param}_DIST"][0]
                if dist == "fixed":
                        paramDict[param] = np.array([float(configDict[param][0])] * numberSims)
                elif dist == "uniform":
                        paramDict[param] = getattr(np.random,dist)(float(configDict[param][0]),float(configDict[param][1]),size=numberSims)
                elif dist == "normal":
                        paramDict[param] = stats.truncnorm.rvs((0.0000000000000001 - float(configDict[param][0]))/float(configDict[param][1]), (np.inf - float(configDict[param][0]))/float(configDict[param][1]), size = int(numberSims), loc = float(configDict[param][0]), scale = float(configDict[param][1]))
                elif dist == "exponential":
                        paramDict[param] = getattr(np.random,dist)(float(configDict[param][0]), size=numberSims)
                else:
                        print(f"Somehow you've gotten to this point with an unacceptable distribution parameter for {param}: {dist}. Please replace with one of 'fixed', 'uniform', 'normal', or 'exponential.'")
                        sys.exit(1)
                if param == "NE":
                        np.clip(paramDict[param],1,None,paramDict[param])
                else:
                        np.clip(paramDict[param],0.00000000000001,1,paramDict[param])

        index = np.array(range(1,int(numberSims)+1))
        neutral = np.column_stack((index, paramDict['NE'], paramDict['MU'], paramDict['RHO']))
        np.savetxt(f"{outputDir}/training_data/neutral_data/array_neutral.txt",neutral,fmt=['%i','%i','%.16F','%.16F'],delimiter='\t',newline='\n',header="#ArrayIndex\tNE\tMU\tRHO", comments='')

        return neutral

def makeSweep(numberSims, configDict, outputDir):
        if "CHANGESIZE" in configDict and "CHANGETIME" in configDict:
                print(f"Making array for {numberSims} sweep simulations with the following parameters:\n \
                                \tNe: {configDict['NE_DIST']} with {configDict['NE']}\n \
                                \tper-bp mutation rate: {configDict['MU_DIST']} with {configDict['MU']}\n \
                                \tper-bp recombination rate: {configDict['RHO_DIST']} with {configDict['RHO']}\n \
                                \tselection strength: {configDict['SEL_DIST']} with {configDict['SEL']}\n \
                                \tselection time, in generations: {configDict['TIME_DIST']} with {configDict['TIME']}\n \
                                \tstarting frequency of selected allele: {configDict['SAF_DIST']} with {configDict['SAF']}\n \
                                \tending frequency of selected allele: {configDict['EAF_DIST']} with {configDict['EAF']}\n \
                                \tpopulation size changes relative to current, going back in time: {configDict['CHANGESIZE']}\n \
                                \ttimes of population size changes, in generations going back in time: {configDict['CHANGETIME']}")
        else:
                print(f"Making array for {numberSims} sweep simulations with the following parameters:\n \
                                \tNe: {configDict['NE_DIST']} with {configDict['NE']}\n \
                                \tper-bp mutation rate: {configDict['MU_DIST']} with {configDict['MU']}\n \
                                \tper-bp recombination rate: {configDict['RHO_DIST']} with {configDict['RHO']}\n \
                                \tselection strength: {configDict['SEL_DIST']} with {configDict['SEL']}\n \
                                \tselection time, in generations: {configDict['TIME_DIST']} with {configDict['TIME']}\n \
                                \tstarting frequency of selected allele: {configDict['SAF_DIST']} with {configDict['SAF']}\n \
                                \tending frequency of selected allele: {configDict['EAF_DIST']} with {configDict['EAF']}\n \
                                \tno demographic changes")
        paramDict = {}
        for param in ["NE", "MU", "RHO", "SEL", "TIME", "SAF", "EAF"]:
                dist = configDict[f"{param}_DIST"][0]
                if dist == "fixed":
                        paramDict[param] = np.array([float(configDict[param][0])] * numberSims)
                elif dist == "uniform":
                        paramDict[param] = getattr(np.random,dist)(float(configDict[param][0]),float(configDict[param][1]),size=numberSims)
                elif dist == "normal":
                        paramDict[param] = stats.truncnorm.rvs((0.0000000000000001 - float(configDict[param][0]))/float(configDict[param][1]), (np.inf - float(configDict[param][0]))/float(configDict[param][1]), size = int(numberSims), loc = float(configDict[param][0]), scale = float(configDict[param][1]))
                elif dist == "exponential":
                        paramDict[param] = getattr(np.random,dist)(float(configDict[param][0]), size=numberSims)
                else:
                        print(f"Somehow you've gotten to this point with an unacceptable distribution parameter for {param}: {dist}. Please replace with one of 'fixed', 'uniform', 'normal', or 'exponential.'")
                        sys.exit(1)
                if param == "NE":
                        np.clip(paramDict[param],1,None,paramDict[param])
                elif param != "TIME":
                        np.clip(paramDict[param],0.00000000000001,1,paramDict[param])
        index = np.array(range(1,int(numberSims)+1))
        sweep = np.column_stack((index, paramDict['NE'], paramDict['MU'], paramDict['RHO'], paramDict['SEL'], paramDict['TIME'], paramDict['SAF'], paramDict['EAF']))
        np.savetxt(f"{outputDir}/training_data/sweep_data/array_sweep.txt",sweep,fmt=['%i','%i','%.16F','%.16F','%.16F','%i','%.16F','%.16F'],delimiter='\t',newline='\n',header="#ArrayIndex\tNE\tMU\tRHO\tSEL\tTIME\tSAF\tEAF", comments='')

        return sweep

def main(numberSims, configFile, outputDir):
        configDict = parseConfig(configFile)
        neutral = makeNeutral(numberSims, configDict, outputDir)
        sweep = makeSweep(numberSims, configDict, outputDir)
        return neutral, sweep, configDict

if __name__=="__main__":
        # run script directly
        import argparse
        parser = argparse.ArgumentParser(
                                prog = "makeTrainingArray",
                                description = "make array of parameters to parameterise simulations for training Flex-sweep")
        parser.add_argument('-t', '--type', action='store', help='neutral or sweep simulations', required=True)
        parser.add_argument('-s', '--numberSims', action='store', type='int', help='number of simulations to run', required=True)
        parser.add_argument('-c', '--configFile', action='store', help='full or relative file path to configuration file', required=True)
        args = parser.parse_args()
        main(args.type, args.numberSims, args.configFile)
