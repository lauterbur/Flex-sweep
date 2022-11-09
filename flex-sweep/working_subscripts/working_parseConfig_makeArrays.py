#!/bin/python3

import numpy as np
import sys
import random
import os
import makeTrainingArray
from makeTrainingArray import makeNeutral, makeSweep

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

        if len(configDict) == 16: # should have all parameters
                for param in params: # CHANGESIZE and CHANGETIME don't currently have distributions
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
                                                sys.exit()
                                elif dist == "uniform" or dist == "normal":
                                        if len(configDict[param]) != 2:
                                                print(f"You have specified a {dist} distribution for {param} but specified more than one value: {param_val}, please modify config file appropriately")
                                                sys.exit()
        return configDict

def trainingArray(argsDict, configDict):
        outputDir = argsDict["outputDir"]
        numberSims = argsDict["numberSims"]
        configFile = argsDict["configFile"]

        if not os.path.exists(f"{outputDir}/training_data/neutral_data/"):
                os.makedirs(f"{outputDir}/training_data/neutral_data/")
        if not os.path.exists(f"{outputDir}/training_data/sweep_data/"):
                os.makedirs(f"{outputDir}/training_data/sweep_data/")
        if argsDict["continue"]:
                if os.path.exists(f"{argsDict['outputDir']}/training_data/neutral/array_neutral.txt") and os.path.getsize(f"{argsDict['outputDir']}/training_data/neutral/array_neutral.txt") > 0:
                        print("--continue flag used, using existing neutral training array")
                        neutral = open(f"{argsDict['outputDir']}/training_data/neutral/array_neutral.txt", "r")
                else:
                        neutral = makeTrainingArray.makeNeutral(numberSims, configDict, outputDir)
                if os.path.exists(f"{argsDict['outputDir']}/training_data/sweep/array_sweep.txt") and os.path.getsize(f"{argsDict['outputDir']}/training_data/sweep/array_sweep.txt") > 0:
                        print("--continue flag used, using existing sweep training array")
                        sweep = open(f"{argsDict['outputDir']}/training_data/sweep/array_sweep.txt", "r")
                else:
                        sweep = makeTrainingArray.makeSweep(numberSims, configDict, outputDir)
        else:
                neutral, sweep = makeTrainingArray(numberSims, configFile, outputDir)

        return neutral, sweep

argsDict={"outputDir": "test", "numberSims": 5, "configFile": "../example/config.txt", "continue": True}
configDict=parseConfig(argsDict["configFile"])
print(configDict)
trainingArray(argsDict,configDict)
