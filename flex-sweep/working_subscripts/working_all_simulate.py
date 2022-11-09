#!/bin/python3

import numpy as np
import sys
import random
import os
import makeTrainingArray
from makeTrainingArray import makeNeutral, makeSweep, parseConfig
import runDiscoal
from multiprocessing import Pool, cpu_count
from joblib import Parallel, delayed

def trainingArray(argsDict):
        outputDir = argsDict["outputDir"]
        numberSims = argsDict["numberSims"]
        configFile = argsDict["configFile"]

        if not os.path.exists(f"{outputDir}/training_data/neutral_data/"):
                os.makedirs(f"{outputDir}/training_data/neutral_data/")
        if not os.path.exists(f"{outputDir}/training_data/sweep_data/"):
                os.makedirs(f"{outputDir}/training_data/sweep_data/")
        if argsDict["continue"]:
                if os.path.exists(f"{argsDict['outputDir']}/training_data/neutral_data/array_neutral.txt") and os.path.getsize(f"{argsDict['outputDir']}/training_data/neutral_data/array_neutral.txt") > 0:
                        print("--continue flag used, using existing neutral training array")
                        neutral = open(f"{argsDict['outputDir']}/training_data/neutral_data/array_neutral.txt", "r")
                else:
                        configDict = parseConfig(configFile)
                        neutral = makeNeutral(numberSims, configDict, outputDir)
                if os.path.exists(f"{argsDict['outputDir']}/training_data/sweep_data/array_sweep.txt") and os.path.getsize(f"{argsDict['outputDir']}/training_data/sweep_data/array_sweep.txt") > 0:
                        print("--continue flag used, using existing sweep training array")
                        configDict = parseConfig(configFile)
                        sweep = open(f"{argsDict['outputDir']}/training_data/sweep_data/array_sweep.txt", "r")
                else:
                        sweep = makeTrainingArray.makeSweep(numberSims, configDict, outputDir)
        else:
                neutral, sweep, configDict = makeTrainingArray.main(numberSims, configFile, outputDir)

        return neutral, sweep, configDict

def simulate(argsDict, configDict, simtype, run, chroms, locus):
        # read array into dictionary with column headers as keys
        import csv
        with open(f"{argsDict['outputDir']}/training_data/{simtype}_data/array_{simtype}.txt", "r") as arrayfile:
                reader = csv.DictReader(arrayfile, skipinitialspace=True, delimiter='\t')
                arrayDict = {name: [] for name in reader.fieldnames}
                for row in reader:
                        for name in reader.fieldnames:
                                arrayDict[name].append(row[name])
        runDict = {}
        runDict["NE"] = int(arrayDict["NE"][run-1]) 
        runDict["MU"] = float(arrayDict["MU"][run-1]) 
        runDict["RHO"] = float(arrayDict["RHO"][run-1]) 
        if simtype == "sweep":
                runDict["SEL"] = float(arrayDict["SEL"][run-1]) 
                runDict["TIME"] = int(arrayDict["TIME"][run-1]) 
                runDict["SAF"] = float(arrayDict["SAF"][run-1]) 
                runDict["EAF"] = float(arrayDict["EAF"][run-1]) 
                if "CHANGETIME" in configDict and "CHANGESIZE" in configDict:
                        runDict["CHANGETIME"] = configDict["CHANGETIME"]
                        runDict["CHANGESIZE"] = configDict["CHANGESIZE"]
        runDiscoal.main(runDict, argsDict['outputDir'], simtype, run, chroms, locus) 

def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

def simulateWrap(argsDict):
        import numpy as np
        runs=list(range(1,argsDict["numberSims"]+1))
        neutral, sweep, configDict = trainingArray(argsDict) # make training arrays
        if argsDict["continue"]:
                print("checking previous simulation runs because you ran with --continue, this can take a little while")
                neutralUnfinished=[]
                sweepUnfinished=[]
                for run in range(1,argsDict["numberSims"]+1):
                      #  print("checking "+str(run))
                        # check if there are missing neutral and sweep simulations
                        neutralCount=0
                        sweepCount=0
                        run_str = str(run)
                        neutralPath=f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_{run_str}.out"
                        sweepPath=f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_{run_str}.out"
                        try:
                                with open(neutralPath, "rb") as fp:
                                        c_generator = _count_generator(fp.raw.read)
                                        neutralCount = sum(buffer.count(b'\n') for buffer in c_generator)
                                        neutralCount += 1
                        except:
                                neutralCount=0
                        try:
                                with open(sweepPath, "rb") as fp:
                                        c_generator = _count_generator(fp.raw.read)
                                        sweepCount = sum(buffer.count(b'\n') for buffer in c_generator)
                                        sweepCount += 1
                        except:
                                sweepCount=0
                        if neutralCount < argsDict["numberChroms"]+6:
                                neutralUnfinished.append(run)
                        if sweepCount < argsDict["numberChroms"]+6:
                                sweepUnfinished.append(run)
                numMissNeut = len(neutralUnfinished)
                numMissSweep = len(sweepUnfinished)
                if neutralUnfinished and sweepUnfinished: # if not all neutral and sweep sims got correctly made from a previous run, try them again
                        print(f"running {numMissNeut} missing neutral simulations")
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(simulate)(argsDict, configDict, "neutral", i, argsDict["numberChroms"], argsDict["locus"]) for i in neutralUnfinished)
                        unFinished(argsDict, "neutral")
                        print(f"running {numMissSweep} missing sweep simulations")
                        unFinished(argsDict, "sweep")
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(simulate)(argsDict, configDict, "sweep", i, argsDict["numberChroms"], argsDict["locus"]) for i in sweepUnfinished)
                elif neutralUnfinished and not sweepUnfinished: # if not all neutral sims got correctly made from a previous run, try them again
                        print(f"running {numMissNeut} missing neutral simulations")
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(simulate)(argsDict, configDict, "neutral", i, argsDict["numberChroms"], argsDict["locus"]) for i in neutralUnfinished)
                        unFinished(argsDict, "neutral")
                elif sweepUnfinished and not neutralUnfinished:
                        print(f"running {numMissSweep} missing sweep simulations")
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(simulate)(argsDict, configDict, "sweep", i, argsDict["numberChroms"], argsDict["locus"]) for i in sweepUnfinished)
                        unFinished(argsDict, "sweep")
                else: # everything made
                        print("all sims complete")
                        pass
        else:  # make training simulations
                print("simulate neutral")
                Parallel(n_jobs=argsDict["numJobs"])(delayed(simulate)(argsDict, configDict, "neutral", i, argsDict["numberChroms"], argsDict["locus"]) for i in range(1,argsDict["numberSims"] + 1))
                unFinished(argsDict, "neutral")
                print("simulate sweep")
                Parallel(n_jobs=argsDict["numJobs"])(delayed(simulate)(argsDict, configDict, "sweep", i, argsDict["numberChroms"], argsDict["locus"]) for i in range(1,argsDict["numberSims"] + 1))
                unFinished(argsDict, "sweep")

def unFinished(argsDict, simtype):
        # check if all sim files >= # samples + 6 lines
        parallelUnfinished=[]
        for run in range(1,argsDict["numberSims"]+1):
                count=0
                run_str = str(run)
                path=f"{argsDict['outputDir']}/training_data/{simtype}_data/{simtype}_data_{run_str}.out"
                try:
                        with open(path, "rb") as fp:
                                c_generator = _count_generator(fp.raw.read)
                                count = sum(buffer.count(b'\n') for buffer in c_generator)
                                count += 1
                except:
                        count=0
                if count < argsDict["numberChroms"]+6:
                        parallelUnfinished.append(run)
        if parallelUnfinished:
                sys.exit(f"Not all {simtype} simulations completed, check runs {parallelUnfinished} then rerun with --continue flag")

argsDict={"outputDir": "test_sif", "numberSims": 5, "numberChroms": 20, "numJobs": cpu_count(), "configFile": "config_testsif.txt", "continue": False, "locus": 1200000}
#configDict=parseConfig(argsDict["configFile"])
simulateWrap(argsDict)
#jobs=cpu_count()
#Parallel(n_jobs=jobs)(delayed(simulate)(argsDict, configDict, "neutral", i, 5, argsDict["locus"]) for i in range(5))
#Parallel(n_jobs=jobs)(delayed(simulate)(argsDict, configDict, "sweep", i, 5, argsDict["locus"]) for i in range(5))
