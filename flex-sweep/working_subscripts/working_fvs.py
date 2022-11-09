#!/bin/python3

import sys, os, warnings, subprocess
from multiprocessing import cpu_count
from joblib import Parallel, delayed
import glob
import shutil
from utils import runNormNeutral
from utils import runNormStats
from utils import runFV

def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

def binsMissing(argsDict):
        # check if neutral bins for normalization exist
        binsUnfinished=set()
        for stat in ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl"]:
                if not os.path.exists(f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins/neutral_data_bins.{stat}") or not os.path.getsize(f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins/neutral_data_bins.{stat}") > 0:
#                        print(f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins/neutral_data_bins.{stat} doesn't exist?")
                        binsUnfinished.add(stat)
        return binsUnfinished

def chromAllStatsMissing(argsDict, sims):
        # check if all chromosome-wide stats exist
        simPaths = [os.path.dirname(x) for x in sims]
        sims = [os.path.basename(x) for x in sims]
        simPre = [os.path.splitext(x)[0] for x in sims]
        simDict = dict(zip(simPre,simPaths))
        center = int(argsDict["locusLength"]/2)
        statsUnfinished=set()
        for sim in simDict:
                if "neutral" in sims[0]:
                        for stat in ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl","H12","HAF"]:
                                if not os.path.exists(f"{simDict[sim]}/stats/{sim}_c{center}.{stat}") or not os.path.getsize(f"{simDict[sim]}/stats/{sim}_c{center}.{stat}") > 0:
                                        statsUnfinished.add(f"{simDict[sim]}/stats/{sim}_c{center}.{stat}")
                else:
                        for stat in ["DIND","hDo","hDs","hf","lf","S","H12","HAF"]:
                                if not os.path.exists(f"{simDict[sim]}/stats/{sim}_c{center}.{stat}") or not os.path.getsize(f"{simDict[sim]}/stats/{sim}_c{center}.{stat}") > 0:
                                        statsUnfinished.add(f"{simDict[sim]}/stats/{sim}_c{center}.{stat}")
        return statsUnfinished

def centerAllStatsMissing(argsDict, simList):
        paths = [os.path.dirname(x) for x in simList]
        simList = [os.path.basename(x) for x in simList]
        pre = [os.path.splitext(x)[0] for x in simList]
        simDict = dict(zip(pre,paths))
        statsUnfinished=set()
        for sim in simDict:
                for stat in ["ihs","iSAFE","nsl"]:
                        for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                                if not os.path.exists(f"{simDict[sim]}/stats/center_{center}/{sim}_c{center}.{stat}") or not os.path.getsize(f"{simDict[sim]}/stats/center_{center}/{sim}_c{center}.{stat}") > 0:
                                        statsUnfinished.add(f"{simDict[sim]}/stats/center_{center}/{sim}_c{center}.{stat}")
        return statsUnfinished

def chromIndStatsMissing(argsDict, simFile):
        # check if neutralnorm stats and exist for simFile
        chromStatsUnfinished=set()
        minCenter = argsDict["minCenter"]
        maxCenter = argsDict["maxCenter"]
        simtype="neutral"
        simFile=simFile.split(".")[0]
        center = argsDict["locusLength"]/2
        if os.path.exists(f"{argsDict['outputDir']}/training_data/neutral_data/stats"): # check if individual stats exist
                for stat in ["HAF","H12","DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl"]:
                        if not os.path.exists(f"{argsDict['outputDir']}/training_data/neutral_data/stats/{simFile}_c{center}.{stat}") or not os.path.getsize(f"{argsDict['outputDir']}/training_data/neutral_data/stats/{simFile}_c{center}.{stat}") > 0:
                                chromStatsUnfinished.add(f"{argsDict['outputDir']}/training_data/neutral_data/stats/{simFile}_c{center}.{stat}")
        return chromStatsUnfinished

def centerIndStatsMissing(argsDict, simFile):
        # check if norm stats and exist for simFile
        centerStatsUnfinished=set()
        minCenter = argsDict["minCenter"]
        maxCenter = argsDict["maxCenter"]
        path = os.path.dirname(simFile)
        sim = os.path.basename(simFile)
        pre = os.path.splitext(sim)[0]

        if os.path.exists(f"{path}/stats"):
                for stat in ["HAF","H12","DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl"]:
                        for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                                if not os.path.exists(f"{path}/stats/center_${center}/{pre}_c{center}.{stat}") or not os.path.getsize(f"{path}/stats/center_{center}/{pre}_c{center}.{stat}") > 0:
                                        centerStatsUnfinished.add(f"{path}/stats/center_{center}/{pre}_c{center}.{stat}")
        return centerStatsUnfinished

def makeNormBins(argsDict):
        stats = ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl","H12","HAF"]
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runNormNeutral.main)(argsDict, i) for i in stats)

def calculateWrap(argsDict):
        neutralSims = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*.out")
        if len(neutralSims) < 1:
                print("No neutral simulations to use.")
                sys.exit(1)
        # check for all simulations
        for simFile in neutralSims:
                checkSims(argsDict, simFile)
        sweepSims = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*.out")
        if len(sweepSims) < 1:
                print("No sweep simulations to use.")
                sys.exit(1)
        # check for all simulations
        for simFile in sweepSims:
                checkSims(argsDict, simFile)

        if argsDict['continue']:
                os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats", exist_ok=True)
                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_{center}", exist_ok=True)
                print("Checking for any unfinished chromosome-wide neutral statistics because of the --continue flag, this may take a while")
                chromNeutralStatsUnfinished = chromAllStatsMissing(argsDict, neutralSims)
                if chromNeutralStatsUnfinished:
#                        print(f"Calculating missing statistics {chromNeutralStatsUnfinished}")
                        print(f"Calculating missing chromosome-wide neutral statistics")
                        chromNeutralStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in chromNeutralStatsUnfinished]
                        chromNeutralStatsUnfinishedSims = set([f"{'_'.join(str.split(x,'_')[0:-1])}.out".replace("stats/","") for x in chromNeutralStatsUnfinishedPartial])
                        # check if neutral chromosome-wide stats completed
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateChromStats)(argsDict, i) for i in chromNeutralStatsUnfinishedSims)
                        print("Making neutral normalization bins")
                        bins = makeNormBins(argsDict)
                binsUnfinished = binsMissing(argsDict)
                if binsUnfinished:
                        bins = makeNormBins(argsDict)
                # check for unfinished neutral sliding center stats
                print("Checking for any unfinished sliding center neutral statistics because of the --continue flag, this may take a while")
                centerNeutralStatsUnfinished = centerAllStatsMissing(argsDict, neutralSims)
                if centerNeutralStatsUnfinished:
#                        print(f"Calculating missing statistics {centerNeutralStatsUnfinished}")
                        print(f"Calculating missing sliding center neutral statistics")
                        centerNeutralStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in centerNeutralStatsUnfinished]
                        centerNeutralStatsUnfinishedSims = [f"{'_'.join(str.split(x,'_')[0:-1])}.out" for x in centerNeutralStatsUnfinishedPartial]
                        centerNeutralStatsUnfinishedPaths = ['/'.join(str.split(x,"/")[0:-3]) for x in centerNeutralStatsUnfinishedSims]
                        centerNeutralStatsUnfinishedSims = [str.split(x,"/")[-1] for x in centerNeutralStatsUnfinishedSims]
                        centerNeutralStatsUnfinishedSimsPaths = set([f"{i}/{j}" for i, j in zip(centerNeutralStatsUnfinishedPaths, centerNeutralStatsUnfinishedSims)])
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateCenterStats)(argsDict, i) for i in centerNeutralStatsUnfinishedSimsPaths)
                if not argsDict['keepSims']:
                        print("Neutral statistics finished and --keepSims is false, removing simulations.")
                        neutralFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*ap*")
                        for i in neutralFiles:
                                os.remove(i)
                os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats", exist_ok=True)
                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}", exist_ok=True)
                print("Checking for any unfinished chromosome-wide sweep statistics because of the --continue flag, this may take a while")
                chromSweepStatsUnfinished = chromAllStatsMissing(argsDict, sweepSims)
                chromSweepStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in chromSweepStatsUnfinished]
                chromSweepStatsUnfinishedSims = set([f"{'_'.join(str.split(x,'_')[0:-1])}.out".replace("stats/","") for x in chromSweepStatsUnfinishedPartial])

                print("Checking for any unfinished sliding center sweep statistics because of the --continue flag, this may take a while")
                centerSweepStatsUnfinished = centerAllStatsMissing(argsDict, sweepSims)
                centerSweepStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in centerSweepStatsUnfinished]
                centerSweepStatsUnfinishedSims = [f"{'_'.join(str.split(x,'_')[0:-1])}.out" for x in centerSweepStatsUnfinishedPartial]
                centerSweepStatsUnfinishedPaths = ['/'.join(str.split(x,"/")[0:-3]) for x in centerSweepStatsUnfinishedSims]
                centerSweepStatsUnfinishedSims = [str.split(x,"/")[-1] for x in centerSweepStatsUnfinishedSims]
                centerSweepStatsUnfinishedSimsPaths = set([f"{i}/{j}" for i, j in zip(centerSweepStatsUnfinishedPaths, centerSweepStatsUnfinishedSims)])

                sweepStatsUnfinished = list(chromSweepStatsUnfinishedSims | centerSweepStatsUnfinishedSimsPaths)

                # check for unfinished sweep stats
                if sweepStatsUnfinished:
#                        print(f"Calculating missing statistics {sweepStatsUnfinished}")
                        print(f"Calculating missing sweep statistics")
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateCenterStats)(argsDict, i) for i in sweepStatsUnfinished)
                if not argsDict['keepSims']:
                        print("Sweep statistics finished and --keepSims is false, removing simulations.")
                        sweepFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*ap*")
                        for i in sweepFiles:
                                os.remove(i)
                        sweepFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*.out")
                        for i in sweepFiles:
                                os.remove(i)

        else:
                if os.path.exists(f"{argsDict['outputDir']}/training_data/neutral_data/stats"):
                        shutil.rmtree(f"{argsDict['outputDir']}/training_data/neutral_data/stats")
                os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats")

                Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateChromStats)(argsDict, i) for i in neutralSims)
                # check if all stats ran
                chromStatsUnfinished = chromAllStatsMissing(argsDict, neutralSims)
                if chromStatsUnfinished:
                        print(f"Not all neutral chromosome-wide statistics were calculated successfully, check {chromStatsUnfinished} and rerun with --continue")
                        sys.exit(1)
                # check if neutral bins finished before continuing calculations
                bins = makeNormBins(argsDict)
                binsUnfinished = binsMissing(argsDict)
                if binsUnfinished:
                        print("Normalization bins from locus-wide neutral statistics have not been created. Please check simulations and run again with --continue.")
                        sys.exit(1)

                # calculate neutral center stats
                os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats", exist_ok=True)
                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_{center}", exist_ok=True)
                Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateCenterStats)(argsDict, i) for i in neutralSims)
                # check if all stats ran
                centerStatsUnfinished = centerAllStatsMissing(argsDict, neutralSims)
                if centerStatsUnfinished:
                        print(f"Not all sliding center statistics for neutral were calculated successfully, check {list(centerStatsUnfinished).pop()} and rerun with --continue")
                        sys.exit(1)
                elif not centerStatsUnfinished and not argsDict['keepSims']:
                        print("Neutral statistics finished and --keepSims is false, removing simulations.")
                        neutralFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*ap*")
                        for i in neutralFiles:
                                os.remove(i)
                        neutralFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*.out")
                        for i in neutralFiles:
                                os.remove(i)

                os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats", exist_ok=True)
                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}", exist_ok=True)
                # calculate sweep stats
                Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateCenterStats)(argsDict, i) for i in sweepSims)
                chromStatsUnfinished = chromAllStatsMissing(argsDict, sweepSims)
                if chromStatsUnfinished:
                        print(f"Not all sweep chromosome-wide statistics were calculated successfully, check {chromStatsUnfinished} and rerun with --continue")
                        sys.exit(1)
                centerStatsUnfinished = centerAllStatsMissing(argsDict, sweepSims)
                if centerStatsUnfinished:
                        print(f"Not all sliding center statistics for sweeps were calculated successfully, check {centerStatsUnfinished} and rerun with --continue")
                        sys.exit(1)
                elif not centerStatsUnfinished and not argsDict['keepSims']:
                        print("Sweep statistics finished and --keepSims is false, removing simulations.")
                        sweepFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*ap*")
                        for i in sweepFiles:
                                os.remove(i)
                        sweepFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*.out")
                        for i in sweepFiles:
                                os.remove(i)

def checkSims(argsDict, simFile):
        neutralUnfinished = []
        neutralCount = 0
        try:
                with open(simFile, "rb") as fp:
                        c_generator = _count_generator(fp.raw.read)
                        neutralCount = sum(buffer.count(b'\n') for buffer in c_generator)
                        neutralCount += 1
        except:
                neutralCount=0
        if neutralCount < argsDict["numberChroms"]+6:
                neutralUnfinished.append(simFile)
        if neutralUnfinished:
                print(f"Neutral training simulations {neutralUnfinished} do not have enough lines to be complete. Rerun Flexsweep_simulate.py using --continue flag if necessary, then restart.")
                sys.exit(1)

def calculateChromStats(argsDict, simFile):
# this function will calculate full chromosome statistics for a single neutral simulation, parallelized in calculateWrap()
        if argsDict["continue"]:
                # check if any chromstats unfinished
                chromIndStatsUnfinished = chromIndStatsMissing(argsDict, simFile)
                if chromIndStatsUnfinished:
                        if os.path.exists(f"{simFile}.statout"):
                                os.remove(f"{simFile}.statout")
                        if os.path.exists(f"{simFile}.err"):
                                os.remove(f"{simFile}.err")
                        if argsDict['rMap']:
                                calcCommandLine=f"calculate_stats/runCalculateNormStats.sh {simFile} {argsDict['locusLength']} {argsDict['rMap']}"
                        else:
                                calcCommandLine=f"calculate_stats/runCalculateNormStats.sh {simFile} {argsDict['locusLength']}"
#                        with open(f"{simFile}.statout","w") as out, open(f"{simFile}.err","w") as err:
#                                calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=out, stderr=err)
                        calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        calcOutput,calcError=calcCommandRun.communicate()
        else:
                print(f"calculating statistics for {simFile}")
                if os.path.exists(f"{simFile}.statout"):
                        os.remove(f"{simFile}.statout")
                if os.path.exists(f"{simFile}.err"):
                        os.remove(f"{simFile}.err")
                if argsDict['rMap']:
                        calcCommandLine=f"calculate_stats/runCalculateNormStats.sh {simFile} {argsDict['locusLength']} {argsDict['rMap']}"
                else:
                        calcCommandLine=f"calculate_stats/runCalculateNormStats.sh {simFile} {argsDict['locusLength']}"
#                with open(f"{simFile}.statout","w") as out, open(f"{simFile}.err","w") as err:
#                        calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=out, stderr=err)
                calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                calcOutput,calcError=calcCommandRun.communicate()

def calculateCenterStats(argsDict, simFile):
# this function will calculate full chromosome statistics for a single neutral simulation, parallelized in calculateWrap()
        if argsDict["continue"]:
                # check if this stat is missing
                centerIndStatsUnfinished = centerIndStatsMissing(argsDict, simFile)
                if centerIndStatsUnfinished:
                        if os.path.exists(f"{simFile}.statout"):
                                os.remove(f"{simFile}.statout")
                        if os.path.exists(f"{simFile}.err"):
                                os.remove(f"{simFile}.err")
                        if argsDict['rMap']:
                                calcCommandLine=f"calculate_stats/runCalculateTrainingStats.sh {simFile} {argsDict['locusLength']} {argsDict['rMap']}"
                        else:
                                calcCommandLine=f"calculate_stats/runCalculateTrainingStats.sh {simFile} {argsDict['locusLength']}"
#                        with open(f"{simFile}.statout","w") as out, open(f"{simFile}.err","w") as err:
#                                calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=out, stderr=err)
                        calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        calcOutput,calcError=calcCommandRun.communicate()
        else:
                #print(f"calculating statistics for {simFile}")
                if os.path.exists(f"{simFile}.statout"):
                        os.remove(f"{simFile}.statout")
                if os.path.exists(f"{simFile}.err"):
                        os.remove(f"{simFile}.err")
                if argsDict['rMap']:
                        calcCommandLine=f"calculate_stats/runCalculateTrainingStats.sh {simFile} {argsDict['locusLength']} {argsDict['rMap']}"
                else:
                        calcCommandLine=f"calculate_stats/runCalculateTrainingStats.sh {simFile} {argsDict['locusLength']}"
#                with open(f"{simFile}.statout","w") as out, open(f"{simFile}.err","w") as err:
#                        calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=out, stderr=err)
                calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                calcOutput,calcError=calcCommandRun.communicate()

def normStats(argsDict):
        # check if neutral bins finished before normalizing
        binsUnfinished = binsMissing(argsDict)
        if binsUnfinished:
                print("Normalization bins from locus-wide neutral statistics have not been created. Please check simulations and run again with --continue.")
                sys.exit(1)
        # probably faster to renormalize everything than to figure out which are missing first, so no need for separate --continue
        simFiles = glob.glob(f"{argsDict['outputDir']}/training_data/*_data/*_data_*.out")
        for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                for window in [50000, 100000, 200000, 500000, 1000000]:
                        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_{center}/window_{window}", exist_ok=True)
                        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_{center}/window_{window}/norm", exist_ok=True)
                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}/window_{window}", exist_ok=True)
                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}/window_{window}/norm", exist_ok=True)
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runNormStats.main)(argsDict, i) for i in simFiles)

def makeFV(argsDict):
        neutralStats = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_*/window_*/norm/neutral_data_*norm.*")
        neutralSims = set(["_".join(os.path.basename(x).split("_")[0:3]) for x in neutralStats])
        runFV.main(argsDict, neutralSims, "neutral")

        sweepStats = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_*/window_*/norm/sweep_data_*norm.*")
        sweepSims = set(["_".join(os.path.basename(x).split("_")[0:3]) for x in sweepStats])
        runFV.main(argsDict, sweepSims, "sweep")

argsDict = {"outputDir": "test", "numberChroms": 20, "numJobs": cpu_count(), "locusLength": 1200000, "distCenters": 10000, "minCenter": 500000, "maxCenter": 700000, "continue": False, "rMap": False, "keepStats": True, "keepSims": True, "rMap": False, "onlyFV": False}
print(argsDict)

if argsDict['onlyFV']:
        makeFV(argsDict)
else:
        # calculate statistics
        print("calculate")
        calculateWrap(argsDict) # this has --continue in it

        # normalize statistics
        print("normalize")
        normStats(argsDict)

        # make feature vectors
        print("fv")
        makeFV(argsDict)

        # remove statistics if not --keepStats
        sys.exit()
#if not argsDict['keepStats']:        
#        simFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/*.{*}", GLOB_BRACE) # get only files (to keep norm bins)
#        simFiles.extend(glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_*
#        for i in neutralFiles:
#                os.remove(i)
