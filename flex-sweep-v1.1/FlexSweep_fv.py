#!/bin/python3

import argparse
import sys
import os
import warnings
import subprocess
from multiprocessing import cpu_count
from joblib import Parallel, delayed
import glob
import shutil
import tarfile
from utils import runNormNeutral, runNormStats, runFV

def parse_arguments(commandline):
        # get command line arguments
        parser = argparse.ArgumentParser(description='calculate statistics and create feature vectors from existing ms-formatted simulated data')
                
        ### this makes the feature vector for training data        
        parser.add_argument('outputDir', help='name of output directory to use in run directory')
        parser.add_argument('numberChroms', type=int, help='the number of chromosomes sampled (eg. 100 diploid individuals = 200 chromosomes')
        parser.add_argument('--rMap', default=False, help='optional, the path to a recombination map')
        parser.add_argument('--numJobs', type=int, default=cpu_count(), help='the number of processors available to run the program')
        parser.add_argument('--locusLength', type=int, default=1200000, help='UNTESTED; length of locus to simulate, in bp')
        parser.add_argument('--distCenters', type=int, default=10000, help='UNTESTED; distance between center points for statistic calculation, in bp')
        parser.add_argument('--minCenter', type=int, default=500000, help='UNTESTED; minimum center point in locus, in bp')
        parser.add_argument('--maxCenter', type=int, default=700000, help='UNTESTEDmaximum center point in locus, in bp')
        parser.add_argument('--continue', action='store_true', help='continue with data from a previous run')
        parser.add_argument('--keepStats', action='store_true', help='keep intermediate statistics files after feature vector has been generated; default is False\nWARNING: this will be large')
        parser.add_argument('--keepSims', action='store_true', help='keep simulation files after statistics have been calculated; default is True')
        parser.add_argument('--onlyFV', default=False, help='only calculate feature vectors, use only when statistics have already been calculated and normalized')

        args = parser.parse_args(commandline)
        argsDict = vars(args)
        return argsDict

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
                        binsUnfinished.add(stat)
        return binsUnfinished

def chromAllStatsMissing(argsDict, sims):
        # check if all chromosome-wide stats exist
        simPaths = [os.path.dirname(x) for x in sims]
        sims = [os.path.basename(x) for x in sims]
        simPre = [os.path.splitext(x)[0] for x in sims]
        simDict = dict(zip(simPre,simPaths))
        statsUnfinished=set()
        for sim in simDict:
                if "neutral" in sims[0]:
                        statFile = pd.read_csv(f"{argsDict['outputDir']}/training_data/neutral_data/stats/{sim}_fulllocus.stats", sep='\s', names=["stat", "position","value","DAF"], engine='python', skiprows=1) # load full statfile
                        for stat in ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl","H12","HAF"]:
                                statVals = statFile.loc[statFile['stat'] == stat]
                                if statVals.empty:
                                        statsUnfinished.add((sim, stat))
                else:
                        for stat in ["DIND","hDo","hDs","hf","lf","S","H12","HAF"]:
                                statVals = statFile.loc[statFile['stat'] == stat]
                                if statVals.empty:
                                        statsUnfinished.add((sim, stat))
        return statsUnfinished

def centerAllStatsMissing(argsDict, simList):
        paths = [os.path.dirname(x) for x in simList]
        simList = [os.path.basename(x) for x in simList]
        pre = [os.path.splitext(x)[0] for x in simList]
        simDict = dict(zip(pre,paths))
        statsUnfinished=set()
        for sim in simDict:
                statFile = pd.read_csv(f"{argsDict['outputDir']}/training_data/neutral_data/stats/{sim}.stats", sep='\s', names=["center", "stat", "position","value","DAF"], engine='python', skiprows=1) # load full statfile
                for stat in ["ihs","iSAFE","nsl"]:
                        for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                                statVals = statFile.loc[(statFile['stat'] == stat) & (statFile['center'] == center)]
                                if statVals.empty:
                                        statsUnfinished.add((sim, center, stat))
        return statsUnfinished

def chromIndStatsMissing(argsDict, simFile):
        # check if neutralnorm stats and exist for simFile
        chromStatsUnfinished=set()
        simtype="neutral"
        simFile=simFile.split(".")[0]
        if os.path.exists((f"{argsDict['outputDir']}/training_data/neutral_data/stats/{sim}_fulllocus.stats"): # check if individual stats exist
                statFile = pd.read_csv(f"{argsDict['outputDir']}/training_data/neutral_data/stats/{sim}_fulllocus.stats", sep='\s', names=["stat", "position","value","DAF"], engine='python', skiprows=1) # load full statfile
                for stat in ["HAF","H12","DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl"]:
                        statVals = statFile.loc[statFile['stat'] == stat]
                        if statVals.empty:
                                chromStatsUnfinished.add((sim, stat))
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
                                if not os.path.exists(f"{path}/stats/center_{center}/{pre}_c{center}.{stat}") or not os.path.getsize(f"{path}/stats/center_{center}/{pre}_c{center}.{stat}") > 0:
                                        centerStatsUnfinished.add(f"{path}/stats/center_{center}/{pre}_c{center}.{stat}")
        return centerStatsUnfinished

def checkSims(argsDict, simFile, simType):
        simsUnfinished = []
        simsCount = 0
        try:
                with open(simFile, "rb") as fp:
                        c_generator = _count_generator(fp.raw.read)
                        simsCount = sum(buffer.count(b'\n') for buffer in c_generator)
                        simsCount += 1
        except:
                simsCount=0
        if simsCount < argsDict["numberChroms"]+6:
                simsUnfinished.append(simFile)
        if simsUnfinished:
                print(f"{simType} training simulations {simsUnfinished} do not have enough lines to be complete. Rerun Flexsweep_simulate.py using --continue flag if necessary, then restart.")
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
                        calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        calcOutput,calcError=calcCommandRun.communicate()
        else:
                if os.path.exists(f"{simFile}.statout"):
                        os.remove(f"{simFile}.statout")
                if os.path.exists(f"{simFile}.err"):
                        os.remove(f"{simFile}.err")
                if argsDict['rMap']:
                        calcCommandLine=f"calculate_stats/runCalculateNormStats.sh {simFile} {argsDict['locusLength']} {argsDict['rMap']}"
                else:
                        calcCommandLine=f"calculate_stats/runCalculateNormStats.sh {simFile} {argsDict['locusLength']}"
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
                        calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                        calcOutput,calcError=calcCommandRun.communicate()
        else:
                if os.path.exists(f"{simFile}.statout"):
                        os.remove(f"{simFile}.statout")
                if os.path.exists(f"{simFile}.err"):
                        os.remove(f"{simFile}.err")
                if argsDict['rMap']:
                        calcCommandLine=f"calculate_stats/runCalculateTrainingStats.sh {simFile} {argsDict['locusLength']} {argsDict['rMap']}"
                else:
                        calcCommandLine=f"calculate_stats/runCalculateTrainingStats.sh {simFile} {argsDict['locusLength']}"
                calcCommandRun=subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                calcOutput,calcError=calcCommandRun.communicate()

def normStats(argsDict):
        # check if neutral bins finished before normalizing
        binsUnfinished = binsMissing(argsDict)
        if binsUnfinished:
                print("Normalization bins from locus-wide neutral statistics have not been created. Please check simulations and run again with --continue.")
                sys.exit(1)
        # probably faster to renormalize everything than to figure out which are missing first, so no need for separate --continue
        simFiles = glob.glob(f"{argsDict['outputDir']}/training_data/*_data/stats/*_data_*.stats")
#        for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
#                for window in [50000, 100000, 200000, 500000, 1000000]:
#                        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_{center}/window_{window}", exist_ok=True)
#                        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_{center}/window_{window}/norm", exist_ok=True)
#                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}/window_{window}", exist_ok=True)
#                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}/window_{window}/norm", exist_ok=True)
        os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats/norm", exist_ok=True)
        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/norm", exist_ok=True)
        binsDir = f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins"
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runNormStats.main)(argsDict, binsDir, i) for i in simFiles)

def makeNormBins(argsDict):
        stats = ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl","H12","HAF"]
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runNormNeutral.main)(argsDict, i) for i in stats)

def calculateWrap(argsDict):
        neutralSims = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*.out")
        if len(neutralSims) < 1:
                print("No neutral simulations to use, are the compressed?")
                sys.exit(1)
        # check for all simulations
        for simFile in neutralSims:
                checkSims(argsDict, simFile, "neutral")
        sweepSims = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*.out")
        if len(sweepSims) < 1:
                print("No sweep simulations to use, are they compressed?")
                sys.exit(1)
        # check for all simulations
        for simFile in sweepSims:
                checkSims(argsDict, simFile, "sweep")

        if argsDict['continue']:
                print("continue")
                os.makedirs(f"{argsDict['outputDir']}/training_data/neutral_data/stats", exist_ok=True)
                print("Checking for any unfinished chromosome-wide neutral statistics because of the --continue flag, this may take a while")
                chromNeutralStatsUnfinished = chromAllStatsMissing(argsDict, neutralSims)
                if chromNeutralStatsUnfinished:
                        print(f"Calculating missing chromosome-wide neutral statistics")
#                        chromNeutralStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in chromNeutralStatsUnfinished]
#                        chromNeutralStatsUnfinishedSims = set([f"{'_'.join(str.split(x,'_')[0:-1])}.out".replace("stats/","") for x in chromNeutralStatsUnfinishedPartial])
                        chromNeutralStatsUnfinishedSims = set(f"{argsDict['outputDir']}/training_data/neutral_data/{x[0]}" for x in chromNeutralStatsUnfinished)
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
                        print(f"Calculating missing sliding center neutral statistics")
#                        centerNeutralStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in centerNeutralStatsUnfinished]
#                        centerNeutralStatsUnfinishedSims = [f"{'_'.join(str.split(x,'_')[0:-1])}.out" for x in centerNeutralStatsUnfinishedPartial]
#                        centerNeutralStatsUnfinishedPaths = ['/'.join(str.split(x,"/")[0:-3]) for x in centerNeutralStatsUnfinishedSims]
#                        centerNeutralStatsUnfinishedSims = [str.split(x,"/")[-1] for x in centerNeutralStatsUnfinishedSims]
#                        centerNeutralStatsUnfinishedSimsPaths = set([f"{i}/{j}" for i, j in zip(centerNeutralStatsUnfinishedPaths, centerNeutralStatsUnfinishedSims)])
                        centerNeutralStatsUnfinishedSims = set(f"{argsDict['outputDir']}/training_data/neutral_data/{x[0]}" for x in centerNeutralStatsUnfinished)
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateCenterStats)(argsDict, i) for i in centerNeutralStatsUnfinishedSims)
                if not argsDict['keepSims']:
                        print("Neutral statistics finished and --keepSims is false, compressing simulations.")
                        neutralFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*.out")
                        tarFile = tarfile.open(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_sims.tgz")
                        for i in neutralFiles:
                                tarFile.add(i)
#                                os.remove(i)
                os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats", exist_ok=True)
#                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
#                        os.makedirs(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_{center}", exist_ok=True)
                print("Checking for any unfinished chromosome-wide sweep statistics because of the --continue flag, this may take a while")
                chromSweepStatsUnfinished = chromAllStatsMissing(argsDict, sweepSims)
                chromSweepStatsUnfinished = s = set(f"{argsDict['outputDir']}/training_data/sweep_data/{x[0]}" for x in chromSweepStatsUnfinished)
#                chromSweepStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in chromSweepStatsUnfinished]
#                chromSweepStatsUnfinishedSims = set([f"{'_'.join(str.split(x,'_')[0:-1])}.out".replace("stats/","") for x in chromSweepStatsUnfinishedPartial])

                print("Checking for any unfinished sliding center sweep statistics because of the --continue flag, this may take a while")
                centerSweepStatsUnfinished = centerAllStatsMissing(argsDict, sweepSims)
#                centerSweepStatsUnfinishedPartial = [os.path.splitext(x)[0] for x in centerSweepStatsUnfinished]
#                centerSweepStatsUnfinishedSims = [f"{'_'.join(str.split(x,'_')[0:-1])}.out" for x in centerSweepStatsUnfinishedPartial]
#                centerSweepStatsUnfinishedPaths = ['/'.join(str.split(x,"/")[0:-3]) for x in centerSweepStatsUnfinishedSims]
#                centerSweepStatsUnfinishedSims = [str.split(x,"/")[-1] for x in centerSweepStatsUnfinishedSims]
#                centerSweepStatsUnfinishedSimsPaths = set([f"{i}/{j}" for i, j in zip(centerSweepStatsUnfinishedPaths, centerSweepStatsUnfinishedSims)])
                 centerSweepStatsUnfinishedSims = set(f"{argsDict['outputDir']}/training_data/sweep_data/{x[0]}" for x in centerSweepStatsUnfinished)

                sweepStatsUnfinished = list(chromSweepStatsUnfinishedSims | centerSweepStatsUnfinishedSimsPaths)

                # check for unfinished sweep stats
                if sweepStatsUnfinished:
                        print(f"Calculating missing sweep statistics")
                        Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateCenterStats)(argsDict, i) for i in sweepStatsUnfinished)
                if not argsDict['keepSims']:
                        print("Sweep statistics finished and --keepSims is false, compressing simulations.")
                        sweepFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*.out")
                        tarFile = tarfile.open(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_sims.tgz")
                        for i in sweepFiles:
                                tarFile.add(i)

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
                        print("Neutral statistics finished and --keepSims is false, compressing simulations.")
                        neutralFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_*.out")
                        tarFile = tarfile.open(f"{argsDict['outputDir']}/training_data/neutral_data/neutral_data_sims.tgz")
                        for i in neutralFiles:
                                tarFile.add(i)
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
                if not argsDict['keepSims']:
                        print("Sweep statistics finished and --keepSims is false, compressing simulations.")
                        sweepFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_*.out")
                        tarFile = tarfile.open(f"{argsDict['outputDir']}/training_data/sweep_data/sweep_data_sims.tgz")
                        for i in sweepFiles:
                                tarFile.add(i)

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
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runNormStats.main)(argsDict, f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins/", i) for i in simFiles)

def makeFV(argsDict):
        neutralStats = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_*/window_*/norm/neutral_data_*norm.*")
        neutralSims = set(["_".join(os.path.basename(x).split("_")[0:3]) for x in neutralStats])
        runFV.main(argsDict, neutralSims, "neutral")

        sweepStats = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_*/window_*/norm/sweep_data_*norm.*")
        sweepSims = set(["_".join(os.path.basename(x).split("_")[0:3]) for x in sweepStats])
        runFV.main(argsDict, sweepSims, "sweep")

def main(commandline):
        argsDict = parse_arguments(commandline)
        print("running with arguments:")
        print(argsDict)

        if not argsDict['outputDir']:
                print("Specify the output directory")
                sys.exit(1)

        if os.path.exists(argsDict["outputDir"]) and argsDict["continue"] == False:
                warnings.warn("output directory "+argsDict["outputDir"]+" already exists, files will be overwritten without --continue flag")
        elif not os.path.exists(argsDict["outputDir"]):
                os.makedirs(argsDict["outputDir"])

        if argsDict['onlyFV']:
                makeFV(argsDict)
        else:
                # calculate statistics
                calculateWrap(argsDict) # this has --continue in it

                # normalize statistics
                normStats(argsDict)

                # make feature vectors
                makeFV(argsDict)

        # remove statistics if not --keepStats
        if not argsDict['keepStats']:        
                simFiles = glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/*.*", GLOB_BRACE) # get only files (to keep norm bins)
                simFiles.extend(glob.glob(f"{argsDict['outputDir']}/training_data/neutral_data/stats/center_*"))
                tarFile = tarfile.open(f"{argsDict['outputDir']}/training_data/neutral_data/stats/neutral_data_stats.tgz")
                for i in neutralFiles:
                        tarFile.add(i)
                simFiles = glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/stats/*.*", GLOB_BRACE) # get only files (to keep norm bins)
                simFiles.extend(glob.glob(f"{argsDict['outputDir']}/training_data/sweep_data/stats/center_*"))
                tarFile = tarfile.open(f"{argsDict['outputDir']}/training_data/sweep_data/stats/sweep_data_stats.tgz")
                for i in sweepFiles:
                        tarFile.add(i)

if __name__=="__main__":
        main(sys.argv[1:])
