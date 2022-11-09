#!/bin/python3.6
import sys, os, warnings, subprocess, shutil
os.environ['TF_XLA_FLAGS'] = '--tf_xla_enable_xla_devices'
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3'

import csv
import glob
from multiprocessing import cpu_count
from joblib import Parallel, delayed
from utils import runNormStats, runClassFV, runClassify

def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

def binsMissing(argsDict):
        # check if neutral bins for normalization exist
        binsUnfinished=set()
        for stat in ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl"]:
                if not os.path.exists(f"{argsDict['normLoc']}/neutral_data_bins.{stat}") or not os.path.getsize(f"{argsDict['normLoc']}/neutral_data_bins.{stat}") > 0:
                        binsUnfinished.add(stat)
        return binsUnfinished

def cutHapMap(argsDict):
        print(f"Cutting data into sliding windows with {argsDict['windowStep']}bp step size") 
        os.makedirs(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows", exist_ok=True)
        windowsFile = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/classification_windows.txt"
        if os.path.exists(windowsFile):
                os.remove(windowsFile)
        cutCommandLine = f"utils/cut_to_1.2Mb.sh {argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows {argsDict['classifyName']} {argsDict['hapmap'][0]} {argsDict['hapmap'][1]} {argsDict['windowStep']}"
        cutCommandRun=subprocess.Popen(cutCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        calcOutput,calcError=cutCommandRun.communicate()
        # make directory for each window
        with open(windowsFile, "r") as wF:
                windows = wF.read().splitlines()
        return windows # list of minimum and maximum value for each window

def convertVCF():
        pass

def calculateStats(argsDict, windowFile):
        windowStart = windowFile.split("_")[1].split("-")[0]
        os.makedirs(f"{windowFile}", exist_ok=True)
        os.makedirs(f"{windowFile}/stats", exist_ok=True)
        for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                os.makedirs(f"{windowFile}/stats/center_{center}", exist_ok=True)
        calcCommandLine = f"calculate_stats/runCalculateClassifyStats.sh {windowFile} {argsDict['locusLength']} {windowStart}" 
        print(calcCommandLine)
        calcCommandRun = subprocess.Popen(calcCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        calcOutput,calcError = calcCommandRun.communicate()

def calculateWrap(argsDict, windows):
        print("Calculating statistics in each sliding window")
        windowFiles = [f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{x}" for x in windows]
        Parallel(n_jobs=argsDict["numJobs"])(delayed(calculateStats)(argsDict, i) for i in windowFiles)

def normStats(argsDict, windows):
        print("Normalizing statistics from data at {argsDict['normLoc']}")
        binsUnfinished = binsMissing(argsDict)
        if binsUnfinished:
                print("Normalization bins from locus-wide neutral statistics not present at expected location {argsDict['normLoc']}.")
                sys.exit(1)
        windowFiles = [f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{x}/{argsDict['classifyName']}_{x}" for x in windows]
        for windowFile in windowFiles:
                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                        for window in [50000, 100000, 200000, 500000, 1000000]:
                                os.makedirs(f"{windowFile}/stats/center_{center}/window_{window}", exist_ok=True)
                                os.makedirs(f"{windowFile}/stats/center_{center}/window_{window}/norm", exist_ok=True)
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runNormStats.main)(argsDict, argsDict['normLoc'], i) for i in windowFiles)

def wrapFV(argsDict, windows):
        print("Creating feature vectors")
        windowFiles = [f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{x}" for x in windows]
        Parallel(n_jobs=argsDict["numJobs"])(delayed(runClassFV.main)(argsDict, i) for i in windowFiles)

def classify(argsDict, windows):
        runWindows = []
        for window in windows:
                if os.path.exists(f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv") and os.path.getsize(f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv") > 0:
                        runWindows.append(window)                
        else:
                print(f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv doesn't exist")
        runClassify.main(argsDict, runWindows)

def classifyWrap(argsDict, windows):
        outputFile = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}_classes.txt"
#        print(outputFile)
#        print(windows)
        if os.path.exists(outputFile):
                os.remove(outputFile)
        classify(argsDict, windows)

def onlyClassify(argsDict, featureVecs):
        outputFile = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}_classes.txt"
        if os.path.exists(outputFile):
                os.remove(outputFile)
        runVecs = []
        for featureVec in featureVecs:
                if os.path.exists(featureVec) and os.path.getsize(featureVec) > 0:
                        lineCount = 0
                        with open(featureVec, "rb") as fp:
                                c_generator = _count_generator(fp.raw.read)
                                lineCount = sum(buffer.count(b'\n') for buffer in c_generator)
                                lineCount += 1
                        if lineCount < 2:
                                print(f"Specified feature vector {featureVec} doesn't exist, is empty, or has only header")
                        else:
                                runVecs.append(featureVec)
        featureVecFiles = [os.path.basename(x) for x in runVecs]
        runVecs = [os.path.splitext(x)[0].split("_")[1] for x in featureVecFiles]
        runClassify.main(argsDict, runVecs)
        
argsDict={"outputDir": "fullTest", "classifyName": "fullTestData", "hapmap": ["test.hap", "test.map"], "keepWindows": True, "keepStats": True, "continue": False, "vcf": False, "modelLoc": False, "normLoc": False, "windowStep": 10000, "numJobs": 1, "locusLength": 1200000, "distCenters": 10000, "minCenter": 500000, "maxCenter": 700000, "threshold": 0.95, "classifyOnly": False}
print("running with arguments:")
print(argsDict)

if not argsDict['outputDir']:
        print("Specify the output directory")
        sys.exit(1)
if argsDict['hapmap']:
        if argsDict['vcf']:
                print("Specify either a single .hap and .map file, or a vcf(.gz) file, not both")
                sys.exit(1)
        elif len(argsDict['hapmap']) != 2:
                print(f"Specify a single .hap and a single .map file with --hapmap, you provided {argsDict['hapmap']}")
                sys.exit(1)
elif argsDict['vcf'] and len(argsDict['vcf']) != 1:
        print(f"Specify a single .vcf file with --vcf, you provided {argsDict['vcf']}")
        sys.exit(1)
elif not argsDict['hapmap'] and not argsDict['vcf']:
        print("Must specify either --hapmap or --vcf")
        sys.exit(1)

if not argsDict['modelLoc']:
        argsDict['modelLoc'] = f"{argsDict['outputDir']}/{argsDict['outputDir']}Model"
if not os.path.exists(argsDict['modelLoc']) and not os.path.getsize(f"{argsDict['modelLoc']}/saved_model.pb") > 0 and \
  not os.path.exists(f"{argsDict['modelLoc']}/assets") and not os.path.exists(f"{argsDict['modelLoc']}/variables"):
        print(f"Model location {argsDict['modelLoc']} either doesn't exist or doesn't contain required saved_model.pb, assets, or variables files")
        sys.exit(1)
if not argsDict['normLoc']:
        argsDict['normLoc'] = f"{argsDict['outputDir']}/training_data/neutral_data/stats/bins"
stats = ["ihs","iSAFE","nsl","DIND","hDo","hDs","hf","lf","S","HAF","H12"]
if not os.path.exists(argsDict['normLoc']):
        print(f"Normalization data location {argsDict['normLoc']} doesn't exist")
        sys.exit(1)
for stat in stats:
        if not os.path.exists(f"{argsDict['normLoc']}/neutral_data_bins.{stat}") and not os.path.getsize(f"{argsDict['normLoc']}/neutral_data_bins.{stat}"):
                print(f"Normalization bin data for {stat} doesn't exist at {argsDict['normLoc']}/neutral_data_bins.{stat}")
                sys.exit(1)

os.makedirs(argsDict["outputDir"], exist_ok=True)
os.makedirs(f"{argsDict['outputDir']}/classification", exist_ok=True)
os.makedirs(f"{argsDict['outputDir']}/classification/fvs", exist_ok=True)

if argsDict['classifyOnly']:
        featureVecs = argsDict['classifyOnly']
        featureVecs = glob.glob(f"{featureVecs}/*")
#        featureVecs = glob.glob(f"{argsDict['outputDir']}/classification/fvs/*")
        print(f"{argsDict['outputDir']}/classification/fvs")
        onlyClassify(argsDict, featureVecs)
        print("Classification complete")

        sys.exit()        
if argsDict['continue']:
        windowsFile = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/classification_windows.txt"
        # check for existing windows
        if os.path.exists(windowsFile) and os.path.getsize(windowsFile) > 0:
                with open(windowsFile, "r") as wF:
                        windows = wF.read().splitlines()        
                # check for missing hap and map files
                missingHapMap = False
                for window in windows:
                        if os.path.exists(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{window}.hap") and os.path.getsize(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{window}.hap") > 0 and \
                        os.path.exists(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{window}.map") and os.path.getsize(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{window}.map") > 0:
                                missingHapMap = False
                        else:
                                missingHapMap = True
                                break
                if missingHapMap:
                        cutHapMap(argsDict)
                # check for missing stats
                missingStats = []
                print("missing stats?")
                for window in windows:
                        break_flag = False
                        for stat in ["ihs","iSAFE","nsl"]:
                                for center in range(argsDict["minCenter"], argsDict["maxCenter"] + argsDict["distCenters"], argsDict["distCenters"]):
                                        file = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{window}/stats/center_{center}/{argsDict['classifyName']}_{window}_c{center}.{stat}"
                                        if os.path.exists(file) and os.path.getsize(file) > 0:
                                                pass
                                        else:
                                                print(f"{file} missing")
                                                missingStats.append(window)
                                                break_flag = True
                                                break
                                if break_flag:
                                        break_flag = False
                                        break
                        for stat in ["DIND","hDo","hDs","hf","lf","S"]:
                                file = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}Windows/{argsDict['classifyName']}_{window}/stats/{argsDict['classifyName']}_{window}_c600000.{stat}"
                                if os.path.exists(file) and os.path.getsize(file) > 0:
                                        pass
                                else:
                                        print(f"{file} missing")
                                        missingStats.append(window)
                                        break_flag = True
                                        break        
                if missingStats:
                        print("yes")
                        sys.exit()
                        missingStats = set(missingStats)
                        calculateWrap(argsDict, missingStats)
                        normStats(argsDict, windows)
                        wrapFV(argsDict, windows)
                        classifyWrap(argsDict, windows)
                # check for missing feature vectors
                else:
                        print("no")
                        print("missing FVs?")
                        missingFVs = []
                        for window in windows:
                                if not os.path.exists(f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv") or os.path.getsize(f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv") < 0:
                                        #print(f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv")
                                        missingFVs.append(window)
                        if missingFVs:
                                print("yes")
#                                normStats(argsDict, windows)
                                wrapFV(argsDict, windows)
                                classifyWrap(argsDict, windows)
                        else:
                                print("no")
                                classifyWrap(argsDict, windows)
                        
        else:
                if argsDict['hapmap']:
                        cutHapMap(argsDict)
                elif argsDict['vcf']:
                        convertVCF(argsDict)
                windows = cutHapMap(argsDict)
                calculateWrap(argsDict, windows)
                normStats(argsDict, windows)
                wrapFV(argsDict, windows)
                classifyWrap(argsDict, windows)

else:
        if argsDict['hapmap']:
                cutHapMap(argsDict)
        elif argsDict['vcf']:
                convertVCF(argsDict)
        windows = cutHapMap(argsDict)
        calculateWrap(argsDict, windows)
        normStats(argsDict, windows)
        wrapFV(argsDict, windows)
        classifyWrap(argsDict, windows)

if not argsDict['keepStats']:
        files = glob.glob(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}/*/")
        for file in files:
                shutil.rmtree(file)

if not argsDict['keepWindows']:
        shutil.rmtree(f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}/")
        
