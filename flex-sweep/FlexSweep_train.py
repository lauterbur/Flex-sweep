#!/bin/python3

import sys, os, warnings, subprocess
import csv
from utils import runTrain
from multiprocessing import cpu_count

def _count_generator(reader):
    b = reader(1024 * 1024)
    while b:
        yield b
        b = reader(1024 * 1024)

def parse_arguments(commandline):
        import argparse, os

        # get command line arguments
        parser = argparse.ArgumentParser(description='train a Flexsweep model using existing feature vectors')

        parser.add_argument('outputDir', help='name of output directory to use in run directory')
        parser.add_argument('--fvLoc', help='optional, the path to the directory containing the feature vector files neutral_train.fv, sweep_train.fv, neutral_predict,  and sweep_predict.fv\nif not used, assumes feature vectors exist in outputDir/fvs from previous run of fv mode')
        parser.add_argument('--fvSplit', nargs=2, type=int, default="10000 80", help='optional, integer values for the number of simulations in the feature vector files to use for training+testing and, the percent of these to use for training. The rest will be used for test predictions.\ndefault = 10000 80')
        parser.add_argument('--numJobs', type=int, nargs=1, default=cpu_count(), help='the number of processors available to run the program')
        parser.add_argument('--locusLength', type=int, default=1200000, help='UNTESTED; length of locus to simulate, in bp')
        parser.add_argument('--distCenters', type=int, default=10000, help='UNTESTED; distance between center points for statistic calculation, in bp')
        parser.add_argument('--minCenter', type=int, default=500000, help='UNTESTED; minimum center point in locus, in bp')
        parser.add_argument('--maxCenter', type=int, default=700000, help='UNTESTEDmaximum center point in locus, in bp')

        args = parser.parse_args(commandline)
        argsDict = vars(args)
        return argsDict

def splitFV(argsDict, type):
        if argsDict['fvLoc']:
                if os.path.exists(f"{argsDict['fvLoc']}/{type}_train.fv") and os.path.getsize(f"{argsDict['fvLoc']}/{type}_train.fv") > 0:
                        lineCount = 0
                        with open(f"{argsDict['fvLoc']}/{type}_train.fv", "rb") as fp:
                                c_generator = _count_generator(fp.raw.read)
                                lineCount = sum(buffer.count(b'\n') for buffer in c_generator)
                                lineCount += 1
                        if lineCount < 2:
                                print(f"Specified feature vector {argsDict['fvLoc']}/{type}_train.fv has no data")
                                sys.exit(1)
                else:
                        print(f"Feature vector {argsDict['fvLoc']}/{type}_train.fv does not exist or is empty")
                        sys.exit(1)
                if os.path.exists(f"{argsDict['fvLoc']}/{type}_predict.fv") and os.path.getsize(f"{argsDict['fvLoc']}/{type}_predict.fv") > 0:
                        lineCount = 0
                        with open(f"{argsDict['fvLoc']}/{type}_predict.fv", "rb") as fp:
                                c_generator = _count_generator(fp.raw.read)
                                lineCount = sum(buffer.count(b'\n') for buffer in c_generator)
                                lineCount += 1
                        if lineCount < 2:
                                print(f"Specified feature vector {argsDict['fvLoc']}/{type}_predict.fv has no data")
                                sys.exit(1)
                else:
                        print(f"Feature vector {argsDict['fvLoc']}/{type}_predict.fv does not exist or is empty")
                        sys.exit(1)
        elif not argsDict['fvLoc']:
                argsDict['fvLoc'] = f"{argsDict['outputDir']}/fvs"
                numTrainTest = int(argsDict["fvSplit"][0])
                allFVs = []
                with open(f"{argsDict['fvLoc']}/{type}.fv", "r") as fvFile:
                        read = csv.reader(fvFile)
                        header = next(read)
                        for row in read:
                                allFVs.append(row)
                if len(allFVs) < numTrainTest:
                        print(f"You specified {numTrainTest} feature vectors for training and testing the model, but only have {len(allFVs)} {type} feature vectors")
                        sys.exit(1)
                numPred = len(allFVs) - numTrainTest
                if numPred < 100:
                        print(f"You are using fewer than 100 {type} feature vectors ({numPred}) for generating accuracy data, recommended to use at least 100.")
                        if numPred == 1:
                                print(f"You only have one {type} feature vector for generating accuracy data, this is not enough!")
                                sys.exit(1)
                        elif len(allFVs) < 32:
                                print(f"You must have at least 32 {type} feature vectors, you only have {numPred}")
                                sys.exit(1)
                trainFVs = [header] + allFVs[0:numTrainTest]
                predFVs = [header] + allFVs[-numPred:]

                with open(f"{argsDict['outputDir']}/fvs/{type}_train.fv", "w") as trainFile:
                        write = csv.writer(trainFile)
                        write.writerows(trainFVs)
                with open(f"{argsDict['outputDir']}/fvs/{type}_predict.fv", "w") as predFile:
                        write = csv.writer(predFile)
                        write.writerows(predFVs)
        else:
                print("Something went wrong with fvLoc")
                sys.exit(1)

def main(commandline):
        argsDict = parse_arguments(commandline)
        print("running with arguments:")
        print(argsDict)

        if not argsDict['outputDir']:
                print("Specify the output directory")
                sys.exit(1)

        splitFV(argsDict, "neutral")
        splitFV(argsDict, "sweep")
        runTrain.main(argsDict)

if __name__=="__main__":
        import sys, subprocess, os, argparse
        main(sys.argv[1:])
