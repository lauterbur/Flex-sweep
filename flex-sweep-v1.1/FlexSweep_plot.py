#!/bin/python3

import sys, os, warnings, subprocess
import csv

def parse_arguments(commandline):
        import argparse, os

        # get command line arguments
        parser = argparse.ArgumentParser(description='simulate data for training or testing, calculate statistics and feature vectors, train, or classify empirical or simulated data')

        parser.add_argument('outputDir', help='path and name of existing output directory to use')
        parser.add_argument('mode', help='mode to plot: train (for plotting accuracy results from training) or classify (for plotting classification results)')
        parser.add_argument('--fvLoc', help='optional, the path to the directory containing the feature vector files neutral.fv and sweep.fv\nif not used, assumes feature vectors exist in outputDir/fvs from previous run of fv mode')
        parser.add_argument('--classifyName', help='simple name of classified data')

        args = parser.parse_args(commandline)
        argsDict = vars(args)
        return argsDict

def plotAcc(argsDict):
        rCommand = f"Rscript --vanilla utils/trainPlots.R {argsDict['fvLoc']}/neutral.preds {argsDict['fvLoc']}/sweep.preds {argsDict['outputDir']}"
        rCommandRun = subprocess.Popen(rCommand.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def plotClass(argsDict):
        rCommand = f"Rscript --vanilla utils/classifyPlot.R {argsDict['outputDir']}/classification/{argsDict['classifyName']}_classes.txt {argsDict['outputDir']} {argsDict['classifyName']}"
        rCommandRun = subprocess.Popen(rCommand.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def main(commandline):

        argsDict = parse_arguments(commandline) # needs outputDir + (optional) fvLoc (train) or classifyname (classify)
        print("running with arguments:")
        print(argsDict)

        os.makedirs(f"{argsDict['outputDir']}/plots", exist_ok=True)

        if argsDict['mode'] == "train":
                if not argsDict['fvLoc']:
                        argsDict['fvLoc'] = f"{argsDict['outputDir']}/fvs"
                plotAcc(argsDict)
        elif argsDict['mode'] == "classify":
                if argsDict['classifyName']:
                        plotClass(argsDict)
                else:
                        print("Must provide --classifyName when mode is classify")
        else:
                print(f"'mode' must be either train or classify, you provided {argsDict['mode']}")
                sys.exit(1)

if __name__=="__main__":
        import sys, subprocess, os, argparse
        main(sys.argv[1:])
