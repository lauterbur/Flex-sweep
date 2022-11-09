#!/bin/python3.6
import sys, os, warnings, subprocess
import csv

def plotAcc(argsDict):
        rCommand = f"Rscript --vanilla utils/trainPlots.R {argsDict['outputDir']}/{argsDict['outputDir']}Model/{argsDict['outputDir']}Predictions/neutral.preds {argsDict['outputDir']}/{argsDict['outputDir']}Model/{argsDict['outputDir']}Predictions/sweep.preds {argsDict['outputDir']}"
        print(rCommand)
        rCommandRun = subprocess.Popen(rCommand.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def plotClass(argsDict):
        rCommand = f"Rscript --vanilla utils/classifyPlots.R {argsDict['fvLoc']}/neutral.preds {argsDict['fvLoc']}/sweep.preds {argsDict['outputDir']}"
        rCommandRun = subprocess.Popen(rCommand.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)

def main(commandline):

        argsDict={"mode": "train", "outputDir": "test", "classifyLoc": False, "classifyName": "testData"}

        if not argsDict['classifyLoc']:
                argsDict['classifyLoc'] = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}_classes.txt"
        print("running with arguments:")
        print(argsDict)

        os.makedirs(f"{argsDict['outputDir']}/plots", exist_ok=True)

        if argsDict['mode'] == "train":
                plotAcc(argsDict)
        elif argsDict['mode'] == "classify":
                plotClass(argsDict)
        else:
                print(f"'mode' must be either train or classify, you provided {argsDict['mode']}")
                sys.exit(1)

if __name__=="__main__":
        import sys, subprocess, os, argparse
        main(sys.argv[1:])
