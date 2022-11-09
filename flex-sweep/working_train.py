import sys, os, warnings, subprocess
import csv
import runTrainWithWorkingTrain

def splitFV(argsDict, type):
        numTrainTest = int(argsDict["fvSplit"].split(" ")[0])
        allFVs = []
        with open(f"{argsDict['fvLoc']}/{type}.fv", "r") as fvFile:
                read = csv.reader(fvFile)
                header = next(read)
                for row in read:
                        allFVs.append(row)
        if len(allFVs) < numTrainTest:
                print(f"You specified {numTrainTest} feature vectors for training and testing the model, but only have {len(allFVs)} feature vectors")
                sys.exit(1)
        numPred = len(allFVs) - numTrainTest
        if numPred < 100:
                print(f"You are using fewer than 100 {type} feature vectors ({numPred}) for generating accuracy data, recommended to use at least 100.")
                if numPred == 1:
                        print(f"You only have one {type} feature vector for generating accuracy data, this is not enough!")
                        sys.exit(1)
        trainFVs = [header] + allFVs[0:numTrainTest]
        predFVs = [header] + allFVs[-numPred:]
        
        with open(f"{argsDict['outputDir']}/fvs/{type}_train.fv", "w") as trainFile:
                write = csv.writer(trainFile)
                write.writerows(trainFVs)
        with open(f"{argsDict['outputDir']}/fvs/{type}_predict.fv", "w") as predFile:
                write = csv.writer(predFile)
                write.writerows(predFVs)

argsDict={"fvLoc": False, "outputDir": "sifTest", "fvSplit": "5 80", "continue": True}
print(argsDict)

if not argsDict['fvLoc']:
        argsDict['fvLoc'] = f"{argsDict['outputDir']}/fvs"
if argsDict['continue']:
        if os.path.exists(f"{argsDict['outputDir']}/fvs/neutral_train.fv") and os.path.getsize(f"{argsDict['outputDir']}/fvs/neutral_train.fv") > 0 and \
            os.path.exists(f"{argsDict['outputDir']}/fvs/neutral_predict.fv") and os.path.getsize(f"{argsDict['outputDir']}/fvs/neutral_predict.fv") > 0:
                print("Neutral train and predict fvs already exist")
        else:
                splitFV(argsDict, "neutral")
        if os.path.exists(f"{argsDict['outputDir']}/fvs/sweep_train.fv") and os.path.getsize(f"{argsDict['outputDir']}/fvs/sweep_train.fv") > 0 and \
            os.path.exists(f"{argsDict['outputDir']}/fvs/sweep_predict.fv") and os.path.getsize(f"{argsDict['outputDir']}/fvs/sweep_predict.fv") > 0:
                print("Sweep train and predict fvs already exist")
        else:
                splitFV(argsDict, "sweep")
        runTrainWithWorkingTrain.main(argsDict)
else:
        splitFV(argsDict, "neutral")
        splitFV(argsDict, "sweep")
        runTrainWithWorkingTrain.main(argsDict)
