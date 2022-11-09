#!/bin/python3.6
import subprocess
import sys

def runNeutral(NE, MU, RHO, CHROMS, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R}"
        #print(discoal)
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runNeutralDemog(NE, MU, RHO, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        CT_SCALED = CHANGETIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -en {CT_SCALED} 0 {CHANGESIZE}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runCompleteHardSweep(NE, MU, RHO, SEL, TIME, CHROMS, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runIncompleteHardSweep(NE, MU, RHO, SEL, TIME, EAF, CHROMS, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 -c {EAF}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runCompleteSoftSweep(NE, MU, RHO, SEL, TIME, SAF, CHROMS, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 -f {SAF}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runIncompleteSoftSweep(NE, MU, RHO, SEL, TIME, SAF, EAF, CHROMS, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 -f {SAF} -c {EAF}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runCompleteHardSweepDemog(NE, MU, RHO, SEL, TIME, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        CT_SCALED = [float(x) / (4 * NE) for x in CHANGETIME]
        timesize = ""
        for i in range(len(CHANGETIME)):
                timesize = f"{timesize}-en {CT_SCALED[i]} 0 {CHANGESIZE[i]} "
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 {timesize}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runIncompleteHardSweepDemog(NE, MU, RHO, SEL, TIME, EAF, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        CT_SCALED = [float(x) / (4 * NE) for x in CHANGETIME]
        timesize = ""
        for i in range(len(CHANGETIME)):
                timesize = f"{timesize}-en {CT_SCALED[i]} 0 {CHANGESIZE[i]} "
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 -c {EAF} {timesize}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runCompleteSoftSweepDemog(NE, MU, RHO, SEL, TIME, SAF, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        CT_SCALED = [float(x) / (4 * NE) for x in CHANGETIME]
        timesize = ""
        for i in range(len(CHANGETIME)):
                timesize = f"{timesize}-en {CT_SCALED[i]} 0 {CHANGESIZE[i]} "
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 -f {SAF} {timesize}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def runIncompleteSoftSweepDemog(NE, MU, RHO, SEL, TIME, SAF, EAF, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile):
        THETA = 4 * LOCUS * NE * MU
        R = 4 * LOCUS * NE * RHO
        CT_SCALED = [float(x) / (4 * NE) for x in CHANGETIME]
        timesize = ""
        for i in range(len(CHANGETIME)):
                timesize = f"{timesize}-en {CT_SCALED[i]} 0 {CHANGESIZE[i]} "
        S = SEL * 2 * NE
        SELTIME = TIME / (4 * NE)
        discoal = f"/discoal/discoal {CHROMS} 1 {LOCUS} -t {THETA} -r {R} -ws {SELTIME} -a {S} -x 0.5 -f {SAF} -c {EAF} {timesize}"
        with open(outFile, "w") as outfile:
                subprocess.run(discoal.split(), stdout=outfile, stderr=subprocess.PIPE)

def main(runDict, outputDir, type, run, CHROMS, LOCUS):
        NE = runDict["NE"]
        MU = runDict["MU"]
        RHO = runDict["RHO"]
        outFile = f"{outputDir}/training_data/{type}_data/{type}_data_{run}.out"

        if type == "neutral":
                if "CHANGETIME" in runDict and "CHANGESIZE" in runDict:
                        CHANGETIME = runDict["CHANGETIME"]
                        CHANGESIZE = runDict["CHANGESIZE"]
                        runNeutralDemog(NE, MU, RHO, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile) 
                else:
                        runNeutral(NE, MU, RHO, CHROMS, LOCUS, outFile) 

        if type == "sweep":
                SEL = runDict["SEL"]
                TIME = runDict["TIME"]
                EAF = runDict["EAF"]
                SAF = runDict["SAF"]

                if "CHANGETIME" in runDict and "CHANGESIZE" in runDict:
                        CHANGETIME = runDict["CHANGETIME"]
                        CHANGESIZE = runDict["CHANGESIZE"]
                        timesize = ""
                        for i in range(len(runDict["CHANGETIME"])):
                                timesize = f"{timesize}-en {runDict['CHANGETIME'][i]} 0 {runDict['CHANGESIZE'][i]} "

                        if runDict["EAF"] == 1 and runDict["SAF"] == 0:
                                runCompleteHardSweepDemog(NE, MU, RHO, SEL, TIME, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile) 
                        elif runDict["EAF"] == 1 and runDict["SAF"] > 0:
                                runCompleteSoftSweepDemog(NE, MU, RHO, SEL, TIME, SAF, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile) 
                        elif runDict["EAF"] < 1 and runDict["SAF"] == 0:
                                runIncompleteHardSweepDemog(NE, MU, RHO, SEL, TIME, EAF, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile) 
                        elif runDict["EAF"] < 1 and runDict["SAF"] > 0:
                                runIncompleteSoftSweepDemog(NE, MU, RHO, SEL, TIME, SAF, EAF, CHROMS, CHANGETIME, CHANGESIZE, LOCUS, outFile) 
                else:
                        if runDict["EAF"] == 1 and runDict["SAF"] == 1:
                                runCompleteHardSweep(NE, MU, RHO, SEL, TIME, CHROMS, LOCUS, outFile) 
                        elif runDict["EAF"] == 1 and runDict["SAF"] > 0:
                                runCompleteSoftSweep(NE, MU, RHO, SEL, TIME, SAF, CHROMS, LOCUS, outFile) 
                        elif runDict["EAF"] < 1 and runDict["SAF"] == 0:
                                runIncompleteHardSweep(NE, MU, RHO, SEL, TIME, EAF, CHROMS, LOCUS, outFile) 
                        elif runDict["EAF"] < 1 and runDict["SAF"] > 0:
                                runIncompleteSoftSweep(NE, MU, RHO, SEL, TIME, SAF, EAF, CHROMS, LOCUS, outFile) 

if __name__=="__main__":
        pass
