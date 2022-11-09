#!/bin/python3.6
# run using singularity container available at TODO
# run as FlexSweep.py outputDir <--continue> <options>
# will detect the number of threads available and run using GNU and python parallel

# TODO: change nodeFile situation to something that makes more sense (# threads?)

###### Sequence to SIMULATE, TRAIN, and CLASSIFY ######
# singularity exec dl_sweeps-16.04-centos-parallel-scikitallel_hapbin.sif python3 FlexSweep_simulate.py outputDir configFile numberSims numberChroms [--locusLength 120] [--continue] [--numJobs 1]
# singularity TODO FlexSweep.py fv outputDir (--continue) (--num_task) (--rmap </path/to/rmap>) (--keep)
# singularity TODO FlexSweep.py train outputDir (--fv_split #to_use_for_train+test %train %test) (--num_task)
# singularity TODO FlexSweep.py classify outputDir (--threshold) (--num_task)

###### Sequence to make FVs using existing data, TRAIN, and CLASSIFY ######
# singularity TODO FlexSweep.py fv outputDir --data_loc </path/to/training_data_dir> (--continue) (--num_task) (--rmap </path/to/rmap>) (--keep)
# singularity TODO FlexSweep.py train outputDir --fv_loc </path/to/feature_vector_dir> (--fv_split #to_use_for_train+test %train %test) (--num_task)
# singularity TODO FlexSweep.py classify outputDir --model_loc </path/to/model_dir> (--threshold) (--num_task)

###### Sequence to TRAIN using existing fvs and CLASSIFY ######
# singularity TODO FlexSweep.py train outputDir --fv_loc </path/to/feature_vector_dir> (--fv_split #to_use_for_train+test %train %test) (--num_task)
# singularity TODO FlexSweep.py classify outputDir --model_loc </path/to/model_dir> --norm_loc </path/to/normalization_data> (--threshold) (--num_task)

###### Sequence to CLASSIFY with existing model (must include the neutral bins for normalizing data to classify, has to be same normalization data as used for training)
# singularity TODO FlexSweep.py classify outputDir --model_loc </path/to/model_dir> --norm_loc </path/to/normalization_data> (--threshold) (--num_task)

###### Directory structure ######
# should create output directory from which to grab everything unless otherwise specified with a flag
# outputDir
#       |
#       -- training_data
#               |
#               -- neutral (includes a file for each simulation)
#                       |
#                       -- stats
#                               |
#                               -- bins
#                               -- ...
#               -- sweep (includes a file for each simulation)
#                       |
#                       -- stats
#                               |
#                               -- ...
#       -- fvs (includes neutral.fv and sweep.fv)
#       -- model (includes model files, history, model test predictions)
#               |
#               -- predictions
#       -- data_windows (includes a separate file and subdirectory for each window, and a fv for each window)
#               |
#               -- window_subdirectories
#                       |
#                       -- stats
#                               |
#                               -- ... 
#       -- classification (will include a predictions file)


# to run from 
def parse_arguments(commandline):
        import argparse, os

        # get command line arguments
        parser = argparse.ArgumentParser(description='simulate data for training or testing, calculate statistics and feature vectors, train, or classify empirical or simulated data')
        parser._positionals.title = 'Modes'
        
        # set up subparsers for individual modes
        subparsers = parser.add_subparsers(dest='subparser_name', help='Mode help')
        parser_simulate = subparsers.add_parser('simulate', help='simulate data for training or testing')
        parser_fv = subparsers.add_parser('fv', help='calculate statistics and create feature vector')
        parser_train = subparsers.add_parser('train', help='train Flex-sweep')
        parser_classify = subparsers.add_parser('classify', help='classify empirical or simulated data')
        
        ## simulate mode
        ### this includes making array then simulating with discoal
        parser_simulate.add_argument('outputDir', help='path and name of output directory to create', default=os.getcwd()+'dlSweepsOutput')
        parser_simulate.add_argument('configFile', type=argparse.FileType("r"), help='configuration file for simulations, see example and template')
        parser_simulate.add_argument('numberSims', type=int, help='the number of simulations of each type (neutral and sweep) to run')
        parser_simulate.add_argument('numberChroms', type=int, help='the number of chromosomes sampled (eg. 100 diploid individuals = 200 chromosomes')
        parser_simulate.add_argument('chromSize', type=int, help='the size of the chromosome, in bp')
        parser_simulate.add_argument('nodeFile', type=str, help='file with list of allocated nodes, eg. the file at $PBS_NODEFILE')
        parser_simulate.add_argument('--continue', action='store_true', help='continue with data from a previous run')
        parser_simulate.add_argument('--num_task', type=int, nargs=1, default=1, help='the number of processors available to run the program')
        
        ## fv mode
        ### this makes the feature vector for training data        
        parser_fv.add_argument('outputDir', help='path and name of output directory to use', default=os.getcwd()+'dlSweepsOutput')
        parser_fv.add_argument('nodeFile', type=str, help='file with list of allocated nodes, eg. the file at $PBS_NODEFILE')
        parser_fv.add_argument('--data_loc', nargs=1, help='optional, the path to directory of user-provided, ms-formatted training data\nif not used, assumes simulated training data exists in outputDir from simulate mode')
        parser_fv.add_argument('--continue', action='store_true', help='continue with data from a previous run')
        parser_fv.add_argument('--num_task', type=int, nargs=1, default=1, help='the number of processors available to run the program')
        parser_fv.add_argument('--rmap', nargs=1, type=argparse.FileType("r"), help='optional, the path to a recombination map file')
        parser_fv.add_argument('--keep_stats', action='store_true', help='keep intermediate statistics files\nWARNING: this will be large')
        parser_fv.add_argument('--keep_sims', action='store_true', help='keep simulation files after feature vector has been generated')

        ## train mode
        parser_train.add_argument('outputDir', help='path and name of output directory to create', default=os.getcwd()+'dlSweepsOutput')
        parser_train.add_argument('nodeFile', type=str, help='file with list of allocated nodes, eg. the file at $PBS_NODEFILE')
        parser_train.add_argument('--fv_loc', nargs=1, help='optional, the path to the directory containing the feature vector files neutral.fv and sweep.fv\nif not used, assumes feature vectors exist in outputDir from previous run of fv mode')
        parser_train.add_argument('--fv_split', nargs=3, type=int, default="10000 80 20", help='optional, integer values for the number of simulations in the feature vector files to use for training+testing, the percent of these to use for training, and the percent to use for testing\ndefault = 10000 80 20')
        parser_train.add_argument('--num_task', type=int, nargs=1, default=1, help='the number of processors available to run the program')
        
        ## classify mode
        parser_classify.add_argument('outputDir', help='path and name of output directory to create', default=os.getcwd()+'dlSweepsOutput')
        parser_classify.add_argument('nodeFile', type=str, help='file with list of allocated nodes, eg. the file at $PBS_NODEFILE')
        parser_classify.add_argument('--model_loc', nargs=1, help='optional, the path to the directory containing the model files\nif not used, assumes model files exist in outputDir from previous run of train mode')
        parser_classify.add_argument('--norm_loc', nargs=1, help='optional, the path to the directory containing the neutral normalization data files\niif not used, assumes normalization data exists in outputDir from previous run of fv mode')
        parser_classify.add_argument('--threshold', nargs=1, type=float, help='optional, the confidence threshold at which to classify a sweep\ndefault = 0.5') 
        parser_classify.add_argument('--rmap', nargs=1, type=argparse.FileType("r"), help='optional, the path to a recombination map file')
        parser_classify.add_argument('--keep_stats', action='store_true', help='keep intermediate statistics files\nWARNING: this will be large')
        parser_classify.add_argument('--continue', action='store_true', help='continue with data from a previous run')
        parser_classify.add_argument('--num_task', type=int, nargs=1, default=1, help='the number of processors available to run the program')

        args = parser.parse_args(commandline)
        argsDict = vars(args)
        return(argsDict)


def parseConfig(configFile):
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
                                        sys.exit()
                        except(IndexError):
                                print(f"{param} does not have expected format in config file. Are you missing an '='?")
                                sys.exit()

        if len(configDict) == 16: # should have all parameters
                for param in params[0:-2]: # CHANGESIZE and CHANGETIME don't currently have distributions
                        if configDict[param+'_DIST'] == "fixed" or configDict[param+'_DIST'] == "exponential":
                                if len(configDict[param]) > 1:
                                        print(f"You have specified a configDict[param+'_DIST'] value for {param} but specified more than one value: {configDict[param]}, please modify config file appropriately")
                                        sys.exit()
                        elif configDict[param+'_DIST'] == "uniform" or configDict[param+'_DIST'] == "normal":
                                if len(configDict[param]) != 2:
                                        print(f"You have specified a {configDict[param+'_DIST']} value for {param} but specified more than one value: {configDict[param]}, please modify config file appropriatel$
                                        sys.exit()
        return(configDict)

def makeArray(configDict, datatype):
        # import MakeTrainingArray

def simulateWrap(argsDict, dlSweepsDir, datatype):
        import numpy as np
        runs=list(range(1,argsDict["numberSims"]+1))
        trainingArray(argsDict,datatype) # make training arrays
        if argsDict["continue"]:
                print("checking previous simulation runs")
                neutralUnfinished=[]
                sweepUnfinished=[]
                for run in range(1,argsDict["numberSims"]+1):
                      #  print("checking "+str(run))
                        # check if there are missing neutral and sweep simulations
                        neutralCount=0
                        sweepCount=0
                        neutralPath=argsDict["outputDir"]+"/"+datatype+"_data/neutral/neutral_data_{}.out".format(run)
                        sweepPath=argsDict["outputDir"]+"/"+datatype+"_data/sweep/sweep_data_{}.out".format(run)
                        try:
                                for line in open(neutralPath):
                                        neutralCount+=1
                        except:
                                neutralCount=0
                        try:
                                for line in open(sweepPath):
                                        sweepCount+=1
                        except:
                                sweepCount=0
                        if neutralCount < argsDict["numberChroms"]+6:
                                neutralUnfinished.append(run)
                        if sweepCount < argsDict["numberChroms"]+6:
                                sweepUnfinished.append(run)
                if neutralUnfinished and sweepUnfinished: # if not all neutral and sweep sims got correctly made from a previous run, try them again
                        print("running {} neutral simulations".format(datatype))
                        simulate(argsDict, dlSweepsDir, "neutral", neutralUnfinished,datatype)
                        print("running {} sweep simulations".format(datatype))
                        simulate(argsDict, dlSweepsDir, "sweep", sweepUnfinished,datatype)
                elif neutralUnfinished and not sweepUnfinished: # if not all neutral sims got correctly made from a previous run, try them again
                        print("running {} neutral simulations".format(datatype))
                        simulate(argsDict, dlSweepsDir, "neutral", neutralUnfinished,datatype)
                elif sweepUnfinished and not neutralUnfinished:
                        print("running {} sweep simulations".format(datatype))
                        simulate(argsDict, dlSweepsDir, "sweep", sweepUnfinished,datatype)
                else: # everything made
                        print("all sims complete")
                        pass    
        else:  # make training simulations
                print("simulate neutral")
                simulate(argsDict, dlSweepsDir, "neutral", runs, datatype)
                print("simulate sweep")
                simulate(argsDict, dlSweepsDir, "sweep", runs, datatype)

def simulate(argsDict, dlSweepsDir, simtype, runs, datatype):
        ## run simulations
        ### make sure .sh files reference scripts correctly, put output in the right place - probably need to write these programattically based on argsDict
#        runs=list(range(1,argsDict["numberSims"]+1))
        nodeLines=open(argsDict["nodeFile"]).readlines()
        numNodes=len(set(nodeLines))
#        numCores=int(len(nodeLines)/numNodes)                        
        numTasks=argsDict["num_task"]
        try:
                os.remove(dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/simulate.joblog")
        except OSError:
                pass
        parallelCommand="parallel --joblog "+dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/simulate.joblog --resume-failed --memfree 42g -u -j "+str(numTasks)+" "+dlSweepsDir+"/"+argsDict["outputDir"]+"/run_training_"+simtype+".sh "+argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/array_"+simtype+".txt {} "+str(argsDict["simParams"][0])+" "+str(argsDict["numberChroms"])+" "+argsDict["outputDir"]+" "+dlSweepsDir+" "+datatype+" "+str(argsDict["chromSize"])+" --slf "+argsDict["nodeFile"]+" --wd "+dlSweepsDir+" ::: "+" ".join(["{:d}".format(x) for x in runs])
        print(parallelCommand)
        parallelCommandRun=subprocess.Popen(parallelCommand.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        while True:
                line=parallelCommandRun.stdout.readline()
                if not line:
                        break
                print(line.strip().decode('ascii'))
                sys.stdout.flush()
        parallelOutput,parallelError=parallelCommandRun.communicate()
        # check if all sim files >= # samples + 7 lines
        parallelUnfinished=[]
        for run in range(1,argsDict["numberSims"]+1):
                count=0
                path=argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/"+simtype+"_data_{}.out".format(run)
                print(path)
                try:
                        for line in open(path):
                                count+=1
                except:
                        count=0
                if count < argsDict["numberChroms"]+6:
                        parallelUnfinished.append(run)
        if parallelUnfinished:
                sys.exit("Not all {} simulations completed, check runs {} then rerun with --continue flag".format(simtype,parallelUnfinished))

def trainingArray(argsDict,datatype):
        import makeTrainingArray
        import numpy as np
        if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data"):
                os.makedirs(argsDict["outputDir"]+"/"+datatype+"_data")
        if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/neutral"):
                os.makedirs(argsDict["outputDir"]+"/"+datatype+"_data/neutral")
        if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/sweep"):
                os.makedirs(argsDict["outputDir"]+"/"+datatype+"_data/sweep")
        if argsDict["continue"]:
                if os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/neutral/array_neutral.txt") and os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/neutral/array_neutral.txt") > 0:
                        print("neutral training array already exists")
                        pass
                else:
                        neutral=makeTrainingArray.makeNeutral(argsDict)
                if os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/sweep/array_sweep.txt") and os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/sweep/array_sweep.txt") > 0:
                        print("sweep training array already exists")
                        pass
                else:
                        sweep=makeTrainingArray.makeSweep(argsDict)
        else:
                ## make training arrays
                neutral=makeTrainingArray.makeNeutral(argsDict,datatype)
                sweep=makeTrainingArray.makeSweep(argsDict,datatype)

def fvMissing(argsDict, simtype, minCenter, maxCenter,datatype): # find what statistics are missing for each run
        runs=list(range(1,argsDict["numberSims"]+1))
        if os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/"+simtype+".fv") and os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/"+simtype+".fv") > 0:
                print("{} feature vector already exists".format(simtype))
                normUnfinished=set()
                nonNormUnfinished=set()
        elif os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats"): # check if individual stats exist
                normUnfinished=normStatsMissing(argsDict, simtype, minCenter, maxCenter,simtype,datatype)
                nonNormUnfinished=statsMissing(argsDict, simtype, minCenter, maxCenter,simtype,datatype)
        else:
                normUnfinished=runs
                nonNormUnfinished=runs
        return normUnfinished, nonNormUnfinished
        
def normStatsMissing(argsDict, simtype, minCenter, maxCenter,datatype):
        runs=list(range(1,argsDict["numberSims"]+1))
        unfinished=set()
        if os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats"): # check if individual stats exist
                for run in runs:
                        for stat in ["HAF","H12"]:
                                if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/"+simtype+"_data_"+str(run)+"."+stat) or not os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/"+simtype+"_data_"+str(run)+"."+stat):
                                        unfinished.add(run) 
                        for center in range(minCenter, maxCenter, argsDict["distCenters"]):
                                for stat in ["ihs","iSAFE","nsl"]:
                                        print('est')
#                                        print(os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/center_"+str(center)+"/norm/"+simtype+"_data_"+str(run)+"c"+str(center)+"_norm."+stat) or not os.path.getsize(argsDict["outputDir"]+"/"+datatype"+_data/"+simtype+"/stats/center_"+str(center)+"/norm/"+simtype+"_data_"+str(run)+"c"+str(center)+"_norm."+stat))
#                                        if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/center_"+str(center)+"/norm/"+simtype+"_data_"+str(run)+"c"+str(center)+"_norm."+stat) or not os.path.getsize(argsDict["outputDir"]+"/"+datatype"+_data/"+simtype+"/stats/center_"+str(center)+"/norm/"+simtype+"_data_"+str(run)+"c"+str(center)+"_norm."+stat) > 0:
#                                                unfinished.add(run)
                                for window in [50000, 100000, 200000, 500000, 1000000]:
                                        for stat in ["DIND","hDo","hDs","hf","lf","S"]:
                                                if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/center_"+str(center)+"/window_"+str(window)+"/norm/"+simtype+"_data_c"+str(center)+"_w"+str(window)+"_norm."+stat) or not os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/center_"+str(center)+"/window_"+str(window)+"/norm/"+simtype+"_data_c"+str(center)+"_w"+str(window)+"_norm."+stat) > 0:
                                                        unfinished.add(run)
        else:
                print("no {} {} stats folder".format(datatype,simtype))
        return unfinished

def statsMissing(argsDict, simtype, minCenter, maxCenter,datatype):
        runs=list(range(1,argsDict["numberSims"]+1))
        unfinished=set()
        if os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats"): # check if individual stats exist
                for run in runs:
                        for stat in ["HAF","H12"]:
                                if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/"+simtype+"_data_"+str(run)+"."+stat) or not os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/"+simtype+"_data_"+str(run)+"."+stat):
                                        unfinished.add(run) 
                        for center in range(minCenter, maxCenter, argsDict["distCenters"]):
                                for stat in ["DIND","hDo","hDs","hf","lf","S","ihs","iSAFE","nsl"]:
                                        if not os.path.exists(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/center_"+str(center)+"/"+simtype+"_data_"+str(run)+"."+stat) or not os.path.getsize(argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats/center_"+str(center)+"/"+simtype+"_data_"+str(run)+"."+stat) > 0:
                                                unfinished.add(run)
        else:
                print("no {} {} stats folder".format(datatype,simtype))
        return unfinished

def chromStatsMissing(argsDict, minCenter, maxCenter):
        # check if neutralnorm stats and norm bins exist
        chromStatsUnfinished=set()
        binsUnfinished=set()
        simtype="neutral"
        runs=list(range(1,argsDict["numberSims"]+1))
        if os.path.exists(argsDict["outputDir"]+"/training_data/"+simtype+"/stats"): # check if individual stats and bins exist
                for run in runs:
                        for stat in ["HAF","H12"]:
                                if not os.path.exists(argsDict["outputDir"]+"/training_data/"+simtype+"/stats/"+simtype+"_data_"+str(run)+"."+stat) or not os.path.getsize(argsDict["outputDir"]+"/training_data/"+simtype+"/stats/"+simtype+"_data_"+str(run)+"."+stat):
                                        chromStatsUnfinished.add(run)
                        for center in range(minCenter, maxCenter, argsDict["distCenters"]):
                                for stat in ["DIND","hDo","hDs","hf","lf","S","ih","iSAFE","nsl"]:
                                        if not os.path.exists(argsDict["outputDir"]+"/training_data/"+simtype+"/stats/center_"+str(center)+"/"+simtype+"_data_"+str(run)+"."+stat) or not os.path.getsize(argsDict["outputDir"]+"/training_data/"+simtype+"/stats/center_"+str(center)+"/"+simtype+"_data_"+str(run)+"."+stat) > 0:
                                                chromStatsUnfinished.add(run)
                for stat in ["DIND","hDo","hDs","hf","lf","S","ih","iSAFE","nsl"]: #"{}/bins/{}_bins.{}".format(path,basename,stat)
                        if not os.path.exists(argsDict["outputDir"]+"/training_data/"+simtype+"/stats/bins/"+simtype+"_data_"+str(run)+"_bins."+stat) or not os.path.getsize(argsDict["outputDir"]+"/training_data/"+simtype+"/stats/bins"+simtype+"_data_"+str(run)+"_bins."+stat) > 0:
                                binsUnfinished.add(stat)
        return chromStatsUnfinished, binsUnfinished

def fv(argsDict,datatype):

# fv mode
        # check which mode
        minCenter=500000
        maxCenter=argsDict['chromSize']-500000
        if argsDict['fv_subparser_name'] == "training": # if training
                if argsDict["continue"]:
                        chromStats, neutralBins = chromStatsMissing(argsDict)
                        neutralNormUnfinished, neutralNonNormUnfinished = fvMissing(argsDict, "neutral",datatype)
                        if neutralNonNormUnfinished and neutralNormUnfinished:
                                pass
                                # calculate stats for neutralNonNormUnfinished
                                # check if neutral bins exist, if missing check if neutralnorm stats exist, rerun if missing
                                # normalize nonNomUnfinished and neutralNormUnfinished (first combine these two in set: unfinished=neutralNonNormUnfinished.union(neutralNormUnfinished)                                       
                        elif neutralNormUnfinished and not neutralNonNormUnfinished:
                                pass
                                # check if neutral bins exist, if missing check if neutralnorm stats exist, rerun if missing
                                # normalize neutralNormUnfinished
                        elif neutralNonNormUnfinished and not neutralNormUnfinished:
                                pass
                                # calculate stats for neutralNormUnfinished because something must have gone weird
                                # check if neutral bins exist, if missing check if neutralnorm stats exist, rerun if missing
                                # normalize neutralNonNormUnfinished
                        sweepNormUnfinished, sweepNonNormUnfinished = fvMissing(argsDict, "sweep",datatype)
                        ## ditto above with sweep        
                else:
                        calculateChromStats(argsDict, dlSweepsDir)
                        calculateStats(argsDict, dlSweepsDir, "neutral", datatype)
                        calculateStats(argsDict, dlSweepsDir, "sweep", datatype)
                        getNormBins(argsDict, dlSweepsDir)                                                
                        # calculate neutralnorm, neutral, and sweep stats
                        # run neutral bins
                        # normalize neutral and sweep stats
        # if empirical:
                # check if ms flag
                # check if continue flag
                        # check for what already exists
                        # run empirical data
        # if all:
                # check if continue flag
                       # check what already exists
                       # run neutral and sim 
# make sure when stats and fv created from user-provided data, gets put into sweep and neutral directories in outputDir, just like it would if from pipeline-simulated data

def getNormBins(argsDict, dlSweepsDir):
        import glob
        # normalize neutral sims chromosome-wide
        nodeLines=open(argsDict["nodeFile"]).readlines()
        numNodes=len(set(nodeLines))
        numCores=int(len(nodeLines)/numNodes)
        numTasks=argsDict["num_task"]
        runs=list(range(1,argsDict["numberSims"]+1))
        # check if neutral training stats exist
        neutralUnfinished=[]
        for run in runs:
                count=0
                path=argsDict["outputDir"]+"/training_data/neutral/stats/neutral_data_{}.*".format(run)
                allFiles=glob.glob(path)
                if len(allFiles) == 11:
                        pass
                elif len(allFiles) < 11:
                        neutralUnfinished.append(run)
                elif len(allFiles) > 11:
                        neutralUnfinished.append(run)
        try:
                os.remove(dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/bins.joblog")
        except OSError:
                pass
        import normalizeBins
        from itertools import repeat
        from multiprocessing import Pool
        base_args=(simtype+"_data",dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/",argsDict['chromSize'],argsDict["numberSims"])
        stats=["ihs","iSAFE","nsl","DIND","hDo","hDs","hf","lf","S"]
        pool=Pool(os.cpu_count()-1) # change this because it will grab number of all cpus, not those allocated
        pool.starmap(normalizeBins.main,zip(repeat(base_args),stats)) # default freq in normalizeNeutral is .02

def calculateChromStats(argsDict, dlSweepsDir):
        nodeLines=open(argsDict["nodeFile"]).readlines()
        numNodes=len(set(nodeLines))
        numCores=int(len(nodeLines)/numNodes)
        numTasks=argsDict["num_task"]
        runs=list(range(1,argsDict["numberSims"]+1))
        neutralUnfinished=[]
        for run in runs:
                count=0
                path=argsDict["outputDir"]+"/training_data/neutral/neutral_data_{}.out".format(run)
                try:
                        for line in open(path):
                                count+=1
                except:
                        count=0
                if count < argsDict["numberChroms"]+6:
                        neutralUnfinished.append(run)
        if argsDict["continue"]:
                if neutralUnfinished:
                        print("running training neutral simulations".format(training))
                        simulate(argsDict, dlSweepsDir, "neutral", neutralUnfinished,"training")
                else:
                        pass # no neutral simulations unfinished
                try:
                        os.remove(dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/neutral/chromstats.joblog")
                except OSError:
                        pass
                normCommandLine="parallel --joblog "+dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/neutral/chromstats.joblog --resume-failed -u -j "+str(numNodes)+" singularity exec --bind "+dlSweepsDir+":/home -e "+dlSweepsDir+"/dl_sweeps-16.04-centos-parallel-scikitallel_hapbin.sif "+dlSweepsDir+"/run_calculate_stats_neutral_norm.sh "+argsDict["outputDir"]+"/training_data/neutral/neutral_data_{}.out "+str(argsDict["chromSize"])+" --slf "+argsDict["nodeFile"]+" --wd "+dlSweepsDir+ "::: "+" ".join(["{:d}".format(x) for x in runs])
                normCommandRun=subprocess.Popen(normCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                normOutput,normError=normCommandRun.communicate()
        else:
                if neutralUnfinished:
                        sys.exit("Neutral training simulations {} are not finished, necessary for normalization. Rerun using --continue flag in simulate mode or provide your own, then restart.".format(neutralUnfinished))                 
                try:
                        os.remove(dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/neutral/chromstats.joblog")
                except OSError:
                        pass
                normCommandLine="parallel --joblog "+dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/neutral/chromstats.joblog --resume-failed -u -j "+str(numNodes)+" singularity exec --bind "+dlSweepsDir+":/home -e "+dlSweepsDir+"/dl_sweeps-16.04-centos-parallel-scikitallel_hapbin.sif "+dlSweepsDir+"/run_calculate_stats_neutral_norm.sh "+argsDict["outputDir"]+"/training_data/neutral/neutral_data_{}.out "+str(argsDict["chromSize"])+" --slf "+argsDict["nodeFile"]+" --wd "+dlSweepsDir+" ::: "+" ".join(["{:d}".format(x) for x in runs])
                normCommandRun=subprocess.Popen(normCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
                normOutput,normError=normCommandRun.communicate()
        
def calculateStats(argsDict, dlSweepsDir, simtype, datatype):
        nodeLines=open(argsDict["nodeFile"]).readlines()
        numNodes=len(set(nodeLines))
        numCores=int(len(nodeLines)/numNodes)
        numTasks=argsDict["num_task"]
        runs=list(range(1,argsDict["numberSims"]+1))
        unfinished=[]
        for run in runs:
                count=0
                path=argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/"+simtype+"_data_{}.out".format(run)
                try:
                        for line in open(path):
                                count+=1
                except:
                        count=0
                if count < argsDict["numberChroms"]+6:
                        unfinished.append(run)
        if unfinished:
                sys.exit("{} {} simulations {} are not finished. Rerun using --continue flag in simulate mode or provide your own, then restart.".format(simtype,datatype,unfinished))
        else:
                dataFiles=argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/"+simtype+"_data_{}.out"
        try:
                os.remove(dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats.joblog")
        except OSError:
                pass
        statsCommandLine="parallel --joblog "+dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/"+simtype+"/stats.joblog --resume-failed -u -j "+str(numNodes)+" singularity exec --bind "+dlSweepsDir+":/home -e "+dlSweepsDir+"/dl_sweeps-16.04-centos-parallel-scikitallel_hapbin.sif "+dlSweepsDir+"/run_calculate_stats.sh "+dataFiles+" --slf "+argsDict["nodeFile"]+" --wd "+dlSweepsDir+" ::: "+" ".join(["{:d}".format(x) for x in runs])
        statsCommandRun=subprocess.Popen(statsCommandLine.split(),stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        statsOutput,statsError=statsCommandRun.communicate()

def normStats(argsDict, dlSweepsDir, simtype, datatype,minCenter,maxCenter):
        # check if bins there
        chromStatsUnfinished,binsUnfinished=chromStatsMissing(argsDict, minCenter, maxCenter)
        if binsUnfinished:
                sys.exit("no neutral training frequency bins with which to normalize stats")
        statsUnfinished=statsMissing(argsDict, simtype, minCenter, maxCenter,simtype,datatype)
        if statsUnfinished:
                sys.exit("{} {} stats incomplete".format(simtype, datatype))
        import normalizeStats
        from multiprocessing import Pool
        runs=list(range(1,argsDict["numberSims"]+1))
        binpath=dlSweepsDir+"/"+argsDict["outputDir"]+"/trainin_data/neutral/neutral_data_bins"
        base_args=(dlSweepsDir+"/"+argsDict["outputDir"]+"/"+datatype+"_data/stats",binpath,argsDict["numberSims"],minCenter,maxCenter,argsDict["distCenters"])
        pool=Pool(os.cpu_count()-1)
        pool.starmap(normalizeStats.main,zip(runs,repeat(base_args)))

#        normCommandLine=
        # run_normalize.py
        # singularity exec --bind /rsgrps/denard/lauterbur/dl_sweeps/:/home -e scripts/dl_sweeps-16.04.sif python3 /rsgrps/denard/lauterbur/dl_sweeps/scripts/calculate_stats/run_normalize_sweep.py $BASENAME     28  ${SIMPATH}/stats /rsgrps/denard/lauterbur/dl_sweeps/training_data/new_neutral/stats/bins/neutral_data_bins
        pass

def cutWindows():
        # check that neutral stats there
        # check that sweep stats there
        # run_cut_window.sh
        # singularity exec --bind /rsgrps/denard/lauterbur/dl_sweeps/:/home -e scripts/dl_sweeps-16.04.sif bash /rsgrps/denard/lauterbur/dl_sweeps/scripts/calculate_stats/run_cut_window.sh ${SIM}
        pass

def compressData():
        pass
        # make sure it only does this in pipeline once all stats are made

def train(argsDict):
        pass

def classify(argsDict):
        pass

def main(commandline):
        import sys, os, warnings, subprocess
        ## check that at least 2 (outputDir mode) arguments provided on command line
        if len(commandline) < 2 and "-h" not in commandline:
                print("make sure to include mode and outputDir")
                sys.exit(1)
        argsDict = parse_arguments(commandline)
        print("running with arguments:")
        print(argsDict)
        if os.path.exists(argsDict["outputDir"]) and argsDict["continue"] == False:
                warnings.warn("output directory "+argsDict["outputDir"]+" already exists, files may be overwritten without --continue flag")
        elif not os.path.exists(argsDict["outputDir"]):
                os.makedirs(argsDict["outputDir"])

        pythonRun = sys.executable
        dlSweepsDir = os.path.dirname(os.path.abspath(__file__))
        #2. if vars(args)['continue'] is True, look to see what files exist (so where it left off previously)

        if argsDict['pseudoempirical']:
                datatype="pseudoempirical"
        else:
                datatype="training"
        if argsDict['subparser_name'] == "pipeline":
                pipeline(argsDict)
        elif argsDict['subparser_name'] == "simulate":
                simulateWrap(argsDict, dlSweepsDir, datatype)
#                simulate(argsDict, dlSweepsDir, datatype)
        elif argsDict['subparser_name'] == "fv":
                runs=list(range(1,argsDict["numberSims"]+1))
#                fv(argsDict,datatype)
                calculateChromStats(argsDict, dlSweepsDir)
        elif argsDict['subparser_name'] == "train":
                train(argsDict)
        elif argsDict['subparser_name'] == "classify":
                classify(argsDict)
        else:
                print("must provide a mode from: pipeline, simulate, fv, train, or classify")

if __name__=="__main__":
        import sys, subprocess, os, argparse
        main(sys.argv[1:])
