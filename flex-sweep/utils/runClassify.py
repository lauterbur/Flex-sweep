import argparse,time,sys,subprocess, csv, os
import numpy as np
import tensorflow as tf
from tensorflow.keras.models import Sequential, Model
from tensorflow.keras import optimizers
#from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input
#from tensorflow.keras.layers import Conv2D, MaxPooling2D, concatenate
from tensorflow.keras import utils
#from sklearn.model_selection import train_test_split
from tensorflow.keras.preprocessing.image import ImageDataGenerator
#from tensorflow.keras.callbacks import EarlyStopping, ModelCheckpoint, History
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.experimental import CosineDecayRestarts
#from tensorflow.keras.metrics import TruePositives, FalsePositives, TrueNegatives, FalseNegatives, BinaryAccuracy, Precision, Recall, AUC
import tensorflow.keras.backend as K
import fnmatch
from multiprocessing import Process, Manager
from joblib import Parallel, delayed

session_conf = tf.compat.v1.ConfigProto(
        intra_op_parallelism_threads=1
        )
#sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)

def classifyData(argsDict, window):
        numberStats = 11
        numberWindowSizes = 5
        numberCenters = int((argsDict['maxCenter'] - argsDict['minCenter'])/argsDict['distCenters']) + 1

        outputModel = argsDict['modelLoc']

        # restore from tensorflow
        model = tf.keras.models.load_model(outputModel)

        lr_decayed_fn = CosineDecayRestarts(initial_learning_rate=0.001, first_decay_steps=300)

        opt_adam=Adam(learning_rate=lr_decayed_fn, epsilon=0.0000001, amsgrad=True)
        model.compile(loss='binary_crossentropy',
                  optimizer=opt_adam)

        weights=model.get_weights()
        print("Loaded model from disk")

        #import data from predictFile
        file = f"{argsDict['outputDir']}/classification/fvs/{argsDict['classifyName']}_{window}.fv"
        testX = np.loadtxt(file,skiprows=1, delimiter=',')
        np.reshape(testX, (1, numberStats, numberWindowSizes * numberCenters))
        testX = testX.reshape(1, numberStats, numberWindowSizes * numberCenters,1)
        validation_gen = ImageDataGenerator(
                featurewise_center=True,
                featurewise_std_normalization=True,
                horizontal_flip=False)
        validation_gen.fit(testX)

        #get predictions
#        preds = model.predict(validation_gen.standardize(testX))
        preds = model(validation_gen.standardize(testX))
        if preds[0][0] >= argsDict['threshold']:
                predictions = 0
        else:
                predictions = 1
        classDict = {1:'neutral',0:'sweep'}

        windowPair = os.path.splitext(os.path.basename(file).split("_")[1])[0]
        startPoint = windowPair.split("-")[0]
        endPoint = windowPair.split("-")[1]

        outputFile = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}_classes.txt"
        if np.isnan(preds[0][0]) or np.isnan(preds[0][1]):
                line = f"{argsDict['classifyName']}\t{startPoint}\t{endPoint}\tnan\t{preds[0][0]}\t{preds[0][1]}\n"                
        else:
                line = f"{argsDict['classifyName']}\t{startPoint}\t{endPoint}\t{classDict[predictions]}\t{preds[0][0]}\t{preds[0][1]}\n"
        return line

def main(argsDict, windows):
        lines = Parallel(n_jobs=argsDict['numJobs'], verbose=100, backend="multiprocessing")(delayed(classifyData)(argsDict, i) for i in windows)

        outputFile = f"{argsDict['outputDir']}/classification/{argsDict['classifyName']}_classes.txt"
        with open(outputFile, 'a') as out:
                out.write(f"data_classified\tstart_point\tend_point\tprediction\tprob(sweep)\tprob(neutral)\n")
                for line in lines:
                        out.write(line)
        print(f"Wrote all classifications to {outputFile}")
if __name__=="__main__":
        main()
