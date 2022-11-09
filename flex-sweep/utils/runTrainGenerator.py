import time, argparse, sys, csv, os, pathlib

import os
import numpy as np
import tensorflow as tf
import tensorflow.keras.backend as K
from tensorflow.keras import optimizers
from tensorflow.keras import utils
from tensorflow.keras.callbacks import EarlyStopping,ModelCheckpoint, History
from tensorflow.keras.experimental import CosineDecayRestarts
from tensorflow.keras.initializers import glorot_uniform
from tensorflow.keras.layers import Dense, Dropout, Activation, Flatten, Input, Conv2D, MaxPooling2D, concatenate, ReLU
from tensorflow.keras.metrics import TruePositives, FalsePositives, TrueNegatives, FalseNegatives, BinaryAccuracy, Precision, Recall, AUC
from tensorflow.keras.models import Sequential, Model, model_from_json
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.preprocessing.image import ImageDataGenerator
from sklearn.model_selection import train_test_split
import fnmatch
#from tensorflow.keras.layers import Conv2D, MaxPooling2D,concatenate
#from tensorflow.keras.layers import LeakyReLU, PReLU, ReLU
#from tensorflow.keras.models import model_from_json

session_conf = tf.compat.v1.ConfigProto(
inter_op_parallelism_threads=1)

sess = tf.compat.v1.Session(graph=tf.compat.v1.get_default_graph(), config=session_conf)

def trainModel(argsDict):
        print("Training model")
        #training data
        numberStats = 11
        numberWindowSizes = 5
        numberCenters = 21

        neutral_data = np.loadtxt(f"{argsDict['fvLoc']}/neutral_train.fv",skiprows=1, delimiter=',')
        neutral_data = np.delete(neutral_data, np.s_[0:4], 1)
        n1 = np.reshape(neutral_data,(neutral_data.shape[0],numberStats,numberWindowSizes * numberCenters))
        
        sweep_data = np.loadtxt(f"{argsDict['fvLoc']}/sweep.fv",skiprows=1, delimiter=',')
        sweep_data = np.delete(sweep_data, np.s_[0:8], 1)
        s1 = np.reshape(sweep_data, (sweep_data.shape[0],numberStats,numberWindowSizes * numberCenters))
        
        both=np.concatenate((s1,n1))
        y=np.concatenate((np.repeat(0,len(s1)),np.repeat(1,len(n1))))
        
        both = both.reshape(both.shape[0],numberStats,numberWindowSizes * numberCenters,1)
        
        trainSplit = int(argsDict["fvSplit"][1])/100
        testSplit = 1 - trainSplit
        
        X_train, X_test, y_train, y_test = train_test_split(both, y, test_size=testSplit)
        
        Y_train = utils.to_categorical(y_train, 2)
        Y_test = utils.to_categorical(y_test, 2)
        X_valid, X_test, Y_valid, Y_test = train_test_split(X_test, Y_test, test_size=0.5)

        datagen = ImageDataGenerator(
                featurewise_center=True,
                featurewise_std_normalization=True,
                horizontal_flip=True)
        
        validation_gen = ImageDataGenerator(
                featurewise_center=True,
                featurewise_std_normalization=True,
                horizontal_flip=False)
        
        test_gen = ImageDataGenerator(
                featurewise_center=True,
                featurewise_std_normalization=True,
                horizontal_flip=False)
        
        # put model together
        model_in = Input(X_train.shape[1:])
        
        # 3x3 layer
        layer1 = Conv2D(64, 3, padding="same", name='convlayer1_1', kernel_initializer='glorot_uniform')(model_in)
        layer1 = ReLU(negative_slope=0)(layer1)
        layer1 = Conv2D(128, 3, padding="same", name='convlayer1_2', kernel_initializer='glorot_uniform')(layer1)
        layer1 = ReLU(negative_slope=0)(layer1)
        layer1 = Conv2D(256, 3, padding="same", name='convlayer1_3')(layer1)
        layer1 = ReLU(negative_slope=0)(layer1)
        layer1 = MaxPooling2D(pool_size=3, name='poollayer1', padding="same")(layer1)
        layer1 = Dropout(0.15, name='droplayer1')(layer1)
        layer1 = Flatten(name='flatlayer1')(layer1)
        
        # 2x2 layer with 1x3 dilation
        layer2 = Conv2D(64, 2,dilation_rate=[1,3], padding="same", name='convlayer2_1', kernel_initializer='glorot_uniform')(model_in)
        layer2 = ReLU(negative_slope=0)(layer2)
        layer2 = Conv2D(128, 2,dilation_rate=[1,3], padding="same", name='convlayer2_2', kernel_initializer='glorot_uniform')(layer2)
        layer2 = ReLU(negative_slope=0)(layer2)
        layer2 = Conv2D(256, 2,dilation_rate=[1,3], padding="same", name='convlayer2_3')(layer2)
        layer2 = ReLU(negative_slope=0)(layer2)
        layer2 = MaxPooling2D(pool_size=2, name='poollayer2')(layer2)
        layer2 = Dropout(0.15, name='droplayer2')(layer2)
        layer2 = Flatten(name='flatlayer2')(layer2)
        
        # 2x2 with 1x5 dilation
        layer3 = Conv2D(64, 2,dilation_rate=[1,5], padding="same", name='convlayer4_1', kernel_initializer='glorot_uniform')(model_in)
        layer3 = ReLU(negative_slope=0)(layer3)
        layer3 = Conv2D(128, 2,dilation_rate=[1,5], padding="same", name='convlayer4_2', kernel_initializer='glorot_uniform')(layer3)
        layer3 = ReLU(negative_slope=0)(layer3)
        layer3 = Conv2D(256, 2,dilation_rate=[1,5], padding="same", name='convlayer4_3')(layer3)
        layer3 = ReLU(negative_slope=0)(layer3)
        layer3 = MaxPooling2D(pool_size=2, name='poollayer3')(layer3)
        layer3 = Dropout(0.15, name='droplayer3')(layer3)
        layer3 = Flatten(name='flatlayer3')(layer3)
        
        # concatenate convolution layers
        concat =  concatenate([layer1, layer2, layer3])
        concat = Dense(512,name="512dense", activation='relu')(concat)
        concat = Dropout(0.2, name='dropconcat1')(concat)
        concat = Dense(128,name="last_dense", activation='relu')(concat)
        concat = Dropout(0.2/2, name='dropconcat2')(concat)
        output = Dense(2, name="out_dense", activation='sigmoid', kernel_initializer='glorot_uniform')(concat) # because binary classification
        model = Model(inputs=[model_in], outputs=[output])
        
        METRICS = [
          TruePositives(name='tp'),
          FalsePositives(name='fp'),
          TrueNegatives(name='tn'),
          FalseNegatives(name='fn'), 
          BinaryAccuracy(name='accuracy'),
          Precision(name='precision'),
          Recall(name='recall'),
          AUC(name='auc'),
          AUC(name='prc', curve='PR'), # precision-recall curve
        ]
        
        lr_decayed_fn = CosineDecayRestarts(initial_learning_rate=0.001, first_decay_steps=300)
        
        opt_adam = Adam(learning_rate=lr_decayed_fn, epsilon=0.0000001, amsgrad=True)
        model.compile(loss='binary_crossentropy',
          optimizer=opt_adam,
          metrics=METRICS)
        
        earlystop = EarlyStopping(monitor='val_accuracy', min_delta=0.001, patience=5, 
                                  verbose=1, mode='max',  restore_best_weights = True)
        
        outputModel = f"{argsDict['outputDir']}/{argsDict['outputDir']}Model/{argsDict['outputDir']}"
        modelPath=f"{argsDict['outputDir']}/{argsDict['outputDir']}Model"
        weightsPath = f"{outputModel}.weights.hdf5"
        
        if not argsDict["continue"]:
                if os.path.exists(modelPath):
                        print(f"Directory for model already exists at {modelPath}. Either rerun with --continue flag or delete directory")
                        sys.exit(1)
                else:
                        os.makedirs(f"{modelPath}")
        else:
                os.makedirs(f"{modelPath}", exist_ok=True)
        
        # train model
        checkpoint = ModelCheckpoint(modelPath, monitor='val_accuracy', verbose=1, save_best_only=True, mode='max')
        callbacks_list = [checkpoint, earlystop]
        
        datagen.fit(X_train)
        
        validation_gen.fit(X_valid)
        test_gen.fit(X_test)
        start = time.time()
        
        history = model.fit_generator(datagen.flow(X_train, Y_train, batch_size=32), 
                                      steps_per_epoch=len(X_train)/32, epochs=100, verbose=1, 
                                      callbacks=callbacks_list, 
                                      validation_data=validation_gen.flow(X_valid,Y_valid, batch_size=32), 
                                      validation_steps=len(X_valid)/32)
        val_score = model.evaluate_generator(validation_gen.flow(X_valid, Y_valid, batch_size=32), len(Y_valid)/32)
        test_score = model.evaluate_generator(test_gen.flow(X_test, Y_test, batch_size=32), len(Y_test)/32)
        
        train_score = model.evaluate_generator(datagen.flow(X_train, Y_train, batch_size=32), len(Y_train)/32)
        print(f"Training and testing model took {time.time()-start} seconds")
        
        weights = model.get_weights()
        with open(f"{weightsPath}.npy", 'w') as outfile:
                for slice in weights:
                        for data_slice in slice:
                                np.savetxt(outfile, data_slice.reshape((-1,)), fmt='%-7.12f')
                                outfile.write('# New slice\n')
        return modelPath
        
def testPredict(argsDict, modelPath, type):
        print("Running prediction tests with model")
        # test model by predicting on held back data
        # restore from tensorflow
        model = tf.keras.models.load_model(modelPath)
        
        METRICS = [
          TruePositives(name='tp'),
          FalsePositives(name='fp'),
          TrueNegatives(name='tn'),
          FalseNegatives(name='fn'), 
          BinaryAccuracy(name='accuracy'),
          Precision(name='precision'),
          Recall(name='recall'),
          AUC(name='auc'),
          AUC(name='prc', curve='PR'), # precision-recall curve
        ]
        
        lr_decayed_fn = CosineDecayRestarts(initial_learning_rate=0.001, first_decay_steps=300)
        opt_adam=Adam(learning_rate=lr_decayed_fn, epsilon=0.0000001, amsgrad=True)
        model.compile(loss='binary_crossentropy',
                optimizer=opt_adam,
                metrics=METRICS)
        
        # import data to predict
        numberStats = 11
        numberWindowSizes = 5
        numberCenters = 21

        testXDataParams = np.loadtxt(f"{argsDict['fvLoc']}/{type}_predict.fv", skiprows=1, delimiter=',')
        numRows, numCols = testXDataParams.shape
        if type == "neutral" and numCols == 1159:
                print("predicting neutral")
                testX = np.delete(testXDataParams, np.s_[0:4], 1) # get rid of first four columns of simulation parameters
                testXParams = testXDataParams[:,0:4]
        elif type == "sweep" and numCols == 1163: # sweep, extra params are iter, SAF, SEL, TIME, EAF, THETA, RHO
                print("predicting sweep")
                testX=np.delete(testXDataParams, np.s_[0:8], 1) # get rid of first eight columns of simulation parameters
                testXParams=testXDataParams[:,0:8]
        else:
                print(f"Predict fv at {argsDict['fvLoc']}/{type}_predict.fv has {numCols} columns, expected 1159 (neutral) or 1163 (sweep)")
                sys.exit(1)
        
        y = np.repeat(0, len(testX))
        testX = np.reshape(testX, (testX.shape[0], numberStats, numberWindowSizes * numberCenters))
        testX = testX.reshape(testX.shape[0], numberStats, numberWindowSizes * numberCenters, 1) # batch size, image width, image height, number of channels
        testY = utils.to_categorical(y, 2)
        
        validation_gen = ImageDataGenerator(
                featurewise_center=True,
                featurewise_std_normalization=True,
                horizontal_flip=False)
        
        validation_gen.fit(testX)
        
        # make predictions
        preds = model.predict(validation_gen.standardize(testX))
        predictions = np.argmax(preds,axis=1)
        
        classDict = {1:'neutral',0:'sweep'}
        
        #output the predictions
        os.makedirs(f"{modelPath}/testPredictions", exist_ok=True)
        with open(f"{modelPath}/testPredictions/{type}.preds",'w') as outFile:
                if numCols == 1155+4:
                        outFile.write('iter\tNE\tMU\tRHO\tpredicted_class\tprob(sweep)\tprob(neutral)\n')
                        for row in range(0,len(preds)):
                                outFile.write(f"{testXParams[row,0]}\t{testXParams[row,1]}\t{testXParams[row,2]}\t{testXParams[row,3]}\t{classDict[predictions[row]]}\t{preds[row][0]}\t{preds[row][1]}\n")
                elif numCols == 1155+8:
                        outFile.write('iter\tNE\tMU\tRHO\t\tSEL\tTIME\tSAF\tEAF\tpredicted_class\tprob(sweep)\tprob(neutral)\n')
                        for row in range(0,len(preds)):
                                outFile.write(f"{testXParams[row,0]}\t{testXParams[row,1]}\t{testXParams[row,2]}\t{testXParams[row,3]}\t{testXParams[row,4]}\t{testXParams[row,5]}\t{testXParams[row,6]}\t{testXParams[row,7]}\t{classDict[predictions[row]]}\t{preds[row][0]}\t{preds[row][1]}\n")

def main(argsDict):
        if argsDict['continue']:
                weightsPath = f"{argsDict['outputDir']}/{argsDict['outputDir']}Model/{argsDict['outputDir']}.weights.hdf5.npy"
                if os.path.exists(weightsPath) and os.path.getsize(weightsPath) > 0:
                        modelPath = f"{argsDict['outputDir']}/{argsDict['outputDir']}Model"
                        testPredict(argsDict, modelPath, "neutral")
                        testPredict(argsDict, modelPath, "sweep")
                else:
                        modelPath = trainModel(argsDict)
                        testPredict(argsDict, modelPath, "neutral")
                        testPredict(argsDict, modelPath, "sweep")                
        else:
                modelPath = trainModel(argsDict)
                testPredict(argsDict, modelPath, "neutral")
                testPredict(argsDict, modelPath, "sweep")
