# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 17:33:39 2017
    
@author: chlian
"""

import numpy as np
import SimpleITK as sitk
import scipy.io as sio
import os


datapath ="/shenlab/lab_stor5/cflian/PVS/train_data/"

mskpath = "/shenlab/lab_stor5/cflian/PVS/masks/"

savepath = "saved_model/"
if not os.path.exists(savepath):
    os.makedirs(savepath)

resultpath = "/shenlab/lab_stor4/cflian/PVS/saved_results/" 


# get model
numofseg = 1
sizeofpatch = 96
imagesize = [sizeofpatch, sizeofpatch, sizeofpatch]


traininglist = [1,8,236,237,238,239,240,241]
validatlist = [7]

batchsize = 4

samplesinepoch = 2400
samplesofvalidate = 300
numofepoch = 40

from M2EDN import get_m2edn
from M2EDN import fmeasure_loss
from generator_multi_channel_patches_prepro import Generator_MultiChannel_Patches_Balance
#from Generate_Test_Chunks import Generator_MultiChannel_Test_Chunks_prepro, Chunks_Back_To_Image
import keras


numofchannel = 1 # image + enhance_img


###########
itr = 0
###########

pathofchannels = []
pathofchannels.append(datapath)
#pathofchannels.append(resultpath+'Prob_Itr{0}_PVS'.format(itr-1))


#####
Model = get_m2edn(numofseg,numofchannel,imagesize)
#Model.load_weights(savepath+'m2edn_nmsk_iter{0}.best.hd5'.format(itr))

Trainingdata = Generator_MultiChannel_Patches_Balance(pathofchannels, mskpath,
                                                      traininglist,
                                                      batchsize, sizeofpatch,
                                                      numofseg, numofchannel)
Validatdata = Generator_MultiChannel_Patches_Balance(pathofchannels, mskpath,
                                                     validatlist, batchsize,
                                                     sizeofpatch, 
                                                     numofseg, numofchannel)

checkpoint = keras.callbacks.ModelCheckpoint(
                      filepath=savepath+'m2edn_nmsk_iter{0}.best.hd5'.format(itr),
                      save_best_only=True, mode='auto')
history = keras.callbacks.History()


Model.fit_generator(generator=Trainingdata,
                    steps_per_epoch=samplesinepoch/batchsize,
                    epochs=numofepoch,
                    validation_data=Validatdata,
                    validation_steps=samplesofvalidate/batchsize,
                    callbacks=[checkpoint,history])

sio.savemat(savepath+'trn_history_iter{0}.mat'.format(itr),
            {'train_history': history.history})

