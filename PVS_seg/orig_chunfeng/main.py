#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 17:33:39 2017
    
@author: chlian
"""

import numpy as np
import SimpleITK as sitk
import os
import keras

from M2EDN import get_m2edn
from Generate_Test_Chunks import Generator_MultiChannel_Test_Chunks_cl
from Generate_Test_Chunks import Chunks_Back_To_Image_cl
from keras import backend as K

datapath = '/nas/longleaf/home/zongx/pine/EveMask_2T2Carb/T2Carb/' #path of input data
sample_list = ['T2_PVS46_CARBOGEN.nii.gz'] #sample index

resultpath = 'saved_result/'
modelpath = 'saved_model/'

###########
itr = 0 #for auto-context iteration 0
#itr = 1 #for auto-context iteration 1
numofchannel = 1 #for auto-context iteration 0
#numofchannel = 2 #for auto-context iteration 1
###########

numofseg = 1

################################## 
pathoftest = []
pathoftest.append(datapath)
#pathoftest.append(resultpath+'Prob_Itr{0}_'.format(itr-1)) #for auto-context interation 1
###################################

sizeofchunk, sizeofchunk_expand = 96, 168
imagesize_predict = [sizeofchunk_expand, sizeofchunk_expand, sizeofchunk_expand]

if os.path.exists(modelpath+'m2edn_nmsk_iter{0}.best.hd5'.format(itr)):
    print 'loading weights ...'
    Model = get_m2edn(numofseg, numofchannel, imagesize_predict)
    Model.load_weights(modelpath+'m2edn_nmsk_iter{0}.best.hd5'.format(itr))


for idx_subject in sample_list:
    
    chunk_batch, nb_chunks, idx_xyz, sizeofimage = Generator_MultiChannel_Test_Chunks_cl(pathoftest, idx_subject,
                                        sizeofchunk, sizeofchunk_expand, numofchannel)
    seg_batch = Model.predict(chunk_batch, batch_size=1)
    prob_image = Chunks_Back_To_Image_cl(seg_batch, nb_chunks, sizeofchunk,
                                      sizeofchunk_expand, idx_xyz, sizeofimage)
    
    print 'Finish segmenting ' + idx_subject
    
    I = sitk.ReadImage(datapath+idx_subject)

    Heat_image = sitk.GetImageFromArray(prob_image, isVector=False)
    Heat_image.SetSpacing(I.GetSpacing())
    sitk.WriteImage(Heat_image,resultpath+'Prob_Itr{0}_'.format(itr)+idx_subject) 


