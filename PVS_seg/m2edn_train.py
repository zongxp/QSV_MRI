#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 17:33:39 2017
    
@author: chlian
"""
import time
import numpy as np
import SimpleITK as sitk
import os
import keras


os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID";
# The GPU id to use, usually either "0" or "1";
os.environ["CUDA_VISIBLE_DEVICES"]="1";  
from M2EDN import get_m2edn
from Generate_Test_Chunks import Generator_MultiChannel_Test_Chunks_cl
from Generate_Test_Chunks import Chunks_Back_To_Image_cl
from Generate_Train_Chunks import Generator_MultiChannel_Train_Chunks_cl
from keras import backend as K

datapath = 't2_train/' #path of input data
segpath= 'seg_train/'

id=[1,2,3,4,5,7,9,10,11,12,14,15,16,17,18,20,21,22,23,24,25,27,28,30,32,33,34,35,36,37,38,39,40,41,42,43,45,46]
modelpath = 'saved_model/'

###########
itr = 0 #for auto-context iteration 0
#itr = 1 #for auto-context iteration 1
numofchannel = 1 #for auto-context iteration 0
#numofchannel = 2 #for auto-context iteration 1
###########

numofseg = 1

################################## 
patht2 = []
patht2.append(datapath)

pathseg = []
pathseg.append(segpath)

batch_size=4;
sizeofchunk, sizeofchunk_expand = 96, 96
imagesize_predict = [sizeofchunk_expand, sizeofchunk_expand, sizeofchunk_expand]

if os.path.exists(modelpath+'m2edn_nmsk_iter{0}.best.hd5'.format(itr)):
    print 'loading weights ...'
    Model = get_m2edn(numofseg, numofchannel, imagesize_predict)
    #Model.load_weights(modelpath+'m2edn_nmsk_iter{0}.best.hd5'.format(itr))
    Model.load_weights(modelpath+'m2edn_retrain_fromModel_newcost_alpha10.hd5')

t2_name=[]
seg_name=[]
for idx in id:
    t2_name.append('T2_PVS%02d.nii.gz'%idx)
    seg_name.append('PVSLength_Prob_Itr1_mask_PVS%02d_koji2_Curv_Dv_stat.nii.gz'%idx)

gen= Generator_MultiChannel_Train_Chunks_cl(datapath, t2_name,segpath,seg_name,sizeofchunk, batch_size)

a=time.time()
Model.fit_generator(gen,steps_per_epoch=1026,epochs=2)
#Model.fit_generator(gen,steps_per_epoch=27,epochs=20,max_queue_size=1)

Model.save('m2edn_retrain_fromModel_newcost_alpha10_4epoch.hd5')
b=time.time()
print ('total time = %d s'%int(b-a))




