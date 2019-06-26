#!/usr/bin/python
# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 17:33:39 2017
    
@author: chlian
"""
import sys
import numpy as np
import SimpleITK as sitk
import os
import fnmatch
import keras
import time
os.environ["CUDA_DEVICE_ORDER"]="PCI_BUS_ID";
# The GPU id to use, usually either "0" or "1";
os.environ["CUDA_VISIBLE_DEVICES"]="1";  

from M2EDN import get_m2edn
from Generate_Test_Chunks import Generator_MultiChannel_Test_Chunks_cl
from Generate_Test_Chunks import Chunks_Back_To_Image_cl
from keras import backend as K



if sys.argv[1][-4:]=='.hd5':
   model=sys.argv[1]
   sample_list=sys.argv[2:]
else:
   model = '/home/zongx/pvs_m2edn/saved_model/m2edn_nmsk_iter0.best.hd5'
   sample_list=sys.argv[1:]


resultpath = '/home/zongx/pvs_m2edn/saved_result/'


itr = 0 #for auto-context iteration 0
numofchannel = 1 #for auto-context iteration 0

numofseg = 1


sizeofchunk, sizeofchunk_expand = 96, 168
imagesize_predict = [sizeofchunk_expand, sizeofchunk_expand, sizeofchunk_expand]

drive,name_model=os.path.split(model)
name_model=name_model.split('.')
print 'loading weights from ' + model
Model = get_m2edn(numofseg, numofchannel, imagesize_predict)

Model.load_weights(model)

for idx_subject in sample_list:
    start=time.time()
    print 'Start segmenting ' + idx_subject
    chunk_batch, nb_chunks, idx_xyz, sizeofimage = Generator_MultiChannel_Test_Chunks_cl(idx_subject,
                                        sizeofchunk, sizeofchunk_expand, numofchannel)
    seg_batch = Model.predict(chunk_batch, batch_size=1)
    prob_image = Chunks_Back_To_Image_cl(seg_batch, nb_chunks, sizeofchunk,
                                      sizeofchunk_expand, idx_xyz, sizeofimage)
    print 'Finish segmenting ' + idx_subject
    
    I = sitk.ReadImage(idx_subject)
    Heat_image = sitk.GetImageFromArray(prob_image, isVector=False)
    Heat_image.SetSpacing(I.GetSpacing())
    Heat_image.SetOrigin(I.GetOrigin())
    Heat_image.SetDirection(I.GetDirection())
    drive,prefix=os.path.split(idx_subject)
    sitk.WriteImage(Heat_image,os.path.join(resultpath,'Prob_'+name_model[0]+'_'+prefix)) 
    end=time.time()
    print 'Total time = %d s\n'%int(end-start)

