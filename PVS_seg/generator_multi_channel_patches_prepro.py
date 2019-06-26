# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 17:42:08 2017

@author: chlian
"""

import numpy as np
import SimpleITK as sitk
import random

def Generator_MultiChannel_Patches_Balance(pathofchannels,pathofreference,samplelist,batchsize,sizeofpatch,numofseg,numofchannel,ave_im=0,std_im=1.0):
    
    while True:
        random.shuffle(samplelist)

        patch_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofchannel),dtype='float32')
        reference_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofseg),dtype='uint8')
        i_batch = 0

        for i_subject in samplelist:
            
            inputs = []
            for i_channel in range(numofchannel):
                if i_channel == 0:
                    image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject))))
                    if ave_im == 0:
                        ave_im, std_im = np.mean(image.astype('float32')), np.std(image.astype('float32'))
                    inputs.append((image.astype('float32')-ave_im)/std_im)
                else:
                    inputs.append(np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject)))))
            reference = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofreference+'Manual_{0}_alex.hdr'.format(i_subject))))
            
            where_pvs = np.where(reference > 0)
            cor_pvs = np.array([np.min(np.array(where_pvs),axis=1),np.max(np.array(where_pvs),axis=1)])            
            
            # random part
            sizeofimage = np.shape(reference)
            shiftrange_random = np.array(sizeofimage) - sizeofpatch
            cor_shiftrange_random = []
            for i_cor in range(3):
                if shiftrange_random[i_cor] <= sizeofpatch:
                    cor_shiftrange_random.append(range(np.random.randint(shiftrange_random[i_cor]),shiftrange_random[i_cor],sizeofpatch/2))
                else:
                    cor_shiftrange_random.append(range(np.random.randint(shiftrange_random[i_cor]%sizeofpatch+1),shiftrange_random[i_cor],sizeofpatch/2))
                    
            # balance part
            shiftrange_balance = (cor_pvs[1,:] - cor_pvs[0,:]) - sizeofpatch
            cor_shiftrange_balance = []
            for i_cor in range(3):
                if (cor_pvs[1,i_cor] - cor_pvs[0,i_cor]) <= sizeofpatch:
                    step_sampling = 10
                    cor_shiftrange_balance.append(range(np.random.randint(cor_pvs[0,i_cor]-step_sampling,cor_pvs[0,i_cor]),np.random.randint(cor_pvs[1,i_cor],cor_pvs[1,i_cor]+step_sampling),sizeofpatch))
                elif shiftrange_balance[i_cor] <= sizeofpatch:
                    cor_shiftrange_balance.append(range(cor_pvs[0,i_cor]+np.random.randint(shiftrange_balance[i_cor]),cor_pvs[0,i_cor]+shiftrange_balance[i_cor],sizeofpatch/6))
                else:
                    cor_shiftrange_balance.append(range(cor_pvs[0,i_cor]+np.random.randint(shiftrange_balance[i_cor]%sizeofpatch+1),cor_pvs[0,i_cor]+shiftrange_balance[i_cor],sizeofpatch/4))
           
            # generator 
            for x_begin in cor_shiftrange_balance[0]:
                for y_begin in cor_shiftrange_balance[1]:
                    for z_begin in cor_shiftrange_balance[2]:
                        
                        for i_channel in range(numofchannel):
                            patch_batch[i_batch,...,i_channel] = inputs[i_channel][x_begin:x_begin+sizeofpatch,y_begin:y_begin+sizeofpatch,z_begin:z_begin+sizeofpatch]
                        reference_batch[i_batch,...,0] = reference[x_begin:x_begin+sizeofpatch,y_begin:y_begin+sizeofpatch,z_begin:z_begin+sizeofpatch]
                       
                        i_batch += 1
                       
                        if i_batch == batchsize:
                             yield (patch_batch, reference_batch)
                             patch_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofchannel),dtype='float32')
                             reference_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofseg),dtype='uint8')
                             i_batch = 0
                        
                        x_random = cor_shiftrange_random[0][random.randrange(len(cor_shiftrange_random[0]))]
                        y_random = cor_shiftrange_random[1][random.randrange(len(cor_shiftrange_random[1]))]
                        z_random = cor_shiftrange_random[2][random.randrange(len(cor_shiftrange_random[2]))]
                       
                        for i_channel in range(numofchannel):
                            patch_batch[i_batch,...,i_channel] = inputs[i_channel][x_random:x_random+sizeofpatch,y_random:y_random+sizeofpatch,z_random:z_random+sizeofpatch]
                        reference_batch[i_batch,...,0] = reference[x_random:x_random+sizeofpatch,y_random:y_random+sizeofpatch,z_random:z_random+sizeofpatch]
                       
                        i_batch += 1
                       
                        if i_batch == batchsize:
                             yield (patch_batch, reference_batch)
                             patch_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofchannel),dtype='float32')
                             reference_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofseg),dtype='uint8')
                             i_batch = 0



def Generator_MultiChannel_Patches_Uniform(pathofchannels,pathofreference,samplelist,batchsize,sizeofchunk,numofseg,numofchannel):

    while True:
        random.shuffle(samplelist)
        
        patch_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofchannel),dtype='float32')
        reference_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofseg),dtype='uint8')
        i_batch = 0
        
        for i_subject in samplelist:
            inputs = []
            for i_channel in range(numofchannel):
                image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject))))
                ave_im, std_im = np.mean(image.astype('float32')), np.std(image.astype('float32'))
                inputs.append((image.astype('float32')-ave_im)/std_im)
            reference = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofreference+'Manual_{0}.hdr'.format(i_subject))))

            sizeofimage = np.shape(reference)    
            nb_chunks = (np.ceil(np.array(sizeofimage)/float(sizeofchunk))).astype(int)

            pad_inputs = []
            for i_channel in range(numofchannel):
                pad_inputs.append(np.zeros(nb_chunks*sizeofchunk, dtype='float32'))
                pad_inputs[i_channel][:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = inputs[i_channel]
            pad_ref = np.zeros(nb_chunks*sizeofchunk, dtype='float32')
            pad_ref[:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = reference

            # generator
            for x_idx in range(nb_chunks[0]):
                for y_idx in range(nb_chunks[1]):
                    for z_idx in range(nb_chunks[2]):
                        for i_channel in range(numofchannel):
                            patch_batch[i_batch,...,i_channel] = \
                                pad_inputs[i_channel][x_idx*sizeofchunk:x_idx*sizeofchunk+sizeofchunk,
                                                      y_idx*sizeofchunk:y_idx*sizeofchunk+sizeofchunk,
                                                      z_idx*sizeofchunk:z_idx*sizeofchunk+sizeofchunk]
                        reference_batch[i_batch,...,0] = pad_ref[x_idx*sizeofchunk:x_idx*sizeofchunk+sizeofchunk,
                                                                 y_idx*sizeofchunk:y_idx*sizeofchunk+sizeofchunk,
                                                                 z_idx*sizeofchunk:z_idx*sizeofchunk+sizeofchunk]
                        i_batch += 1
                        if i_batch == batchsize:
                             yield (patch_batch, reference_batch)
                             patch_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofchannel),dtype='float32')
                             reference_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofseg),dtype='uint8')
                             i_batch = 0



def Generator_MultiChannel_Patches_Uniform2(pathofchannels,pathofreference,samplelist,batchsize,sizeofchunk,numofseg,numofchannel):

    while True:
        random.shuffle(samplelist)
        
        patch_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofchannel),dtype='float32')
        reference_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofseg),dtype='uint8')
        i_batch = 0

        for i_subject in samplelist:
            inputs = []
            for i_channel in range(numofchannel):
                if i_channel == 0:
                    image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject))))
                    ave_im, std_im = np.mean(image.astype('float32')), np.std(image.astype('float32'))
                    inputs.append((image.astype('float32')-ave_im)/std_im)
                else:
                    inputs.append(np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject)))))
            reference = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofreference+'Manual_{0}.hdr'.format(i_subject))))

            sizeofimage = np.shape(reference)    
            nb_chunks = (np.ceil(np.array(sizeofimage)/float(sizeofchunk))).astype(int)

            pad_inputs = []
            for i_channel in range(numofchannel):
                pad_inputs.append(np.zeros(nb_chunks*sizeofchunk, dtype='float32'))
            for i_channel in range(numofchannel):
                pad_inputs[i_channel][:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = inputs[i_channel]
            pad_ref = np.zeros(nb_chunks*sizeofchunk, dtype='uint8')
            pad_ref[:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = reference

            # generator
            for x_idx in range(nb_chunks[0]):
                for y_idx in range(nb_chunks[1]):
                    for z_idx in range(nb_chunks[2]):
                        for i_channel in range(numofchannel):
                            patch_batch[i_batch,...,i_channel] = \
                                pad_inputs[i_channel][x_idx*sizeofchunk:(x_idx+1)*sizeofchunk,
                                                      y_idx*sizeofchunk:(y_idx+1)*sizeofchunk,
                                                      z_idx*sizeofchunk:(z_idx+1)*sizeofchunk]
                        reference_batch[i_batch,...,0] = pad_ref[x_idx*sizeofchunk:(x_idx+1)*sizeofchunk,
                                                                 y_idx*sizeofchunk:(y_idx+1)*sizeofchunk,
                                                                 z_idx*sizeofchunk:(z_idx+1)*sizeofchunk]
                        i_batch += 1
                        if i_batch == batchsize:
                            yield (patch_batch, reference_batch)
                            patch_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofchannel),dtype='float32')
                            reference_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofseg),dtype='uint8')
                            i_batch = 0



def Generator_MultiChannel_Patches_Random(pathofchannels,pathofreference,samplelist,batchsize,sizeofpatch,numofseg,numofchannel,ave_im,std_im):
    
    while True:
        random.shuffle(samplelist)
        
        patch_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofchannel),dtype='float32')
        reference_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofseg),dtype='uint8')
        i_batch = 0
        
        for i_subject in samplelist:
            
            inputs = []
            for i_channel in range(numofchannel):
                if i_channel == 0:
                    image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject))))
                    if ave_im == 0:
                        ave_im, std_im = np.mean(image.astype('float32')), np.std(image.astype('float32'))
                    inputs.append((image.astype('float32')-ave_im)/std_im)
                else:
                    inputs.append(np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels[i_channel]+'{0}.hdr'.format(i_subject)))))
            reference = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofreference+'Manual_{0}.hdr'.format(i_subject))))

            sizeofimage = np.shape(reference)
            shiftrange_random = np.array(sizeofimage) - sizeofpatch
            cor_shiftrange_random = []
            for i_cor in range(3):
                if shiftrange_random[i_cor] <= sizeofpatch:
                    cor_shiftrange_random.append(range(np.random.randint(shiftrange_random[i_cor]),shiftrange_random[i_cor],sizeofpatch/6))
                else:
                    cor_shiftrange_random.append(range(np.random.randint(shiftrange_random[i_cor]%sizeofpatch+1),shiftrange_random[i_cor],sizeofpatch/4))

    
            # generator
            for x_begin in cor_shiftrange_random[0]:
                for y_begin in cor_shiftrange_random[1]:
                    for z_begin in cor_shiftrange_random[2]:
                
                        for i_channel in range(numofchannel):
                            patch_batch[i_batch,...,i_channel] = inputs[i_channel][x_begin:x_begin+sizeofpatch,y_begin:y_begin+sizeofpatch,z_begin:z_begin+sizeofpatch]
                        reference_batch[i_batch,...,0] = reference[x_begin:x_begin+sizeofpatch,y_begin:y_begin+sizeofpatch,z_begin:z_begin+sizeofpatch]
                        
                        i_batch += 1
                        
                        if i_batch == batchsize:
                            yield (patch_batch, reference_batch)
                            patch_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofchannel),dtype='float32')
                            reference_batch = np.zeros((batchsize,sizeofpatch,sizeofpatch,sizeofpatch,numofseg),dtype='uint8')
                            i_batch = 0

