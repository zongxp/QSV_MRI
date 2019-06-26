# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:53:09 2017

@author: chlian
"""

          
import SimpleITK as sitk
import numpy as np
import os


def Generator_MultiChannel_Test_Chunks_cl(idx_subject, sizeofchunk, 
                                          sizeofchunk_expand, numofchannel, ave_im=0, std_im=1.0):
 while True:   
    inputs = []
    for i_channel in range(numofchannel):
        if i_channel == 0:
            image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(idx_subject)))
            ave_im, std_im = np.mean(image.astype('float32')), np.std(image.astype('float32'))
            inputs.append((image.astype('float32')-ave_im)/std_im)
        else:
            image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(idx_subject)))
            inputs.append(image)
    
    sizeofimage = np.shape(inputs[0])    
    nb_chunks = (np.ceil(np.array(sizeofimage)/float(sizeofchunk))).astype(int)
    
    pad_inputs = []
    for i_channel in range(numofchannel):
        pad_inputs.append(np.zeros(nb_chunks*sizeofchunk, dtype='float32'))
        pad_inputs[i_channel][:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = inputs[i_channel]
        
    width = int(np.ceil((sizeofchunk_expand-sizeofchunk)/2.0))
    size_expand_im = np.array(np.shape(pad_inputs[0])) + 2*width

    expand_inputs = []
    for i_channel in range(numofchannel):
        expand_inputs.append(np.zeros(size_expand_im, dtype='float32'))
        expand_inputs[i_channel][width:-width, width:-width, width:-width] = pad_inputs[i_channel]

    batchsize = np.prod(nb_chunks)
    idx_chunk = 0
    chunk_batch = np.zeros((batchsize,sizeofchunk_expand,sizeofchunk_expand,sizeofchunk_expand,numofchannel),dtype='float32')

    
    idx_xyz = np.zeros((batchsize,3),dtype='uint16')
    for x_idx in range(nb_chunks[0]):
        for y_idx in range(nb_chunks[1]):
            for z_idx in range(nb_chunks[2]):
                
                idx_xyz[idx_chunk,:] = [x_idx,y_idx,z_idx]
                
                for i_channel in range(numofchannel):
                    chunk_batch[idx_chunk,...,i_channel] = expand_inputs[i_channel][x_idx*sizeofchunk:x_idx*sizeofchunk+sizeofchunk_expand,\
                      y_idx*sizeofchunk:y_idx*sizeofchunk+sizeofchunk_expand,\
                      z_idx*sizeofchunk:z_idx*sizeofchunk+sizeofchunk_expand]

                idx_chunk += 1
                    
    return chunk_batch, nb_chunks, idx_xyz, sizeofimage


def Chunks_Back_To_Image_cl(segment_chunks, nb_chunks, sizeofchunk, sizeofchunk_expand, idx_xyz, sizeofimage):
    
    batchsize = np.size(segment_chunks,0)
    width = int(np.ceil((sizeofchunk_expand-sizeofchunk)/2.0))
    segment_image = np.zeros((nb_chunks[0]*sizeofchunk,nb_chunks[1]*sizeofchunk,nb_chunks[2]*sizeofchunk))
    
    for idx_chunk in range(batchsize):
        
        idx_low = idx_xyz[idx_chunk,:] * sizeofchunk
        idx_upp = (idx_xyz[idx_chunk,:]+1) * sizeofchunk
        
        segment_image[idx_low[0]:idx_upp[0],idx_low[1]:idx_upp[1],idx_low[2]:idx_upp[2]] = \
        segment_chunks[idx_chunk,...,0][width:-width,width:-width,width:-width]
        
    segment_image = segment_image[:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]]
    return segment_image

