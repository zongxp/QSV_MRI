# -*- coding: utf-8 -*-
"""
Created on Thu Mar 23 15:53:09 2017

@author: chlian
"""

          
import SimpleITK as sitk
import numpy as np



def Generator_MultiChannel_Train_Chunks_cl(pathofchannels, T2_name, pathofseg,seg_name,sizeofchunk, 
                                            batchsize,ave_im=0, std_im=1.0):
    
 while True: 
    for ifile in range(len(T2_name)):
        image = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofchannels+T2_name[ifile])))
        ave_im, std_im = np.mean(image.astype('float32')), np.std(image.astype('float32'))
        inputs=((image.astype('float32')-ave_im)/std_im)
        image2 = np.array(sitk.GetArrayFromImage(sitk.ReadImage(pathofseg+seg_name[ifile])))
        inputs2=image2>0


        sizeofimage = np.shape(inputs)    
        nb_chunks = (np.ceil(np.array(sizeofimage)/float(sizeofchunk))).astype(int)


        pad_inputs=np.zeros(nb_chunks*sizeofchunk, dtype='float32')
        pad_inputs[:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = inputs
        pad_inputs2=np.zeros(nb_chunks*sizeofchunk, dtype='float32')
        pad_inputs2[:sizeofimage[0], :sizeofimage[1], :sizeofimage[2]] = inputs2

        numofchannel=1;
        i_channel=0;
        chunk_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofchannel),dtype='float32')
        seg_batch = np.zeros((batchsize,sizeofchunk,sizeofchunk,sizeofchunk,numofchannel),dtype='float32')
        ibatch=0
        for x_idx in range(nb_chunks[0]):
            for y_idx in range(nb_chunks[1]):
                for z_idx in range(nb_chunks[2]):
                        chunk_batch[ibatch,...,i_channel] = pad_inputs[x_idx*sizeofchunk:x_idx*sizeofchunk+sizeofchunk,\
                          y_idx*sizeofchunk:y_idx*sizeofchunk+sizeofchunk,\
                          z_idx*sizeofchunk:z_idx*sizeofchunk+sizeofchunk]   
                        seg_batch[ibatch,...,i_channel] = pad_inputs2[x_idx*sizeofchunk:x_idx*sizeofchunk+sizeofchunk,\
                          y_idx*sizeofchunk:y_idx*sizeofchunk+sizeofchunk,\
                          z_idx*sizeofchunk:z_idx*sizeofchunk+sizeofchunk]   
                        if ibatch==batchsize-1:  
                          ibatch=0
                          yield chunk_batch,seg_batch
                        else:
                          ibatch=ibatch+1

                
#   return chunk_batch, nb_chunks, idx_xyz, sizeofimage



