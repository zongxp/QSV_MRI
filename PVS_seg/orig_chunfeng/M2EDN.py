# -*- coding: utf-8 -*-
"""
Created on Fri Aug  4 10:47:39 2017

@author: chlian
"""

from keras.models import Model
from keras import backend as K

from keras.layers import Input
from keras.layers.convolutional import Conv3D, UpSampling3D
from keras.layers.pooling import MaxPooling3D
from keras.layers.merge import Concatenate, Add, concatenate


# F_{\beta}-Measure loss
def fmeasure_loss(y_true, y_pred):

    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    
    alpha = 1.0 # alpha = beta**2

    tp = K.sum(y_pred_f * y_true_f)

    FL = 1- ((1.0 + alpha) * tp + K.epsilon()) / \
         (alpha * K.sum(y_true_f) + K.sum(y_pred_f) + K.epsilon())

    return FL

def sen(y_true, y_pred):
    
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    
    tp = K.sum(y_pred_f * y_true_f)
    re = (tp + K.epsilon()) / (K.sum(y_true_f)+K.epsilon())

    return re

def ppv(y_true, y_pred):
    
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    
    tp = K.sum(y_pred_f * y_true_f)
    re = (tp + K.epsilon()) / (K.sum(y_pred_f) + K.epsilon())
    
    return re

def dsc(y_true, y_pred):
    
    y_true_f = K.flatten(y_true)
    y_pred_f = K.flatten(y_pred)
    
    tp = K.sum(y_pred_f * y_true_f)
    re = (2 * tp + K.epsilon()) / \
         (K.sum(y_true_f) + K.sum(y_pred_f) + K.epsilon())

    return re

def binary_cross_entropy(y_true, y_pred):
    y_true = K.flatten(y_true)
    y_pred = K.flatten(y_pred)
    return K.mean(K.binary_crossentropy(y_true, y_pred))

def mse(y_true, y_pred):
    y_true = K.flatten(y_true)
    y_pred = K.flatten(y_pred)
    return K.mean(K.square(y_pred - y_true))

"""
M$^2$EDN
"""
def get_m2edn(numoflandmarks,inputchannel,imagesize):
    
    inputs = Input((imagesize[0], imagesize[1], imagesize[2], inputchannel))
    
    
    ###
    conv1 = Conv3D(filters=64, kernel_size=[3, 3, 3],
                   activation='relu', padding='same', data_format='channels_last')(inputs)
                   
    pool1 = MaxPooling3D(pool_size=(2, 2, 2), 
                         data_format='channels_last')(conv1)
    
    down1 = MaxPooling3D(pool_size=(2, 2, 2), 
                         data_format='channels_last')(inputs)
    
    conv_down1 = Conv3D(filters=64, kernel_size=[3, 3, 3],
                        activation='relu', padding='same', data_format='channels_last')(down1)
                        
    concat1 = Concatenate(axis=-1)([conv_down1, pool1])
    
    ###
    conv2 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                   activation='relu', padding='same', data_format='channels_last')(concat1)
                   
    pool2 = MaxPooling3D(pool_size=(2, 2, 2), 
                         data_format='channels_last')(conv2)
    
    down2 = MaxPooling3D(pool_size=(4, 4, 4), strides=(4,4,4),
                         data_format='channels_last')(inputs)
    
    conv_down2 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                        activation='relu', padding='same', data_format='channels_last')(down2)
                        
    concat2 = Concatenate(axis=-1)([conv_down2, pool2])

    
    ###
    conv3 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                   activation='relu', padding='same', data_format='channels_last')(concat2)
                   
    pool3 = MaxPooling3D(pool_size=(2, 2, 2), 
                         data_format='channels_last')(conv3)
    
    conv4 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                   activation='relu', padding='same', data_format='channels_last')(pool3)
    
    
    ###
    up1 = UpSampling3D(size=(2, 2, 2), data_format='channels_last')(conv4)

    concat_up1 = Concatenate(axis=-1)([up1,conv3])
    
    conv5 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                   activation='relu', padding='same', data_format='channels_last')(concat_up1)
                   
    
    ###
    up2 = UpSampling3D(size=(2, 2, 2), data_format='channels_last')(conv5)

    concat_up2 = Concatenate(axis=-1)([up2, conv2])
    
    conv6 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                   activation='relu', padding='same', data_format='channels_last')(concat_up2)
                   
    
    ###
    up3 = UpSampling3D(size=(2, 2, 2), data_format='channels_last')(conv6)

    concat_up3 = Concatenate(axis=-1)([up3, conv1])
    
    conv7 = Conv3D(filters=64, kernel_size=[3, 3, 3], 
                   activation='relu', padding='same', data_format='channels_last')(concat_up3)
    
    
    ###
    conv8 = Conv3D(filters=numoflandmarks, kernel_size=[1, 1, 1],
                   activation='sigmoid', data_format='channels_last')(conv7)
    
    
    ###
    model = Model(inputs=inputs, outputs=conv8)
    
    model.compile(optimizer='Adam', 
                  metrics=[dsc, sen, ppv], 
                  loss=fmeasure_loss)
    
    return model


if __name__ == '__main__':
	model = get_unet_deeper(numoflandmarks=1,inputchannel=2,imagesize=[168,168,168])
	model.summary()

    
