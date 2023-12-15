# -*- coding: utf-8 -*-
"""
Created on Mon Jun 19 08:09:40 2023

@author: likun
"""
'''
# A neural network model using Oridinal Encoding. You need to specify the input file on the following line:
line 63: x1,y1=getdate('E:/imm/file5.csv',1) # autoantigen epitopes
line 64: x2,y2=getdate('E:/imm/file3.csv',0) #pathogen epitopes
file3 and file5 can be found in the https://github.com/LKbiosoup/epitope-core 
line 120 plt.savefig('out_picture1.png',dpi=600)
line 131 plt.savefig('out_picture1.png',dpi=600) #Two output image directories that need to be specified
'''

import tensorflow as tf
from tensorflow import keras
import numpy as np
import math
import random
from tensorflow.python.client import device_lib
import matplotlib.pyplot as plt
def cl(seq):
    seq=seq.replace('"','')
    seq=seq.split('+')[0]
    seq=seq.strip()
    return seq
def getdate(file,num):
    am='XGAPVLISTQNMCDERKHFYW'
    am=list(am)
    lisx=[];lisy=[]
   
    with open(file,'r') as f:
            lines=f.readlines()
    f.close()
  
        
    for i in range(len(lines)):
        if i%1==0:#
            strs=cl(lines[i])
            strs=strs.strip()
            try:
                if len(strs) <20:
                    strs=strs+(20-len(strs))*'X'
                    #print(strs)
                    datax=[am.index(i) for i in strs]
                    print(datax)
                    datax2=[]
                 
                    for x in datax:
                        data0=21*[0]
                        data0[x]=1
                        datax2.append(data0)
                        
                    lisx.append(datax)#
                    lisy.append(num)
            except:
                tag=1
    return lisx,lisy

print(tf.__version__)
x1,y1=getdate('E:/imm/file5.csv',1)
x2,y2=getdate('E:/imm/file3.csv',0)

len1=len(x1)
x1.extend(x2[:len1])
y1.extend(y2[:len1])
lisx=x1
lisy=y1

data = [(x, y) for x, y in zip(lisx, lisy)]
random.shuffle(data)
lisx = [d[0] for d in data]
lisy = [d[1] for d in data]

x_train=lisx[:math.ceil(len(lisx)*0.7)]
y_train=lisy[:math.ceil(len(lisy)*0.7)]
x_val=lisx[math.ceil(len(lisx)*0.7):]
y_val=lisy[math.ceil(len(lisy)*0.7):]

vocab_size = 21
input_length = 20
output_dim = 16 # 

model = keras.Sequential()
#model.add(keras.layers.Flatten(input_shape=(20, 21)))
model.add(keras.layers.Embedding(vocab_size, output_dim))
model.add(keras.layers.GlobalAveragePooling1D())#
model.add(keras.layers.Dense(16, activation='relu'))
model.add(keras.layers.Dense(1, activation='sigmoid'))
model.summary()
model.compile(optimizer='adam',
              loss='binary_crossentropy',
              metrics=['accuracy'])

history = model.fit(x_train,
                   y_train,
                    epochs=40,
                    batch_size=512,
                    validation_data=(x_val, y_val),
                    verbose=1)
history_dict = history.history
history_dict.keys()
acc = history_dict['accuracy']
val_acc = history_dict['val_accuracy']
loss = history_dict['loss']
val_loss = history_dict['val_loss']

epochs = range(1, len(acc) + 1)
plt.rc('font', family='Times New Roman')     
plt.plot(epochs, loss, 'bo', label='Training loss')
plt.plot(epochs, val_loss, 'b', label='Validation loss')
plt.title('Training and validation loss')
plt.xlabel('Epochs')
plt.ylabel('Loss')
plt.legend()
plt.savefig('out_picture1.png',dpi=600)
plt.show()

plt.clf()   

plt.plot(epochs, acc, 'bo', label='Training acc')
plt.plot(epochs, val_acc, 'b', label='Validation acc')
plt.title('Training and validation accuracy')
plt.xlabel('Epochs')
plt.ylabel('Accuracy')
plt.legend()
plt.savefig('out_picture1.png',dpi=600)
plt.show()

            
        
    