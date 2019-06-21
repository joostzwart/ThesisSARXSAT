import numpy as np
import Sparsification as sp
#import SAT
#import logging
#import hitting_set_L0
#import time
#from random import randint
import storage
import dataset
import DataGenerator
import random
import parameters
import ReaderData
from PyQt5.QtWidgets import QApplication, QWidget, QPushButton, QVBoxLayout

T=90
ny=1
nu=1
delta=0.5
dw=2
inputs=1
outputs=1
nod=1
nt=0.01
modes=3
blocks=70
seed=random.randint(1,100)
## Preferences
merging=0
splitlargedataset=1
chuncks=5
certtype=2
modelgeneration=2
theta=[]
H={}
pwarx=1
"""
(input, output, inputs, outputs, T, r)=ReaderData.Read(nu,ny,modelgeneration)
switches=[7, 9, 11, 13, 15, 17, 19, 22, 25, 26, 28, 31, 33, 35, 37, 39, 41, 43, 47, 48, 52, 53, 55, 56, 58, 59, 62, 64, 67, 68, 70, 72, 73, 75, 78, 80, 81, 83, 85, 87, 89, 91, 93, 95, 97, 99]
b={}
c={}
a=0
for i, k in enumerate(switches[0:]):
    b[i] = np.take(r, range(a, k), 0)
    c[i] = np.take(output, range(a, k), 1)
    a = k

u2=np.split(input,switches,axis=1)
y2=np.split(output,switches,axis=1)
r2=np.split(r,switches,axis=0)

print("r2",r2)
print("b[0]",b[0])
print("b[1]",b[1])
print("y2shape",output.shape)"""
Nl=[[1,2,3],[],[],[4,5],[],[6,7]]
model=["m1","m2","m3","m4","m5","m6"]
storage.action("create")
