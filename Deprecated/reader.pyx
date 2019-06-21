from __future__ import print_function

import numpy as np
cimport numpy as np
DTYPE = np.int
FTYPE=  np.double
ctypedef np.int_t DTYPE_t
ctypedef np.double_t FTYPE_t

## This file reads in a comma seperated file and simulates the model
cpdef read(int nu,int ny):
    ## Reading in the input
    data_file = open("Data/TestfileInput.txt", "r")
    lines = data_file.read().split("\n")
    cdef list input=[]
    for line in lines:
        L=line.split(",")
        del L[-1]
        input.append([float(i) for i in L])    
    data_file.close()
    input=np.array(input)

    ## Reading in the output
    data_file = open("Data/TestfileOutput.txt", "r")
    lines = data_file.read().split("\n")
    cdef list output=[]
    for line in lines:
        L=line.split(",")
        del L[-1]
        output.append([float(i) for i in L])    
    data_file.close()
    output=np.array(output)
    if output.ndim==1:
        output=output.flatten()

    ## Obtaining number of inputs/outputs
    inputs=len(input)
    outputs=len(output)
    T=len(output[0])

    ## Creating regressor
    r=[]
    for k in range(0,T-ny):
        ux=np.hstack(tuple([input[a][k-nu+ny:k+ny] for a in range(inputs)]))
        yx=np.hstack(tuple([output[a][k:k+ny] for a in range(outputs)]))
        r.append(np.append(yx,ux))
    r=np.array(r)
    if output.ndim==2:
        output=output[:,ny:]
    else:
        output=output[ny:]
    return(input,output,inputs,outputs,T,r)