import numpy as np
cimport numpy as np
cdef class Dataset:
    """Class for keeping track of information in a dataset"""

    cdef public double[:,:] y, u, r
    cdef public list certificate
    cdef public list switches
    cdef public list model
    cdef public list Nl
    cdef public list modelfinal
    cdef public list dictionary
    cdef public list Sigma
    cdef public list n
    cdef public list pwarx_certificate
    cdef public int T

    def __init__(self,double[:,:] output,double[:,:] regressor,double[:,:] input):
        self.r=regressor
        self.y=output
        self.u=input
        self.n=[]
        self.certificate=[]
        self.pwarx_certificate=[]
        self.switches=[]
        self.model=[]
        self.Nl=[]
        self.modelfinal=[]
        self.dictionary=[]
        self.Sigma=[]
        self.T=len(output[0])

    def printData(self):
        print("TESTING")
        print("Regressor = ",self.r)
        print("Output = ",self.y)
        print("Input = ",self.u)