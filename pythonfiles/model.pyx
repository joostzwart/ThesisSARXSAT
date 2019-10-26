from __future__ import print_function
import numpy as np
cimport numpy as np
DTYPE = np.int
FTYPE=  np.double
ctypedef np.int_t DTYPE_t
ctypedef np.double_t FTYPE_t

cpdef simulate(int T, int dw, str example,float nt,int seed,int inputs,int outputs):
    """This function simulates a switched affine system. The types of systems are taken from
    the literature. Check the readme for references (deprecated).
    
    Arguments:
        T {integer} -- Length of dataset
        dw {integer} -- Dwell time
        example {string} -- Type of model
        nt {float} -- Bound on the noise
        seed {integer} -- Seed for the randomness
        inputs {integer} -- Number of inputs
        outputs {integer} -- Number of outputs
    
    Returns:
        y{list} -- Output
        r{list} -- Regressor
        n{list} -- Switching sequence
        u{list} -- Input
    """

    ## If a seed is given use this seed
    if seed !=0:
        np.random.seed(seed)

    ## Initializing order of the input and output
    cdef int ny=2
    cdef int nu=1

    ## Type casting
    cdef int i, x
    cdef list n
    cdef np.ndarray thetaS, ux, yx


    ## Generating uniformly distributed random input
    cdef np.ndarray u=2*np.random.rand(inputs,T)-1
    assert u.dtype == FTYPE

    
    ## Generating uniformly distributed random y0
    cdef np.ndarray y=np.zeros((outputs,T))
    y[:,:ny]=2*np.random.rand(outputs,ny)-1
    assert y.dtype == FTYPE
