from __future__ import print_function
import numpy as np
cimport numpy as np
DTYPE = np.int
FTYPE=  np.double
ctypedef np.int_t DTYPE_t
ctypedef np.double_t FTYPE_t

cpdef simulate(int T, int dw, str example,float nt,int seed,int inputs,int outputs):
    """This function simulates a switched affine system. The types of systems are taken from
    the literature. Check the readme for references.
    
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


    '''
    if example=='Bako':                                                                                         
    ## Simulating model(Bako)

        ## setting parameters
        ny=2
        nu=1
        [y]=np.random.rand(1,ny)
        thetaS=np.array([[-0.4,0.25,-0.15,0.08],[1.55,-0.58,-2.1,0.96],[1,-0.24,-0.65,0.3]])

        ## Choosing active mode 
        n=np.random.randint(1,4,int(np.floor(T/dw)))
        n=np.repeat(n,dw)
        n=np.concatenate((n, np.ones(int(T-np.floor(T/dw)*dw))), axis=0)
        n=[1]*10+[3]*10

        ## Simulate model
        for x in range(ny,T):   
            if n[x]==1:
                y=np.append(y,np.dot(np.append(y[x-ny:],u[x-nu:x]),thetaS[0]))
            elif n[x]==2:
                y=np.append(y,np.dot(np.append(y[x-ny:],u[x-nu:x]),thetaS[1]))
            else:
                y=np.append(y,np.dot(np.append(y[x-ny:],u[x-nu:x]),thetaS[2]))
    
    if example=='Ozay':
    ## Simulating model(Ozay)
        
        ##Setting Parameters
        thetaS=np.array([[0.24,0.2,2],[-0.53,-1.4,1]])
        ny=2
        nu=1
        [y]=np.random.rand(1,ny)
        
        ## Choosing active mode
        n=[1]*25+[2]*25+[1]*25+[2]*25
        


        ## Simulating model
        for x in range(ny,T):
            if n[x]==1:
                y=np.append(y,np.dot(np.append(y[x-ny:],u[x-nu:x]),thetaS[0])+float(nt*(np.random.rand(1)-0.5)))
            else:
                y=np.append(y,np.dot(np.append(y[x-ny:],u[x-nu:x]),thetaS[1])+float(nt*(np.random.rand(1)-0.5)))
    '''
    if example=='OzayMIMO':
    ## Simulating model(Ozay)
               
        ##Setting Parameters

        """SISO ny=2 nu=1"""
        thetaS=np.array([[0.2,0.6,2],[0.1,0.26,1.5],[-0.2,0.6,1]],dtype=FTYPE)

        """MIMO outputs=2 inputs=2 ny=2 nu=1"""
        thetaS=np.array([[[0.3,-0.26,-0.1,0.2,1.5,2],[0.4,-2,-0.2,0.4,1.5,2]],\
                         [[0.1,0.2,0.3,-0.2,2,3],[0.2,0.25,0.3,-0.25,2,3]],\
                        [[0.05,-0.15,-0.05,0.4,1.5,2],[0.2,-0.1,-0.3,0.1,1.5,2]]],dtype=FTYPE)


        ## Choosing active mode
        n=[1]*10+[2]*10+[3]*10+[1]*10+[2]*10+[3]*10+[1]*10+[2]*10+[3]*10+[1]*10#+[2]*30+[3]*30
        r=np.empty((T-ny,nu*inputs+ny*outputs), dtype=FTYPE)

        ## Simulating model
        for x in range(ny,T):
            ## Obtain the input and output values used to calculate the next time step
            ux=np.reshape(u[:,x-nu:x],nu*inputs)
            yx=np.reshape(y[:,x-ny:x],ny*outputs)
            r[x-ny,:]=np.append(yx,ux)

            ## Simulate the model
            y[:,x]=np.dot(thetaS[n[x]-1],np.append(yx,ux))+float(nt*(np.random.rand(1)-0.5))
            

    y=y[:,ny:]
    n=n[ny:]
    return y,r,n,u,thetaS