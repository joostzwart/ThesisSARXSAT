import L0
import parameters
import dataset
import numpy as np
cimport numpy as np
import ReaderData
import DataGenerator
import random
import block_operations
import logging
import storage
import time
import export

def main(self,**kwargs):
    """
    Main function of this application. Can be run independently or used via the GUI
    """
    ## type casting
    cdef int ny, nu, dw, blocks, splitlargedataset, inputs, outputs, modes,\
        T, modelgeneration, NOD, i
    cdef double delta, nt, error, mse
    cdef np.ndarray Yt, Ut, Rt, theta, output, input, r
    cdef list n, block, Y, U, R, Fmodels, sigma, swt

    ## Initialize parameters from gui or use static ones if no gui was used
    #if __name__ == '__main__':
    ny=2
    nu=1
    dw=15
    delta=0.1
    blocks = 2
    splitlargedataset = 1
    inputs = 2
    outputs = 1
    modes = 3
    nt = 0.01
    T=180
    seed=random.randint(1,100)
    seed=95
    modelgeneration=2

    ## Converting kwargs to locals
    if 'ny' in kwargs:
        ny=kwargs['ny']
    if 'nu' in kwargs:
        nu=kwargs['nu']
    if 'dw' in kwargs:
        dw=kwargs['dw']
    if 'delta' in kwargs:
        delta=kwargs['delta']
    if 'blocks' in kwargs:
        blocks=kwargs['blocks']
    if 'split' in kwargs:
        splitlargedataset = kwargs['split']
    if 'models' in kwargs:
        modes=kwargs['models']
    if 'outputs' in kwargs:
        outputs = kwargs['outputs']
    if 'inputs' in kwargs:
        inputs = kwargs['inputs']
    if 'nt' in kwargs:
        nt = kwargs['nt']
    if 'T' in kwargs:
        T = kwargs['T']
    if 'seed' in kwargs:
        seed=kwargs['seed']
    if 'pwarx' in kwargs:
        pwarx=kwargs['pwarx']
    if 'modelgeneration' in kwargs:
        modelgeneration = kwargs['modelgeneration']

    ## Parameters
    nod=1

    ## Preferences
    cdef int merging=0
    cdef int chuncks=1
    cdef int unstuck=2 ## Method of dealing with the L0 minimization when it gets stuck.
                     # 1 for starting with different weights.
                     # 2 for finding a model for the 3 intervals which have the smallest error.

    ## initialization, subfix "t" means total and differentiates between features of the total dataset,
    # and features for specific blocks
    oldmodels=set(())
    cdef list Run=[None] * blocks
    logging.basicConfig(level=10)

    t1=time.time()
    if modelgeneration==1:
        ## Create preference object
        Par=parameters.Parameters(T,ny,nu,delta,dw,nt,inputs,outputs,chuncks,oldmodels,nod)
        Pref=parameters.Preferences(merging,splitlargedataset,modelgeneration,seed,unstuck)
        L0.main(Par, Pref)

    if modelgeneration==2 or modelgeneration == 3:
        ## Simulate data to the file and read it again (seems superfluous but
        # necessary feature for importing textfiles in a later stadium)
        (theta,n)=DataGenerator.generate(T,ny,nu,inputs,outputs,nt,modes,dw,seed,pwarx)
        (input, output, inputs, outputs, T, r)=ReaderData.Read(nu,ny,modelgeneration)

        ## plot input and output when using gui
        if not __name__ == '__main__':
            self.ax.plot(output[0], 'r-')
            self.ax.plot(input[0], 'b-')
            self.ax.legend(['Output', 'Input'])

            for i in range(1,len(output)):
                self.ax.plot(output[i], 'r-')
            for i in range(1,len(input)):
                self.ax.plot(input[i],'b-')
            self.plotOutput.draw()

        ## Create preference object
        Par=parameters.Parameters(T,ny,nu,delta,dw,nt,inputs,outputs,chuncks,oldmodels,nod)
        Pref=parameters.Preferences(merging,splitlargedataset,modelgeneration,seed,unstuck)

        ## Some information for debugging purposes
        logging.debug("The seed that was used = {}".format(seed))
        logging.debug("Theta used for simulation = {}".format(theta))

        ## Split the Data
        Y=np.array_split(output,blocks,axis=1)
        U=np.array_split(input,blocks,axis=1)
        R=np.array_split(r,blocks,axis=0)

        ## Create a list to store the models per block
        block=[None] * blocks
        for i in range(blocks):
            block[i] = dataset.Dataset(Y[i], R[i], U[i])

        ## Run the solver for each data "block"
        for i in range(blocks):
            Run[i]=L0.main(Par, Pref, ProvidedDataset=block[i], theta=theta, n=n)
            if not __name__ == '__main__':
                self.progressbar.setValue(int((100/blocks)*(i+1)))
            for model in Run[i][0].modelfinal:
                Par.oldmodels.add(tuple(model))
            logging.debug("Nl for this run = {}".format(Run[i][0].Nl))

        ## Retrieve final switching sequence
        Yt=np.hstack(([Y[i] for i in range(blocks)]))
        Ut=np.hstack(([U[i] for i in range(blocks)]))
        Rt=np.vstack(([R[i] for i in range(blocks)]))
        Total_ds=dataset.Dataset(Yt,Rt,Ut)

        ## Retrieve models, switching sequence and switches
        (Fmodels,sigmat,swt) = block_operations.merge_blocks(Run,blocks)
        storage.datastorage(Fmodels, Par, 'HVAC', 'Full HVAC dataset', 0, 0, sigmat, 0)
        (error,mse,sigmat)=block_operations.hitting_set_model_reduction(input,output,r,nu,ny,inputs,outputs,swt,Fmodels,delta,theta,n,modelgeneration)
        t2=time.time()
        logging.info("Total time for identification process: {}".format(t2-t1))
        storage.datastorage(Fmodels, Par,'HVAC','Full HVAC dataset, with error correct',error,mse,sigmat,float(t2-t1))
        export.exportToMatlab(Fmodels,sigmat,Ut,Yt,Rt)

if __name__ == '__main__':
    main([])