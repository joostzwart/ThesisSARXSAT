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
import partition
import matplotlib.pyplot as plt
from sys import stdout, stderr
from io import StringIO

def main(self,gui,**kwargs):
    """
    Main function of this application. Can be run independently or used via the GUI
    """
    ## type casting
    cdef int ny, nu, dw, blocks, splitlargedataset, inputs, outputs, modes,\
        T, modelgeneration, NOD, i, pwarx
    cdef double delta, nt, error, mse
    cdef np.ndarray Yt, Ut, Rt, theta, output, input, r
    cdef list n, block, Y, U, R, Fmodels, sigma, swt

    ## Initialize parameters from gui or use static ones if no gui was used
    #if __name__ == '__main__':
    Ls=1
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
    modelgeneration=2
    pwarx=0
    inputfile = "Data/DatasetInput_Ts_20(J).txt"
    outputfile = "Data/Datasetoutput_Ts_20(J).txt"
    inputtype = 1
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
    if 'outputfile' in kwargs:
        outputfile = kwargs['outputfile']
    if 'inputfile' in kwargs:
        inputfile = kwargs['inputfile']
    if 'Ls' in kwargs:
        Ls = kwargs['Ls']
    if 'inputtype' in kwargs:
        inputtype = kwargs['inputtype']


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

    ## Configure logger
    logging.info("initialize logger")
    baseLogger = logging.getLogger('dev')
    baseLogger.setLevel(20)
    strio = StringIO()
    Console = logging.StreamHandler(strio)
    baseLogger.addHandler(Console)
    #logger2=logging.getLogger('mainlogger')
    #logger2.setLevel(10)
    t1=time.time()
    if modelgeneration==1:
        ## Create preference object
        Par=parameters.Parameters(T,ny,nu,delta,dw,nt,inputs,outputs,chuncks,oldmodels,nod,pwarx,Ls)
        Pref=parameters.Preferences(merging,splitlargedataset,modelgeneration,seed,unstuck)
        L0.main(Par, Pref)

    if modelgeneration==2 or modelgeneration == 3:
        ## Simulate data to the file and read it again (seems superfluous but
        # necessary feature for importing textfiles in a later stadium)
        (theta,n,H)=DataGenerator.generate(T,ny,nu,inputs,outputs,nt,modes,dw,seed,pwarx,inputtype)
        baseLogger.info("Hyperplanes used: {}".format(H))
        (input, output, inputs, outputs, T, r)=ReaderData.Read(nu,ny,modelgeneration,input=inputfile,output=outputfile)

        ## extend regressor in case of pwarx
        if pwarx==1:
            extendedregressor=np.ones((r.shape[0],r.shape[1]+1))
            extendedregressor[:,:-1]=r
            r=extendedregressor

        ## plot input and output when using gui
        if not __name__ == '__main__' and gui and modelgeneration==3:
            self.ax.clear()
            self.ax.plot(output[0], 'r-')
            self.ax.plot(input[0][max(ny,nu):], 'b-')
            self.ax.legend(['Output', 'Input'])#

            for i in range(1,len(output)):
                self.ax.plot(output[i], 'r-')
            for i in range(1,len(input)):
                self.ax.plot(input[i][max(ny,nu):],'b-')
            self.plotOutput.draw()

        ## Create preference object
        Par=parameters.Parameters(T,ny,nu,delta,dw,nt,inputs,outputs,chuncks,oldmodels,nod,pwarx,Ls)
        Pref=parameters.Preferences(merging,splitlargedataset,modelgeneration,seed,unstuck)

        ## Some information for debugging purposes
        baseLogger.debug("The seed that was used = {}".format(seed))
        baseLogger.debug("Theta used for simulation = {}".format(theta))
        baseLogger.debug("The switching sequence used: {}".format(n))
        #logger2.debug("The seed that was used = {}".format(seed))
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
            if not __name__ == '__main__' and gui:
                self.progressbar.setValue(int((100/blocks)*(i+1)))
            for model in Run[i][0].modelfinal:
                Par.oldmodels.add(tuple(model))
            baseLogger.debug("Nl for this run = {}".format(Run[i][0].Nl))
            if gui:
                block_operations.updateLog(self,Console)

        ## Retrieve final dataset
        Yt=np.hstack(([Y[i] for i in range(blocks)]))
        Ut=np.hstack(([U[i] for i in range(blocks)]))
        Rt=np.vstack(([R[i] for i in range(blocks)]))
        Total_ds=dataset.Dataset(Yt,Rt,Ut)

        ## Retrieve models, switching sequence and switches
        (Fmodels,sigmat,swt) = block_operations.merge_blocks(Run,blocks)
        storage.datastorage(Fmodels, Par, 'HVAC', 'Full HVAC dataset', 0, 0, sigmat, 0)
        (error,mse,sigmat)=block_operations.hitting_set_model_reduction(self,Console,input,output,r,nu,ny,inputs,outputs,swt,Fmodels,delta,theta,n,modelgeneration,gui)
        t2=time.time()
        baseLogger.info("Total time for identification process: {}".format(t2-t1))
        if gui:
            block_operations.updateLog(self,Console)
        storage.datastorage(Fmodels, Par,'HVAC','Full HVAC dataset, with error correct',error,mse,sigmat,float(t2-t1))


        ## plot switching sequence in case of using the gui
        if not __name__ == '__main__' and gui:
            cm = plt.cm.get_cmap('tab20')
            self.ax_switching.clear()
            scatterplot=self.ax_switching.scatter(np.array(list(range(1,len(sigmat)+1))),sigmat,c=sigmat,s=200,marker='x')
            legend1=self.ax_switching.legend(*scatterplot.legend_elements(),
                    loc="upper right", title="Modes")
            self.ax_switching.add_artist(legend1)
            self.ax_switching.grid()
            self.ax_switching.set_xlabel('[k]')
            self.ax_switching.set_ylabel('mode')
            self.plot_switching.draw()

        ## in case of a PWARX model, partition the regressorspace
        if pwarx==1:
            Ht=partition.partition(Rt,sigmat,max(sigmat)+1)
            baseLogger.info("Coefficients of separating hyperplanes: {}".format(Ht))
        else:
            Ht=[]
        export.exportToMatlab(Fmodels,sigmat,Ut,Yt,Rt,Ht,pwarx)
        return(n,sigmat,error,mse,t2-t1,Fmodels)
if __name__ == '__main__':
    main([])