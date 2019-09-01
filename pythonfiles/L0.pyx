
""" L0.pyx takes care of running one iteration of the identification scheme.
    It schedules the SAT solver and the Theory solver.

    Returns:
        ds {list} -- Dataset
"""
import numpy as np
cimport numpy as np
import SAT
import TheorySolver as TS
import logging
import models
import Sparsification as sp
import SimplifyCertificate as SC

def main(UserParameters,UserPreferences,**kwargs):
    import model as model
    import dataset

    ## initializing variables
    cdef int T=UserParameters.T                                                                 #Length of Dataset
    cdef int ny=UserParameters.ny                                                               #Order of the output
    cdef int nu=UserParameters.nu                                                               #Order of the input
    cdef float delta=UserParameters.delta                                                       #Bound on the error
    cdef int dw=UserParameters.dw                                                               #Dwell time
    cdef float nt=UserParameters.nt                                                             #Bound on the noise
    cdef int NOD= UserParameters.nod                                                            #Number of datasets used.
    cdef int inputs=UserParameters.inputs                                                       #Number of inputs
    cdef int outputs=UserParameters.outputs                                                     #Number of outputs
    cdef int pwarx=UserParameters.pwarx
    cdef int Ls=UserParameters.Ls
    ## User choises
    cdef int seed=UserPreferences.seed                              #Seed of the randomness. 0: for a random seed
    cdef int modelgeneration=UserPreferences.modelgeneration        #Method of generating models. 1: use python to generate a model.
    cdef int info=10                                                #Extend of Logging. DEBUG: No info. INFO: Process information. WARNING: Errors.
    cdef int Merging=UserPreferences.merging                        #Option to use more datasets. 0: for using one dataset. 1: for using more datasets
    cdef int SplitLargeDataset=UserPreferences.splitlargedataset    #Option used for large datasets.
    cdef int chunks=UserParameters.chunks                           #Number of chunks used in the dataset
    cdef set oldmodels=UserParameters.oldmodels
    cdef int unstuck = UserPreferences.unstuck                      #Method of dealing with the L0 minimization getting stuck

    ## Setting log level
    logger=logging.getLogger('dev')
    logger.setLevel(info)
    np.set_printoptions(formatter={'float': '{: 0.10f}'.format})

    ## typecasting
    cdef int m, i, status, N, j, k
    cdef double[:,:] y, r, u,
    cdef str satisfied
    cdef list n, ds,certificate, cert, TotalSwitches, Nl, SigmaT, model_ident, modelfinal, pwarx_certificate
    cdef list dsnew=[]
    cdef np.ndarray theta

    ## Simulating model
    if modelgeneration==1:
        (y,r,n,u,theta)= model.simulate(T, dw, 'OzayMIMO', nt, seed, inputs, outputs)
        ds=[dataset.Dataset(y, r, u)]

        ## Simulate more models when using multiple datasets
        if Merging==1:
            for m in range(NOD):
                (y, r, n, u, theta)= model.simulate(T, dw, 'OzayMIMO', nt, seed, inputs, outputs)
                ds.append(dataset.Dataset(y, r, u))

    #Loading in dataset
    if modelgeneration==3 or modelgeneration==2:
        if 'ProvidedDataset' in kwargs:
            ds=[kwargs.get("ProvidedDataset")]
        else:
            logger.warning("no dataset provided. exiting program....")
            exit()

    # adding model parameters for calculating error
    if modelgeneration==3:
        theta=kwargs.get("theta")
        n=kwargs.get("n")

    ## Split large Datasets
    if SplitLargeDataset==1:

        ## Split up the dataset
        ysplitted=np.array_split(ds[0].y,chunks,axis=1)
        rsplitted=np.array_split(ds[0].r,chunks)
        usplitted=np.array_split(ds[0].u,chunks,axis=1)

        ## Generate datasets from one large dataset
        for i in range(chunks):
            dsnew.append(
                dataset.Dataset(ysplitted[i], rsplitted[i], usplitted[i]))
        ds=dsnew

        ## Proceed as if more datasets were presented to the program
        Merging=1
        NOD=chunks

    ## Keep track of feasible intervals
    cdef list dictionary=[]
    T=T-max(ny,nu)

    # Run solver till a feasible model is found
    for N in range(NOD):
        certificate=[]

        ## Generating initial switching seqence
        (ds[N].switches,satisfied)= SAT.satisfying( ds[N].T, dw, ny, ds[N].certificate, pwarx, ds[N].pwarx_certificate,Ls)

        ## Run the iterative identification process till a feasible model is found or this is impossible
        status=0
        while status!=1:
            (cert,ds[N].model,ds[N].Nl,status,ds[N].dictionary,pwarx_cert)=TS.Theory(ds[N],
            delta,T,ds[N].switches,ds[N].dictionary,ny,nu,inputs,outputs,oldmodels,unstuck,pwarx)
            if status==1:
                break

            ## Add and simplify certificate
            if cert:
                for i in range(len(cert)):
                    ds[N].certificate.append(cert[i])
                ds[N].certificate=SC.simplify(ds[N].certificate)

            ## Add pwarx certificate
            if pwarx_cert:
                for i in range(len(pwarx_cert)):
                    ds[N].pwarx_certificate.append(pwarx_cert[i])

            ## Add total switching sequence when no certificate is generated
            else:
                certificate.append(ds[N].switches)
            (ds[N].switches,satisfied)= SAT.satisfying( ds[N].T, dw, ny, ds[N].certificate, pwarx, ds[N].pwarx_certificate,Ls)

            ## Sanity check whether a switching sequence was found
            if not ds[N].switches:
                logger.warning("No Satisfying switching sequence found. Exiting program...")
                exit()

        ## Removing models not corresponding to an interval
        (ds[N].model,ds[N].Nl)=models.RemovingUnusedModels(ds[N].model,ds[N].Nl)

        ## Merging models that appear to be the same
        (ds[N].modelfinal,ds[N].Nl)= models.mergingModels(ds[N].model, ds[N].Nl)
        ds[N].model=ds[N].modelfinal

        ## Ordering models based on largest first parameter
        (ds[N].Nl,ds[N].modelfinal)= models.orderModels(ds[N].Nl, ds[N].model, ds[N].modelfinal)

        ## Recover switching sequence
        ds[N].switches=[0]+ds[N].switches+[ds[N].T]
        ds[N].Sigma=[0]*(ds[N].T)
        logger.debug("NL {}".format(ds[N].Nl))
        for j in range(len(ds[N].Nl)):
            for k in ds[N].Nl[j]:
                ds[N].Sigma[ds[N].switches[k]:ds[N].switches[k+1]]=[j+1]*(ds[N].switches[k+1]-ds[N].switches[k])

    ## L1 minimization on whole dataset to find the final model
    if Merging==1:
        regressor=[]
        output=[]

        ## Group the regressor and output based on to which interval they belong
        for dataset in ds:
            for interval in range(len(dataset.switches)-1):
                regressor.append(np.take(dataset.r,range(dataset.switches[interval],dataset.switches[interval+1]),0))
                output.append(np.take(dataset.y,range(dataset.switches[interval],dataset.switches[interval+1]),1))
        (model_ident,Nl,status)=sp.L1(regressor, ny, nu, delta, output, 10, 0.01, [], inputs, outputs,oldmodels,unstuck)

        ## Group switching sequences
        TotalSwitches=[]
        for dataset in ds:
            TotalSwitches.append(dataset.switches)
        ds[N].switches=TotalSwitches
        temp=ds[N].switches[0]

        temp=[]
        offset=0
        for i in range(0,len(ds[0].switches)):
            dx=[x+offset for x in ds[N].switches[i]]
            if i!=0:
                del dx[0]
            temp=temp+dx
            offset=offset+ds[N].switches[i][-1]

        ds[N].switches=temp
        ## Removing models not corresponding to an interval
        (model_ident,Nl)=models.RemovingUnusedModels(model_ident,Nl)

        ## Merging models that appear to be the same
        (modelfinal,Nl)= models.mergingModels(model_ident, Nl)

        ## Ordering models based on largest first parameter
        (ds[N].Nl,ds[N].modelfinal)= models.orderModels(Nl, model_ident, ds[N].modelfinal)

        ## Recover Switching sequence
        ds[N].Sigma=models.RecoverMergedSwitching(TotalSwitches,ds[N].Nl)

    ## Printing found parameters
    #if inputs>1 or outputs>1:
    #    (A,B)= models.recoverParameters(ds[N].modelfinal, nu, ny, inputs, outputs)
    #    logging.info("A = {}".format(A))
    #    logging.info("B = {}".format(B))

    ## Print information about found model
    logger.info("satisfied")
    logger.debug("certificate {}".format(ds[N].certificate))
    if modelgeneration==3:
            logger.info("normalized error: {}".format(models.CalculateError(ds[N].modelfinal,theta,n,ds[N].Sigma,ny,nu)))

    ## Print information in case of dataset split in chuncks
    if Merging==1:
        logger.info("Switching sequence {}".format(ds[N].Sigma))
        logger.info("Final Model after L1 minimization: {}".format(ds[N].modelfinal))
        logger.info("Final interval separetion {}".format(ds[N].Nl))
        logger.info("Switches present at {}".format(TotalSwitches))

    else:
        logger.info("Relation between models and intervals = {}".format(ds[N].Nl))
        logger.info("Identified switches {}".format(ds[N].switches))
        logger.info("Identified models {}".format(ds[N].modelfinal))
        logger.info("Identified switching sequence {}".format(ds[N].Sigma))

    ## Writing information to a comma seperated text file for further analysis
    with open("DataSolved.txt", "wb") as f:
        np.savetxt(f, [ds[N].Sigma], fmt='%i', delimiter=",")
        np.savetxt(f, ds[N].modelfinal, fmt='%f6', delimiter=",")
    with open("InputOutput.txt", "wb") as f2:
        for i in range(outputs):
            np.savetxt(f2, ds[N].y[i], fmt='%f6', delimiter=",")
        for i in range(inputs):
            np.savetxt(f2, ds[N].u[i], fmt='%f6', delimiter=",")
    return(ds)

if __name__ == "__main__":
    main()