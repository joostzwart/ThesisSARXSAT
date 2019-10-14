from __future__ import print_function
import Sparsification as sp
import numpy as np
import LinProg as LP
import logging
import checkPoints

cpdef Theory(dataset, double delta, int NOPV, list switches, list dictionary, int ny, int nu,
             int inputs, int outputs, set oldmodels, int unstuck, int pwarx):
    """Calling the function that runs cplex
    
    Arguments:
        dataset {dataset} -- A dataset
        delta {float} -- Bound on the noise
        NOPV {integer} -- Number of parameter vectors
        switches {list} -- Switching sequence
        dictionary {list} -- Feasible intervals
        ny {integer} -- Order of the output
        nu {integer} -- Order of the input
        inputs {integer} -- Number of inputs
        outputs {integer} -- Number of outputs
    
    Returns:
        certificate {list} -- Certificate for the SAT solver
        model {list} -- Identified models
        Nl {list} -- Datapoints/intervals per model
        status {int} -- status of Cplex
        dictionary {list} -- Feasible intervals
    """
    ## Initialization
    cdef int a=0
    cdef dict b={}
    cdef dict c={}
    cdef list regressor=[]
    cdef list output=[]
    switches.sort()
    cdef list models = []
    cdef int count=len(switches)+1
    cdef list feasible=[]
    cdef list certificate=[]
    cdef list pwarx_certificate=[]
    cdef int i, j
    cdef int status
    slack=False

    ## Setting log level
    logger=logging.getLogger('dev')

    ## Refine the switching sequence
    switches=[0]+switches + [len(dataset.r)]

    ## Split the regressor according to the proposed switching sequence
    for i, k in enumerate(switches[1:]):
        b[i]=np.take(dataset.r,range(a,k),0)
        c[i]=np.take(dataset.y,range(a,k),1)
        a=k

    ## Print switching sequence to give an indication of progression
    logger.debug("switches {}".format(switches))
    
    ## Check the intervals between switching for feasibility
    for i in range(count):
        if [switches[i],switches[i+1]] not in dictionary:
            (_,status)=LP.populatebyrow(c[i],b[i],delta,ny,nu,inputs,outputs,slack)

            ## Create an certificate whenever an infeasible switching sequence is detected
            if status > 1:
                certificate.append(switches[i:i+2])
                feasible.append(status)

            ## Add a dictionary entry whenever a mode is feasible
            elif status == 1:
                dictionary.append([switches[i],switches[i+1]])

    ## Return to the SAT solver if any interval is deemed infeasible
    if 2 in feasible or 3 in feasible or 0 in feasible:
        return(certificate,[],[],3,dictionary,[])

    feasible=[]
    ## In case of an PWARX system perform a check to
    #  determine whether the datapoints from adjecent intervals are linearly seperable
    if pwarx==1:
        for i in range(count-1):
            try:
                status=checkPoints.intersect(10000*b[i],10000*b[i+1],delta)
            except:
                status=1
            feasible.append(status)

            ## Add certificate when a switch is not correct
            if status > 1:
                pwarx_certificate.append(switches[i:i+3])

    ## Return to the SAT solver if any interval is deemed infeasible
    if 2 in feasible or 3 in feasible:
        return(certificate,[],[],3,dictionary,pwarx_certificate)

    ## Create regressor and output
    for j in b:
        regressor.append(b[j])
        output.append(c[j])

    ## Return found models from L1 relaxed optimization problem
    (model,Nl,status)=sp.L1(regressor,ny,nu,delta,output,10,0.01,models,inputs,outputs,oldmodels,unstuck)
    if status == 1:
        return([],model,Nl,status,dictionary,[])
            

 