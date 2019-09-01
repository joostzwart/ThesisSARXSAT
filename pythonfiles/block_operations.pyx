import numpy as np
cimport numpy as np
import logging
import models
import SAT
import LinProg as LP

cpdef merge_blocks(Run, int blocks):
    """ 
    This function merges the blocks of the original dataset back to one dataset again
    
    Arguments:
        Run {dict} -- Information on all the runs
        blocks {int} -- Number of blocks in which the dataset was split up
    Returns:
        Fmodels {list} -- List of models
        sigmat {list} -- switching sequence
        swt {list} -- switches 
    """
    cdef list Nlt=[]
    cdef list Fmodels=[]
    cdef list modelst=[]
    cdef list swt=[]
    cdef list sw, sigma2t, sigmat
    cdef int i, counter,
    cdef int offset=0

    ## Set logger
    logger=logging.getLogger('dev')

    ## Retrieve list of intervals
    Nlt=Nlt+Run[0][0].Nl

    for i in range(1,blocks):
        offset=offset+max([max(a) for a in Run[i-1][0].Nl])+1
        Nlt=Nlt+[[x+offset for x in A] for A in Run[i][0].Nl]

    ## Retrieve list of models
    for i in range(0,blocks):
        modelst=modelst+Run[i][0].modelfinal

    ## Some information for debug purposes
    logger.info("NL for the whole dataset = {}".format(Nlt))
    logger.debug("Total number of models identified = {}".format(len(Nlt)))
    logger.debug("Total models identified with duplicates still present = {}".format(modelst))
    logger.debug("Number of models identified with duplicates still present = {}".format(len(modelst)))

    ## Retrieve final switchting sequence
    offset=0
    for i in range(0,blocks):
        sw=[x+offset for x in Run[i][0].switches]
        del sw[0]
        swt=swt+sw
        offset=offset+Run[i][0].switches[-1]
    logger.info("final identified switches present at {}".format(swt))

    ## Retrieve which model is active at which interval
    sigma2t=[0]*(max([max(a) for a in Nlt])+1)
    offset=0
    for counter,value in enumerate(Nlt):
        for x in value:
            sigma2t[x]=counter
    logger.debug("sigma2t {}".format(sigma2t))
    sigmat=sigma(swt,sigma2t)

    ## Remove duplicate models
    for sublist in modelst:
        if sublist not in Fmodels:
            Fmodels.append(sublist)
    logger.debug("models without duplicates {}".format(Fmodels))
    logger.debug("number of models left {}".format(len(Fmodels)))
    return(Fmodels,sigmat,swt)

cpdef sigma(list swt, sigma2t):
    """ 
    This function merges the switching sequences
    
    Arguments:
        swt {list} -- switches
        sigma2t {int} -- switching sequence
        
    Returns:
        sigmat {list} -- switching sequence
    """
    cdef list sigmat=[]
    cdef int offset=0
    cdef int counter

    for counter,value in enumerate(swt):
        sigmat=sigmat+[sigma2t[counter]]*(value-offset)
        offset=value
    return(sigmat)

cpdef hitting_set_model_reduction(self,Console,np.ndarray input, np.ndarray output, np.ndarray r, int nu, int ny,
                                  int inputs, int outputs, swt,list Fmodels, double delta, np.ndarray theta, list n,
                                  int modelgeneration,gui):
    """ 
    This function solves the minimum hitting set problem by SAT solving
    
    Arguments:
        input {array} -- Input
        output {array} -- Output
        r {array} -- Regressor
        nu {int} -- Order of the input
        ny {int} -- Order of the output
        inputs {int} -- Number of inputs
        outputs {int} -- Number of outputs
        swt {list} -- switches
        Fmodels {list} -- List of old models
        delta {double} -- Bound on the error
        n {list} -- True switching sequence
        modelgeneration {int} -- Synthetic or real system
    
    Returns:
        datafiterror {double} -- List of models
        mse {double} -- Relates intervals to models
        sigma4 {list} -- Switching sequence    
    """

    ## Set logger
    logger=logging.getLogger('dev')

    ## Try to reduce number of models further by fitting models to intervals
    u2=np.split(input,swt,axis=1)
    y2=np.split(output,swt,axis=1)
    r2=np.split(r,swt,axis=0)
    if swt[-1]==len(r):
        del u2[-1]
        del y2[-1]
        del r2[-1]

    ## Typecasting
    cdef list Fittingmatrix=[]
    cdef int interval, i
    cdef list Ft, St, m3x
    cdef np.ndarray y3x, r3x, min_mod
    cdef double datafiterror = 0
    cdef double mse = 0

    ## Create fittingmatrix
    for interval in range(len(u2)):
        Ft=[]
        i=0
        for m in Fmodels:
            try:
                if np.linalg.norm(y2[interval].flatten() - np.kron(np.eye(outputs), r2[interval]).dot(np.array(list(m))), np.inf)<delta+delta/100:
                    Ft.append(i)
                i = i + 1
            except:
                print("Empty array was added.")
        Fittingmatrix.append(Ft)
    logger.debug("Fittingmatrix {}".format(Fittingmatrix))

    ## Hitting set problem
    cdef str satisfiability="unsat"
    k=1
    while satisfiability!="sat":
        (satisfiability,min_mod)=SAT.hitting_set(Fittingmatrix,list(range(len(Fmodels))),k)
        k=k+1

    logger.info("Reduced number of models = {}".format(min_mod))
    logger.info("cardinality of models = {}".format(len(set(min_mod))))

    ## merge fitting datasets
    cdef list S=[]
    cdef list model_mapping=[]
    for m in set(min_mod):
        St = []
        for interval in range(len(y2)):
            if int(min_mod[interval]) == int(m):
                St.append(interval)
        S.append(St)
        model_mapping.append(m)

    ## Split the input, output and the regressor into parts
    cdef list y3=[]
    cdef list r3=[]
    cdef list m3=[]
    for i in range(len(S)):
        r3x = np.empty((0,len(r[0])))
        y3x = np.empty((outputs,0))
        for j in S[i]:
            y3x = np.hstack((y3x, y2[j]))
            r3x = np.vstack((r3x, r2[j]))
        y3.append(y3x)
        r3.append(r3x)

    ## Identify models one final time
    for i in range(len(y3)):
        (m3x,status)=LP.populatebyrow(y3[i],r3[i],delta,ny,nu,inputs,outputs,True)
        m3.append(m3x)
    cdef list sigma3=sigma(swt,min_mod)
    cdef list sigma4=[0]*len(sigma3)
    for i in range(len(sigma3)):
        sigma4[i]=model_mapping.index(sigma3[i])

    ## Print some information about final models
    logger.info("Switching sequence = {}".format(sigma4))
    logger.info("Final Models: {}".format(m3))
    if modelgeneration==3:
        mse=models.CalculateError(m3,theta,n,[x+1 for x in sigma4],ny,nu)
        logger.info("MSE = {}".format(mse))
    datafiterror=models.CalculateDatafitError(m3,[x+1 for x in sigma4],output,r,outputs)
    logger.info("Data fit error = {}".format(datafiterror))
    if gui:
        updateLog(self,Console)
    return(datafiterror,mse,sigma4)

def updateLog(self,console):
        if not __name__ == '__main__':
            newtext=console.stream.getvalue()
            newtext=newtext.splitlines()
            if len(newtext)>0:
                newtext=newtext[-1]
            else:
                newtext="---"
            self.logLabel3.setText(self.logLabel2.text())
            self.logLabel2.setText(self.logLabel1.text())
            self.logLabel1.setText(newtext)
