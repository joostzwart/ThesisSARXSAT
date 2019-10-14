from __future__ import print_function

import cplex
import numpy as np
cimport numpy as np
import logging
import itertools
import random
DTYPE = np.int
FTYPE=  np.double
ctypedef np.int_t DTYPE_t
ctypedef np.double_t FTYPE_t

cpdef minimizeL1(regressor, int ny, int nu, double delta, output,
                weights, int inputs, int outputs):
    """Calls Cplex
    
    Arguments:
        regressor {list} -- Set of regressor data
        ny {integer} -- Order of the output
        nu {integer} -- Order of the input
        delta {float} -- Bound of error
        output {list} -- Set of output data
        weights {list} -- Weights for the L1 minimization
        inputs {integer} -- Number of inputs
        outputs {integer} -- Number of outputs
    
    Returns:
        values {list} -- Optimized variables
        status {list} -- Status
    """
    ##typecasting
    cdef int i, j, k, n, z
    cdef np.ndarray[dtype=FTYPE_t,ndim=2] c1, c2, c, cons, q, qmimo, qt, nt
    ## initializing the problem and turning of printing to screen
    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)

    ## defining parameters
    cdef int NOP=(len(regressor[0][0]))*outputs                                                              #number of parameters
    cdef int NOPV=len(output)                                                                                #NOPV = number of parameter vectors

    ## defining cplex specific parameters
    cdef list my_obj = weights+[0.0]*(NOP)*(NOPV+1)                                                           #Objective function
    cdef list my_ub = [3.0]*(NOP*(1+NOPV)+NOPV)                                                              #Upperbound on variabels
    cdef list my_lb = [-3.0]*(NOP*(1+NOPV)+NOPV)                                                             #Lowerbound on variabels
    cdef list my_colnames =   ["z({})".format(i) for i in range(1,NOPV+1)]+\
                    ["pb{0}".format(i) for i in range(1,NOP+1)]+\
                    ["p{}({})".format(i,j) for i in range(1,NOPV+1) for j in range(1,NOP+1)]        #Name of variables
    cdef list rows=[]                                                                                         #Constraints
                                                                        
    ## Creating right side of constraints
    cdef list rhs2=[]
    for k in range(outputs):
        for j in range(len(output)):
            rhs2=rhs2+(delta-output[j][k]).tolist()+(delta+output[j][k]).tolist()
    
    cdef list my_rhs = [0.0]*NOP*2*NOPV+rhs2
    cdef int size=len(my_rhs)
    cdef int rlen=len(regressor[0][0])
    cdef list my_rownames = ["c{0}".format(i) for i in range(1,size+1)]
    cdef list my_sense = ["L"] * size

    ## Creating constraint matrix
    c1=np.hstack((np.zeros((NOP,NOPV)),-np.identity(NOP),np.zeros((NOP,NOP*NOPV))))
    c2=np.hstack((np.zeros((NOP,NOPV)),np.identity(NOP),np.zeros((NOP,NOP*NOPV))))
    c=np.vstack((c1,c2))
    cons=np.tile(c,(NOPV,1))
    

    for k in range (0,NOPV):
        n=k*2
        cons[n*NOP:n*NOP+NOP,(k+1)*NOP+NOPV:(k+2)*NOP+NOPV]=np.identity(NOP)
        cons[(n+1)*NOP:(n+1)*NOP+NOP,(k+1)*NOP+NOPV:(k+2)*NOP+NOPV]=-np.identity(NOP)

        cons[n*NOP:n*NOP+NOP,k]=-1
        cons[(n+1)*NOP:(n+1)*NOP+NOP,k]=-1

    
    ## Creating constraints for the noise description
    qmimo=np.empty((0,(1+NOPV)*NOP+NOPV), float)
    for i in range(outputs):
        q= np.empty((0,(1+NOPV)*NOP+NOPV), float)
        for j in range(len(regressor)):
            nt=np.zeros((len(regressor[j]),NOP))
            nt[:,i*rlen:(i+1)*rlen]=regressor[j]
            qt=np.zeros((2*len(regressor[j]),(1+NOPV)*NOP+NOPV))
            qt[:len(regressor[j]),(1+j)*NOP+NOPV:(2+j)*NOP+NOPV]=-nt
            qt[len(regressor[j]):,(1+j)*NOP+NOPV:(2+j)*NOP+NOPV]=nt
            q=np.vstack((q,qt))
        qmimo=np.vstack((qmimo,q))

    cdef list cons_list=np.vstack((cons,qmimo)).tolist()

    ## Creating the total contraint matrix (rows)
    for z in range(len(cons_list)):
        rows.append([my_colnames,cons_list[z]])

    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, names=my_colnames)
    prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                                rhs=my_rhs, names=my_rownames)

    ## Choose Solver (optional)                          
    #prob.parameters.lpmethod.set(prob.parameters.lpmethod.values.dual)

    ##Solve the problem and return the solution
    prob.solve()
    if prob.solution.get_status()==1:
        return(prob.solution.get_values(),prob.solution.get_status())
    else:
        return([],prob.solution.get_status())




cpdef L1(list regressor, int ny, int nu, double delta, list output, int NOI,
       double regc, list models, int inputs, int outputs, set oldmodels, int unstuck):
    """ 
    This function is the outer shell of the L1 minimization. 
    It groups datapoints to found models and identifies models
    till all datapoints are accounted for.
    
    Arguments:
        regressor {array} -- Regressor
        ny {integer} -- Order of the output
        nu {integer} -- Order of the input
        delta {double} -- Bound on the error
        output {array} -- Output
        NOPV {int} -- Number of intervals (Number of parameter vectors)
        NOI {integer} -- Number of iterations
        regc {double} -- Regulization constant
        models {list} -- List of old models
        inputs {integer} -- Number of inputs
        outputs {integer} -- Number of outputs
        oldmodels {set} -- Set of old models used for other blocks of the dataset
    
    Returns:
        model{list} -- List of models
        Nl{list} -- Relates intervals to models
        status{int} -- Status of optimization program    
    """
    
    ## Initialize variables
    cdef int Difweights=0
    cdef double treshhold=delta/100
    cdef list model=[]
    cdef list Errors, NError
    cdef list Nl=[range(len(regressor))]
    cdef list weights, pb, NT
    cdef tuple oldmodel
    cdef int status, i, x, sorted, q, em
    cdef double err, acceptederror
    cdef int errorcounter=0
    r=regressor
    y=output

    ## Add oldmodels from different blocks
    for oldmodel in oldmodels:
        model.append(list(oldmodel))
    stuck=0
    ## Run the solver till all datapoints are accounted for
    while Nl[-1]:

        ## Identify a model fitting as much intervals as possible
        if Difweights==0:
            weights=[1]*len(y)
        else:
            weights=[random.random() for i in range(len(y))]

        ## Find the next model
        if stuck==0:
            (pb,status)=weightedL1(r,ny,nu,delta,y,NOI,regc,weights,inputs,outputs)
            model.append(pb)
        else:
            model=modelForSmallestErrors(Errors,regressor,ny,nu,delta,output,NOI,regc,weights,inputs,outputs,model,NError,1)
            stuck=0

        ## return if no feasible model can be found
        if status==0:
            logging.debug("No correct model can be identified")
            return([],[],0)

        ## Reset Error matrices
        Errors=[]
        NError=[]

        ## Classify intervals to models
        NT=[[] for i in range(0,len(model)+1)]
        for x in range(len(output)):
            sorted=-1
            err=0
            for i in range(len(model)):
                err=np.linalg.norm(output[x].flatten()-np.kron(np.eye(outputs),regressor[x]).dot(model[i]),np.inf)
                if err<(delta+treshhold) and sorted==-1:
                    NT[i].append(x)
                    sorted=i
                    acceptederror=err

                ## Calculate errors for last model and keep track of the corresponding intervals
                elif i==len(model)-1 and sorted==-1:
                    Errors.append(err)
                    NError.append(x)

                ## When a model fits better, attribute the interval to this new model
                if sorted>-1 and err<acceptederror:
                    for j in range(len(NT)):
                        if x in NT[j]:
                            NT[j].remove(x)
                    NT[i].append(x)

            ## Add interval to the "to be sorted" list if no corresponding fitting model has been found
            if sorted==-1:
                NT[i+1].append(x)

        ## If no convergence is obtained, initialize new random weights
        if not NT[-2] and unstuck==1:
            sorted=-1
            Difweights=1
        else:
            Difweights=0

        ##Second method when no convergence is obtained. This time search for model for the intervals with lowest erros.
        if not NT[-2] and unstuck==2:
            stuck=1

        ## Create regressor of datapoints that are not yet assigned to a model
        r=np.array([regressor[k] for k in NT[-1]])
        y=[output[j] for j in NT[-1]]

        ## Break if all intervals are accounted for
        if not NT[-1]:
            del NT[-1]
            Nl=NT
            break
        Nl=NT

    if len(model)>len(Nl):
        for em in range(len(model)-len(Nl)):
            Nl.append([])

    return(model,Nl,status)

cpdef weightedL1(regressor, int ny, int nu, double delta, output,
                 int NOI, double regc, list W, int inputs, int outputs):
    """
    Weighted L1 optimization.
    
     Arguments:
        regressor {array} -- Regressor
        ny {integer} -- Order of the output
        nu {integer} -- Order of the input
        delta {double} -- Bound on the error
        output {array} -- Output
        NOI {integer} -- Number of iterations
        regc {double} -- regulization constant
        W {list} -- Weights
        inputs {integer} -- Number of inputs
        outputs {integer} -- Number of outputs
    
    Returns:
        model{list} -- model
        status{int} -- Status of optimization program
    """
    
    ##Perform the weighted L1 relaxation
    cdef int NOPV=len(output)
    cdef list vari, z, pb
    cdef int status

    for _ in itertools.repeat(None, NOI): 
        (vari,status)=minimizeL1(regressor,ny,nu,delta,output,W,inputs,outputs)
        if status != 1:
            return([],0)
        
        ## Substract the important parameters and update the weights
        z=vari[0:NOPV]
        model=vari[NOPV:NOPV+(len(regressor[0][0]))*outputs]
        W=[1/(x+regc) for x in z]
    return(model,status)

cpdef modelForSmallestErrors(list Errors, r, int ny, int nu, double delta, y, int NOI, double regc, list weights,
                             int inputs, int outputs, model, list NError, k):
    """
    Calculate new model for intervals with the lowest error
    
     Arguments:
        Errors {list} -- List of Errors
        r {array} -- Regressor
        ny {integer} -- Order of the output
        nu {integer} -- Order of the input
        delta {double} -- Bound on the error
        y {array} -- Output
        NOI {integer} -- Number of iterations
        regc {double} -- regulization constant
        weights {list} -- Weights
        inputs {integer} -- Number of inputs
        outputs {integer} -- Number of outputs
        model {list} -- List of models
        Nerrors {list} -- List that connects intervals to the respective errors
    
    Returns:
        model{list} -- List of models with new model included
    """

    E=np.array(Errors)
    id=np.argpartition(E, k)
    id=id[:k]
    ix=[NError[a] for a in id]
    r=np.array([r[k] for k in ix])
    y=[y[j] for j in ix]
    (pb,status)=weightedL1(r,ny,nu,delta,y,NOI,regc,[1]*len(ix),inputs,outputs)
    model.append(pb)
    return model

