import numpy as np
cimport numpy as np

FTYPE=  np.double
ctypedef np.double_t FTYPE_t

#
#
### This class contains functions that work on models
#
cpdef mergingModels(list model,list Nl):
    """This function merges models that appear to be the same.

    Arguments:
        model {list} -- List of models.
        Nl {list} -- List of intervals belonging to a model.

    Returns:
        modelfinal {list} -- Merged list of models.
        Nl {list} -- Merged list of intervals belonging to a model.
    """
    ##typecasting
    cdef int i, j, same
    cdef list modelfinal=[]

    ## Merge models
    for j in range(len(model)):
        same=0
        for i in range(j+1,len(model)):
            if np.linalg.norm(np.array(model[i])-np.array(model[j]))<0.00001:
                same=1
                Nl[i]=Nl[i]+Nl[j]
                Nl[j]=[]
        if same==0:
            modelfinal.append(model[j])
    Nl = [k for k in Nl if k != []]
    return(modelfinal,Nl)

cpdef RemovingUnusedModels(list model, list Nl):
    """This function removes superfluous models
    Arguments:
        model {list} -- List of models.
        Nl {list} -- List of intervals belonging to a model.

    Returns:
        model {list} -- Merged list of models.
        Nl {list} -- Merged list of intervals belonging to a model.
    """
    ## typecasting
    cdef int i

    for i,m in reversed(list(enumerate(Nl))):
        if not m:
            del model[i]
            del Nl[i]

    return(model,Nl)

cpdef orderModels(list Nl, list model, list modelfinal):
    """This function orders models based on the first parameter.

    Arguments:
        model {list} -- List of models.
        Nl {list} -- List of intervals belonging to a model.
        modelfinal{list} -- List of models.

    Returns:
        modelfinal {list} -- Ordered list of models.
        Nl {list} -- Ordered list of intervals belonging to a model.
    """
    ## Typecasting
    cdef int i,j

    ## Order the models
    for j in range(len(modelfinal)):
        for i in range(len(modelfinal)-1):  
            if model[i][0]<model[i+1][0]:
                model[i], model[i+1] = model[i+1],  model[i]
                Nl[i],Nl[i+1]=Nl[i+1], Nl[i]
    return (Nl,model)

cpdef recoverParameters(list model, int nu, int ny, int inputs, int outputs):
    """This function Recovers the system matrices

    Arguments:
        model {list} -- List of models.
        nu {integer} -- Order of the input.
        ny {integer} -- Order of the output.
        inputs {integer} -- Number of inputs.
        outputs {integer} -- Number of outputs

    Returns:
        A {array} -- System matrix A.
        B {array} -- System matrix B.
    """

    ##typecasting
    cdef list A, At, B, Bt, M, submodel
    cdef int i
    ## Creating A
    A=[]
    for submodel in model:
        At=[]
        for y in range(ny):
            M=[submodel[y+i*(ny*outputs+nu*inputs)+j*ny] for i in range(outputs) for j in range(ny)]
            M_np=np.array(M).reshape((outputs,ny))
            At.append(M_np)
        A.append(At)


    ## Creating B
    B=[]
    for submodel in model:
        Bt=[]
        for u in range(nu):
            M=[submodel[ny*outputs+u+i*(ny*outputs+nu*inputs)+j*nu] for i in range(outputs) for j in range(nu)]
            M_np=np.array(M).reshape((outputs,nu))
            Bt.append(M_np)
        B.append(Bt)        

    return (A,B)

cpdef CalculateError(list model, np.ndarray theta, list n, list switchingsequence, int ny, int nu):
    """This function Recovers the merged switching sequence

    Arguments:
        model {list} -- List of models.
        theta {list} -- True parameters used to simulate model.
        n     {list} -- True switching sequence.
        switchingsequence {list} -- Identified switching sequence.

    Returns:
        error {float} -- normalized error
    """
    n=n[max(ny,nu):]
    cdef int i
    cdef double error=0

    ## Creating arrays from lists for the necessary computations
    theta=np.array(theta,dtype=FTYPE).reshape(len(theta),len(model[0]))
    for i in range(len(switchingsequence)):
        error=error+np.sum(np.abs(theta[n[i]-1,:]-model[switchingsequence[i]-1]))
    return error/(len(n)*len(theta[0,:]))

cpdef CalculateDatafitError(list model, list switchingsequence, np.ndarray output, np.ndarray regressor, int outputs):
    """This function Recovers the merged switching sequence

    Arguments:
        model {list} -- List of models.
        switchingsequence {list} -- Identified switching sequence.

    Returns:
        error {float} -- normalized error of the datafit
    """
    cdef int i
    cdef double error = 0
    for i in range(len(switchingsequence)):
        error=error+np.linalg.norm(output[:,i] - np.kron(np.eye(outputs), regressor[i]).dot(model[switchingsequence[i]-1]), 1)
    return error/(len(switchingsequence)*outputs)



cpdef RecoverMergedSwitching(list switches,list Nl):
    """This function Recovers the merged switching sequence

    Arguments:
        Nl {list} -- List of intervals attributed to a given model.
        switches {list} -- List of switches.

    Returns:
        sequence {list} -- Final switching sequence after merging of multiple datasets.
    """
    ##typecasting
    cdef int i, j

    ##Obtaining which model is active at which interval
    cdef list intervals=[0]*sum([len(x) for x in Nl])
    for i in range(len(Nl)):
        for j in Nl[i]:
            intervals[j]=i+1


    ## Translating this to which model is active for which datapoint
    cdef list sequence=[]
    cdef int k=0
    for chunck in switches:
        for i in range(1,len(chunck)):
            sequence=sequence+[intervals[k]]*(chunck[i]-chunck[i-1])
            k=k+1
    return sequence
