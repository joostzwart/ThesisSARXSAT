#!/usr/bin/python
from z3 import * 
import numpy as np
cimport numpy as np

cpdef hitting_set(list S, list M,int v):
    """This function solves the hitting set problem
    
    Arguments:
        S {list} -- Sets that need at least one element in common with M  
        M {list} -- Set that needs to have minimal cardinality
        v {integer} -- Cardinality of M
    """
    ## type casting
    cdef np.ndarray T=np.zeros((len(S),len(M)),dtype=bool)
    cdef int i, j, counter, N, l
    cdef list b=[]
    cdef list aux = []
    cdef list countvars = []
    cdef list F = []
    cdef list Ft
    cdef int K=0
    cdef list c = []
    cdef list ct = []
    cdef list Q
    cdef np.ndarray H

    ## Creating Boolean Truth matrix and required Booleans
    for i in range(len(S)):
        for j in range(len(M)):
            if M[j] in S[i]:
                T[i,j]=True
            b.append(Bool('x{} - {}'.format(i,j)))
    for i in range(len(M)):
        aux.append(Bool('aux{}'.format(i)))
    for i in range(K):
        countvars.append(Bool('c{}'.format(i)))

    ## Create solver
    s = Solver()

    ## Constraint at least one item from M in each set s
    for i in range(len(S)):
        F.append(Or([b[i*len(M)+counter] for counter,y in enumerate(T[i,:]) if y]))
    Ft=[And([f for f in F])]
    s.add(Ft)


    ## Constraint auxillairy variables
    F=[]
    for i in range(len(M)):
        F.append(And([Or(Not(b[i+counter*len(M)]),aux[i]) for counter,e in enumerate(T[:,i])]))
    s.add(F)

    ## Create list of extra variables for the cardinality constraint
    N = len(aux)

    for i in range(1,N+1):
        for j in range(i+1):
            ct.append(Bool('Layer {}, Node {}'.format(i,j)))
        c.append(ct)
        ct=[]

    ##add constraint on first node
    s.add(Or(Not(aux[0]),c[0][1]))
    s.add(Or(aux[0],c[0][0]))

    for l in range(N-1):
        s.add([Or(Not(And(c[l][i],aux[l+1])),c[l+1][i+1]) for i in range(l+2)])
        s.add([Or(Not(And(c[l][i],Not(aux[l+1]))),c[l+1][i]) for i in range(l+2)])

    ## Constrain cardinality by setting all upper layer nodes to False except the Vth node
    cdef list P=list(range(N+1))
    P.remove(v)
    s.add([Not(c[-1][j]) for j in P])

    ##Recover parameters if a satisfiable solution was found
    satisfied=str(s.check())
    if satisfied=="sat":
        sol=s.model()
        Q=[sol.evaluate(b[i],model_completion=True) for i in range(len(b))]
        Q=[str(z)=="True" for z in Q]
        H=np.array(Q).reshape((-1,len(M)))
        H2=[]
        for i in range(H.shape[0]):
            for j in range(len(M)):
                if H[i,j]:
                    H2.append(j)
                    break
        return(satisfied,np.array(H2))
    else:
        return(satisfied,np.zeros((2,1)))

cpdef satisfying(int Td, int dwell, int ny, list certificate, int pwarx, list pwarx_certificate,int Ls):
    """This function handles the Boolean logic and proposes a switching sequence
    
    Arguments:
        Td {integer} -- [Length of dataset]
        dwell {integer} -- [Specified dwell time]
        ny {integer} -- [Order of the output]
        certificate {list} -- [List of certificates]
        pwarx {int} -- 1 for PWARX, 0 for SARX
    """

    ##typecasting
    cdef int i, j, k, val
    ## Creating Booleans
    cdef list b=[]
    cdef list a=[]
    for j in range(Td):
        a.append([])
    for i in range(Td):
        b.append(Bool('b{}'.format(i)))
    cdef int T=len(b)
    if pwarx==1:
        for i in range(Td):
            b.append(Bool('b{}'.format(i)))
            for j in range(Ls):
                a[i].append(Bool('a{}{}'.format(i,j)))

    s = Solver()

    if pwarx!=1:
        ## Adding Dwell time constraint
        for k in range(0,len(b)):
            if (k+dwell-1<len(b)):
                s.add(Or(Not(b[k]),And([Not(b[k+i]) for i in range(1,dwell)])))
            else:
                s.add(Or(Not(b[k]),And([Not(b[k+i]) for i in range(1,T-k)])))

    ## adding certificate to prune the searchspace        
    for j in range(len(certificate)):
        if len(certificate[j])==1:
            certificate[j]=[0]+certificate[j]
        ix=range(certificate[j][0]+1,certificate[j][1])
        s.add([Or([b[val] for val in ix])])

    ## add certificates in case of PWARX model
    if pwarx==1:
        for j in range(len(pwarx_certificate)):
            ix=list(range(pwarx_certificate[j][0]+1,pwarx_certificate[j][1]))
            ix2=list(range(pwarx_certificate[j][1]+1,pwarx_certificate[j][2]))
            s.add(Or([b[val] for val in ix]+[Not(b[pwarx_certificate[j][1]])] + [b[val] for val in ix2]))

        ## Adding switching limit constraint
        #constraint for first variable
        s.add(Or(Not(b[0]),a[0][0]))
        for i in range(1,Ls):
            s.add(Not(a[0][i]))

        # constraint that enforces that at least as much booleans are true in the sequence
        # for the next one as the previous one
        for i in range(1,Td):
            for j in range(1,Ls):
                s.add(Or(Not(a[i-1][j]),a[i][j]))
        # add one to counter
        for i in range(1,Td):
            for j in range(1,Ls):
                s.add(Or(Not(And(b[i],a[i-1][j-1])),a[i][j]))

        #when count is reached set rest to zero
        for i in range(Td-1):
            s.add(Or(Not(a[i][Ls-1]),Not(b[i+1])))

        ## First row constraint
        for i in range(1,Td):
            s.add(Or(Not(b[i]),a[i][0]))
            s.add(Or(Not(a[i-1][0]),a[i][0]))

    ##Translating SAT problem to index   
    cdef str satisfied=str(s.check())
    cdef list index=[]
    cdef str bool
    if satisfied=="sat":
        sol=s.model()
        for id in sol.decls():
            if is_true(sol[id]):
                bool=str(id)
                if bool[0]=='b':
                    index.append(int(bool[1:]))
    else:
        print("UNSATISFIED")
    return(index,satisfied)

cpdef satisfying_limit_switching(int Td, int Ls, int ny, list certificate, list pwarx_certificate):
    """This function handles the Boolean logic and proposes a switching sequence
    
    Arguments:
        Td {integer} -- [Length of dataset]
        dwell {integer} -- [Specified dwell time]
        ny {integer} -- [Order of the output]
        certificate {list} -- [List of certificates]
        pwarx {int} -- 1 for PWARX, 0 for SARX
    """
    ##typecasting
    cdef int i, j, k, val
    ## Creating Booleans
    cdef list b=[]
    cdef list a=[]
    for j in range(Td):
        a.append([])
    for i in range(Td):
        b.append(Bool('b{}'.format(i)))
        for j in range(Ls):
            a[i].append(Bool('a{}{}'.format(i,j)))


    cdef int T=len(b)
    s = Solver()


    ## Adding switching limit constraint
    #constraint for first variable
    s.add(Or(Not(b[0]),a[0][0]))
    for i in range(1,Ls):
        s.add(Not(a[0][i]))

    # constraint that enforces that at least as much booleans are true in the sequence
    # for the next one as the previous one
    for i in range(1,Td):
        for j in range(1,Ls):
            s.add(Or(Not(a[i-1][j]),a[i][j]))
    # add one to counter
    for i in range(1,Td):
        for j in range(1,Ls):
            s.add(Or(Not(And(b[i],a[i-1][j-1])),a[i][j]))

    #when count is reached set rest to zero
    for i in range(Td-1):
        s.add(Or(Not(a[i][Ls-1]),Not(b[i+1])))

    ## First row constraint
    for i in range(1,Td):
        s.add(Or(Not(b[i]),a[i][0]))
        s.add(Or(Not(a[i-1][0]),a[i][0]))

    ## adding certificate to prune the searchspace
    for j in range(len(certificate)):
        if len(certificate[j])==1:
            certificate[j]=[0]+certificate[j]
        ix=range(certificate[j][0]+1,certificate[j][1])
        s.add([Or([b[val] for val in ix])])

    ## add certificates for PWARX
    for j in range(len(pwarx_certificate)):
        ix=list(range(pwarx_certificate[j][0]+1,pwarx_certificate[j][1]))
        ix2=list(range(pwarx_certificate[j][1]+1,pwarx_certificate[j][2]))
        s.add(Or([b[val] for val in ix]+[Not(b[pwarx_certificate[j][1]])] + [b[val] for val in ix2]))

    ##Translating SAT problem to index
    cdef str satisfied=str(s.check())
    cdef list index=[]
    cdef str bool
    if satisfied=="sat":
        sol=s.model()
        for id in sol.decls():
            if is_true(sol[id]):
                bool=str(id)
                if bool[0]=='b':
                    print(bool)
                    index.append(int(bool[1:]))
    else:
        print("UNSATISFIED")
    return(index,satisfied)

