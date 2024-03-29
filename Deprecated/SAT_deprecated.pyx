#!/usr/bin/python
from z3 import * 
import numpy as np
cimport numpy as np
'''
cpdef constraint_cardinality(int N,int K):
    cdef list b = []
    cdef list c = []
    cdef list ct = []

    for i in range(N):
        b.append(Bool('Node {}'.format(i)))
    for i in range(1,N+1):
        for j in range(i+1):
            ct.append(Bool('Layer {}, Node {}'.format(i,j)))
        c.append(ct)
        ct=[]

    s = Solver()
    s.add(Or(Not(b[0]),c[0][1]))
    s.add(Or(b[0],c[0][0]))
    for l in range(N-1):
        s.add([Or(Not(And(c[l][i],b[l+1])),c[l+1][i+1]) for i in range(l+2)])
        s.add([Or(Not(And(c[l][i],Not(b[l+1]))),c[l+1][i]) for i in range(l+2)])

    x=list(range(N))
    x.remove(K)
    s.add([Not(c[-1][j]) for j in x])

    satisfied=str(s.check())
    if satisfied=="sat":
        sol=s.model()
        print("sol",sol)
        Q=[sol.evaluate(b[i],model_completion=True) for i in range(N)]
        Q2=[str(x)=="True" for x in Q]
        print("Q",Q2)
    else:
        print("NOT satisfied")
'''
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

##source: https://pdfs.semanticscholar.org/3500/d44ab62291421c5be82949dfd4279e5a95eb.pdf


cpdef satisfying(int Td,int dwell,int ny,list certificate):
    """This function handles the Boolean logic and proposes a switching sequence
    
    Arguments:
        Td {integer} -- [Length of dataset]
        dwell {integer} -- [Specified dwell time]
        ny {integer} -- [Order of the output]
        certificate {list} -- [List of certificates]
    """

    ##typecasting
    cdef int i, j, k, val

    ## Creating Booleans
    cdef list b=[]
    for i in range(Td):
        b.append(Bool('b{}'.format(i)))
    cdef int T=len(b)
    s = Solver()

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
    
    ##Translating SAT problem to index   
    cdef str satisfied=str(s.check())
    cdef list index=[]
    cdef str bool
    if satisfied=="sat":
        sol=s.model()
        for id in sol.decls():
            if is_true(sol[id]):
                bool=str(id)
                index.append(int(bool[1:]))
    else:
        print("UNSATISFIED")
    return(index,satisfied)