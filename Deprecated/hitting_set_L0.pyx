import cplex
import numpy as np
cimport numpy as np
import itertools

cpdef solve(np.ndarray S,list K,list weights):
    """
        Deprecated file. Can be switched for solving the hitting set problem via L0 minimization. It is quicker than SAT solving 
        but does often not result in the true minimum hitting set. 
    """
    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)

    cdef int i, z
    cdef int L=len(S)
    cdef int l=len(K)

    cdef list my_obj=weights+[0.0]*(l+L*l)
    cdef list my_ub=[l]*(L+l+L*l)
    cdef list my_lb=[0]*(L+l+L*l)
    cdef list my_colnames = ["z({})".format(i) for i in range(1, L + 1)] + \
                  ["xb{}".format(i) for i in range(1, l + 1)] + \
                  ["x{}({})".format(j, i) for i in range(1, L + 1) for j in range(1, l + 1)]
    cdef list rows = []
    cdef list my_rhs=[0]*(2*l*L)+[-1.0]*(L)+[-1.0]*(L)

    cdef int size = len(my_rhs)
    cdef list my_rownames = ["c{0}".format(i) for i in range(1, size + 1)]
    cdef list my_sense = ["L"] * size

    cdef np.ndarray cons=np.zeros((size,L+l+L*l))
    for i in range(L):
        ## zt
        cons[i * l:(i + 1) * l, i] = -1
        cons[L*l+(i * l):L*l+((i + 1) * l), i] = -1

        ## xt
        cons[i * l:(i + 1) * l, (L + l) + i * l:(L + l) + (i + 1) * l] = np.eye(l)
        cons[(L*l)+i * l:(L*l)+((i + 1) * l), (L + l) + i * l:(L + l) + (i + 1) * l] = -np.eye(l)

        ##mt
        cons[2 * l * L + i,L+l+(i*l):L+l+((i+1)*l)] = -S[i, :]
        cons[2 * l * L + L+i,L:L+l] = -S[i, :]
    ##xbar
    cons[0:L*l,L:L+l]=-np.tile(np.eye(l), (L, 1))
    cons[L*l:2*L * l, L:L + l] = np.tile(np.eye(l), (L, 1))

    cdef list cons_list=cons.tolist()

    ## Creating the total contraint matrix (rows)
    for z in range(len(cons_list)):
        rows.append([my_colnames,cons_list[z]])

    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, names=my_colnames)
    prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                                    rhs=my_rhs, names=my_rownames)
    prob.write("lpex1.lp")
    ##Solve the problem and return the solution
    prob.solve()
    if prob.solution.get_status()==1:
        return(prob.solution.get_values(),prob.solution.get_status())
    else:
        return([],prob.solution.get_status())

cpdef weighting(np.ndarray S,list K):

    ##Perform the weighted L1 relaxation
    cdef int status
    cdef list vari, z, pb
    cdef int L=len(S)
    cdef int l=len(K)
    cdef list weights = [1.0] * L
    for _ in itertools.repeat(None, 5):
        (vari, status) = solve(S,K,weights)
        if status != 1:
            return ([], 0)

        ## Substract the important parameters and update the weights
        z = vari[0:L]
        pb = vari[L:L+l]
        weights = [1 / (x + 0.001) for x in z]
    return (pb, status)