import cplex
import numpy as np
cimport numpy as np
DTYPE = np.int
FTYPE=  np.double
ctypedef np.int_t DTYPE_t
ctypedef np.double_t FTYPE_t
from cplex.exceptions import CplexError


cpdef populatebyrow(np.ndarray[FTYPE_t, ndim=2] output, np.ndarray[FTYPE_t, ndim=2] regressor, double delta, int ny,
                    int nu, int inputs, int outputs):
    """Linear programming for checking feasible models
    
    Arguments:
        output {list} -- Set of outputs datapoints belonging to the specific interval that is being checked for feasiblity
        regressor {list} -- Set of regressor datapoints belonging to the specific interval that is being checked for feasiblity
        delta {double} -- Bound on the error
        ny {int} -- Order of the regressor
        nu {int} -- Order of the input
        inputs {int} -- Number of inputs
        outputs {int} -- Number of outputs

    Returns:
        variables {list} -- Variables of the minimization problem
        status {list} -- Status of the minimization problem
    """

    ## Initializing the problem and turning of printing to screen
    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)

    ## Type casting
    cdef int i, j, z
    cdef double x

    ## Defining parameters
    cdef list output_list=output.tolist()
    cdef int size=(len(output_list[0]))*2*outputs
    cdef int s=len(output_list[0])
    cdef int l=len(regressor[0])
    cdef int L=len(regressor[0])*outputs
    cdef list my_obj = [0.0]*L
    cdef list my_ub = [3.0]*L
    cdef list my_lb = [-3.0]*L
    cdef list my_colnames = ["x{0}".format(i) for i in range(1,L+1)]
    cdef list my_rhs = [x for j in range(outputs) for x in output_list[j]]+[-x for j in range(outputs) for x in output_list[j] ] + s*[1.01*delta]
    cdef list my_rownames = ["c{0}".format(i) for i in range(1,size+1+s)]
    cdef list my_sense = ["L"] * (size+s)

    ## Turn on L1 minimization
    my_colnames=my_colnames+["s{0}".format(i) for i in range(1,s+1)]
    my_obj = my_obj+[0.0]*s
    my_ub = my_ub+[1000.0]*s
    my_lb = my_lb+[0.0]*s

    ## Creating left side of constraints
    cdef np.ndarray R=np.zeros((size,L),dtype=FTYPE)
    for i in range(outputs):
        R[i*s:(i+1)*s,i*l:(i+1)*l]=-regressor
        R[(i+outputs)*s:(i+outputs+1)*s,i*l:(i+1)*l]=regressor
    cdef list cons = R.tolist()

    R=np.hstack((R,np.tile(-np.identity(s),(2*outputs,1))))

    ## constraint on s
    R2=np.zeros((s,L),dtype=FTYPE)
    R2=np.hstack((R2,np.identity(s)))
    R3=np.vstack((R,R2))
    cons=R3.tolist()

    ## Adding names to the constraints (rows)
    cdef list rows=[]
    for z in range(size+s):
        rows.append([my_colnames,cons[z]])

    ## creating qmatrix
    qmat=[]
    for i in range(L):
        qmat.append([[], []])
    for i in range(s):
        qmat.append([[L+i],[1.0]])


    ## Add the required information to Cplex and solve the problem
    prob.objective.set_sense(prob.objective.sense.minimize)
    prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, names=my_colnames)
    prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                                rhs=my_rhs, names=my_rownames)
    prob.objective.set_quadratic(qmat)
    prob.solve()

    ## Return the model and the status
    if prob.solution.get_status()==1:
        sol=prob.solution.get_values()
        sol=sol[0:L]
        sol=[-x for x in sol]
        return sol, prob.solution.get_status()
    else:
        return([],3)

