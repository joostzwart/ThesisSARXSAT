
import cplex
import numpy as np




def minimizeL1(regressor,ny,nu,delta,output,weights):
    """Calls Cplex

    Arguments:
        regressor {list} -- Set of regressor data
        ny {integer} -- Order of the output
        nu {integer} -- Order of the input
        delta {float} -- Bound of error
        output {list} -- Set of output data
        weights {list} -- Weights for the L1 minimization

    Returns:
        values {list} -- Optimized variables
        status {list} -- Status
    """

    ## initializing the problem and turning of printing to screen
    prob = cplex.Cplex()
    prob.set_log_stream(None)
    prob.set_error_stream(None)
    prob.set_warning_stream(None)
    prob.set_results_stream(None)

    ## defining parameters
    NOP=ny+nu                                                                                       #number of parameters
    NOPV=len(output)                                                                                #NOPV = number of parameter vectors

    ## defining cplex specific parameters
    my_obj = weights+[0.0]*(NOP)*(NOPV+1)                                                           #Objective function
    my_ub = [3.0]*(NOP)*(1+2*NOPV)                                                                  #Upperbound on variabels
    my_lb = [-3.0]*(NOP)*(1+2*NOPV)                                                                 #Lowerbound on variabels
    my_colnames =   ["z{}({})".format(i,j) for i in range(1,NOPV+1) for j in range(1,NOP+1)]+\
                ["pb{0}".format(i) for i in range(1,NOP+1)]+\
                ["p{}({})".format(i,j) for i in range(1,NOPV+1) for j in range(1,NOP+1)]             #Name of variables
    rows=[]                                                                                         #Constraints
    ## Creating right side of constraints
    rhs2=(delta-output).tolist()+(delta+output).tolist()
    my_rhs = [0.0]*(ny+nu)*2*NOPV+rhs2
    size=len(my_rhs)
    my_rownames = ["c{0}".format(i) for i in range(1,size+1)]
    my_sense = ["L"] * size


    ## Creating constraint matrix
    c1=np.hstack((np.zeros((NOP,NOP*NOPV)),-np.identity(NOP),np.zeros((NOP,NOP*NOPV))))
    c2=np.hstack((np.zeros((NOP,NOP*NOPV)),np.identity(NOP),np.zeros((NOP,NOP*NOPV))))
    c=np.vstack((c1,c2))
    cons=np.tile(c,(NOPV,1))

    for k in range (0,NOPV):
        n=k*2
        cons[n*NOP:n*NOP+NOP,(k+NOPV+1)*NOP:(k+NOPV+2)*NOP]=np.identity(NOP)
        cons[(n+1)*NOP:(n+1)*NOP+NOP,(k+NOPV+1)*NOP:(k+NOPV+2)*NOP]=-np.identity(NOP)

        cons[n*NOP:n*NOP+NOP,k*NOP:(k+1)*NOP]=-np.identity(NOP)
        cons[(n+1)*NOP:(n+1)*NOP+NOP,k*NOP:(k+1)*NOP]=-np.identity(NOP)

        q= np.zeros((2*NOPV,(1+NOPV+NOPV)*NOP))

    for j in range(0,NOPV):
        q[j,(1+NOPV+j)*NOP:(2+NOPV+j)*NOP]=-regressor[j]
        q[j+NOPV,(1+NOPV+j)*NOP:(2+NOPV+j)*NOP]=regressor[j]
        cons=np.vstack((cons,q)).tolist()



    ## Creating the total contraint matrix (rows)
    for z in range(size):
        rows.append([my_colnames,cons[z]])

        prob.objective.set_sense(prob.objective.sense.minimize)
        prob.variables.add(obj=my_obj, lb=my_lb, ub=my_ub, names=my_colnames)
        prob.linear_constraints.add(lin_expr=rows, senses=my_sense,
                                rhs=my_rhs, names=my_rownames)

    ##Solve the problem and return the solution
    prob.solve()
    if prob.solution.get_status()==1:
        return(prob.solution.get_values(),prob.solution.get_status())
    else:
        return(0,prob.solution.get_status())

def weightedL1(regressor,ny,nu,delta,output,NOI,regc,W):
    NOPV=len(output)
    for i in range(0,NOI):
        (vari,status)=minimizeL1(regressor,ny,nu,delta,output,W)
        if status != 1:
            return(0,0)
        z=vari[:(ny+nu)*NOPV]
        pb=vari[NOPV*(ny+nu):(NOPV+1)*(ny+nu)]
        W=[1/(x+regc) for x in z]
        return(pb,status,vari)

