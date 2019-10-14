import numpy as np
cimport numpy as np
import random

cpdef generate(int T,int ny,int nu,int inputs,int outputs,double nt,int modes,int dw,seed,int pwarx,int inputtype):
    """Generate random SARX/PWARX datasets

    Arguments:
        T {int} -- Length of the dataset
        ny {int} -- Order of the regressor
        nu {int} -- Order of the input
        inputs {int} -- Number of inputs
        outputs {int} -- Number of outputs
        nt {double} -- Bound on the noise
        modes {int} -- Number of modes
        dw {int} -- Minimal dwell time
        seed {int} -- Seed for the randomness
        pwarx {int} -- True if a PWARX system is desired
    Returns:
        thetaR {array} -- Parametervector of the system
        n{list} --   Switching sequence
    """

    if seed != 0:
        random.seed(seed)
        np.random.seed(seed)
    else:
        random.seed(random.randint(1,100000))
        np.random.seed(random.randint(1,100000))

    ## generate the parameter vectors
    cdef np.ndarray thetaY=(2.0/ny*outputs)*(np.random.rand(modes,outputs,ny*outputs)-0.5)
    cdef np.ndarray thetaU=4*np.random.rand(modes,outputs,nu*inputs)-2
    cdef np.ndarray thetaR
    cdef np.ndarray thetaAffine = np.random.rand(modes, outputs, 1)
    cdef list u = []
    cdef int m=0
    cdef int i
    cdef np.ndarray newmode, oldmodes

    ## generate the regressor vector
    if pwarx==1:
        thetaR = 0.5*np.dstack((thetaY,thetaU,thetaAffine))
    else:
        thetaR = 0.5*np.dstack((thetaY,thetaU))
    H={}
    ## Generating uniformly distributed random input
    for i in range(inputs):
        if inputtype==2:
            time=np.arange(0, T, 1)
            u_i=(5*random.random()-2.5)*np.sin(0.5*random.random()*time)
        else:
            [u_i] = (5 * np.random.rand(1, T) - 2.5)
        u.append(u_i)

    ## Generating uniformly distributed random y0
    cdef np.ndarray y = np.zeros((outputs, T))
    for i in range(outputs):
        [y_i] = 2 * np.random.rand(1, ny) - 1
        y[i,:ny] = y_i

    cdef list n=[1]*random.randint(1,dw)
    ## Generating switching sequence
    while len(n)<T:
        chance=random.random()
        mode=random.randint(1, modes)
        n=n+[mode]*dw
        while chance>0.5:
            n=n+[mode]
            chance=random.random()

    ## Generate partition of regressorspace
    if pwarx==1:
        n=[modes+1 for x in n]
        for i in range(modes-1):
            newmode = 10*np.random.rand(1, nu * inputs + ny * outputs+1)-5
            if i > 0:
                H[i] = np.vstack((newmode, oldmodes))
                oldmodes = np.vstack((oldmodes, -newmode))
            else:
                H[i] = newmode
                oldmodes = -newmode
        H[modes-1]=oldmodes
    for i in range(max(ny,nu),T):
        ## Obtain the input and output values used to calculate the next time step
        ux = np.hstack(tuple([u[a][i - nu:i] for a in range(inputs)]))
        uy = np.hstack(tuple([y[a][i - ny:i] for a in range(outputs)]))

        ## Simulate the model
        if pwarx==1:
            for j in range(modes):
                er=np.dot(H[j], np.append(np.append(uy, ux), 1))
                if (er<0).all():
                    y[:, i] = np.dot(thetaR[j], np.append(uy, np.append(ux,1))) + float(nt * (np.random.rand(1) - 0.5))
                    n[i]=j+1
        else:
            y[:, i] = np.dot(thetaR[n[i] - 1], np.append(uy, ux)) + float(nt * (np.random.rand(1) - 0.5))

    with open("Data/TestfileInput.txt", "wb") as f:
        np.savetxt(f, u, fmt='%f9', delimiter=",")

    with open("Data/TestfileOutput.txt", "wb") as f:
        np.savetxt(f, y, fmt='%f9', delimiter=",")

    return (thetaR,n,H)
