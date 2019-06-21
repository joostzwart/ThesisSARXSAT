import numpy as np
cimport numpy as np
import random

cpdef generate(int T,int ny,int nu,int inputs,int outputs,double nt,int modes,int dw,seed,int pwarx):
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

    ## Generate random stable discrete systems
    cdef np.ndarray thetaY=(2/ny*outputs)*(np.random.rand(modes,outputs,ny*outputs)-0.5)
    cdef np.ndarray thetaU=2*np.random.rand(modes,outputs,nu*inputs)-1
    cdef np.ndarray thetaR=0.5*np.dstack((thetaY,thetaU))
    cdef list u = []
    H={}
    cdef int m=0
    cdef int i
    cdef np.ndarray Temp, Temp2
    ## Generating uniformly distributed random input
    for i in range(inputs):
        [u_i] = 10*(2 * np.random.rand(1, T) - 1)
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
            Temp = 10*np.random.rand(1, nu * inputs + ny * outputs+1)-5
            if i > 0:
                H[i] = np.vstack((Temp, -Temp2))
                Temp2 = np.vstack((Temp, H[i - 1]))
            else:
                H[i] = Temp
                Temp2 = Temp
        H[i+1]=-Temp2
    for i in range(max(ny,nu),T):
        ## Obtain the input and output values used to calculate the next time step
        ux = np.hstack(tuple([u[a][i - nu:i] for a in range(inputs)]))
        uy = np.hstack(tuple([y[a][i - ny:i] for a in range(outputs)]))

        ## Simulate the model
        if pwarx==1:
            for j in range(modes):
                er=np.dot(H[j], np.append(np.append(uy, ux), 1))
                print("err",er)
                if (er<0).all():
                    y[:, i] = np.dot(thetaR[j], np.append(uy, ux)) + float(nt * (np.random.rand(1) - 0.5))
                    n[i]=j+1
        else:
            y[:, i] = np.dot(thetaR[n[i] - 1], np.append(uy, ux)) + float(nt * (np.random.rand(1) - 0.5))

    with open("Data/TestfileInput.txt", "wb") as f:
        np.savetxt(f, u, fmt='%f9', delimiter=",")

    with open("Data/TestfileOutput.txt", "wb") as f:
        np.savetxt(f, y, fmt='%f9', delimiter=",")

    return (thetaR,n)
