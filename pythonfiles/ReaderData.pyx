import numpy as np
cimport numpy as np
def Read(int nu,int ny,int modelgeneration,**kwargs):

    cdef np.ndarray input, output,  ux, yx
    cdef int T, inputs, outputs, maxorder

    ## Checking for user defined dataset
    if 'input' in kwargs:
        input_path=kwargs.get('input')
    else:
        input_path="Data/DatasetInput_Ts_20(J).txt"

    if 'output' in kwargs:
        output_path=kwargs.get('output')
    else:
        output_path="Data/DatasetOutput_Ts_20(J).txt"

    ## Reading textfiles
    if modelgeneration==2:
        output=np.genfromtxt(output_path, delimiter=',')
        input=np.genfromtxt(input_path, delimiter=',')
    else:
        output = np.genfromtxt("Data/TestfileOutput.txt", delimiter=',')
        input = np.genfromtxt("Data/TestfileInput.txt", delimiter=',')
    if input.ndim==1:
        input=np.array([input])

    if output.ndim==1:
        output=np.array([output])
    ## Obtaining number of inputs/outputs
    inputs = len(input)
    outputs = len(output)
    T = len(output[0])
    maxorder=max(ny, nu)
    ## Creating regressor
    r=[]

    for k in range(T - maxorder):
        ux = np.hstack(tuple([input[a][k - nu + maxorder:k + maxorder] for a in range(inputs)]))
        yx = np.hstack(tuple([output[a][k-ny+maxorder:k + maxorder] for a in range(outputs)]))
        r.append(np.append(yx, ux))

    r = np.array(r)

    output=output[:,maxorder:]
    return (input, output, inputs, outputs, T, r)