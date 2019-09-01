cdef class Parameters:
    """Class that keeps track of user parameters"""
    cdef public int T, ny, nu, dw, inputs, outputs, nod, pwarx, Ls
    cdef public int chunks, merging, splitlargedataset, modelgeneration, unstuck
    cdef public double delta, nt
    cdef public set oldmodels

    def __init__(self,int T,int ny,int nu,double delta,int dw,double nt,int inputs,
                 int outputs,int chuncks,oldmodels,int nod, int pwarx,int Ls):
        self.T=T
        self.ny=ny
        self.nu=nu
        self.delta=delta
        self.dw=dw
        self.nt=nt
        self.inputs=inputs
        self.outputs=outputs
        self.chunks=chuncks
        self.oldmodels=oldmodels
        self.nod=nod
        self.pwarx=pwarx
        self.Ls=Ls

class Preferences:
    def __init__(self,int merging,int splitlargedataset,int modelgeneration,seed,int unstuck):
        self.merging = merging
        self.splitlargedataset = splitlargedataset
        self.modelgeneration=modelgeneration
        self.seed=seed
        self.unstuck=unstuck


