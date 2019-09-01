from datetime import datetime
import scipy.io
import numpy as np
cimport numpy as np
import os

def exportToMatlab(list Fmodels,list sigmat,np.ndarray Ut,np.ndarray Yt,np.ndarray Rt,H,pwarx):
    """
    This function exports information to .mat files for further analysis in matlab.

    Arguments:
        Fmodels {list} -- Final models after identification
        sigmat {list} -- Switching sequence
        Ut {array} -- Input
        Yt {array} -- Output
        Rt {array} -- regressor
    """
    time = datetime.now()
    path = "Data/IdentifiedDatasets/" + time.strftime('%m_%d_%Y_%H_%M_%S')
    try:
        os.mkdir(path)
    except OSError:
        print ("Creation of the directory %s failed" % path)
    else:
        print ("Successfully created the directory %s " % path)

    scipy.io.savemat(path + "/Fmodels.mat", mdict={'Fmodels': Fmodels})
    scipy.io.savemat(path + "/Sigmat.mat", mdict={'sigmat': sigmat})
    scipy.io.savemat(path + "/Ut.mat", mdict={'Ut': Ut})
    scipy.io.savemat(path + "/Yt.mat", mdict={'Yt': Yt})
    scipy.io.savemat(path + "/Rt.mat", mdict={'Rt': Rt})
    if pwarx==1:
        scipy.io.savemat(path + "/H.mat", mdict={'H': H})

