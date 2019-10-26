from sklearn import svm
import numpy as np
cimport numpy as np
import matplotlib.pyplot as plt
import math


cpdef mesh(np.ndarray X):
    cdef np.ndarray x=np.arange(X[:,0].min()-1, X[:,0].max()+1, 0.05)
    cdef np.ndarray y=np.arange(X[:,1].min()-1, X[:,1].max()+1, 0.05)
    xx, yy = np.meshgrid(x, y, sparse=False)
    return xx,yy

def heightmap(xx,yy,H,int modes):
    cpdef np.ndarray Z=np.zeros((len(xx[:,0]),len(xx[0,:])))
    for i in range(len(xx[:,0])):
        for j in range(len(xx[0,:])):
            for k in range(modes):
                if (np.dot(H[k],np.array((xx[i,j],yy[i,j],1)))>=0).all():
                    Z[i,j]=k+1
    return(Z)

cpdef constructH(H,int modes):

    Hs = []
    for i in range(modes):
        Hs.append(np.empty((0,len(H[0]))))
    for j in range(modes):
        for i in range(modes):
            if j != i:
                Hi = H[j] - H[i]
                Hs[j]=np.vstack((Hs[j],Hi))
        Hs[j]=-Hs[j]
    return Hs

def createLabels(Y):
    L=[]
    for i in Y:
        L.append("Mode {}".format(i))
    return L
def partition(X,Y,modes):
    ## check to determine the number of modes is larger than 1
    if modes==1:
        return []
    X=X[:,:-1]

    ##Initialize weights
    Yw=np.array(Y)
    for i in range(max(Y)+1):
        w1=(Yw == i).sum()
        w2=float(w1)
        weight={i: len(Y)/math.sqrt(math.sqrt(w2))}
    ## Initialize a SVM
    clf = svm.LinearSVC(max_iter = 10000)
    clf.fit(X,Y)

    ## Construct proper H
    H = np.hstack((np.array(clf.coef_),np.transpose(np.array([clf.intercept_]))))
    if modes==2:
        H=[H[0],-H[0]]

    H2 = constructH(H, modes)
    ## Plot contour map
    if len(X[0, :]) == 2:
        xx, yy = mesh(X)
        Z = heightmap(xx, yy, H2, modes)
        plt.figure(1)
        plt.contourf(xx, yy, Z, alpha=0.6, cmap='gist_rainbow')
        plt.scatter(X[:, 0], X[:, 1], c=Y,s=50, cmap='plasma')
        ax=plt.gca()
        ax.set_xlim([X[:,0].min()-1, X[:,0].max()+0.9])
        ax.set_ylim([X[:,1].min()-1, X[:,1].max()+0.9])
        plt.show()
    return H2