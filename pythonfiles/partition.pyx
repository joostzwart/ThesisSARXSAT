from sklearn import svm
import numpy as np
cimport numpy as np
import matplotlib.pyplot as plt


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
    Ht=[]
    cdef int i
    cdef int count1=0
    cdef int count2=1
    cdef int count3=1
    for i in range(modes):
        Ht.append(np.empty((0,len(H[0])),dtype=float))


    for i in range(len(H[:,0])):
        Ht[count1]=np.append(Ht[count1],np.expand_dims(H[i,:],axis=0), axis=0)
        Ht[count2]=np.append(Ht[count2],np.expand_dims(-H[i, :],axis=0),axis=0)
        if count2<modes-1:
            count2=count2+1
        else:
            count1=count1+1
            count2=count3+1
            count3=count3+1
    return Ht

def createLabels(Y):
    L=[]
    for i in Y:
        L.append("Mode {}".format(i))
    return L
def partition(X,Y,modes):
    ## check to determine the number of modes is larger than 1
    if max(Y)==1:
        return []
    X=X[:,:-1]
    Hs = {}
    ys = {}
    ## Initialize a SVM
    clf = svm.SVC(kernel ='linear', C = 1)
    clf.fit(X,Y)


    H = np.hstack((np.array(clf.coef_),np.transpose(np.array([clf.intercept_]))))
    ## Construct proper H
    H2 = constructH(H, modes)

    ## Plot contour map
    if len(X[0, :]) == 2:
        xx, yy = mesh(X)
        Z = heightmap(xx, yy, H2, modes)
        plt.contourf(xx, yy, Z, alpha=0.6, cmap='gist_rainbow')

        ## Plot lines
        x=np.linspace(X[:,0].min()-1,X[:,0].max()+1,1000)
        for i in range(len(H)):
            ys[i]= -H[i,0] / H[i,1] * x - H[i,2] / H[i,1]
            plt.plot(x, ys[i], color='green')
        L=createLabels(Y)
        plt.scatter(X[:, 0], X[:, 1], c=Y,s=50, cmap='plasma')
        ax=plt.gca()
        ax.set_xlim([X[:,0].min()-1, X[:,0].max()+0.9])
        ax.set_ylim([X[:,1].min()-1, X[:,1].max()+0.9])
        plt.show()
    return H2