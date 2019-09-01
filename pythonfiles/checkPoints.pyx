from ppl import Variable, Generator_System, C_Polyhedron, point
from scipy.spatial import ConvexHull
import numpy as np
cimport numpy as np


def erosion(P,epsilon):
    generator=P.generators()
    D=[]
    dimensions = P.space_dimension()
    for i in range(dimensions):
        D.append(Variable(i))

    for i in range(dimensions):
        p_shifted=C_Polyhedron(generator)
        p_shifted.affine_image(D[i], D[i] - epsilon)
        P.intersection_assign(p_shifted)
        p_shifted=C_Polyhedron(generator)
        p_shifted.affine_image(D[i], D[i] + epsilon)
        P.intersection_assign(p_shifted)

    P.minimized_constraints()
    return P

cpdef intersect(np.ndarray group_1,np.ndarray group_2,delta):
    """Checks if two groups of points overlap

    Arguments:
        P {array} -- First group of points
        Q {array} -- Second group of points

    Returns:
        status {int} -- status: 1 if the groups do not overlap
    """
    group_1=group_1[:,:-1]
    group_2=group_2[:,:-1]
    cpdef int status=1
    cpdef list hull
    cpdef list P
    cpdef list Q
    if group_1.shape[0]>3:
        hull=list(ConvexHull(group_1).vertices)
        P=group_1[hull].tolist()
    else:

        P=group_1.tolist()
    if group_2.shape[0]>3:
        hull=list(ConvexHull(group_2).vertices)
        Q=group_2[hull].tolist()
    else:
        Q=group_2.tolist()
    cpdef int dimensions = len(P[0])
    cpdef list D=[]

    ## Create one variable per dimension
    for i in range(dimensions):
        D.append(Variable(i))

    ## Create first polyhedron
    gsP = Generator_System()
    for i in range(len(P)):
        point_i=P[i][0]*D[0]
        for j in range(1,dimensions):
            point_i=point_i+P[i][j]*D[j]
        gsP.insert(point(point_i))
    polyP = C_Polyhedron(gsP)

    ## erode first polyhedron
    polyP=erosion(polyP,delta)

    ## Create second polyhedron
    gsQ = Generator_System()
    for i in range(len(Q)):
        point_i=Q[i][0]*D[0]
        for j in range(1,dimensions):
            point_i=point_i+Q[i][j]*D[j]
        gsQ.insert(point(point_i))
    polyQ = C_Polyhedron(gsQ)

    ## erode second polyhedron
    polyQ=erosion(polyQ,delta)
    ## Check for intersection
    polyP.intersection_assign(polyQ)
    if polyP!=C_Polyhedron(dimensions, 'empty'):
        status=2
    return(status)

