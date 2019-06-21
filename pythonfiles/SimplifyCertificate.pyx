cpdef simplify(list cert):
    """
    Removes redundant certificates
    
    Arguments:
        cert {list} -- List of certificates.
    
    Returns:
        cert {list} -- Reduced list of certificates.
    """

    cdef int i
    cdef list C1, C2, a
    for C1 in cert:
        for i,C2 in enumerate(cert):
            if C1 is not C2:
                if C2[0]<=C1[0] and C2[1]>=C1[1]:
                    cert[i]=[0,0]
    
    cert=[a for a in cert if a != [0,0]]
    return cert