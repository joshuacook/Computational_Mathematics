def KTBC(type='K',n=10,sparse=False):
    # KTBC Create finite difference model matrix.
    # K=KTBC(TYPE,N,SPARSE) creates model matrix TYPE of size N-by-N.
    # TYPE is one of the characters 'K', 'T', 'B', or 'C'.
    # The command K = KTBC('K', 100, 1) gives a sparse representation
    # K=KTBC uses the defaults TYPE='K', n=10, and SPARSE=False.
    # Change the 3rd argument from 1 to 0 for dense representation!
    # If no 3rd argument is given, the default is dense
    # If no argument at all, KTBC will give 10 by 10 matrix K
    e = np.ones(n)
    e_off = np.ones(n-1)
    K = sprs.diags([-e_off,2*e,-e_off],[-1,0,1])

    if type == 'K': 
        K = K
    if type == 'T':
        K = sprs.csr_matrix(K)
        K[0,0] = 1
    if type == 'B': 
        K = sprs.csr_matrix(K)
        K[0,0] = 1
        K[n-1,n-1] = 1
    if type == 'C':
        K = sprs.csr_matrix(K)
        K[0,n-1] = -1
        K[n-1,0] = -1


    if sparse == False:
        return K.todense()
    else:
        return K