###PROV###
def random_upper_tri(FF,n):
    """
    This function returns a random upper triangular matrix of size (n x n) in the finite field FF.
    """
    M = matrix(FF,n,n)
    
    for i in range(n):
        for j in range(i,n):
            M[i,j] = FF.random_element()
            
    return M

def Upper(FF,M):
    """
    This function take as input a matrix representing a quadratic form P and returns an upper triangular matrix Upper_M representing the same quadratic form P.
    """
    row = M.nrows()
    col = M.ncols()
    Upper_M = matrix(FF,row,col)
    for i in range(row):
        Upper_M[i,i] = M[i,i]
        for j in range(i+1,col):
            Upper_M[i,j] = M[i,j] + M[j,i]
    
    return Upper_M

def KeyGen(q,m,delta,v, verbose = False):
    """
    This function generates a key pair (O,P) for PROV parameters q,m,delta,v where m is the number of public key polynomials, m+delta the dimension of the secret linear subspace O 
    where the public key polynomial vanish and v the number of vinegar variable. This algorithm follows how the private key are generated in PROV specification.
    Important notice: This code is a demonstration tool and should not used be for applications where security matters. 
    """
    n = m+delta+v
    FF = GF(q)
    
    # generate matrix O
    O = random_matrix(FF,v,m+delta)
    col_O = O.ncols()
    Id = identity_matrix(FF,col_O)
    Oe = block_matrix([[O],[Id]])

    P = []
    for e in range(m):
        P1 = random_upper_tri(FF,v)
        P2 = random_matrix(FF,v,m+delta)
        P3 = O.transpose() * P1 * O + O.transpose()*P2
        P3 = Upper(FF,P3)
        P4 = matrix(FF,m+delta,v)
        Pe = block_matrix([[P1,P2],[P4,P3]])
        P.append(Pe)
        
    return (O,Oe), P


def oracle(O) :
    """
    This function takes as input the secret key and returns a vector o of the secret subspace represented by the basis O.
    """
    FF = O.base_ring()
    row_O = O.nrows()
    col_O = O.ncols()
    n = row_O + col_O
    Id = identity_matrix(FF,col_O)
    Oe = block_matrix([[O],[Id]])
    o = matrix(FF,n,1)
    while(equal_zero(o)):
        x = random_matrix(FF,col_O,1)
        o = Oe*x

    return o

def oracle2(Oe) :
    """
    This function takes as input the secret key and returns a vector o of the secret subspace represented by the basis Oe.
    """
    FF = Oe.base_ring()
    col_O = O.ncols()
    
    o = matrix(FF,n,1)
    while(equal_zero(o)):
        x = random_matrix(FF,col_O,1)
        o = Oe*x
    
    return o



def verify_subspace(P,Oe):
    """
    This function verifies if the multivariate quadratic map P vanish on the linear subspace O defined by the basis Oe.
    """
    row = P[0].nrows()
    col = Oe.ncols()

    FF = Oe.base_ring()

    P_bis = []
    Zeros = []
    for p in P:
        Zeros.append(matrix(FF,col,col))
        P_bis.append(Oe.transpose()*p*Oe)
        
    b = equality_poly(FF,P_bis,Zeros)

    return True

def equality_poly(FF,P1,P2):
    """
    This function verifies with a naïve algorithm that two multivariate quadratic maps are equal.
    """
    l1 = len(P1)
    l2 = len(P2)
    
    if l1 != l2:
        return False

    n = P1[0].nrows()
    
    for i in range(10):
        x = random_matrix(FF,n,1)
        for j in range(l1):
            res1 = x.transpose()*P1[j]*x
            res2 = x.transpose()*P2[j]*x
            if res1 != res2:
                return False

    return True

def equal_zero(M):
    """
    This function verifies if the matrix M equals zero and returns True if this is the case and False otherwise.
    """
    FF = M.base_ring()
    row = M.nrows()
    col = M.ncols()
    Zero = matrix(FF,row,col)

    return Zero == M