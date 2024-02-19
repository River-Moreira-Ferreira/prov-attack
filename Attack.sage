"""
This file contains an implementation of the attack described in the paper "Polynomial-Time Key-Recovery Attack on the PROV Specification", by River Moreira Ferreira and Ludovic Perret in 2024.
This code is inspired from the implementation provided with the paper "One vector to rule them all: Key recovery from one vector in UOV schemes" by Pierre Pébereau in 2023.
"""

def attack_prov(P,o,delta) :
    """ 
    This function given a PROV public and a vector o in the secret linear subspace O return a basis of O.
    """
    m = len(P)
    
    P_o = []
    for p in P:   
        P_o.append([o.transpose()*(p+p.transpose())])
    J = block_matrix(P_o)
    B = matrix(J.right_kernel().basis())
    
    #print("Upper bound of the rank of the restriced public polynomial (ghat) = ",2*(n-2*m-delta))
    charac = P[0].base_ring().characteristic()
    B2 = []
    for g in P :
        ghat = B*g*B.transpose()  #restriction of G to the kernel of J
        if charac == 2:
            ghat = ghat + ghat.transpose()
        for b in ghat.kernel().basis() :
            if len(B2) == 0 or b not in span(B2) :
                B2.append(b)
        if len(B2) == m + delta :
            break
    B3 = matrix(B2)
    C = B3*B
    return C


def int_converter(str):
    """
    This function converts a string into a integer and if the string does not represent an integer returns 1.
    """
    try:
        return int(str)
    except ValueError:
        return 1

def get_security_level(argv):
    """
        This function takes as the command line arguments and returns the security level choosen by the user or 1 by default.
    """
    sec_level = 1
    if(len(argv)==2):
        tmp = int_converter(argv[1])
        if(tmp == 3):
            sec_level = tmp
        if(tmp == 5):
            sec_level = tmp

    return sec_level


load("PROV.sage")
import sys

sec_level = get_security_level(sys.argv)
q, m, delta, v= 256, 46, 8, 82 # Level 1
if(sec_level == 3):
    q, m, delta, v= 256, 70, 8, 122 # Level 3
if(sec_level == 5):
    q, m, delta, v= 256, 96, 8, 160 # Level 5

n = m + delta + v

print("Security level:",sec_level)
print('Parameters : q=',q,' m=',m,' delta=',delta,' v=',v,' n=',n)

print('Key Generation')
(O,Oe), P = KeyGen(q, m, delta, v)

b = verify_subspace(P,Oe)
if b == False:
    print("The public polynomials does not vanish on the linear subspace defined by the private key.\n")
else:
    print("The public polynomials vanish on the linear subspace defined by the private key.\n")

print("Run the attack")
o = oracle(O)

import time 
a = time.time()
B = attack_prov(P,o,delta)
b = time.time()
print("Attack performed in",round(b-a,2),"seconds")

print('Dimension of the subspace we found ',B.dimensions()[0], ' vs dimension of O', m+delta)

b = verify_subspace(P,B.transpose())
if b == False:
    print("The public polynomials does not vanish on the linear subspace we found.")
else:
    print("The public polynomials vanish on the linear subspace we found.")