import numpy as np
from scipy.linalg import blas as FB
#INDEX - Use ctrl+F to browse faster
#(1) - BINARY OPERATIONS
#(2) - SPIN OPERATIONS
#(3) - OUT-OF-TIME-ORDERED CORRELATORS
##(3.1) - WITH TEMPERATURE CHOICES (IF T=INFTY USE INFTY TEMPERATURE OTOCS FOR EXTRA SPEED)
##(3.2) - INFTY TEMPERATURE OTOCS

#----------------  (1) BINARY OPERATIONS  ----------------#
#Translations from fortran --> python

#btest (https://gnu.huihoo.org/gcc/gcc-7.1.0/gfortran/BTEST.html)
def btest(i, pos): 
    return bool(i & (1 << pos))

#ibclr/ibset (https://gcc.gnu.org/onlinedocs/gcc-4.6.1/gfortran/IBCLR.html / http://www.lahey.com/docs/lfpro78help/F95ARIBSETFn.htm)
def set_bit(v, index, x): #if x=0 --> ibclr / if x=1 --> ibset 
  """Set the index:th bit of v to 1 if x is truthy, else to 0, and return the new value."""
  mask = 1 << index   
  v &= ~mask
  if x:
    v |= mask
  return v

#----------------  (2) SPIN OPERATIONS  ----------------#

# Pauli at site operators

def S_xi(pos_i,sites):
    dim = 2**sites
    S = np.zeros((dim,dim),dtype=complex)
    estados2 = np.zeros(dim,dtype=np.int64)
    for i in range(dim):
        if btest(i,pos_i) == True:
            estados2[i] = set_bit(i,pos_i,0)
        else:
            estados2[i] = set_bit(i,pos_i,1)
    for i in range(dim):
        for j in range(dim):
            if i == estados2[j]:
                S[i,j] = S[i,j]+1
    return S

def S_yi(pos_i,sites):
    dim = 2**sites
    S = np.zeros((dim,dim),dtype=complex)
    estados2 = np.zeros(dim,dtype=np.int64)
    a = np.zeros(dim,dtype=complex)
    for i in range(dim):
        if btest(i,pos_i) == True:
            estados2[i] = set_bit(i,pos_i,0)
            a[i] = 1j
        else:
            estados2[i] = set_bit(i,pos_i,1)
            a[i] = -1j
    for i in range(dim):
        for j in range(dim):
            if i == estados2[j]:
                S[i,j] = S[i,j]+a[j]
    return S
    
def S_zi(pos_i,sites):
    dim = 2**sites
    S = np.zeros((dim,dim),dtype=complex)
    for i in range(dim):
        if btest(i,pos_i) == True:
            S[i,i] = 1
        else:
            S[i,i] = -1
    return S

# Neighbor i-j site interactions

def S_xxij(pos_i,pos_j,sites):
    dim = 2**sites
    S = np.zeros((dim,dim),dtype=complex)
    estados2 = np.zeros(dim,dtype=np.int64)
    for i in range(dim):
        if btest(i,pos_i) == True:
            estados2[i] = set_bit(i,pos_i,0)
        else:
            estados2[i] = set_bit(i,pos_i,1)
    for i in range(dim):
        if btest(estados2[i],pos_j) == True:
            estados2[i] = set_bit(estados2[i],pos_j,0)
        else:
            estados2[i] = set_bit(estados2[i],pos_j,1)
    for i in range(dim):
        for j in range(dim):
            if i == estados2[j]:
                S[i,j] = S[i,j]+1
    return S

def S_yyij(pos_i,pos_j,sites):
    dim = 2**sites
    S = np.zeros((dim,dim),dtype=complex)
    estados2 = np.zeros(dim,dtype=np.int64)
    a = np.zeros(dim,dtype=complex)
    for i in range(dim):
        if btest(i,pos_i) == True:
            estados2[i] = set_bit(i,pos_i,0)
            a[i] = 1j
        else:
            estados2[i] = set_bit(i,pos_i,1)
            a[i] = -1j
    #estados2 = estados2.astype(np.int64)
    for i in range(dim):
        if btest(estados2[i],pos_j) == True:
            estados2[i] = set_bit(estados2[i],pos_j,0)
            a[i] = a[i]*1j
        else:
            estados2[i] = set_bit(estados2[i],pos_j,1)
            a[i] = -a[i]*1j
    for i in range(dim):
        for j in range(dim):
            if i == estados2[j]:
                S[i,j] = S[i,j]+a[j]
    return S

def S_zzij(pos_i,pos_j,sites):
    dim = 2**sites
    S = np.zeros((dim,dim),dtype=complex)
    for i in range(dim):
        if btest(i,pos_i) == True:
            if btest(i,pos_j) == True:
                S[i,i] = 1
            else:
                S[i,i] = -1
        else:
            if btest(i,pos_j) == True:
                S[i,i] = -1
            else:
                S[i,i] = 1
    return S
    
# Entire spin direction operations

def S_x(sites): # Entire S_x operator
    """ Defines the operator S_x = sum S_x_i where i is the position between (0,N-1) and N = number of sites."""
    dimension = 2 ** sites 
    H = np.zeros((dimension,dimension),dtype=complex)
    for l in range(1,dimension+1):
        for j in range(0,sites):
            if btest(l-1,j) == True:
                hh = set_bit(l-1,j,0) + 1
            if btest(l-1,j) == False:
                hh = set_bit(l-1,j,1) + 1
            H[l-1,hh-1]=H[l-1,hh-1]+1
    return H

def S_y(sites): # Entire S_y operator
    """ Defines the operator S_y = sum S_y_i where i is the position between (0,N-1) and N = number of sites."""
    dimension = 2 ** sites 
    H = np.zeros((dimension,dimension),dtype=complex)
    for l in range(1,dimension+1):
        for u in range(0,sites):
            if btest(l-1,u) == True:
                hh = set_bit(l-1,u,0) + 1
                H[l-1,hh-1]=H[l-1,hh-1] + 1j
            if btest(l-1,u) == False:
                hh = set_bit(l-1,u,1) + 1
                H[l-1,hh-1]=H[l-1,hh-1] - 1j
    return H

def S_z(sites): # Entire S_z operator
    """ Defines the operator S_z = sum S_z_i where i is the position between (0,N-1) and N = number of sites."""
    dimension = 2 ** sites 
    H = np.zeros((dimension,dimension),dtype=complex)
    for l in range(1,dimension+1):
        for j in range(0,sites):
            if btest(l-1,j) == True:
                hh = set_bit(l-1,j,1) + 1
                H[l-1,hh-1]=H[l-1,hh-1] - 1
            if btest(l-1,j) == False:
                hh = set_bit(l-1,j,0) + 1
                H[l-1,hh-1]=H[l-1,hh-1] + 1
    return H
    
def Parity(sites): # Returns the parity operator in the S_z base #### HAY QUE REHACERLO!!!!!! ESTA MAL! (PERO NO MUY MAL)
    dim = 2 ** sites
    identity = np.zeros((dim,dim),dtype=complex)
    zeros = np.zeros((4,4),dtype=complex)
    for n in range(dim):
        identity[n,n] = 1
    Dx = np.matrix([[0,1],[1,0]],dtype=complex)
    Dy = np.matrix([[0,-1j],[1j,0]],dtype=complex)
    Dz = np.matrix([[1,0],[0,-1]],dtype=complex)
    Parity = identity
    if sites % 2 == 0:
        tt = sites/2
    else:
        tt = (sites-1)/2
    for u in range(tt):
        P = 0.5* (identity + np.matmul(S_x_i(u,sites),S_x_i(sites-1-u,sites)) + np.matmul(S_y_i(u,sites),S_y_i(sites-1-u,sites)) + np.matmul(S_z_i(u,sites),S_z_i(sites-1-u,sites)) )
        Parity = np.matmul(Parity,P)
    return Parity
    
#----------------  (3) OUT-OF-TIME-ORDERED CORRELATORS  ----------------#

### (3.1) WITH TEMPERATURE CHOICES (IF T=INFTY USE INFTY TEMPERATURE OTOCS FOR EXTRA SPEED)

#OTOC
def OTOC(V,W,ener,basis,N,dt,t0,beta,ortho):
    dim = len(ener)
    if ortho == True:
        basist = np.transpose(basis)
    else:
        basist = np.linalg.inv(basis)
    ope0 = np.matmul(basist,V)
    ope0 = np.matmul(ope0,basis)
    ope = np.matmul(basist,W)
    ope = np.matmul(ope,basis)
    otoc=[]
    mm = np.zeros((dim,dim),dtype=complex)
    mm1 = np.zeros((dim,dim),dtype=complex)
    ensemblemedio = dim
    ensemble = np.identity(dim,dtype=complex)
    if beta =! 0:
        for i in range(0,dim):
            ensemble[i][i] = np.exp(-beta*ener[i])
        ensemblemedio = np.abs(np.matrix.trace(ensemble))
    U = np.zeros((dim,dim),dtype=complex) # U evolucion temporal
    Udagger = np.zeros((dim,dim),dtype=complex) # U*   
    for count0 in range(0,N):
        tie = count0*dt+t0
        tiempo.append(tie)
        for count1 in range(0,dim):
            U[count1][count1] = np.exp(-1j*tie*ener[count1])
            Udagger[count1][count1] = np.exp(1j*tie*ener[count1])
        mm = np.matmul(ope0,U)
        mm1= np.matmul(Udagger,mm)
        mm = np.matmul(mm1,ope) - np.matmul(ope,mm1)
        mm = np.matmul(np.transpose(np.conjugate(mm)),mm)
        mm = np.matmul(ensemble,mm)
        otoc.append(0.5*np.abs(np.matrix.trace(mm)) / ensemblemedio )
    return otoc

# OTOC 1 - Re(F) - USE ONLY IF V,W ARE BOTH UNITARY & HERMITIAN
def OTOCF(V,W,ener,basis,N,dt,t0,beta,ortho):
    dim = len(ener)
    if ortho == True:
        basist = np.transpose(basis)
    else:
        basist = np.linalg.inv(basis)
    ope0 = np.matmul(basist,V)
    ope0 = np.matmul(ope0,basis)
    ope = np.matmul(basist,W)
    ope = np.matmul(ope,basis)
    otoc=[]
    mm = np.zeros((dim,dim),dtype=complex)
    mm1 = np.zeros((dim,dim),dtype=complex)
    ensemblemedio = dim
    ensemble = np.identity(dim,dtype=complex)
    if beta =! 0:
        for i in range(0,dim):
            ensemble[i][i] = np.exp(-beta*ener[i])
        ensemblemedio = np.abs(np.matrix.trace(ensemble))
    U = np.zeros((dim,dim),dtype=complex) # U evolucion temporal
    Udagger = np.zeros((dim,dim),dtype=complex) # U*   
    for count0 in range(0,N):
        tie = count0*dt+t0
        tiempo.append(tie)
        for count1 in range(0,dim):
            U[count1][count1] = np.exp(-1j*tie*ener[count1])
            Udagger[count1][count1] = np.exp(1j*tie*ener[count1])
        mm = np.matmul(ope0,U)
        mm1= np.matmul(Udagger,mm)
        mm = np.matmul(mm1,ope)
        mm = np.matmul(mm,mm)
        mm = np.matmul(ensemble,mm)
        otoc.append(1 - (np.matrix.trace(mm)/ensemblemedio).real )
    return otoc

### (3.2) INFTY TEMPERATURE OTOCS

# OTOC OPTMIZED SPEED
def OTOC_opt_infty(V,W,ener,basis,N,dt,t0,ortho):
    basis = np.complex64(basis, order = 'F')
    V = np.complex64(V, order = 'F')
    W = np.complex64(W, order = 'F')
    dim = len(ener)
    if ortho == True:
        basist = np.transpose(basis)
    else:
        basist = np.linalg.inv(basis)
    S0 = FB.cgemm(1,V,basis)
    S0 = FB.cgemm(1,basist,S0)
    S = FB.cgemm(1,W,basis)
    S = FB.cgemm(1,basist,S)
    mm = np.zeros((dim,dim),dtype="complex64",order='F')
    mm1 = np.zeros((dim,dim),dtype="complex64",order='F')
    otok = np.zeros(N,dtype = "float32")
    U = np.zeros((dim,dim),dtype="complex64",order='F') # U evolucion temporal
    Udagger = np.zeros((dim,dim),dtype="complex64",order='F') # U*
    for ti in range(0,N):
        for c1 in range(0,dim):
            U[c1][c1] = np.exp(-1j*tiempo[ti]*ener[c1])
            Udagger[c1][c1] = np.exp(1j*tiempo[ti]*ener[c1])
        mm = FB.cgemm(1,S0,U)
        mm1= FB.cgemm(1,Udagger,mm)
        mm = FB.cgemm(1,mm1,S) - FB.cgemm(1,S,mm1)
        mm = FB.cgemm(1,np.transpose(np.conjugate(mm)),mm)
        otok[ti] = 0.5*(np.abs(np.matrix.trace(mm))/dim)
    otok = np.array(otok)
    return otok
    
# OTOC 1-Re(F) OPTMIZED SPEED - USE ONLY IF V,W ARE BOTH UNITARY & HERMITIAN
def OTOCF_opt_infty(V,W,ener,basis,N,dt,t0,ortho):
    basis = np.complex64(basis, order = 'F')
    V = np.complex64(V, order = 'F')
    W = np.complex64(W, order = 'F')
    dim = len(ener)
    if ortho == True:
        basist = np.transpose(basis)
    else:
        basist = np.linalg.inv(basis)
    S0 = FB.cgemm(1,V,basis)
    S0 = FB.cgemm(1,basist,S0)
    S = FB.cgemm(1,W,basis)
    S = FB.cgemm(1,basist,S)
    mm = np.zeros((dim,dim),dtype="complex64",order='F')
    mm1 = np.zeros((dim,dim),dtype="complex64",order='F')
    otok = np.zeros(N,dtype = "float32")
    U = np.zeros((dim,dim),dtype="complex64",order='F') # U evolucion temporal
    Udagger = np.zeros((dim,dim),dtype="complex64",order='F') # U*
    for ti in range(0,N):
        for c1 in range(0,dim):
            U[c1][c1] = np.exp(-1j*tiempo[ti]*ener[c1])
            Udagger[c1][c1] = np.exp(1j*tiempo[ti]*ener[c1])
        mm = FB.cgemm(1,S0,U)
        mm1= FB.cgemm(1,Udagger,mm)
        mm = FB.cgemm(1,mm1,S)
        mm = FB.cgemm(1,mm,mm)
        otok[ti] = 1 - (np.matrix.trace(mm)/len(ener)).real 
    otok = np.array(otok)
    return otok
    
# OTOC SLOW (IF HIGHER PRECISION REQUIRED)
def OTOC_infty(V,W,ener,basis,N,dt,t0,ortho):
    dim = len(ener)
    if ortho == True:
        basist = np.transpose(basis)
    else:
        basist = np.linalg.inv(basis)
    S0 = np.matmul(V,basis)
    S0 = np.matmul(basist,S0)
    S = np.matmul(W,basis)
    S = np.matmul(basist,S)
    mm = np.zeros((dim,dim),dtype=complex)
    mm1 = np.zeros((dim,dim),dtype=complex)
    otok = np.zeros(N,dtype = float)
    U = np.zeros((dim,dim),dtype=complex) # U evolucion temporal
    Udagger = np.zeros((dim,dim),dtype=complex) # U*
    for ti in range(0,N):
        for c1 in range(0,dim):
            U[c1][c1] = np.exp(-1j*tiempo[ti]*ener[c1])
            Udagger[c1][c1] = np.exp(1j*tiempo[ti]*ener[c1])
        mm = np.matmul(S0,U)
        mm1= np.matmul(Udagger,mm)
        mm = np.matmul(mm1,S) - np.matmul(S,mm1)
        mm = np.matmul(np.transpose(np.conjugate(mm)),mm)
        otok[ti] = 0.5*(np.abs(np.matrix.trace(mm))/dim)
    otok = np.array(otok)
    return otok

# OTOC 1-Re(F) SLOW (IF HIGHER PRECISION REQUIRED) - USE ONLY IF V,W ARE BOTH UNITARY & HERMITIAN
def OTOCF_infty(V,W,ener,basis,N,dt,t0,ortho):
    dim = len(ener)
    if ortho == True:
        basist = np.transpose(basis)
    else:
        basist = np.linalg.inv(basis)
    S0 = np.matmul(V,basis)
    S0 = np.matmul(basist,S0)
    S = np.matmul(W,basis)
    S = np.matmul(basist,S)
    mm = np.zeros((dim,dim),dtype=complex)
    mm1 = np.zeros((dim,dim),dtype=complex)
    otok = np.zeros(N,dtype = float)
    U = np.zeros((dim,dim),dtype=complex) # U evolucion temporal
    Udagger = np.zeros((dim,dim),dtype=complex) # U*
    for ti in range(0,N):
        for c1 in range(0,dim):
            U[c1][c1] = np.exp(-1j*tiempo[ti]*ener[c1])
            Udagger[c1][c1] = np.exp(1j*tiempo[ti]*ener[c1])
        mm = np.matmul(S0,U)
        mm1= np.matmul(Udagger,mm)
        mm = np.matmul(mm1,S)
        mm = np.matmul(mm,mm)
        otok[ti] = 1 - (np.matrix.trace(mm)/len(ener)).real 
    otok = np.array(otok)
    return otok

