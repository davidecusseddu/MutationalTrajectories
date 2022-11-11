# Code for PDE model
# below the code for thomas algorithm and  crout algorithm for solving the tridiagonal system and calculation of the integrals
import numpy as np
def thomas(A,f): # Thomas algorithm solves tridiagonal systems
    M,N = np.shape(A)
    
    beta = np.zeros(M-1)  # L lower diagonal
    c = np.diagonal(A, 1) # U upper diagonal
    alpha = np.zeros(M)   # U diagonal
    
    alpha[0]=A[0,0]
    y = np.zeros(M)
    y[0] = f[0]
    
    for i in range(1,M):
        beta[i-1] = A[i,i-1]/alpha[i-1] 
        alpha[i] = A[i,i] - beta[i-1]*c[i-1]
        y[i] = f[i]-beta[i-1]*y[i-1] # solution of Ly=f
    
    x = np.zeros(M)
    x[M-1] = y[M-1]/alpha[M-1] 
    for j in range(1,M):
        i = (M-1)-j
        x[i] = (y[i] - c[i]*x[i+1])/alpha[i] 
    return x

def crout(A,b): # Crout factorisation algorithm for tridiagonal systems Ref Burden Faires 2010
    [m,n] = np.shape(A)
    L = np.zeros([n,n])
    U = np.zeros([n,n])
    z = np.zeros(n)
    
    # Step 1
    L[0,0] = A[0,0]
    U[0,1] = A[0,1]/L[0,0]
    z[0] = b[0]/L[0,0]
    
    # Step 2
    for i in range(1,n-1):
        L[i,i-1] = A[i,i-1]
        L[i,i] = A[i,i]-L[i,i-1]*U[i-1,i]
        U[i,i+1] = A[i,i+1]/L[i,i]
        z[i] = (b[i] - L[i,i-1]*z[i-1])/L[i,i]

    # Step 3    
    L[n-1,n-2] = A[n-1,n-2]
    L[n-1,n-1] = A[n-1,n-1] - L[n-1,n-2]*U[n-2,n-1]
    z[n-1] = (b[n-1] - L[n-1,n-2]*z[n-2])/L[n-1,n-1]
    
    # Step 4
    x = np.zeros(n)
    x[n-1] = z[n-1]

    # Step 5
    for j in range(2,n+1):
        i = n-j
        x[i] = z[i]-U[i,i+1]*x[i+1]

    # Step 6
    return x
    
    
def integral(y,dx): 
    # It calculates the integral of a function y given as a vector
    # using Trapezoid method
    # it works for uniform discretisation steps dx
    L = len(y)
    integ=0.5*(y[0]+y[L-1])
    for i in range(1,L-1):
        integ += y[i]
    integ = dx*integ
    return integ

def bc_homo_L(b):
    b[0] = 0.0
    return b

def bc_homo_R(b):
    b[-1] = 0.0
    return b
