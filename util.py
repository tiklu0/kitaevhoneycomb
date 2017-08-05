# -*- coding: utf-8 -*-
"""
Created on Fri May  5 00:47:40 2017

@author: user
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import math as math




def h1(n,alpha,vz,mu,t,delta):
    i=0
    h=np.zeros((4*n,4*n))
    k=np.zeros((4*n,4*n))
    #first the diagonal submatrices are initialized
    for i in range(n):
        h[4*i][4*i]=-mu
        h[4*i+2][4*i]=delta     #column then row, NOT [row][column]
        h[4*i+3][4*i+1]=delta
        h[4*i+1][4*i+1]=-mu
        h[4*i+2][4*i+2]=mu   #everything else besides rashba and hopping
        h[4*i][4*i+2]=delta
        h[4*i+3][4*i+3]=mu
        h[4*i+1][4*i+3]=delta
        h[4*i][4*i+1]=vz
        h[4*i+1][4*i]=vz
        h[4*i+2][4*i+3]=vz
        h[4*i+3][4*i+2]=vz
    #now for the off diagonal sub-matrices
    for i in range(n):
        if(i>0):
            k[4*(i-1)][4*i]=-t
            k[4*(i-1)+1][4*i+1]=-t
            k[4*(i-1)+2][4*i+2]=t
            k[4*(i-1)+3][4*i+3]=t            #rashba and hopping
            k[4*(i-1)+1][4*i]=alpha
            k[4*(i-1)][4*i+1]=-alpha
            k[4*(i-1)+3][4*i+2]=-alpha
            k[4*(i-1)+2][4*i+3]=alpha
    h=h+k.conj().T+k #conjugate transpose is needed for the hopping strength and rashba part
    return h

H1=h1(60,0.7,2,2,1,1)


def h2(N,n,alpha,vz,mu,tx,delta,ty):
    i=0
    h=np.zeros((4*n*N,4*n*N))
    k=np.zeros((4*n*N,4*n*N))
    #first the diagonal submatrices are initialized
    for i in range(n*N):
        h[4*i][4*i]=-mu
        h[4*i+2][4*i]=delta     #column then row, NOT [row][column]
        h[4*i+3][4*i+1]=delta
        h[4*i+1][4*i+1]=-mu
        h[4*i+2][4*i+2]=mu   #everything else besides rashba and hopping
        h[4*i][4*i+2]=delta
        h[4*i+3][4*i+3]=mu
        h[4*i+1][4*i+3]=delta
        h[4*i][4*i+1]=vz
        h[4*i+1][4*i]=vz
        h[4*i+2][4*i+3]=vz
        h[4*i+3][4*i+2]=vz
    #now for the off diagonal sub-matrices
    for m in range(N): 
        for i in range(n):
            if(i>0):
                k[m*4*n+4*(i-1)][m*4*n+4*i]=-tx
                k[m*4*n+4*(i-1)+1][m*4*n+4*i+1]=-tx
                k[m*4*n+4*(i-1)+2][m*4*n+4*i+2]=tx
                k[m*4*n+4*(i-1)+3][m*4*n+4*i+3]=tx            #rashba and hopping
                k[m*4*n+4*(i-1)+1][m*4*n+4*i]=-alpha
                k[m*4*n+4*(i-1)][m*4*n+4*i+1]=alpha
                k[m*4*n+4*(i-1)+3][m*4*n+4*i+2]=alpha
                k[m*4*n+4*(i-1)+2][m*4*n+4*i+3]=-alpha
            if(m>0):
                k[(m-1)*4*n+4*i][(m)*4*n+4*i]=-ty
                k[(m-1)*4*n+4*i+1][(m)*4*n+4*i+1]=-ty
                k[(m-1)*4*n+4*i+2][(m)*4*n+4*i+2]=ty
                k[(m-1)*4*n+4*i+3][(m)*4*n+4*i+3]=ty
    h=h+k.conj().T+k #conjugate transpose is needed for the hopping strength and rashba part
    return h
 

H2=h2(2,30,0.7,2,2,1,1,0.1)

H=H2-H1
print(H[120][0])

