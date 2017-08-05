# -*- coding: utf-8 -*-
"""
Created on Tue May 23 23:03:07 2017

@author: Nagaosa
"""

import scipy.sparse as sparse
import matplotlib.pyplot as plt
from scipy import interpolate
import math as math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D as ax
from pylab import savefig
"""
N---- number of chains
n ---- number of sites in each chain
"""
n=60
N=4
fu=0
def g0rf(omega):  #o is frequency, NOT label
    d=0.01
    w=(omega+complex(0,d))*np.identity(2*n*N,dtype=complex)
    g0=np.linalg.inv(w-H)
    return g0

def g0af(omega):  #o is frequency, NOT label
    d=0.01
    w=(omega+complex(0,d))*np.identity(2*n*N,dtype=complex)
    g0=np.linalg.inv(w-H)
    return g0.conj().T


def hamiltonian(N,n,mu,dx,tx,dy,ty):
    i=0
    h=np.zeros((2*n*N,2*n*N))
    k=np.zeros((2*n*N,2*n*N))
    #first the diagonal submatrices are initialized
    for i in range(n*N):
        h[2*i][2*i]=-mu/2
        h[2*i+1][2*i+1]=mu/2
    #now for the off diagonal sub-matrices
    for m in range(N): 
        for i in range(n):
            if(i>0):
                k[m*2*n+2*(i-1)][m*2*n+2*i]=-tx/2
                k[m*2*n+2*(i-1)+1][m*2*n+2*i+1]=tx/2            
                k[m*2*n+2*(i-1)][m*2*n+2*i+1]=-dx/2
                k[m*2*n+2*(i-1)+1][m*2*n+2*i]=dx/2   
            if(m>0):
                k[(m-1)*2*n+2*i][(m)*2*n+2*i]=-ty/2
                k[(m-1)*2*n+2*i+1][(m)*2*n+2*i+1]=ty/2
                k[(m-1)*2*n+2*i][(m)*2*n+2*i+1]=-dy/2
                k[(m-1)*2*n+2*i+1][(m)*2*n+2*i]=dy/2
                
    h=h+k.conj().T+k #conjugate transpose is needed for the hopping strength and rashba part
    return h




def sigmarf(x,gamma,n,N):   #no w dependency here
    sigmar=np.zeros((2*n*N,2*n*N),dtype=np.complex_)            #N-----no. of chains
    i=0                                                             #n---- no. of sites
    for i in range(N):
        for j in range(N):
            sigmar[2*n*i+2*x][2*n*j+2*x]=gamma*1J
            sigmar[2*n*i+2*x+1][2*n*j+2*x+1]=gamma*1J

    return sigmar

i=0
def sigmaAf(x,gamma,n,N):   #no w dependency here
    sigmaA=np.zeros((2*n*N,2*n*N),dtype=np.complex_)
    for i in range(N):
        for j in range(N):
            sigmaA[2*n*i+2*x][2*n*j+2*x]=-gamma*1J
            sigmaA[2*n*i+2*x+1][2*n*j+2*x+1]=-gamma*1J
    return sigmaA
    
    
    
    
    
def sigmalf(o,v,beta,x,n,N):
    sigmalk=np.zeros((2*n*N,2*n*N),dtype=np.complex_)
    for i in range(N):
            for j in range(N):
                sigmalk[2*i*n+2*x][2*j*n+2*x]=-2*gamma*(1J)*(1/(1+math.exp(beta*(omega[o]-v))))
                sigmalk[2*i*n+2*x+1][2*j*n+2*x+1]=-2*gamma*(1J)*(1/(1+math.exp(beta*(omega[o]+v))))
                
    return sigmalk

def sigmagf(o,v,beta,x,n,N):
    sigmagk=np.zeros((2*n*N,2*n*N),dtype=np.complex_)
    for i in range(N):
            for j in range(N):
                sigmagk[2*i*n+2*x][2*j*n+2*x]=2*gamma*(1J)*(1-1/(1+math.exp(beta*(omega[o]-v))))
                sigmagk[2*i*n+2*x+1][2*j*n+2*x+1]=2*gamma*(1J)*(1-1/(1+math.exp(beta*(omega[o]+v))))
    return sigmagk


    
    
def grf(x,o,sigmar):
    g0=g0rf(omega[o])
    gr=np.linalg.inv(np.linalg.inv(g0)-sigmar)
    return gr
    
def gAf(x,o,sigmaA):
    g0=g0af(omega[o])    
    gA=np.linalg.inv(np.linalg.inv(g0)-sigmaA)
    return gA




################### spectrum ###########################################################
"""
for i in range(600):
    ty=0.01*(i+1)
    H=hamiltonian(N,60,0.7,0.5,1,0.3*ty,ty) #hamiltonian(N,n,mu,dx,tx,dy,ty)
    w=np.linalg.eigvalsh(H)
    w=w.real
    u=ty*np.ones(2*4*60)
    plt.plot(u,w,'r. ',ms=1)
    
plt.title("$N=3,n=60,\\mu=V_{z}=2,\\Delta=t_x=1,\\alpha=0.5$")
plt.xlabel("$t_y$")
plt.ylabel("Energy")
axes = plt.gca()
axes.set_ylim([-1,1])
plt.show()
"""

##############################333    

   

################### Differential Conductance - normal vs temperature ##########################################3
 
vec=[0.1,1.26,1.875,3.3812,5.5]
conductance=np.zeros(60)
for l in range(5):
    omega1=np.arange(-0.03,-0.000001,0.00005)  #@ccuracy greater than 1/80
    omega2=np.arange(0.000001,0.03,0.00005)
    omega=np.append(omega1,omega2)
    N=4
    n=60
    bane=len(omega)
    #yut=np.zeros(bane)
    domega=0.00005#/(l+1) 
    #x=np.zeros(n)
    print(l)
    gamma=0.5
    
    ty=vec[l]
    H=hamiltonian(N,60,0.7,0.5,1,0.3*ty,ty) #hamiltonian(N,n,mu,dx,tx,dy,ty)
    
   
    #temp=np.linspace(10,200,10)
    #conductance=np.zeros(5)     
    for t in range(1):
        
        for i in range(1):   #site position
            sigmaA=sigmaAf(i,gamma,n,N)
            sigmar=sigmaA.conjugate().T
            I=np.zeros(5)
            V=np.zeros(5)
            j=0
            #x[i]=i
            for j in range(2):
                o=0
                #m=0*sparse.eye(4*n*N,dtype=complex,format="csr")
                #print(m)
                V[j]=-0.005+(j*0.01)
                kor=np.zeros(2*N)
                for o in range(bane):
                    sigmal=sigmalf(o,V[j],500 ,i,60,4)       #sigmakf(o,v,beta,x,n,N)
                    sigmag=sigmagf(o,V[j],500 ,i,60,4)
                    
                    gr=grf(i,o,sigmar)                   #grf(x,o)
                    #gk=gkf(i,o,sigmak,gr)       #gkf(i,o,sigmak,gr)
                    
                    sigmal=sparse.csr_matrix(sigmal)
                    sigmag=sparse.csr_matrix(sigmag)
                    gr=sparse.csr_matrix(gr)
                    gl=(gr*sigmal)*(gr.H)
                    gg=(gr*sigmag)*(gr.H)
                    plo=np.zeros(2*N)
                    ind=0
                    for z in range(2):
                        for u in range(N):
               
                            for f in range(N):
                                plo[ind]=plo[ind]+gg[2*n*f+2*(i)+z,2*n*u+2*(i)+z]*sigmal[2*n*u+2*(i)+z,2*n*f+2*(i)+z]-gl[2*n*f+2*(i)+z,2*n*u+2*(i)+z]*sigmag[2*n*u+2*(i)+z,2*n*f+2*(i)+z]                                    
                            kor[ind]=kor[ind]+domega*plo[ind]
                            ind=ind+1
                    
                    #s=plo[0]+plo[1]+plo[2]+plo[3]-plo[4]-plo[5]-plo[6]-plo[7]
                    #s=plo[0]+plo[1]-plo[2]-plo[3]
                    #plt.plot(omega[o],s,'r. ')
                    """  
                    if(l==0):
                        plt.plot(omega[o],s,'r. ')        
                    if(l==1):
                        plt.plot(omega[o],s,'b. ')    
                    if(l==2):
                        plt.plot(omega[o],s,'g. ')    
                    if(l==3):
                        plt.plot(omega[o],s,'c. ')    
                    if(l==4):
                        plt.plot(omega[o],s,'k. ') 
                    """ 
                 #I[j]=kor[0]+kor[1]-kor[2]-kor[3]
                tu=0
                for q in range(N):
                    if(q<2): #N=4
                      I[j]=I[j]+kor[tu]+kor[tu+1]
                    else:
                      I[j]=I[j]-kor[tu]-kor[tu+1]
                      tu=tu+2       
                
            conductance[l]=(I[1]-I[0])/0.02
            print(conductance[l])
                #plt.plot(omega,yut, 'r. ')
            
 
 