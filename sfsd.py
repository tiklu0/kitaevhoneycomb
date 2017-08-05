# -*- coding: utf-8 -*-
"""
Created on Thu Jul 27 23:07:57 2017

@author: user
"""

import numpy as np

import scipy.sparse as sparse
import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy import interpolate
import math as math
from mpl_toolkits.mplot3d import Axes3D as ax
from pylab import savefig
import matplotlib.patches as mpatches
"""
N---- number of chains
n ---- number of sites in each chain
"""
print('fml')

n=60
N=1

#np.set_printoptions(threshold=np.inf)
#


def g0rf(omega):  #o is frequency, NOT label
    d=0.0001
    w=(omega+complex(0,d))*np.identity(4*n*N,dtype=complex)
    g0=np.linalg.inv(w-H)
    return g0

def g0af(omega):  #o is frequency, NOT label
    d=0.0001
    w=(omega+complex(0,d))*np.identity(4*n*N,dtype=complex)
    g0=np.linalg.inv(w-H)
    return g0.conjugate().T


def func(x,y,g0):
    h=-(g0[4*n*y+4*x][4*n*y+4*x]+g0[4*n*y+4*x+1][4*n*y+4*x+1]+g0[4*n*y+4*x+2][4*n*y+4*x+2]+g0[4*n*y+4*x+3][4*n*y+4*x+3])/np.pi
    return np.imag(h)

def hamiltonian(N,n,alpha,vz,mu,tx,delta,ty):
    i=0
    h=np.zeros((4*n*N,4*n*N))
    k=np.zeros((4*n*N,4*n*N))
    #first the diagonal submatrices are initialized
    for i in range(n*N):
        h[4*i][4*i]=-mu/2
        h[4*i+2][4*i]=delta/2     #column then row, NOT [row][column]
        h[4*i+3][4*i+1]=delta/2
        h[4*i+1][4*i+1]=-mu/2
        h[4*i+2][4*i+2]=mu/2  #everything else besides rashba and hopping
        h[4*i][4*i+2]=delta/2
        h[4*i+3][4*i+3]=mu/2
        h[4*i+1][4*i+3]=delta/2
        h[4*i][4*i+1]=vz/2
        h[4*i+1][4*i]=vz/2
        h[4*i+2][4*i+3]=vz/2
        h[4*i+3][4*i+2]=vz/2
    #now for the off diagonal sub-matrices
    for m in range(N): 
        for i in range(n):
            if(i>0):
                k[m*4*n+4*(i-1)][m*4*n+4*i]=-tx/2
                k[m*4*n+4*(i-1)+1][m*4*n+4*i+1]=-tx/2
                k[m*4*n+4*(i-1)+2][m*4*n+4*i+2]=tx/2
                k[m*4*n+4*(i-1)+3][m*4*n+4*i+3]=tx/2           #rashba and hopping
                k[m*4*n+4*(i-1)+1][m*4*n+4*i]=-alpha/2
                k[m*4*n+4*(i-1)][m*4*n+4*i+1]=alpha/2
                k[m*4*n+4*(i-1)+3][m*4*n+4*i+2]=alpha/2
                k[m*4*n+4*(i-1)+2][m*4*n+4*i+3]=-alpha/2
            if(m>0):
                k[(m-1)*4*n+4*i][(m)*4*n+4*i]=-ty/2
                k[(m-1)*4*n+4*i+1][(m)*4*n+4*i+1]=-ty/2
                k[(m-1)*4*n+4*i+2][(m)*4*n+4*i+2]=ty/2
                k[(m-1)*4*n+4*i+3][(m)*4*n+4*i+3]=ty/2
    h=h+k.conj().T+k #conjugate transpose is needed for the hopping strength and rashba part
    return h
#print(hamiltonian(1,2,0.5,2,2,1,1,0.1)) #(N,n,alpha,vz,mu,tx,delta,ty)
    

########################### Normal function definitions ####################################
def xra(w):
    if(np.abs(w)<(delta/10)):
        u=-w/(np.sqrt(-w*w+(delta/10)*(delta/10)))
    else:
        u=-complex(0,1)*np.abs(w)/(np.sqrt(w*w-(delta/10)*(delta/10)))
    return u



def xrr(w):
    if(np.abs(w)<delta/10):
        u=-w/(np.sqrt(-w*w+(delta/10)*(delta/10)))
    else:
        u=complex(0,1)*np.abs(w)/(np.sqrt(w*w-(delta/10)*(delta/10)))
    return u


def xrf(w,beta):
    if(np.abs(w)>(delta/10)):
        u=np.abs(w)/(np.sqrt(w*w-(delta/10)*(delta/10)))
    else:
        u=0
    return u



i=0
def Ssigmarf(o,v,x,gamma,n,N):   #no w dependency here
    sigmar=np.zeros((4*n*N,4*n*N),dtype=np.complex_)            #N-----no. of chains
                                                                 #n---- no. of sites
    for i in range(N):
        for l in range(N):
            sigmar[4*n*i+4*x][4*n*l+4*x]=gamma*xrr(omega[o]-v)
            sigmar[4*n*i+4*x+1][4*n*l+4*x+1]=gamma*xrr(omega[o]-v)
            sigmar[4*n*i+4*x+2][4*n*l+4*x+2]=gamma*xrr(omega[o]+v)
            sigmar[4*n*i+4*x+3][4*n*l+4*x+3]=gamma*xrr(omega[o]+v)
    return sigmar
i=0

def SsigmaAf(o,v,x,gamma,n,N):   #no w dependency here
    sigmaA=np.zeros((4*n*N,4*n*N),dtype=np.complex_)
    for i in range(N):
        for l in range(N):
            sigmaA[4*n*i+4*x][4*n*l+4*x]=gamma*xra(omega[o]-v)
            sigmaA[4*n*i+4*x+1][4*n*l+4*x+1]=gamma*xra(omega[o]-v)
            sigmaA[4*n*i+4*x+2][4*n*l+4*x+2]=gamma*xra(omega[o]+v)
            sigmaA[4*n*i+4*x+3][4*n*l+4*x+3]=gamma*xra(omega[o]+v)
    return sigmaA
i=0
    
    
    
def Ssigmalf(o,v,beta,x,n,N,xr):
    sigmalk=np.zeros((4*n*N,4*n*N),dtype=np.complex_)
    for i in range(N):
            for j in range(N):
                sigmalk[4*i*n+4*x][4*j*n+4*x]=-2*gamma*(1J)*xr*(1/(1+math.exp(beta*(omega[o]-v))))
                sigmalk[4*i*n+4*x+1][4*j*n+4*x+1]=-2*gamma*(1J)*xr*(1/(1+math.exp(beta*(omega[o]-v))))
                sigmalk[4*i*n+4*x+2][4*j*n+4*x+2]=-2*gamma*(1J)*xr*(1/(1+math.exp(beta*(omega[o]+v))))
                sigmalk[4*i*n+4*x+3][4*j*n+4*x+3]=-2*gamma*(1J)*xr*(1/(1+math.exp(beta*(omega[o]+v))))
    return sigmalk

def Ssigmagf(o,v,beta,x,n,N,xr):
    sigmagk=np.zeros((4*n*N,4*n*N),dtype=np.complex_)
    for i in range(N):
            for j in range(N):
                sigmagk[4*i*n+4*x][4*j*n+4*x]=2*gamma*(1J)*xr*(1-1/(1+math.exp(beta*(omega[o]-v))))
                sigmagk[4*i*n+4*x+1][4*j*n+4*x+1]=2*gamma*(1J)*xr*(1-1/(1+math.exp(beta*(omega[o]-v))))
                sigmagk[4*i*n+4*x+2][4*j*n+4*x+2]=2*gamma*(1J)*xr*(1-1/(1+math.exp(beta*(omega[o]+v))))
                sigmagk[4*i*n+4*x+3][4*j*n+4*x+3]=2*gamma*(1J)*xr*(1-1/(1+math.exp(beta*(omega[o]+v))))
    return sigmagk


    
    
def grf(x,o,Ssigmar):
    g0=g0rf(omega[o])
    gr=np.linalg.inv(np.linalg.inv(g0)-Ssigmar)
    return gr
    
def gAf(x,o,sigmaA):
    g0=g0af(omega[o])    
    gA=np.linalg.inv(np.linalg.inv(g0)-sigmaA)
    return gA
    


###################util###########################################333
 
vec=[0.1,1.23,2,6,10]
conductance=np.zeros(5)
for l in range(1):
  
    N=1
    n=60
   
    #yut=np.zeros(bane)
    #/(l+1) 
    #x=np.zeros(n)
    print(l)
    gamma=0.5
    delta=1
    H=hamiltonian(1,60,-0.5,2,2,1,1,0.1) #hamiltonian(N,n,alpha,vz,mu,tx,delta,ty)
    
   
    #temp=np.linspace(10,200,10)
    #conductance=np.zeros(5)     
    for t in range(1):
        
        for i in range(1):   #site position
           
            I=np.zeros(5)
            V=np.zeros(5)
            j=0
            #x[i]=i
            for j in range(1):
                o=0
                
                #m=0*sparse.eye(4*n*N,dtype=complex,format="csr")
                #print(m)
                V[j]=0.1-0.005+(0.01)
                omega1=np.arange(-0.2,-0.1-V[j]-0.00001,0.0001)  #@ccuracy greater than 1/80
                omega2=np.arange(-0.1-V[j]+0.001,-0.005-0.001,0.0001)
                #omega3=np.arange(-0.005+0.001,0.005-0.001,0.001)
                #omega4=np.arange(0.005+0.001,0.1+V[j]-0.001,0.001)
                #omega5=np.arange(0.1+V[j]+0.00001,0.3,0.0001)
                omega=np.append(omega1,omega2)
                #omega=np.append(omega,omega4)
                #omega=np.append(omega,omega4)
                #omega=np.append(omega,omega5)
                bane=len(omega)
                
                kor=np.zeros(4*N)
                domega=0.0001
                for o in range(bane):
                    
                    xr=xrf(omega[o],200)
                    sigmaA=SsigmaAf(o,V[j],i,gamma,n,N)        #SsigmaAf(o,v,x,gamma,n,N)
                    sigmar=sigmaA.conjugate().T
                    sigmal=Ssigmalf(o,V[j],200,i,60,N,xr)       
                    sigmag=Ssigmagf(o,V[j],200,i,60,N,xr)
                    
                    gr=grf(i,o,sigmar)                   #grf(x,o)
                    #gk=gkf(i,o,sigmak,gr)       #gkf(i,o,sigmak,gr)
                    
                    sigmal=sparse.csr_matrix(sigmal)
                    sigmag=sparse.csr_matrix(sigmag)
                    gr=sparse.csr_matrix(gr)
                    gl=(gr*sigmal)*(gr.H)
                    gg=(gr*sigmag)*(gr.H)
                    plo=np.zeros(4*N)
                    ind=0
                    for z in range(4):
                        for u in range(N):
               
                            for f in range(N):
                                plo[ind]=plo[ind]+gg[4*n*f+4*i+z,4*n*u+4*i+z]*sigmal[4*n*u+4*i+z,4*n*f+4*i+z]-gl[4*n*f+4*i+z,4*n*u+4*i+z]*sigmag[4*n*u+4*i+z,4*n*f+4*i+z]                                    
                            kor[ind]=kor[ind]+domega*plo[ind]
                            ind=ind+1
                    #s=plo[0]+plo[1]+plo[2]+plo[3]+plo[4]+plo[5]+plo[6]+plo[7]-plo[8]-plo[9]-plo[10]-plo[11]-plo[12]-plo[13]-plo[14]-plo[15]
                    s=plo[0]+plo[1]-plo[2]-plo[3]
                    plt.plot(omega[o],s,'r. ')        
                    
            """    
                tu=0
                for q in range(N):
                        if(q<2): #N=4
                          I[j]=I[j]+kor[tu]+kor[tu+1]+kor[tu+2]+kor[tu+3]
                        else:
                          I[j]=I[j]-kor[tu]-kor[tu+1]-kor[tu+2]-kor[tu+3]
                        tu=tu+4       
                
            conductance[l]=(I[1]-I[0])/(0.04*(4-np.pi))
            print(conductance[l])
                #plt.plot(omega,yut, 'r. ')

            """
     
#plt.plot(vec,conductance,'r. ')
#plt.savefig('NDCg.png') 
   
#######################################33333333333333333333##############################################
"""  
vec=[0.1, 2.0, 6.0, 3.5,10]
g=[0.02,0.2, 0.5]

for l in range(3):
    
    n=60
    N=4
    conductance=np.zeros(10)
    #x=np.zeros(n)
    print(l)
    gamma=g[l]
    delta=1
    H=hamiltonian(4,60,0.5,2,2,1,1,0.1) #hamiltonian(N,n,alpha,vz,mu,tx,delta,ty)
    
    temp=np.linspace(10,200,10)
    #conductance=np.zeros(5)     
    for t in range(10):
        omega=np.arange(-2,-0.000001,0.01)  #@ccuracy greater than 1/80
        omega2=np.arange(0.000001 ,2,0.01)
        omega=np.append(omega,omega2)
        bane=len(omega)
        domega=0.01#/(l+1) 
        
        for i in range(1):   #site position
            sigmaA=sigmaAf(i,gamma,n,N)
            #sigmaA=sparse.csr_matrix(sigmaA)
            I=np.zeros(5)
            V=np.zeros(5)
            j=0
            #x[i]=i
            for j in range(2):
              
                o=0
                #m=0*sparse.eye(4*n*N,dtype=complex,format="csr")
                #print(m)
                V[j]=-0.005+(j*0.01)
                kor=np.zeros(4*N)
                for o in range(bane):
                    sigmal=sigmalf(o,V[j],temp[t],i,60,4)       #sigmakf(o,v,beta,x,n,N)
                    sigmag=sigmagf(o,V[j],temp[t],i,60,4)
                    gr=grf(i,o)                   #grf(x,o)
                    #gk=gkf(i,o,sigmak,gr)       #gkf(i,o,sigmak,gr)
                    gl=np.dot(gr,sigmal)
                    gl=np.dot(gl,gr.conjugate().T)
                    gg=np.dot(gr,sigmag)
                    gg=np.dot(gg,gr.conjugate().T)
                    plo=np.zeros(4*N)
                    ind=0
                    for z in range(4):
                        for u in range(N):
               
                            for f in range(N):
                                plo[ind]=plo[ind]+gg[4*n*f+4*i+z][4*n*u+4*i+z]*sigmal[4*n*u+4*i+z][4*n*f+4*i+z]-gl[4*n*f+4*i+z][4*n*u+4*i+z]*sigmag[4*n*u+4*i+z][4*n*f+4*i+z]                                    
                            #kor[ind]=kor[ind]+domega*plo[ind]
                            ind=ind+1  
                    #m=m+domega*((gr*sigmak)+(gk*sigmaA))
                            
                    
                #m=tzf(i)*m
                #I[j]=0.5*m.diagonal().sum()
                 
                tu=0
                for q in range(N):
                    if(q<2): #N=4
                      I[j]=I[j]+kor[tu]+kor[tu+1]+kor[tu+2]+kor[tu+3]
                    else:
                      I[j]=I[j]-kor[tu]-kor[tu+1]-kor[tu+2]-kor[tu+3]
                    tu=tu+4
                 
                
            
            #tck = interpolate.splrep(V,I,s=0)
            #yder = interpolate.splev(0, tck, der=1)
            #conductance[t]=yder/2
            conductance[t]=(I[1]-I[0])/0.04
            print(conductance[t])
                           
        #xnew=np.linspace(x[0],x[n-1],num=900,endpoint=True)
        #spl = interpolate.spline(x,-conductance,xnew,3,'smoothest')
        
 
    #plt.plot(vec,-conductance,'r* ',label='$\\Gamma=0.5$')
    
    if(l==0):
      plt.plot(temp,conductance,'r.-',label='$\\Gamma=0.02$')
        
    if(l==1):
        plt.plot(temp,conductance,'b.-',label='$\\Gamma=0.2$')
    if(l==2):
        plt.plot(temp,conductance,'g.-',label='$\\Gamma=0.5$')
 
#plt.plot(te, conductance,'r. ')   
 
plt.title("$\\mu=V_{z}=2,\\Delta=t_{x}=1,\\alpha=0.2,t_{y}=0.1, N=4$")
plt.xlabel("$\\beta$")
plt.ylabel("Differential Conductance in units of $\\frac{2e^2}{h}$")  
plt.grid(True)
plt.legend() 
     
plt.savefig('NDCtempsimple.png', format='png', dpi=1000)
"""
