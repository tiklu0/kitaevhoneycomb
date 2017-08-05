import numpy as np
import matplotlib.pyplot as plt
from scipy import interpolate
import math as math




def hamiltonian(n,alpha,vz,mu,t,delta):
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
            k[4*(i-1)+1][4*i]=-alpha
            k[4*(i-1)][4*i+1]=alpha
            k[4*(i-1)+3][4*i+2]=alpha
            k[4*(i-1)+2][4*i+3]=-alpha
    h=h+k.conj().T+k #conjugate transpose is needed for the hopping strength and rashba part
    return h

H=hamiltonian(60,-0.7,2,2,1,1)

w=np.linalg.eigvalsh(H)
w=w.real
hist, bins = np.histogram(w, bins=400)
width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2
plt.figure(1)
plt.bar(center, hist, align='center', width=width)
plt.title("$\mu=V_{z}=2,\Delta=t=1,Rashba=0.2$")
plt.xlabel("Energy")
plt.ylabel("DOS")
plt.show()





#H=hamiltonian(n=2,alpha=0.2,vz=2,mu=2,t=1,delta=1)
#print(H[][])





#the nanowire greens function


def g0rf(omega):  #o is frequency, NOT label
    d=0.01
    w=(omega+complex(0,d))*np.identity(4*n*N,dtype=complex)
    g0=np.linalg.inv(w-H)
    return g0

def g0af(omega):  #o is frequency, NOT label
    d=0.01
    w=(omega+complex(0,d))*np.identity(4*n*N,dtype=complex)
    g0=np.linalg.inv(w-H)
    return g0.conj().T













"""
gamma=0.05   #the tunnelling parameter
n=60 #number of sites



x=np.linspace(1,60,60)
V=np.linspace(-0.005,0.005,6)
omega=np.arange(-0.05,0.05,0.001)  #@ccuracy greater than 1/80
bane=len(omega)
domega=0.001



g0list=[]       
i=0       
for i in range(bane):
    g0list.append(g0f(omega[i]))
    
    
#LDOS
def LDOS(n,alpha,vz,mu,t,delta):
    g0=g0f(0)
    l=np.zeros(n,dtype=np.complex_);
    site=np.ones(n);
    for i in range(n):
        l[i]=-((g0[4*i][4*i]+g0[4*i+1][4*i+1]+g0[4*i+2][4*i+2]+g0[4*i+3][4*i+3]))/np.pi
        site[i]=i+1
    l=l.imag
    print( l)
    xnew=np.linspace(site[0],site[n-1],num=500,endpoint=True)
    spl = interpolate.spline(site,l,xnew,3,'smoothest')
    plt.plot(xnew,spl)
    plt.title("($\mu=V_{z}=2,\Delta=t=1,Rashba =0.2$)")
    plt.xlabel("Site")
    plt.ylabel("LDOS")
    plt.show()
    
LDOS(60,0.2,2,2,1,1)

#normal lead   x is lattice site

#g0=g0r=g0a



def sigmarf(x,gamma,n):   #no w dependency here
    sigmar=np.zeros((4*n,4*n),dtype=np.complex_)
    sigmar[4*x][4*x]=complex(0,gamma)
    sigmar[4*x+1][4*x+1]=complex(0,gamma)
    sigmar[4*x+2][4*x+2]=complex(0,gamma)
    sigmar[4*x+3][4*x+3]=complex(0,gamma)
    return sigmar

def sigmaAf(x,gamma,n):   #no w dependency here
    sigmaA=np.zeros((4*n,4*n),dtype=np.complex_)
    sigmaA[4*x][4*x]=complex(0,-gamma)
    sigmaA[4*x+1][4*x+1]=complex(0,-gamma)
    sigmaA[4*x+2][4*x+2]=complex(0,-gamma)
    sigmaA[4*x+3][4*x+3]=complex(0,-gamma)
    return sigmaA



#STORING THE GAMMA_R/A Matrices which do not depend in omega

i=0
sigmarlist=[]
sigmaAlist=[]

for i in range(n):
    sigmarlist.append(sigmarf(i,gamma,n))
    sigmaAlist.append(sigmaAf(i,gamma,n))



def sigmakf(o,j,beta,x,n):
       sigmak=np.zeros((4*n,4*n),dtype=np.complex_)
       sigmak[4*x][4*x]=gamma*(-2J)*math.tanh(beta*(omega[o]-V[j])/2)
       sigmak[4*x+1][4*x+1]=gamma*(-2J)*math.tanh(beta*(omega[o]+V[j])/2)
       sigmak[4*x+2][4*x+2]=gamma*(-2J)*math.tanh(beta*(omega[o]-V[j])/2)
       sigmak[4*x+3][4*x+3]=gamma*(-2J)*math.tanh(beta*(omega[o]+V[j])/2)
       return sigmak


def grf(x,o):
    g0=g0list[o]
    gr=np.linalg.inv(np.linalg.inv(g0)-sigmarlist[x])
    return gr

def gAf(x,o):
    g0=g0list[o]
    gA=np.linalg.inv(np.linalg.inv(g0)-sigmaAlist[x])
    return gA


def gkf(i,o,beta,n,j,sigmak,grf):
    im=np.dot( grf, sigmak)
    return np.dot(im,gAf(i,o))



def tzf(x):
    tz=np.zeros((4*n,4*n))
    tz[4*x][4*x]=1
    tz[4*x+1][4*x+1]=-1
    tz[4*x+2][4*x+2]=1
    tz[4*x+3][4*x+3]=-1
    return tz


conductance=np.zeros(60)
gamma=0.05
i=0
tz=[]
for i in range(60):
    tz.append(tzf(i))
i=0
j=0
o=0

for i in range(60):
    I=np.zeros(6)
    j=0
    for j in range(6):
        o=0
        m=np.zeros((4*n,4*n))
        for o in range(bane):
            sigmak=sigmakf(o,j,200,i,60)
            gr=grf(i,o)
            m=m+domega*(np.dot(gr,sigmak)+np.dot(gkf(i,o,200,60,j,sigmak,gr),sigmaAlist[i]))
        m=np.dot(tz[i],m)
        I[j]=0.5*np.trace(m)

    tck = interpolate.splrep(V,I,s=0)
    yder = interpolate.splev(0, tck, der=1)
    conductance[i]=yder


xnew=np.linspace(x[0],x[n-1],num=500,endpoint=True)
spl = interpolate.spline(x,--conductance,xnew,3,'smoothest')
plt.plot(xnew,spl,'r')
plt.title("($\mu=V_{z}=2,\Delta=t=1,Rashba =0.2,\Gamma=0.05$)")
plt.xlabel("Site")
plt.ylabel("Differential $Conductance$ in units of $\frac{e}{h}$")
plt.show()

"""
##################  Differential Conductance------Superconducting function definitions         #####################################

def xrr(w):
    if(np.abs(w)<delta):
        u=-w/(np.sqrt(-w*w+delta*delta))
    else:
        u=complex(0,1)*np.abs(w)/(np.sqrt(w*w-delta*delta))
    return u

def xra(w):
    if(np.abs(w)<delta):
        u=-w/(np.sqrt(-w*w+delta*delta))
    else:
        u=-complex(0,1)*np.abs(w)/(np.sqrt(w*w-delta*delta))
    return u

def xrk(w,beta):
    if(np.abs(w)>delta):
        u=-2*complex(0,1)*np.abs(w)*math.tanh(beta*w/2)/(np.sqrt(w*w-delta*delta))
    else:
        u=0
    return u



i=0
def Ssigmarf(o,v,x,gamma,n,N):   #no w dependency here
    sigmar=np.zeros((4*n*N,4*n*N),dtype=np.complex_)            #N-----no. of chains
                                                                 #n---- no. of sites
    for i in range(N):
        sigmar[4*n*i+4*x][4*n*i+4*x]=gamma*xrr(omega[o]-v)
        sigmar[4*n*i+4*x+1][4*n*i+4*x+1]=gamma*xrr(omega[o]+v)
        sigmar[4*n*i+4*x+2][4*n*i+4*x+2]=gamma*xrr(omega[o]-v)
        sigmar[4*n*i+4*x+3][4*n*i+4*x+3]=gamma*xrr(omega[o]+v)
    return sigmar
i=0

def SsigmaAf(o,v,x,gamma,n,N):   #no w dependency here
    sigmaA=np.zeros((4*n*N,4*n*N),dtype=np.complex_)
    for i in range(N):
        sigmaA[4*n*i+4*x][4*n*i+4*x]=gamma*xra(omega[o]-v)
        sigmaA[4*n*i+4*x+1][4*n*i+4*x+1]=gamma*xra(omega[o]+v)
        sigmaA[4*n*i+4*x+2][4*n*i+4*x+2]=gamma*xra(omega[o]-v)
        sigmaA[4*n*i+4*x+3][4*n*i+4*x+3]=gamma*xra(omega[o]+v)
    return sigmaA
i=0

def Ssigmakf(o,v,x,n,N,beta):   #no w dependency here
    sigmak=np.zeros((4*n*N,4*n*N),dtype=np.complex_)
    for i in range(N):
        sigmak[4*n*i+4*x][4*n*i+4*x]=gamma*xrk(omega[o]-v,beta)
        sigmak[4*n*i+4*x+1][4*n*i+4*x+1]=gamma*xrk(omega[o]+v,beta)
        sigmak[4*n*i+4*x+2][4*n*i+4*x+2]=gamma*xrk(omega[o]-v,beta)
        sigmak[4*n*i+4*x+3][4*n*i+4*x+3]=gamma*xrk(omega[o]+v,beta)
    return sigmak


def grf(x,o,v):
    g0=g0rf(omega[o])
    gr=np.linalg.inv(np.linalg.inv(g0)-Ssigmarf(o,v,x,gamma,n,N))
    return gr
    
def gAf(x,o,v):
    g0=g0af(omega[o])
    
    gA=np.linalg.inv(np.linalg.inv(g0)-SsigmaAf(o,v,x,gamma,n,N))
    return gA
    
    
def gkf(i,o,Ssigmak,gr,v):
    Ssigmak=sparse.csr_matrix(Ssigmak)
    gr=sparse.csr_matrix(gr)
    im=gr*Ssigmak
    io=sparse.csr_matrix(gAf(i,o,v))
    ret=im*io
    return ret.toarray()
    
i=0
    
def tzf(x):
    n=60
    tz=0*sparse.eye(4*n*N,dtype=complex,format="csr")
    for i in range(N):
        tz[4*n*i+4*x,4*n*i+4*x]=1
        tz[4*n*i+4*x+1,4*n*i+4*x+1]=1
        tz[4*n*i+4*x+2,4*n*i+4*x+2]=-1
        tz[4*n*i+4*x+3,4*n*i+4*x+3]=-1
    return tz
    



##################  Differential Conductance------superconducting vs temperature        #####################################    


i=0
i=0
j=0
o=0
delta=1
gamma=0.05
vec=[0.1, 1.25, 2.0, 3.2, 3.6, 6, 11]
conductance=np.zeros(n)

for l in range(1):
    #V=np.linspace(-0.05,0.05,5)

    gamma=0.05  
    n=60
    N=1
   
    H=hamiltonian(4,60,0.5,2,2,1,1,0.1)
     
    #conductance=np.zeros(4)
    #temp=np.linspace(10,300,15)
    x=np.zeros(n)    
    for t in range(1):   #temperature
        for i in range(n):   #site position
            I=np.zeros(2)
            j=0
            x[i]=i
            for j in range(2):
                o=0
                m=0*sparse.eye(4*n*N,dtype=complex,format="csr")
                
                omega=np.arange(-delta-0.005+0.001,-1,0.00011/(l+5))  #@ccuracy greater than 1/80
                bane=len(omega)
                domega=0.00011/(l+5)
                for o in range(bane):
                        
                    Ssigmak=Ssigmakf(o,-0.005+j*(0.01),i,60,4,200)
                    gr=grf(i,o,-0.005+j*(0.01))
                    gk=gkf(i,o,Ssigmak,gr,-0.005+j*(0.01))
                    SsigmaA=SsigmaAf(o,-0.005+j*(0.01),i,gamma,n,N)
                        
                    gr=sparse.csr_matrix(gr)
                    Ssigmak=sparse.csr_matrix(Ssigmak)
                    gk=sparse.csr_matrix(gk)
                    SsigmaA=sparse.csr_matrix(SsigmaA)
                    m=m+domega*((gr*Ssigmak)+(gk*SsigmaA))
                    print(i)
                
                omega=np.arange(1,delta+0.005-0.001,0.00011/(l+5))  #@ccuracy greater than 1/80
                bane=len(omega)
                domega=0.00011/(l+5)
                for o in range(bane):
                        
                    Ssigmak=Ssigmakf(o,-0.005+j*(0.01),i,60,4,200)
                    gr=grf(i,o,-0.005+j*(0.01))
                    gk=gkf(i,o,Ssigmak,gr,-0.005+j*(0.01))
                    SsigmaA=SsigmaAf(o,-0.005+j*(0.01),i,gamma,n,N)
                        
                    gr=sparse.csr_matrix(gr)
                    Ssigmak=sparse.csr_matrix(Ssigmak)
                    gk=sparse.csr_matrix(gk)
                    SsigmaA=sparse.csr_matrix(SsigmaA)
                    m=m+domega*((gr*Ssigmak)+(gk*SsigmaA))
                    print(i) 
                
                
                
                omega=np.arange(+delta+0.005+0.001,10,0.15/(l+5))  #@ccuracy greater than 1/80
                bane=len(omega)
                domega=0.15/(l+5)
                for o in range(bane):
                        
                    Ssigmak=Ssigmakf(o,-0.005+j*(0.01),i,60,4,200)
                    gr=grf(i,o,-0.005+j*(0.01))
                    gk=gkf(i,o,Ssigmak,gr,-0.005+j*(0.01))
                    SsigmaA=SsigmaAf(o,-0.005+j*(0.01),i,gamma,n,N)
                        
                    gr=sparse.csr_matrix(gr)
                    Ssigmak=sparse.csr_matrix(Ssigmak)
                    gk=sparse.csr_matrix(gk)
                    SsigmaA=sparse.csr_matrix(SsigmaA)
                    m=m+domega*((gr*Ssigmak)+(gk*SsigmaA))
                    print(i)
                        #print(t)
                         
                omega=np.arange(-10,-delta-0.005-0.001,0.15/(l+5))  #@ccuracy greater than 1/80
                bane=len(omega)
                domega=0.15/(l+5)
                for o in range(bane):
                        
                    Ssigmak=Ssigmakf(o,-0.005+j*(0.01),i,60,4,200)
                    gr=grf(i,o,-0.005+j*(0.01))
                    gk=gkf(i,o,Ssigmak,gr,-0.005+j*(0.01))
                    SsigmaA=SsigmaAf(o,-0.005+j*(0.01),i,gamma,n,N)
                        
                    gr=sparse.csr_matrix(gr)
                    sigmak=sparse.csr_matrix(Ssigmak)
                    gk=sparse.csr_matrix(gk)
                    SsigmaA=sparse.csr_matrix(SsigmaA)
                    m=m+domega*((gr*Ssigmak)+(gk*SsigmaA))
                    print(i)
               
                m=tzf(i)*m
                I[j]=0.5*m.diagonal().sum()
                
                
            
            
            conductance[i]=(I[1]-I[0])/(0.02*(4-np.pi))            #tck = interpolate.splrep(V,I,s=0)
            #yder = interpolate.splev(0, tck, der=1)
            #conductance[l]=yder/(2*(4-np.pi))
    
    xnew=np.linspace(x[0],x[n-1],num=900,endpoint=True)
    spl = interpolate.spline(x,-conductance,xnew,3,'smoothest')
    plt.plot(xnew,spl,'r. ')
    plt.title("($\\mu=V_{z}=2,\\Delta=t_{x}=1,t_{y}=0.1, \\alpha=0.2,\\Gamma=0.05,\\beta=200$")
    plt.xlabel("Site")
    plt.ylabel("Differential Conductance in units of $\\frac{2e^2(4-\\pi)}{h}$")             
    #plt.plot(temp,conductance,'r. ')
    #plt.title("($\mu=V_{z}=2,\Delta=t=1,Rashba =0.2,\Gamma=0.05$)")
    
"""
print(conductance)
H=hamiltonian(4,60,0.5,2,2,1,1,0.1)
n=60
N=4     
    #conductance=np.zeros(15)
temp=np.linspace(10,300,15)

omega=np.arange(0,100,1)
print(omega)
for o in range(100):
    Ssigmak=Ssigmakf(o,0.005,i,60,4,temp[12])
    gr=grf(i,o,0.005)
    gk=gkf(i,o,Ssigmak,gr,0.005)
    SsigmaA=SsigmaAf(o,0.005,i,gamma,n,N)
                    
    gr=sparse.csr_matrix(gr)
    Ssigmak=sparse.csr_matrix(Ssigmak)
    gk=sparse.csr_matrix(gk)
    SsigmaA=sparse.csr_matrix(SsigmaA)
    m=((gr*Ssigmak)+(gk*SsigmaA))
    fr=m.toarray()
    print(np.abs(fr).max())
























