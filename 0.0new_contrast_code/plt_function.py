#if type==si/si02: si/si02/hbn/mose2/hbn
#if type==ag/al2o3: ag/al2o3/hbn/mose2/hbn
import numpy as np
import matplotlib.pyplot as plt
from refrac_index_si import get_n,get_k
from refrac_index_sio2 import get_n_sio2
from refrac_index_Ag import get_n_ag,get_k_ag
from refrac_index_al2o3 import get_n_al2o3
#include the thickness of hbn into the plot function
#wavelength:570 to 620 nm,  oscillator at 603nm 
ev=1.602176634*10**(-19)
h=6.626*10**(-34)
c=3*10**(8)
def lam2ev(lam):#lambda:nm, E:ev
    lam=lam*10**(-9)
    E=h*c/lam
    return E/ev

def delta(ni,di,lam):#lam dependent
    return 2*np.pi*ni*di/lam
def tmatrix(delta,n):#lam dependent
    tm=np.zeros((2,2),dtype=np.complex_)#dtype
    tm[0,0]=np.cos(delta)
    tm[0,1]=np.sin(delta)*1j/n
    tm[1,0]=np.sin(delta)*1j*n
    tm[1,1]=np.cos(delta)
    return tm
#print(tmatrix(np.pi/2,3+2*1j))#test
n1=1.85#HBN
n3=1.85#HBN
def chi(lam,E0,f,gamma):
    E=lam2ev(lam)
    return f/(E0**2-E**2-1j*E*gamma)
def epsilon(E,E0,f,gamma):
    chi=f/(E0**2-E**2-1j*E*gamma)#E0:resonance energy?
    return 1+chi
#print(epsilon(2.0,2.2,3.5,2.4))#test

def n2(epsilon):#mose2
    epabs=abs(epsilon)
    epreal=epsilon.real
    n=(((epabs+epreal))/2)**(1/2)
    k=(((epabs-epreal))/2)**(1/2)
    return n+1j*k 


#print(vec1(550))#test




# d1=10:up hbn
# d2=0.65:mose2
# d3=10:down hbn
# d4=90:sio2
def plotfunc(type,wavelengtharr,E0,f,gamma,d1,d2,d3,d4):###
    # if type=="sio2/si":
    def n4(lam):#sio2
        return get_n_sio2(lam)
    def n5(lam):#si
        return get_n(lam)+1j*get_k(lam)
    def n4_(lam):#al2o3
        return get_n_al2o3(lam)
    def n5_(lam):#ag
        return get_n_ag(lam)+1j*get_k_ag(lam)
    def vec1(lam):
        v1=np.zeros((2,1),dtype=np.complex_)#dtype
        v1[0,0]=1
        if type=="sio2/si":
            v1[1,0]=n5(lam)
        elif type=="ag/al2o3":
            v1[1,0]=n5_(lam)
        else:
            print("no such type")
        return v1


    def tmatrix4(lam):
        return tmatrix(delta(n4(lam),d4,lam),n4(lam))
    def tmatrix4_(lam):
        return tmatrix(delta(n4_(lam),d4,lam),n4_(lam))
    #print(tmatrix4(550))#test
    
    def tmatrix3(lam):
        return tmatrix(delta(n3,d3,lam),n3)
    
    def tmatrix2(lam,E0,f,gamma):#########################
        E=lam2ev(lam)
        return tmatrix(delta(n2(lam),d2,lam),n2(epsilon(E,E0,f,gamma)))
    # def tmatrix2s(lam):#########################
    #     return tmatrix(delta(1,d2*10**(-10),lam),1)---isn't used anymore
    def tmatrix1(lam):
        return tmatrix(delta(n1,d1,lam),n1)
    def contrast(x,xs):
        return (x-xs)/xs

    size=np.size(wavelengtharr)
    Barr=np.zeros(size,dtype=np.complex_)#whole
    Carr=np.zeros(size,dtype=np.complex_)

    Barrs=np.zeros(size,dtype=np.complex_)#substrate
    Carrs=np.zeros(size,dtype=np.complex_)

    chiarr_real=np.zeros(size,dtype=np.complex_)
    chiarr_imag=np.zeros(size,dtype=np.complex_)

    for i in range(np.size(wavelengtharr)):
        lam=wavelengtharr[i]
        
        vecBC=np.matmul(tmatrix1(lam),vec1(lam))
        vecBC=np.matmul(tmatrix2(lam,E0,f,gamma),vecBC)
        vecBC=np.matmul(tmatrix3(lam),vecBC)
        if type=="sio2/si":
            vecBC=np.matmul(tmatrix4(lam),vecBC)
        elif type=="ag/al2o3":
            vecBC=np.matmul(tmatrix4_(lam),vecBC)

        Barr[i]=(vecBC[0,0])
        Carr[i]=(vecBC[1,0])
        
        #substrate
        vecBCs=np.matmul(tmatrix1(lam),vec1(lam))
        # vecBCs=np.matmul(tmatrix2s(lam),vecBCs)
        vecBCs=np.matmul(tmatrix3(lam),vecBCs)
        if type=="sio2/si":
            vecBC=np.matmul(tmatrix4(lam),vecBC)
        elif type=="ag/al2o3":
            vecBC=np.matmul(tmatrix4_(lam),vecBC)
        
        Barrs[i]=(vecBCs[0,0])
        Carrs[i]=(vecBCs[1,0])
        chiarr_real[i]=(chi(lam,E0,f,gamma).real)
        chiarr_imag[i]=(chi(lam,E0,f,gamma).imag)
    
    #print(vecBC)#test
    Yarr=np.zeros_like(Barr,dtype=np.complex_)
    Yarrs=np.zeros_like(Barrs,dtype=np.complex_)
    for i in range(np.size(Yarr)):
        Yarr[i]=Carr[i]/Barr[i]
        Yarrs[i]=Carrs[i]/Barrs[i]
    rarr=np.zeros_like(Yarr,dtype=np.complex_)
    rarrs=np.zeros_like(Yarrs,dtype=np.complex_)
    Rarr=np.zeros_like(Yarr,dtype=np.complex_)#reflected intensity
    Rarrs=np.zeros_like(Yarrs,dtype=np.complex_)#reflected intensity
    contrastarr=np.zeros_like(Yarrs,dtype=np.complex_)

    for i in range(np.size(Yarr)):
        rarr[i]=(1-Yarr[i])/(1+Yarr[i])
        Rarr[i]=abs(rarr[i])**2
        rarrs[i]=(1-Yarrs[i])/(1+Yarrs[i])
        Rarrs[i]=abs(rarrs[i])**2
        contrastarr[i]=contrast(Rarr[i],Rarrs[i])
#imag part is 0
#turn it into a real array
    x=[]
    y=[]
    ctr=[]
    z=[]
    a=[]
    rs_prime=[]
    for i in range(np.size(Rarr)):
        x.append(Rarr[i].real)
        y.append(Rarrs[i].real)
        ctr.append(contrastarr[i].real)
        z.append(chiarr_real[i].real)
        a.append(chiarr_imag[i].real)
        rs_prime.append(rarrs[i].real)#

# x:Rarr
# y:Rarrs
# ctr:contrastarr
# z:chiarr_real
# a:chiarr_imag
# rs_prime:frenel coefficient of the substrate
    return [wavelengtharr,x,y,ctr,z,a,rs_prime]#

