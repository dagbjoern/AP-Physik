import numpy as np
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)

T1, T2 =np.genfromtxt('c).txt', unpack=True)
c_k , f1 , f2 = np.genfromtxt('b).txt', unpack=True)
l=23.954e-3
c=0.7932e-9
nM=unp.uarray([1,1,1,3,3,4,5,6],[1,1,1,1,1,1,1,1])
C_k=unp.uarray(c_k,c_k*0.2)
C=unp.uarray(c,c*0.05)
L=unp.uarray(l,l*0.05)
print('C_K=',C_k)
print('C',C)
print('L',L)

t1=unp.uarray(T1,5)


t2=unp.uarray(T2,5)
fstart=unp.uarray(19230,5000)
fend=unp.uarray(97660,5000)
csp=0.028e-9
Csp=unp.uarray(csp,csp*0.05)
print('Csp',Csp)
v1=1/(2*np.pi*unp.sqrt(L*(C+Csp)))
v2=1/(2*np.pi*unp.sqrt(L*((1/C+2/C_k)**-1+Csp)))
print(C, v1 , L , Csp)

fc1=fstart+(fend-fstart)*t1*0.001
fc2=fstart+(fend-fstart)*t2*0.001
nM=nM*2
n=(v1+v2)/(2*(v2-v1))


print('höchste relative',(80.7-56.2)/56.2)
print('niedrigsten relative',(39.3-38.5)/38.5)


a=(noms(nM)-noms(n))/(noms(n))

print(fc1)
print(fc2)
print(a)
#print(nM)
