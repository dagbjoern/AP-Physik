import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const



#Berechnung
def sigma():
    ep=1.295
    re=2.82e-15
    eins=(1+ep)/(ep**2)
    zwei=(2*(1-ep)/(1+2*ep))
    drei=(1/ep)*np.log(1+2*ep)
    vier=(1/(2*ep))*np.log(1+2*ep)
    fünf=(1+3*ep)/((1+2*ep)**2)
    return(2*np.pi*(re**2)*(eins*(zwei-drei)+vier-fünf))



def n(z,V_M):
    N_L=2.6867774e25
    return((z*N_L)/(V_M))
#n für Kupfer
z_cu=29
V_M_cu=7.11e-06
print('\nKufer z=',z_cu,'\n Molares Volumen=',7.11e-06)
print('n_cu=',n(z_cu,V_M_cu))
#n für Blei
z_pb=82
V_M_pb=18.26*10e-6
print('\n Blei z=',z_pb,'\n Molares Volumen=',7.11e-06)
print('n_pb=',n(z_pb,V_M_pb))



sig=sigma()
print('\nsigma',sigma())
print('\nKufper theorie\n mu=',n(z_cu,V_M_cu)*sig)

print('\nBlei theorie\n mu=',n(z_pb,V_M_pb)*sig)




def N_D(D,N_0,mu):
    return(N_0*np.exp(-mu*D))


T_cu ,n_cu ,d_cu =np.genfromtxt('messdatenGcu.txt',unpack=True)
d_cu=d_cu*10**(-2)
n_cu=unp.uarray(n_cu,np.sqrt(n_cu))
print(n_cu)


N_cu=n_cu/T_cu
#print(N_cu)

params_cu, covariance_cu =curve_fit(N_D,d_cu,noms(N_cu))
N_0_cu=params_cu[0]
mu_cu=params_cu[1]
print('\n N_0_cu=',N_0_cu)
print('mu_cu=',mu_cu)


#*params
x=np.linspace(0,0.70)
plt.figure(1)
plt.errorbar(d_cu,noms(N_cu) ,xerr=0,yerr=stds(N_cu), fmt='cx')
plt.plot(d_cu,noms(N_cu),'rx',label=r'$Messwerte$')
plt.plot(x,N_D(x,*params_cu),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Dicke \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Zählrate \ N \ in \ 1/s}$')
plt.savefig('a)cu.pdf')


T_pb ,n_pb ,d_pb =np.genfromtxt('messdatenGpb.txt',unpack=True)
d_pb=d_pb*10**(-2)
n_pb=unp.uarray(n_pb,np.sqrt(n_pb))
print(n_pb)
N_pb=n_pb/T_pb
#print(N_pb)


params_pb, covariance_pb =curve_fit(N_D,d_pb,noms(N_pb))

N_0_pb=params_pb[0]
mu_pb=params_pb[1]

print('N_0_pb=',N_0_pb)
print('mu_pb=',mu_pb)

plt.figure(2)
plt.errorbar(d_pb,noms(N_pb) ,xerr=0,yerr=stds(N_pb), fmt='cx')
plt.plot(d_pb,noms(N_pb),'rx',label=r'$Messwerte$')
plt.plot(x,N_D(x,*params_pb),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Dicke \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Zählrate \ N \ in \ 1/s}$')
plt.savefig('a)pb.pdf')

T_b, n_b ,d_b =np.genfromtxt('messdatenB.txt',unpack=True)
n_b=unp.uarray(n_b,np.sqrt(n_b))
N_b=n_b/T_b
#print(N_b)

#plt.figure(1)
#plt.plot(v,U_A,'rx',label=r'$Messwerte$')
#plt.plot(x,gaus2(x,1.14,35300,290),'k-',label=r'$Messwerte$')
#plt.plot(x,gaus2(x,*params),'b-',label=r'$Ausgleichsfunktion$')
#plt.legend(loc='best')
#plt.xlabel(r'$Frequenz \ \nu /hz$')
#plt.ylabel(r'$frac{U_A}{U_E} $')
#plt.savefig('a).pdf')
