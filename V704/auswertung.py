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
    re=2.82*10**(-15)
    eins=(1+ep)/(ep**2)
    zwei=(2*(1+ep))/(1+2*ep)
    drei=(1/ep)*np.log(1+2*ep)
    vier=(1/(2*ep))*np.log(1+2*ep)
    fünf=(1+3*ep)/((1+2*ep)**2)
    print('thomschonscher Wirkungsquerschnitt',(8/3)*np.pi*(re**2))
    return(2*np.pi*(re**2)*(eins*(zwei-drei)+vier-fünf))



def n(z,V_M):
    N_L=2.6867774e25
    N_L=const.N_A
    return((z*N_L)/(V_M))
#n für Kupfer
z_cu=29
V_M_cu=7.11e-06
print('\nKufer z=',z_cu,'\n Molares Volumen=',7.11*10**(-6))
print('n_cu=',n(z_cu,V_M_cu))
#n für Blei
z_pb=82
V_M_pb=18.26*1e-6
print('\n Blei z=',z_pb,'\n Molares Volumen=',7.11e-06)
print('n_pb=',n(z_pb,V_M_pb))



sig=sigma()
print('\nsigma',sig)
print('\nKufper theorie\n mu=',n(z_cu,V_M_cu)*sig)

print('\nBlei theorie\n mu=',n(z_pb,V_M_pb)*sig)
re=2.82*10**(-15)
print(' thomschonscher schitt mu=',n(z_pb,V_M_pb)*(8/3)*np.pi*(re**2) )



def N_D(D,N_0,mu):
    return(N_0*np.exp(-mu*D))


T_cu ,n_cu ,d_cu =np.genfromtxt('messdatenGcu.txt',unpack=True)
d_cu=d_cu*10**(-3)
D_cu=unp.uarray(d_cu,0.0002)
n_cu=unp.uarray(n_cu,np.sqrt(n_cu))
print(n_cu)

Zy_0=unp.uarray(935,np.sqrt(935))
Ny_0=Zy_0/900
print('Nullmessung yZ',Zy_0,'\n YNs',Ny_0)

N_cu=n_cu/T_cu
print(N_cu)
print('Abgezogen von N_0',N_cu-Ny_0)
N_cu=N_cu-Ny_0
params_cu, covariance_cu =curve_fit(N_D,d_cu,noms(N_cu))
errors_cu = np.sqrt(np.diag(covariance_cu))


N_0_cu=unp.uarray(params_cu[0],errors_cu[0])
mu_cu=unp.uarray(params_cu[1],errors_cu[1])
print('\n N_0_cu=',N_0_cu)
print('mu_cu=',mu_cu)


#*params
x=np.linspace(0,0.070)
plt.figure(1)
plt.errorbar(d_cu,noms(N_cu) ,xerr=stds(D_cu),yerr=stds(N_cu), fmt='cx')
plt.plot(d_cu,noms(N_cu),'rx',label=r'$Messwerte$')
plt.plot(x,N_D(x,*params_cu),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Dicke \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Zählrate \ N \ in \ 1/s}$')
plt.savefig('a)cu.pdf')


T_pb ,n_pb ,d_pb =np.genfromtxt('messdatenGpb.txt',unpack=True)
d_pb=d_pb*10**(-3)
D_pb=unp.uarray(d_pb,0.0002)
n_pb=unp.uarray(n_pb,np.sqrt(n_pb))
print(n_pb)
N_pb=n_pb/T_pb
print(N_pb)
N_pb=N_pb-Ny_0
print('abgezogen ny_0 pb',N_pb)
params_pb, covariance_pb =curve_fit(N_D,d_pb,noms(N_pb))
errors_pb = np.sqrt(np.diag(covariance_pb))


N_0_pb=unp.uarray(params_pb[0],errors_pb[0])
mu_pb=unp.uarray(params_pb[1],errors_pb[1])



print('N_0_pb=',N_0_pb)
print('mu_pb=',mu_pb)

plt.figure(2)
plt.errorbar(d_pb,noms(N_pb) ,xerr=stds(D_pb),yerr=stds(N_pb), fmt='cx')
plt.plot(d_pb,noms(N_pb),'rx',label=r'$Messwerte$')
plt.plot(x,N_D(x,*params_pb),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Dicke \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Zählrate \ N \ in \ 1/s}$')
plt.savefig('a)pb.pdf')




T_b, n_b ,d_b =np.genfromtxt('messdatenB.txt',unpack=True)
n_b=unp.uarray(n_b,np.sqrt(n_b))
df=np.array([1,1,0.5,1,1,1,1,5,1,2,1])
df=df*10**(-6)
D_b=unp.uarray(d_b,df)
D_b=D_b*10**3

R=2.7*D_b


Zb_0=unp.uarray(292,np.sqrt(292))
Nb_0=Zb_0/900
print(' \n Nullmessung yZ',Zb_0,'\n YNs',Nb_0)




print(D_b)
print(n_b)
N_b=n_b/T_b
N_b=N_b-Nb_0
print('Impulsrate beta',N_b)

Zb_0=unp.uarray(292,np.sqrt(292))
Nb_0=Zb_0/900




def fehler_b(std_m,x):
    return(std_m*np.sqrt(np.mean(x**2)))


bereich1x=R[:5]
bereich1y=N_b[:5]
A1 , B1 , r ,p ,stdA1 =stats.linregress(noms(bereich1x),np.log(noms(bereich1y)))

A1=unp.uarray(A1,stdA1)
B1=unp.uarray(B1,fehler_b(stdA1,noms(bereich1x)))
print('\nA1',A1,'B1',B1)

bereich2x=R[5:]
bereich2y=N_b[5:]

A2 , B2 , r ,p ,stdA2 =stats.linregress(noms(bereich2x),np.log(noms(bereich2y)))

A2=unp.uarray(A2,stdA2)
B2=unp.uarray(B2,fehler_b(stdA2,noms(bereich2x)))

print('\nA2',A2,'\nB2',B2)

R_max=(B2-B1)/(A1-A2)
print('\n R_max',(B2-B1)/(A1-A2))

E_max=1.92*unp.sqrt(R_max**2+0.22*R_max)

print('E_max',E_max)


def g(x,A,B):
 return A*x+B
x1=np.linspace(0,0.675)
x=np.linspace(0.675,1.4)
plt.figure(3)
plt.errorbar(noms(R),noms(N_b) ,xerr=stds(R),yerr=stds(N_b), fmt='cx')
plt.plot(noms(bereich1x),noms(bereich1y),'gx',label=r'$\mathrm{Messwerte\ für\ erste\ Gerade}$')
plt.plot(noms(bereich2x),noms(bereich2y),'bx',label=r'$\mathrm{Messwerte\ für\ zweite\ Gerade}$')
plt.plot(x1,np.exp(g(x1,noms(A1),noms(B1))),'g-',label=r'$\mathrm{Ausgleichsgerade\ für\ ersten\ Bereich}$')
plt.plot(x,np.exp(g(x,noms(A2),noms(B2))),'b-',label=r'$\mathrm{Ausgleichsgerade\ für\ zweite\ Bereich}$')
plt.yscale('log')
plt.ylim(0,100)
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Massenbelegung \ R \ in \ g/cm^2\ }$')
plt.ylabel(r'$\mathrm{Zählrate \ N \ in \ 1/s}$')
plt.savefig('b).pdf')

#plt.plot(x,N_D(x,*params_pb),'b-',label=r'$Ausgleichsfunktion$')
