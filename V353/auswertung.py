import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
#a ,b= np.genfromtxt('text.txt', unpack=True)
t, u_c = np.genfromtxt('a.txt', unpack=True)
U_c=unp.uarray(u_c,0.05)
T=unp.uarray(t,0.05)
U_0=unp.uarray(19.4,0.05)
y=1-(U_c/U_0)
#phi=(a/b)*2*np.pi
#plt.figure(1)
#plt.plot(1/(b*10**-3),phi,'rx')
#plt.savefig('plot.pdf')
print(T)
index= [12,13]
Tt=np.delete(T, index)
print(Tt)
yy=np.delete(y, index)

m , b , r ,p ,std =stats.linregress(noms(Tt),np.log(noms(yy)))
print('m', m)
plt.figure(1)
plt.errorbar(t ,noms(y),xerr=stds(T),yerr=stds(y), fmt='rx')
plt.plot(noms(T),noms(y),'kx',label=r'$Messwerte$')
x=np.linspace(0,5)
plt.plot(t, np.exp(m*t+b),label=r'$Ausgleichsfunktion$')
plt.yscale('log')
plt.legend(loc='best')
plt.xlabel(r'$t/s$')
plt.ylabel(r'$ln(1-(U_c/U_0)) $')
plt.savefig('a.pdf')
M=unp.uarray(m,std)
print('RC ist = ', 1/M, M)
print('Delta Rc',std/m**2)

w , a_w =np.genfromtxt('b.txt',unpack=True)
# RC=0.001388
A_w=unp.uarray(a_w,0.0005)
U_0=unp.uarray(6.8, 0.05)
W=unp.uarray(w,0.5)
AU=A_w/U_0
def F(x,  RC ):
    return 1/((1+((x)**2)*RC**2)**(1/2))

params, pcov = curve_fit(F, noms(W), noms(AU))
plt.figure(2)
x=np.logspace(1,5,1000)
plt.plot(x,F(x,params),label=r'$Ausgleichsfunktion$')
plt.errorbar(w,noms(AU),xerr=stds(W),yerr=stds(AU), fmt='rx')
plt.plot(w,noms(AU),'kx',label=r'$Messwerte$')
plt.legend(loc='best')
plt.xscale('log')
plt.xlabel(r'$ f/s^{-1}$')
plt.ylabel(r'$ U_a/U_0$')
plt.savefig('b.pdf')
print('RC b)=',params, 'fehler =' ,pcov)

RC=params


w, a ,b =np.genfromtxt('c.txt',unpack=True)
W=unp.uarray(w,0.5)
A=unp.uarray(a,0.0005)
B=unp.uarray(b,0.0005)
phi=A/B*2*np.pi
def D(x, RC):
    return np.arctan(-x*RC)
params, pcov = curve_fit(D, noms(W),noms(phi))
x=np.logspace(1,4,1000)
plt.figure(3)
plt.plot(x,D(x,params),label=r'$Ausgleichsfunktion$')
plt.errorbar(w, noms(phi),xerr=stds(W),yerr=stds(phi),fmt='rx')
plt.plot(noms(w),noms(phi),'kx',label=r'$Messwerte$')
plt.xscale('log')
plt.xlabel(r'$f/s^{-1}$')
plt.legend(loc='best')
plt.ylabel(r'$\phi/rad$')
plt.savefig('c.pdf')
print('RbC c)=',params, 'fehler =' ,pcov)
print(A/B*2*np.pi)



w, a ,b =np.genfromtxt('c.txt',unpack=True)
w, a_w  =np.genfromtxt('d.txt',unpack=True)

A_w=unp.uarray(a_w,0.0005)
U_0=unp.uarray(6.8, 0.05)

phi=a/b*2 *np.pi
def A(x,RC):
    return (-np.sin(np.arctan(-x*RC))/(x*RC))
AU=A_w/U_0
plt.figure(4)
plt.polar(noms(phi),noms(AU),'rx',label=r'$Messwerte$')
x=np.linspace(0.0001,2000,10000)
plt.polar(D(x,RC),A(x,RC),label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.savefig('d.pdf')
print (phi)
