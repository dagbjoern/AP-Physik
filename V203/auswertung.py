import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit


Tf , Tg , p1 =np.genfromtxt('1messung.txt',unpack=True)
T2 , p2 =np.genfromtxt('2messung.txt',unpack=True)
R=8.3144598

p1=p1*1e-3
Tg=Tg+273.15
T2=T2+273.15
p2=8.14+p2



m, bk, r, p, stdsm=stats.linregress(1/Tg,np.log(p1))

print

def F(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d

params, covariance =curve_fit(F,T2,p2)

print('linrgress')
print('a =', bk)
print('b =', m)



a = params[0]
b = params[1]
c = params[2]
d = params[3]

print('Polynom')
print('a =', params[0])
print('b =', params[1])
print('c =', params[2])
print('d =', params[3])




x1=np.linspace(1/(20+273.15),1/(100+273.15))
plt.figure(1)
plt.plot(1/Tg, np.log(p1),'rx',label=r'$Messwerte$')
plt.plot(x1,m*(x1)+bk,'b-',label=r'$Ausgleichsfunktion$')
plt.xlabel(r'$Kehrwert \ der \ Temperatur \ \frac{K}{T}$')
plt.ylabel(r'$Logarithmus \ des \  Druckes \  \log\left(\frac{p}{bar}\right)$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot1.pdf')


x2=np.linspace(20+273.15,200+273.15)
plt.figure(2)
plt.plot(T2, p2,'rx',label=r'$Messwerte$')
plt.plot(x2,F(x2,*params),'b-',label=r'$Ausgleichsfunktion$')
plt.xlabel(r'$Temperatur \  \frac{T}{K}$')
plt.ylabel(r'$Druck \  \frac{p}{bar}$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot2.pdf')

def p_t(T):
    return a*T**3+b*T**2+c*T+d

def p_t_abl(T):
    return 3*a*T**2+2*b*T+c

def L_max(T):
    return T*((R*T/(2*p_t(T)))+np.sqrt((R*T/(2*p_t(T)))**2-(0.9/p_t(T))))*p_t_abl(T)


def L_min(T):
    return T*((R*T/(2*p_t(T)))-np.sqrt((R*T/(2*p_t(T)))**2-(0.9/p_t(T))))*p_t_abl(T)

#plt.figure(1)
#plt.errorbar(noms(Tb**2) ,noms(A**2),xerr=stds(Tb**2),yerr=stds(A**2), fmt='rx')
#plt.plot(noms(Tb**2),noms(A**2),'kx',label=r'$Messwerte$')
#x=np.linspace(0,16)
#plt.plot(x,m*x+b,label=r'$Ausgleichsfunktion$')
#plt.legend(loc='best')
#plt.xlabel(r'$\frac{T^2}{s^2}  $')
#plt.ylabel(r'$\frac{a^2}{m^2} $')
#plt.savefig('b.pdf')


plt.figure(3)
plt.plot(x2,L_min(x2),'b-',label=r'$Verdampfungswärmekurve$')
plt.xlabel(r'$Temperatur \ \frac{T}{K}$')
plt.ylabel(r'$Verdampfungswärme \ \frac{L}{J \cdot mol^{-1}}$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot3.pdf')

plt.figure(4)
plt.plot(x2,L_max(x2),'b-',label=r'$Verdampfungswärmekurve$')
plt.xlabel(r'$Temperatur \ \frac{T}{K}$')
plt.ylabel(r'$Verdampfungswärme \ \frac{L}{J \cdot mol^{-1}}$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot4.pdf')


np.savetxt('1tabelle.txt',np.column_stack((Tg,p1)),fmt='%r',delimiter=' & ')
np.savetxt('2tabelle.txt',np.column_stack((T2,p2)),fmt='%r',delimiter=' & ')
