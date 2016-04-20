import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit


eich_phi=305.5
lam_H, phi_a =np.genfromtxt('a).txt', unpack=True)
#lam_H=lam_H*10e-9
phi_H=eich_phi-phi_a
phi_H=2*np.pi*(phi_H/360)


m , b , r ,p ,std =stats.linregress(np.sin(phi_H),lam_H)
x=np.linspace(0.4,0.7)
plt.figure(1)
#plt.errorbar(t ,noms(y),xerr=stds(T),yerr=stds(y), fmt='rx')
plt.plot(np.sin(phi_H),lam_H,'kx',label=r'$Messwerte$')
plt.plot(x,m*x+b,label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$sin(\phi)$')
plt.ylabel(r'$Wellenlänge \ \ \frac{\lambda}{\mathrm{nm}} $')
plt.savefig('a).pdf')
m=unp.uarray(m,std)
g=1/m
print(g)
