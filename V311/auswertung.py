import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const

#formeln

#Wiederstand
def URI(I,U):
    return U/I

#ladungsträger pro Volumen

#Ladungsträger pro Atom

#mittlere Flugzeit

#mittlere freie Wellenlänge

#mittlere driftgeschwindigkeit

#Beweglichkeit

#Totalgeschwindigkeit



#zink
I_z, U_z =np.genfromtxt('messungRzink.txt',unpack=True)
I_hall_z, U_hall_z=np.genfromtxt('hallzink.txt',unpack=True)
b=0.024#breite
l=0.034#länge
U_z=U_z*10**(-3)
U_hall_z=U_hall_z*10**(-3)
R_z=URI(I_z,U_z)
print('R_Z=',R_z)



#wolfram
I_w, U_w =np.genfromtxt('messungRwolfram.txt',unpack=True)
I_hall_w, U_hall_w=np.genfromtxt('hallwolfram.txt',unpack=True)
b=0.026#breite
l=0.04#länge
U_w=U_w*10**(-3)
U_hall_w=U_hall_w*10**(-3)
R_w=URI(I_w,U_w)
print('R_w=',R_w)




#hysterese:

I_hys, B_hys=np.genfromtxt('hysterese.txt',unpack=True)



plt.figure(1)
plt.plot(I_hys,B_hys,'rx',label=r'$Messwerte$')
plt.legend(loc='best')
plt.xlabel(r'$Strom I/A$')
plt.ylabel(r'$Flussdichte B/T $')
plt.savefig('hysterese.pdf')
