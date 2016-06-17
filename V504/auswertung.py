import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const




# #Kennlinen
U, I1 , I2 , I3 , I4 , I5 =np.genfromtxt('messung1.txt',unpack=True)
# I_anlauf=np.genfromtxt('messung2.txt',unpack=True)

#U=np.array([[0:5:60],[70:10:250]]
U1=6.5
I1_leucht=2.5
U2=6
I2_leucht=2.4
U3=5.5
I3_leucht=2.3
U4=5.2
I4_leucht=2.2
U5=5
I1_leucht=2.1

I1=I1*1e-6
I2=I2*1e-6
I3=I3*1e-6
I4=I4*1e-6
I5=I5*1e-6
print(U)



plt.figure(1)
plt.plot(U,I1,'rx',label=r'$\mathrm{Messwerte für I=2,5}$')
plt.plot(U,I2,'rx',label=r'$\mathrm{Messwerte für I=2,4}$')
plt.plot(U,I3,'rx',label=r'$\mathrm{Messwerte für I=2,3}$')
plt.plot(U,I4,'rx',label=r'$\mathrm{Messwerte für I=2,2}$')
plt.plot(U,I5,'rx',label=r'$\mathrm{Messwerte für I=2,1}$')

#plt.plot([100,100],[0,100],'b--',label=r'$\mathrm{Theorie}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{spannung U}$')
plt.ylabel(r'$\mathrm{Bildweite \ b/mm}$')
plt.savefig('plot1.pdf')


#plt.plot(x,N_D(x,*params_pb),'b-',label=r'$Ausgleichsfunktion$')

#
#
# def fehler_b(std_m,x):
#     return(std_m*np.sqrt(np.mean(x**2)))
#
#
# Ab , Bb , rb ,pb ,stdb =stats.linregress(1+V,b_abbe)
