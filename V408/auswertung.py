import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const



#Berechnung

def f(b,g):
    D=(1/b)+(1/g)
    #print('Brechkraft',D)
    return 1/D

#Erste Messung
g1, b1 , B1 =np.genfromtxt('messung1.txt',unpack=True)
g2, b2 =np.genfromtxt('messung2.txt',unpack=True)


B1=B1*1e1
g1=g1*1e1
b1=b1*1e1
g2=g2*1e1
b2=b2*1e1

print('verhältniss:')
print('G/B',30/B1)
print('g/b',g1/b1)

print('über die Linsen Gleichung:')
print('Messung1: \n',f(b1,g1))
print('mittelwert: \n', np.mean(f(b1,g1)),'+-',np.std(f(b1,g1)))
print('Messung2: \n',f(b2,g2))
print('mittelwert: \n', np.mean(f(b2,g2)),'+-',np.std(f(b2,g2)))



plt.figure(1)
for n in range(1,10,1):
    #plt.errorbar(noms(R),noms(N_b) ,xerr=stds(R),yerr=stds(N_b), fmt='cx')
    plt.plot([0,g1[n]],[b1[n],0],'r-')
plt.plot(0,0,'r-',label=r'$\mathrm{Messwerte}$')
plt.plot([100,100],[0,100],'b--',label=r'$\mathrm{Theorie}$')
plt.plot([0,100],[100,100],'b--')
plt.plot(100,100,'bx')
#plt.errorbar(95.8,95.8 ,xerr=3.4,yerr=3.4, fmt='cx')
#plt.plot([95.8,95.8],[0,95.8],'c--',label=r'$\mathrm{Linsengleichung}$')


#plt.ylim(0,100)
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Gegenstandweite \ g/mm}$')
plt.ylabel(r'$\mathrm{Bildweite \ b/mm}$')
plt.savefig('plot1.pdf')

#plt.plot(x,N_D(x,*params_pb),'b-',label=r'$Ausgleichsfunktion$')


plt.figure(2)
#plt.errorbar(noms(R),noms(N_b) ,xerr=stds(R),yerr=stds(N_b), fmt='cx')
for n in range(1,10,1):
    plt.plot([0,g2[n]],[b2[n],0],'r-')
plt.plot(0,0,'r-',label=r'$\mathrm{Messwerte}$')
plt.plot([150,150],[0,150],'b--',label=r'$\mathrm{Theorie}$')
plt.plot([0,150],[150,150],'b--')
plt.plot(150,150,'bx')
#plt.ylim(0,100)
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Gegenstandweite \ g/mm}$')
plt.ylabel(r'$\mathrm{Bildweite \ b/mm}$')
plt.savefig('plot2.pdf')


#Methode von Bessel

def f_Bessel(d, e):
    return (e**2-d**2)/(4*e)




b1b  ,  g1b   , b2b ,    g2b   , eb= np.genfromtxt('bessel.txt',unpack=True)

d1=np.abs(b1b-g1b)
d2=np.abs(b2b-g2b)
print(d1,d2,eb)


print('f1_bessel',f_Bessel(d1,eb))
print('f2_bessel',f_Bessel(d2,eb))
