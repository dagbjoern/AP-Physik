import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const
from uncertainties import ufloat


U , n_10 , I=np.genfromtxt('a.txt',unpack=True)
U_lin , n_lin , I_lin=np.genfromtxt('a.txt',skip_header=3,skip_footer=2,unpack=True)

n=n_10/10 #Zählrate pro secunde
n_lin=n_lin/10
I=unp.uarray(I*10**(-6),0.000002)
N=unp.uarray(n,np.sqrt(n))
print(N)
def fehler_b(std_m,x):
    return(std_m*np.sqrt(np.mean(x**2)))

m , b , r ,p ,std =stats.linregress(U_lin,n_lin)
x=np.linspace(375,670)
M=unp.uarray(m,std)
B=unp.uarray(b,fehler_b(std,U_lin))
#Abgeschätzter fehler  für Zähler 10

plt.figure(1)
plt.errorbar(U ,n,xerr=0,yerr=np.sqrt(n), fmt='cx')
plt.plot(U,n,'gx',label=r'$\mathrm{Messwert \ außerhalb \ des \ Plateaus}$')
plt.plot(U_lin,n_lin,'bx',label=r'$\mathrm{Messwerte \ innerhalb \ des \ Plateaus}$')
plt.plot(x,m*x+b,'b-',label=r'$\mathrm{Ausgleichsfunktion}$')
#plt.plot(300,0,'r.')
plt.legend(loc='best')
plt.xlabel(r'$Spannung \ \ U/V$')
plt.ylabel(r'$Teilchenzahl   \ \ N/s $')
plt.savefig('a).pdf')

print('y-Achsenabschnitt B',B)
print('steigung',M)
print('plateau-steigung:',(M*400+B)/(M*300+B),)
print('Fehler:',(m*100)/299.3)
print('zeitlicher Abstand zwischen primär und nachentladungsimpulsen',
'\n bei 700V in micro s',105,
'\n bei 350V in micro s',120)
print('Totzeit bei 700V mit ozill: in micro s',80 )

n1=ufloat(195.2,np.sqrt(195.5))
n12=ufloat(417.9,np.sqrt(417.9))
n2=ufloat(236.0,np.sqrt(236.0))
print('n1=','{:.5u}'.format(n1))
print('n2=','{:.5u}'.format(n2))
print('n12=','{:.5u}'.format(n12))
def T(N1,N12,N2):
    return (N1+N2-N12)/(2*N1*N2)

print('totzeit bei 2.methode', T(n1,n12,n2))


def delta_Q(I,N):
    return(I/N)

deltaQ=delta_Q(I,N)

#print('\n deltaQ \n',deltaQ)
print('\n Q in -Q/e \n',deltaQ/const.e)
