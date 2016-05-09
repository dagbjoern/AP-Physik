import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const

c15, d15 ,c30 ,d30 , c60 , d60 =np.genfromtxt('messung1.txt',unpack=True)
tiefe,vmax,streumax,vmin,streumin =np.genfromtxt('messung2.txt',unpack=True)
c15=c15*10**(-3)
c30=c30*10**(-3)
c60=c60*10**(-3)

x=np.linspace(0,0.16)

plt.figure(1)
plt.plot(c15,d15/np.cos(2*np.pi*80.064/360),'rx',label=r'$Messwerte für \alpha=15°$')
plt.plot(x,2*2000000*x/343,'r-',label=r'$theoriekurve für 15°$')
plt.plot(c30,d30/np.cos(2*np.pi*70.529/360),'bx',label=r'$Messwerte für \alpha=30°$')
plt.plot(c60,d60/np.cos(2*np.pi*54.736/360),'gx',label=r'$Messwerte für \alpha=60°$')

plt.ylabel(r'$\Delta\nu/\cos \alpha$')
plt.xlabel(r'$strömungsgeschwindigkeit \ c$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot1.pdf')

plt.figure(2)
plt.plot(tiefe,vmax,'rx',label=r'$Messwerte \ für \ 60\% \ Strömungsleistung$')
plt.plot(tiefe,vmin,'bx',label=r'$Messwerte \ für \ 40\% \ Strömungsleistung$')
plt.xlabel(r'$Messtiefe  \ s/mm} $')
plt.ylabel(r'$Geschwindigkeit \ v/\mathrm{m}$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot2.pdf')

plt.figure(3)
plt.plot(tiefe,streumax,'rx',label=r'$Messwerte \ für\ 60\% \ Strömungsleistung$')
plt.plot(tiefe,streumin,'bx',label=r'$Messwerte \ für\ 40\% \ Strömungsleistung$')
plt.xlabel(r'$Messtiefe \ s/mm} $')
plt.ylabel(r'$Streuintensität \ \beta/1000\mathrm{V}^2\mathrm{s}^{-1}$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot3.pdf')
