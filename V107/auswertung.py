import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit

k, g  = np.genfromtxt('werte.txt',unpack=True)
dg, dk = np.genfromtxt('wertekugel.txt',unpack=True)
t=np.array([21,27.5,30.5,34,37,42,47,52,57.5,61])
Temp=unp.uarray(t,2)
Temp=Temp+273.15
t1, t2, t3, t4, t5, t6, t7, t8, t9 = np.genfromtxt('temperatur.txt', unpack=True)
mk= 4.46
mg= 4.96
Kkl=0.07640e-3
rhow=0.99
print('Temperatur',Temp)


Dk=unp.uarray(np.average(dk),np.std(dk))
Dg=unp.uarray(np.average(dg),np.std(dg))

print('Dk',Dk)
print('Dg',Dg)


Vok=(4/3)*np.pi*(Dk/2)**3
Vog=(4/3)*np.pi*(Dg/2)**3
rhok=mk/Vok
rhog=mg/Vog


print('rhok',rhok)
print('rhog',rhog)
T= unp.uarray([np.average(t1),np.average(t2),np.average(t3),np.average(t4),np.average(t5),
np.average(t6),np.average(t7),np.average(t8),np.average(t9)],
[np.std(t1),np.std(t2),np.std(t3),np.std(t4),np.std(t5),
np.std(t6),np.std(t7),np.std(t8),np.std(t9)])
Tk=unp.uarray([np.average(k)],[np.std(k)])
Tg=unp.uarray([np.average(g)],[np.std(g)])
print('Tg',Tg)
print('Tg*2',Tg*2)
V= 0.1/T
Vk=0.1/Tk
Vg=0.05/Tg

Tg_10 = 0.1/Vg
print('tg',Tg_10)
V=np.append(Vg,V)
print('Volumen Groß',Vog)
print('Volumen klein',Vok)

T=np.append(Tg_10,T)
n=Kkl*(rhok-rhow)*Tk
print('Viskosität',n)
Kgr= n/((rhog-rhow)*Tg_10)
print('Kgr',Kgr)
ngr=Kgr*(rhog-rhow)*T

m, b, r, p, stdsm=stats.linregress(noms(1/Temp),np.log(noms(Kgr*(rhog-1)*T)))
x=np.linspace(20+273.15,65+273.15)

print('Geschwindigkeit klein=',Vk,'Zeitmittelwert=',Tk)
plt.figure(1)
plt.errorbar(noms(1/Temp),noms(unp.log(Kgr*(rhog-1)*T)),xerr=stds(1/Temp),yerr=stds(unp.log(Kgr*(rhog-1)*T)),fmt='rx')
plt.plot(noms(1/Temp),noms(unp.log(Kgr*(rhog-1)*T)),'xk',label=r'$Messwerte$')
plt.plot(1/x,m*(1/x)+b,'-b',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.grid(True)
plt.xlabel(r'$\frac{1}{T}/{K^{-1}}$')
plt.ylabel(r'$\ln(\eta / Pa \cdot s)$')
plt.savefig('plot.pdf')

def Fehlera(x,sa):
    b=sa*(sum(x**2)/len(x))**(1/2)
    return b
sa1=Fehlera(noms(1/Temp),stdsm)
lnA=unp.uarray(b,sa1)
Reg=(rhow*Dg)/(ngr*T)
Rek=(rhow*Dg)/(n*Tk)
A=unp.exp(lnA)
print('Reg',Reg)
print('Rek',Rek)
#print(Temp)
#np.savetxt('v.txt',np.column_stack((Temp,T,ngr)),fmt='%r',delimiter=' & ')
#print('B',m,stdsm)
#print('lnA',lnA)
#print('A',A)
