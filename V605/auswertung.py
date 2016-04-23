import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const




eich_phi=305.5
lam_H, phi_a =np.genfromtxt('a).txt', unpack=True)
#lam_H=lam_H*10e-9
#phi_H=phi_a
phi_H=eich_phi-phi_a
phi_H=2*np.pi*(phi_H/360)

print(np.sin(phi_H))
print(phi_H,lam_H)
a=np.delete(phi_H,8)
b=np.delete(lam_H,8)
m , b , r ,p ,std =stats.linregress(np.sin(a),b)
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
print(m)
g=1/m
print(g)

stdA=std*(np.mean(np.sin(a)**2))**(1/2)
b=unp.uarray(b,stdA)#kackwert
print(b)
#b)
#berechnung der Eichgröße



def eichgr(lam1,lam2,delta_t,phi1,phi2):
 phimitt = unp.uarray([np.average([phi1,phi2])],[np.std([phi1,phi2])])
 print('phi12',phimitt)
 print('delta t',delta_t)
 return((lam1-lam2)/(delta_t*unp.cos(phimitt)))

grünstarkschwach=eichgr(lam_H[3],lam_H[4],120,phi_H[3],phi_H[4])

print('eichgröße von grün stark und schwach',eichgr(lam_H[3],lam_H[4],120,phi_H[3],phi_H[4]))
print('eichgröße von grün stark und blaugrün',eichgr(lam_H[4],lam_H[5],365,phi_H[4],phi_H[5]))
print('eichgröße gemittelt',(grünstarkschwach+eichgr(lam_H[4],lam_H[5],365,phi_H[4],phi_H[5]))/2)
eichzahl=(grünstarkschwach+eichgr(lam_H[4],lam_H[5],365,phi_H[4],phi_H[5]))/2

#c)




def delta_lambda(phi1,phi2,s):
    phi_mittel=unp.uarray([np.average([rad(phi1),rad(phi2)])],[np.std([rad(phi1),rad(phi2)])])
    return unp.cos(phi_mittel)*eichzahl*s

def rad(phi):
    return (phi/360)*2*np.pi


def lambdamittel(phi1,phi2):
   phi_mittel=unp.uarray([np.average([rad(phi1),rad(phi2)])],[np.std([rad(phi1),rad(phi2)])])
   return(unp.sin(phi_mittel)*g)

def delta_E(delta_l,l):
    return(const.h*const.c*delta_l/l)

def o2(delta_E,n,l,z):
    return(-((delta_E*(n**3)*l*(l+1)/(const.Rydberg*const.alpha**2))**(1/4))+z)


#werte für natrum n

#o2_natrium=32
lNmitr=lambdamittel(268,268)
d_lNr=delta_lambda(268,268,27)
d_ENr=delta_E(d_lNr,lNmitr)
print('rot\n,''\n winkelmittelwert=',rad(268),',\n lambdamittel=',lambdamittel(268,268),
',\n delta_lambda=', delta_lambda(268,268,27),',\n Delta E=',d_ENr,
'\n o2=',o2(d_ENr,3,1,11))


def alles(
element,
farbe,
z, #kernladung,
delta_s,#abstand zwischen den zwei dubletten
phi1,
phi2,
l,# quantendrehimpuls
n,#quantenzahl
):
    lXmitx=lambdamittel(phi1,phi2)
    d_lXx=delta_lambda(phi1,phi2,delta_s)
    d_ENx=delta_E(d_lXx,lXmitx)
    print(element,'\n',farbe,'\n winkelmittelwert=',(rad(-phi1+eich_phi)+rad(-phi2+eich_phi))/2
    ,',\n lambdamittelwert=',lambdamittel(phi1,phi1)
    ,',\n delta_lambda=',delta_lambda(phi1,phi2,delta_s)
    ,',\n Delta E=',delta_E(d_lXx,lXmitx)/const.e
    ,',\n o2=',o2(d_ENx,l,n,z))
    return()

alles('natrium','rot',11,27,268,268,3,1)
