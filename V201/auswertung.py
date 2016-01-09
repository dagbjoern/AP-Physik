import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit

#Messung Zinn
m_z,m_zw,Tzw,Tzm,Tz=np.genfromtxt('zinn.txt', unpack=True)


#Messung Graphit
m_g=107.67
m_gw=465.18#
Tgw=293.40
Tgm=296.63
Tg=352.19


#konstanten für zinn und Graphit
rhoz=7.28/1e-6
rhog=2.25/1e-6
Molz=118.7
Molg=12.0
az=27.0e-6
ag=8.0e-6
kz=55.0e9
kg=33.0e9
cw=4.18


#nullmessung
Mwk=240.06
Mww=327.11
Twk=294.64
Tww=358.4
T0m=327.82
print(Mwk,Mww,Twk,Tww,T0m)

cgmg=(cw*Mww*(Tww-T0m)-cw*Mwk*(T0m-Twk))/(T0m-Twk)
print('cgmg',cgmg)


def ck(mw,mk,Tm,Tw,Tk):
 return ((cw*mw+cgmg)*(Tm-Tw))/(mk*(Tk-Tm))

def C(c,M,a,k,rho,Tm):         #molwärme
    return c*M-9*(a**2)*k*(M/rho)*Tm

def ab(am,at):
    return (am-at)/at


m_g=107.67
m_gw=465.18#
Tgw=293.40
Tgm=296.63
Tg=352.19

ckz=ck(m_zw,m_z,Tzm,Tzw,Tz)
ckg=ck(m_gw,m_g,Tgm,Tgw,Tg)
#ckg=ck(465.18,107.67,295.63,293.40,352.19)

ckz=unp.uarray(np.average(ckz),np.std(ckz))
Tzm=np.average(Tzm)

print('ckz',ckz)
print('ckg',ckg)
Cz=C(ckz,Molz,az,kz,rhoz,Tzm)
Cg=C(ckg,Molg,ag,kg,rhog,Tgm)


print('Cz',Cz)
print('Cg',Cg)
print('3R',3*8.3144598)

print('Abweichung von 3R z',ab(noms(Cz),3*8.3144598))
print('Abweichung von 3R g',ab(Cg,3*8.3144598))
#print(ckz)
#print('Graphit ck',ckg)
#print('Zinn ck',np.average(ckz))
