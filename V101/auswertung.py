import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit

phi,  f=np.genfromtxt('a.txt', unpack=True)
abst, tb=np.genfromtxt('b.txt', unpack=True)
db, hb=np.genfromtxt('zylinderb.txt', unpack=True)
hz, dz, tz= np.genfromtxt('Zylinder.txt', unpack=True)
dk, tk=np.genfromtxt('kugel.txt', unpack=True)
kd=np.genfromtxt('puppekopf.txt', unpack=True)
ad, bd, td, tg, ta=np.genfromtxt('puppe.txt', unpack=True)



Phi=unp.uarray(phi,5)
Phi=2*np.pi*(Phi/360)
print(Phi)
F=unp.uarray(f-0.05,0.03)
Abst=unp.uarray(abst*1e-3,0.0005)
Tb=unp.uarray(tb,0.5)
Tb=Tb/5
Mb=unp.uarray(0.22250,0.000005)
Db=unp.uarray(db*1e-3,0.0005)
Hb=unp.uarray(hb*1e-3,0.0005)
Hz=unp.uarray(hz*1e-3,0.0005)
Dz=unp.uarray(dz*1e-3,0.0005)
Tz=unp.uarray(np.average(tz),np.std(tz))
Tz=Tz/5
Dk=unp.uarray(dk*1e-3,0.0005)
Tk=unp.uarray(np.average(tk),np.std(tk))
Tk=Tk/5
Kd=unp.uarray(np.average(kd*1e-3),np.std(kd*1e-3))
Ad=unp.uarray(np.average(ad*1e-3),np.std(ad*1e-3))
Bd=unp.uarray(np.average(bd*1e-3),np.std(bd*1e-3))
Td=unp.uarray(np.average(td*1e-3),np.std(td*1e-3))
Tg=unp.uarray(np.average(tg),np.std(tg))
Tg=Tg/5
Ta=unp.uarray(np.average(ta),np.std(ta))
Ta=Ta/5
#Aufgabe 1


dbmit=np.average(noms(Db))
dbstd=np.std(noms(Db))
Dbmit=unp.uarray(dbmit,dbstd)

Rbmit=Dbmit/2


hbmit=np.average(noms(Hb))
hbstd=np.std(noms(Hb))
Hbmit=unp.uarray(hbmit,hbstd)

print(Dbmit,Hbmit)
A=Abst+(Hbmit/2)

D_feder=(F*137.5*1e-3/Phi)
print(np.average(noms(F/Phi)))

D_federmit=unp.uarray(np.average(noms(D_feder)),np.std(noms(D_feder)))
def FehlerB(x,sa):
    b=sa*(sum(x**2)/11)**(1/2)
    return b

print('D_feder Mittelwert=', D_federmit)
print('oder',np.average(noms(F/Phi))*137.5*1e-3)




m , b , r ,p , std =stats.linregress((noms(Tb)**2),(noms(A)**2))
mdyn=unp.uarray(m,std)
Ddyn=mdyn*(8*(np.pi**2)*Mb)
print('Ddyn=',Ddyn)

plt.figure(1)
plt.errorbar(noms(Tb**2) ,noms(A**2),xerr=stds(Tb**2),yerr=stds(A**2), fmt='rx')
plt.plot(noms(Tb**2),noms(A**2),'kx',label=r'$Messwerte$')
x=np.linspace(0,16)
plt.plot(x,m*x+b,label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\frac{T^2}{s^2}  $')
plt.ylabel(r'$\frac{a^2}{m^2} $')
plt.savefig('b.pdf')
#np.savetxt('tabelleA.txt',np.column_stack((Phi,F,D_feder)),fmt='%r' ,delimiter="   &   ")



I_D=-(b*2*Mb)-2*Mb*((Rbmit**2)/4 + (Hbmit**2)/12)



print('trägheitsmoment drill achse=',I_D)
print(b)

#Aufgabe 2
#Zylinder
dzmit=np.average(noms(Dz))
dzstd=np.std(noms(Dz))
Dzmit=unp.uarray(dzmit,dzstd)
Rzmit=Dbmit/2

print('Dz',Dzmit)
#Iz=Tz**2*D_federmit/(4*np.pi**2)#Gemssenes Trägheitsmoment
Iz=Tz**2*Ddyn/(4*np.pi**2)#Gemssenes Trägheitsmoment

Iz=Iz#abziehen der Drehachse

Mz=unp.uarray(1.9739,0.00005)#masse Zylinder

Iztheo=(Mz*(Rzmit**2))/2#theorie wert

print('I_Zylindermittelwert=',Iz)
print('I_Zylindermittelwerttheorie=',Iztheo)
print('rel',(Iz-Iztheo)/Iztheo)
#Kugel
Ik=(Tk**2)*Ddyn/(4*np.pi**2)#gemessenes Trägheitsmonent
#Ik=Ik-I_D#abziehen der Drehachse
dkmit=np.average(noms(Dk))
dkstd=np.std(noms(Dk))
Dkmit=unp.uarray(dkmit,dkstd)

print('Dk',Dkmit)
Rzmit=Dkmit/2

#Ik=Tk**2*D_federmit/(4*np.pi**2)
#Ik=Ik-I_D#abziehen der Drehachse

Mk=unp.uarray(0.8125,0.00005)

Iktheo=(2/5)*Mk*Rzmit**2



print('I_Kugelmittelwert=',Ik)
print('I_Kugeltheorie',Iktheo)
print('rel',(Ik-Iktheo)/(Iktheo))

#puppe
#Kd, Ad, Bd, Td, Tg, Ta

Mp=unp.uarray(0.3407,0.00005)
Kh=unp.uarray(72.2*1e-3,0.0005)
Ah=unp.uarray(np.average((0.183,0.1786)),np.std((0.183,0.1786)))
Bh=unp.uarray(np.average((0.2324,0.2329)),np.std((0.2324,0.2329)))
Th=unp.uarray(np.average((0.12384,0.1204,0.122)),np.std((0.12384,0.1204,0.122)))

def Vz(d,h):
    return np.pi*(d/2)**2*h


Vges=Vz(Kd,Kh)+2*Vz(Ad,Ah)+2*Vz(Bd,Bh)+Vz(Td,Th)

rho=Mp/Vges #dichte
def Izp (d,h,a):#trägheitsmoment für zylinder parrallel zu drehachse
    return ((rho*Vz(d,h)*(d/2)**2)/2)+rho*Vz(d,h)*a**2


def Izs (d,h,a): #trägheitsmoment für zylinder senktrecht zu drehachse
    return rho*Vz(d,h)*((((d/2)**2)/4)+(h**2)/12)+rho*Vz(d,h)*a**2


Igemessen_g=(Tg**2)*Ddyn/(4*np.pi**2)
#print(Igemessen_g)
#Igemessen_g=Igemessen_g-I_D#abziehen der Drehachse
#print(Igemessen_g)
Igemessen_a=(Ta**2)*Ddyn/(4*np.pi**2)
#print(Igemessen_a)
# Igemessen_a=Igemessen_a-I_D#abziehen der Drehachse
#print(Igemessen_a)

Itheorie_g=Izp(Kd,Kh,0)+Izp(Td,Th,0)+2*Izp(Bd,Bh,Bd/2)+2*Izp(Ad,Ah,(Td/2)+(Ad/2))
Itheorie_a=Izp(Kd,Kh,0)+Izp(Td,Th,0)+2*Izp(Bd,Bh,Bd/2)+2*Izp(Ad,Ah,(Td/2)+(Ah/2))

print('I_puppe  g =',Igemessen_g)
print('I_puppetheorie g=',Itheorie_g)
print('rel',(Igemessen_g-Itheorie_g)/(Itheorie_g))


print('I_puppe  a=',Igemessen_a)
print('I_puppetheorie a=',Itheorie_a)
print('rel',(Igemessen_a-Itheorie_a)/(Itheorie_a))



#Itheorie_g=(Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g,Itheorie_g)

#Itheorie_a=(Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a,Itheorie_a)

#np.savetxt('tabelleCg.txt',np.column_stack((Tg,Igemessen_g,Itheorie_g)),fmt='%r',delimiter='  &  ')
#np.savetxt('tabelleCa.txt',np.column_stack((Ta,Igemessen_a,Itheorie_a)),fmt='%r',delimiter='  &  ')
