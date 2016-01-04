import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit

#Gänge
G=np.array([6,12,18,24,30,36,42,48,54,60])
#geschwindigkeiten
_V6,_V12,_V18,_V24,_V30,_V36,_V42,_V48,_V54,_V60=np.genfromtxt('Vwagen.txt', unpack=True)
#frequenzen
_6v,_6r,_12v,_12r,_18v,_18r,_24v,_24r,_30v,_30r=np.genfromtxt('Deffekt.txt',unpack=True)
_36v,_36r,_42v,_42r,_48v,_48r,_54v,_54r,_60v,_60r=np.genfromtxt('Deffekt2.txt',unpack=True)
#wellenlängen
w1, w2 =np.genfromtxt('Wellen.txt', unpack=True)
#schwebung
Sg6,Sg12,Sg18,Sg24,Sg30,Sg36,Sg42,Sg48,Sg54,Sg60=np.genfromtxt('schwebung1.txt',unpack=True)
#länge Stecke
L=unp.uarray(0.45,0.005)
#grundfequenz
v_0=unp.uarray(20.7416,0)
#wellenlänge
Lw=np.array([0.0,0.0,0.0,0.0])


Lw[0]=w1[1]-w1[0]
Lw[1]=w1[2]-w1[1]
Lw[2]=w2[1]-w2[0]
Lw[3]=w2[2]-w2[1]

Ld=unp.uarray(np.average(Lw),np.std(Lw))
print('Wellenlänge',Ld)



c=v_0*Ld#Schallgeschwindigkeit

print('1/lamda',1/Ld)
print('v0/c',v_0/c)
print('Schallgeschwindigkeit',c)

#gemittelte Schwebung
Schw=unp.uarray([np.average(Sg6),np.average(Sg12),np.average(Sg18),np.average(Sg24),
np.average(Sg30),np.average(Sg36),
np.average(Sg42),np.average(Sg48),
np.average(Sg54),np.average(Sg60)],[np.std(Sg6),np.std(Sg12),np.std(Sg18),np.std(Sg24),np.std(Sg30),np.std(Sg36),np.std(Sg42),np.std(Sg48),np.std(Sg54),np.std(Sg60)])
Schw=Schw*1e-3
#gemittelte Zeit
T=unp.uarray([np.average(_V6),np.average(_V12),np.average(_V18),np.average(_V24),
np.average(_V30),np.average(_V36),
np.average(_V42),np.average(_V48),
np.average(_V54),np.average(_V60)],[np.std(_V6),np.std(_V12),np.std(_V18),np.std(_V24),np.std(_V30),np.std(_V36),np.std(_V42),np.std(_V48),np.std(_V54),np.std(_V60)])
#geschwindigkeiten der einzelnen Gänge
Gesch=(L/T)

#Frequenz für messung c)
Vgesv=unp.uarray([np.average(_6v),np.average(_12v),np.average(_18v),np.average(_24v),
np.average(_30v),np.average(_36v),np.average(_42v),np.average(_48v),np.average(_54v),
np.average(_60v)],[np.std(_6v),np.std(_12v),np.std(_18v),np.std(_24v),np.std(_30v),
np.std(_36v),np.std(_42v),np.std(_48v),np.std(_54v),np.std(_60v)])#frequenz vorwärts

Vgesr=unp.uarray([np.average(_6r),np.average(_12r),np.average(_18r),np.average(_24r),
np.average(_30r),np.average(_36r),np.average(_42r),np.average(_48r),np.average(_54r),
np.average(_60r)],[np.std(_6r),np.std(_12r),np.std(_18r),np.std(_24r),np.std(_30r),
np.std(_36r),np.std(_42r),np.std(_48r),np.std(_54r),np.std(_60r)])#frequenz rückwärts



ve=v_0+(2*Gesch[9]/Ld)
vq=v_0+(1/(1-(2*Gesch[9]/c)))
print('ve',ve)
print('vq',vq)#unterschied zwischen gleichung 2 und 5
mv , bv , rv ,pv , stdv =stats.linregress(noms(Gesch),noms(Vgesv-v_0))
mr , br , rr ,pr , stdr =stats.linregress(noms(-Gesch),noms(Vgesr-v_0))


ms, bs , rs ,ps , std_s =stats.linregress(noms(Gesch*2),noms(Schw))
print(ms)

print(v_0-Vgesr[0])
print(v_0*(Gesch[0]/c))
print('ruhefrequenz',v_0)
#ist noch nicht ganz richtig
plt.figure(1)
x=np.linspace(-0.6,0.6)
plt.plot(x,mv*x+bv,'--r',label=r'$Ausgleichsfunktion Vorwärts$')
plt.plot(x,mr*x+br,'--b',label=r'$Ausgleichsfunktion Rückwärts$')
plt.plot(x,noms(v_0*(x/c)),'--k',label=r'$Aus Formel 4$')
plt.errorbar(noms(-Gesch),noms(Vgesr-v_0),xerr=stds(-Gesch),yerr=stds(Vgesr-v_0),fmt='bx')
plt.errorbar(noms(Gesch),noms(Vgesv-v_0),xerr=stds(Gesch),yerr=stds(Vgesv-v_0),fmt='rx')
plt.plot(noms(Gesch),noms(Vgesv-v_0),'rx',label=r'$Messwerte Vorwärts$')
plt.plot(noms(-Gesch),noms(Vgesr-v_0),'bx',label=r'$Messwerte Rückwärts$' )
plt.legend(loc='best')
plt.xlabel(r'$ \frac{v}{m\cdot s^{-1}} $')
plt.ylabel(r'$\frac{\Delta \nu}{s^{-1}}$')
plt.savefig('graph.pdf')

print('propfaktor Vorwärts',mv,stdv)
print('propfaktor Rückwärts',mr,stdr)
print('propschwebung',ms,std_s)

plt.figure(2)
x=np.linspace(0,1.2)
plt.errorbar(noms(2*Gesch),noms(Schw),xerr=stds(2*Gesch),yerr=stds(Schw),fmt='cx')
plt.plot(x,noms(v_0*(x/c)),'--k',label=r'$Aus \ Formel \ 4$')
plt.plot(noms(Gesch*2),noms(Schw),'cx',label=r'$Messwerte \ aus \ Schwebung$')
plt.plot(x,ms*x+bs,'--b',label=r'$Ausgleichsfunktion \ Schwebung$')
plt.legend(loc='best')
plt.xlabel(r'$ \frac{V}{m\cdot s^{-1}} $')
plt.ylabel(r'$\frac{\Delta \nu}{s^{-1}}$')
plt.savefig('graph2.pdf')


np.savetxt('tabellegeschwindigkeit.txt',np.column_stack((G,Gesch)),fmt='%r',delimiter=' & ')
np.savetxt('Frequenzänderungvorwärts.txt',np.column_stack((G,Gesch,Vgesv)),fmt='%r',delimiter=' & ')
np.savetxt('Frequenzänderungrückwärts.txt',np.column_stack((G,-Gesch,Vgesr)),fmt='%r',delimiter=' & ')
np.savetxt('Schwebung.txt',np.column_stack((G,Gesch*2,Schw)),fmt='%r',delimiter=' & ')
