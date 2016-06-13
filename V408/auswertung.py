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
#print(d1,d2,eb)

f_bessel=[f_Bessel(d1,eb),f_Bessel(d2,eb)]

#print('f1_bessel',f_Bessel(d1,eb))
#print('f2_bessel',f_Bessel(d2,eb))

print('fbessel mittelwert:',np.mean(f_bessel),'+-',np.std(f_bessel))


#rot
b1rot  ,  g1rot   , b2rot ,    g2rot   , erot = np.genfromtxt('rot.txt',unpack=True)


d1rot=np.abs(b1rot-g1rot)
d2rot=np.abs(b2rot-g2rot)
#print(d1,d2,eb)

f_bessel_rot=[f_Bessel(d1rot,erot),f_Bessel(d2rot,erot)]

#print('f1_bessel rot' ,f_Bessel(d1rot,erot))
#print('f2_bessel rot' ,f_Bessel(d2rot,erot))

#print('fbessel rot mittelwert:',np.mean(f_bessel_rot),'+-',np.std(f_bessel_rot))


b1blau  ,  g1blau   , b2blau ,    g2blau   , eblau= np.genfromtxt('blau.txt',unpack=True)

d1blau=np.abs(b1blau-g1blau)
d2blau=np.abs(b2blau-g2blau)
#print(d1,d2,eb)

f_bessel_blau=[f_Bessel(d1blau,eblau),f_Bessel(d2blau,eblau)]

#print('f1_bessel blau',f_Bessel(d1blau,eblau))
#print('f2_bessel blau',f_Bessel(d2blau,eblau))

#print('fbessel mittelwert blau:',np.mean(f_bessel_blau),'+-',np.std(f_bessel_blau))


#Methode von abbe

def f_abbe(f1,f2):
    D_abbe=(1/f1)+(1/f2)-(61/(f1*f2))
    return 1/D_abbe


print('abbe theorie',f_abbe(-100,100))



B_abbe  ,  b_abbe,  g_abbe   , e_abbbe = np.genfromtxt('abbe.txt',unpack=True)

V=B_abbe/30
print('Abbildungsmaßstab',V)



def fehler_b(std_m,x):
    return(std_m*np.sqrt(np.mean(x**2)))


Ab , Bb , rb ,pb ,stdb =stats.linregress(1+V,b_abbe)
Bb_fehler=fehler_b(stdb,1+V)


x=np.linspace(0,6)
plt.figure(3)
#plt.errorbar(noms(R),noms(N_b) ,xerr=stds(R),yerr=stds(N_b), fmt='cx')
plt.plot(1+V,b_abbe,'rx',label=r'$\mathrm{Messwerte}$')
plt.plot(x,Ab*x+Bb,'b-',label=r'$\mathrm{Ausgleichsfunktion}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{1+V}$')
plt.ylabel(r'$\mathrm{Bildweite \ b/mm}$')
plt.savefig('plotabbeb.pdf')



Ag , Bg , rg ,pg ,stdg =stats.linregress(1+(1/V),g_abbe)

Bg_fehler=fehler_b(stdg,1+(1/V))

plt.figure(4)
#plt.errorbar(noms(R),noms(N_b) ,xerr=stds(R),yerr=stds(N_b), fmt='cx')
plt.plot(1+(1/V),g_abbe,'rx',label=r'$\mathrm{Messwerte}$')
plt.plot(x,Ag*x+Bg,'b-',label=r'$\mathrm{Ausgleichsfunktion}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{1+\frac{1}{V}}$')
plt.ylabel(r'$\mathrm{Gegenstandweite \ g/mm}$')
plt.savefig('plotabbeg.pdf')

print('\n f_abbe 1+V',Ab,'+-',stdb)
print('h`',Bg,'+-',Bg_fehler)


print('\n f_abbe 1+1/V',Ag,'+-',stdg)
print('h',Bb,'+-',Bb_fehler)

f_g=unp.uarray(Ag,stdg)
f_b=unp.uarray(Ab,stdb)

print((f_g+f_b)/2)

print(np.mean([f_g,f_b]))
