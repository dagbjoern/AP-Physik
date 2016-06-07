import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const


def einzel(phi,b,A_0):
    lamda=633e-9
    teil1=(A_0**2)*(b**2)*((lamda/(np.pi*b*np.sin(phi)))**2)
    teil2=np.sin((np.pi*b*np.sin(phi))/lamda)**2
    return(teil1*teil2)


def doppel(phi,b,A_0,s):
    lamda=633e-9
    teil1=(A_0**2)*np.cos((np.pi*s*np.sin(phi)/lamda))**2
    teil2=(lamda/(np.pi*b*np.sin(phi)))**2
    teil3=np.sin((np.pi*b*np.sin(phi))/lamda)**2
    return(teil1*teil2*teil3)


def a(ag,at):
    return(np.absolute(ag-at)/(at))

L=0.9 #abstand schirm spalt
I_dunkelstrom=0.225*1e-9#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

d_1, I_1= np.genfromtxt('spalt1.txt',unpack=True)
d_2, I_2= np.genfromtxt('spalt2.txt',unpack=True)
d_3, I_3= np.genfromtxt('spalt3.txt',unpack=True)
d_dopp , I_dopp =np.genfromtxt('doppelspalt.txt',unpack=True)

print(I_1,1e-3)


print(d_1)
d_1=d_1*1e-3
I_1=I_1*1e-6-I_dunkelstrom
phi_1=d_1/L

print(I_1)

d_2=d_2*1e-3
I_2=I_2*1e-6-I_dunkelstrom
phi_2=d_2/L

d_3=d_3*1e-3
I_3=I_3*1e-6-I_dunkelstrom
phi_3=d_3/L

I_dopp=I_dopp*1e-6-I_dunkelstrom
#np.savetxt('doppelspaltab.txt',np.column_stack((d_dopp[:49],I_dopp[:49],d_dopp[49:],I_dopp[49:])),fmt='%r',delimiter=' & ')
d_dopp=d_dopp*1e-3
phi_dopp=d_dopp/L

phi_1f=phi_1[1:]
phi_2f=phi_2[1:]
phi_3f=phi_3[1:]
phi_doppf=phi_dopp[1:]

I_1f=I_1[1:]
I_2f=I_2[1:]
I_3f=I_3[1:]
I_doppf=I_dopp[1:]

print('phi_dopp',phi_doppf)


params_1, covariance_1 =curve_fit(einzel,phi_1f,I_1f,p0=[0.000075,1])
errors_1 = np.sqrt(np.diag(covariance_1))
print('\n1')
print('breite des Spaltes:', params_1[0],'+-',errors_1[0])
print('Abweichung',a(params_1[0],0.000075))
print('A_0:',params_1[1],'+-',errors_1[1])
#def einzel(phi,b,lamda,A_0)

x=np.linspace(-0.025,0.025,10000)

plt.figure(1)
#plt.errorbar(d_cu,noms(N_cu) ,xerr=stds(D_cu),yerr=stds(N_cu), fmt='cx')
plt.plot(phi_1,I_1,'rx',label=r'$Messwerte$')
plt.plot(x,einzel(x,*params_1),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Abstand \ von \ 0.-Hauptmaxima \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Intensität \ I/A}$')
plt.savefig('spalt1.pdf')

params_2, covariance_2 =curve_fit(einzel,phi_2f,I_2f,p0=[0.00015,5])
errors_2 = np.sqrt(np.diag(covariance_2))
print('\n2')
print('breite des Spaltes:', params_2[0],'+-',errors_2[0])
print('Abweichung',a(params_2[0],0.00015))
print('A_0:',params_2[1],'+-',errors_2[1])

plt.figure(2)
#plt.errorbar(d_cu,noms(N_cu) ,xerr=stds(D_cu),yerr=stds(N_cu), fmt='cx')
plt.plot(phi_2,I_2,'rx',label=r'$Messwerte$')
plt.plot(x,einzel(x,*params_2),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Abstand von 0.-Hauptmaxima \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Intensität \ I/A}$')
plt.savefig('spalt2.pdf')


params_3, covariance_3 =curve_fit(einzel,phi_3f,I_3f,p0=[0.0004,1])
errors_3 = np.sqrt(np.diag(covariance_3))

print('\n3')
print('breite des Spaltes:', params_3[0],'+-',errors_3[0])
print('Abweichung',a(params_3[0],0.0004))
print('A_0:',params_3[1],'+-',errors_3[1])


plt.figure(3)
#plt.errorbar(d_cu,noms(N_cu) ,xerr=stds(D_cu),yerr=stds(N_cu), fmt='cx')
plt.plot(phi_3,I_3,'rx',label=r'$Messwerte$')
plt.plot(x,einzel(x,*params_3),'b-',label=r'$Ausgleichsfunktion$')
#plt.plot(x,einzel(x,0.0004,3),'g-',label=r'$Ausgleichsfunktion$')

plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Abstand von 0.-Hauptmaxima \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Intensität \ I/A}$')
plt.savefig('spalt3.pdf')

#def doppel(phi,b,lamda,A_0,s):

params_dopp, covariance_dopp =curve_fit(doppel,phi_doppf,I_doppf,p0=[0.0001,1,0.0004])
errors_dopp = np.sqrt(np.diag(covariance_dopp))

print('\n Doppel')
print('breite des Spaltes:', params_dopp[0],'+-',errors_dopp[0])
print('Abweichung',a(params_dopp[0],0.0001))
print('A_0:',params_dopp[1],'+-',errors_dopp[1])
print('abstand',params_dopp[2],'+-',errors_dopp[2])
print('Abweichung',a(params_dopp[2],0.0004))




plt.figure(4)
#plt.errorbar(d_cu,noms(N_cu) ,xerr=stds(D_cu),yerr=stds(N_cu), fmt='cx')
plt.plot(d_dopp,I_dopp,'rx',label=r'$Messwerte$')
plt.plot(x,doppel(x,*params_dopp),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{winkel von 0.-Hauptmaxima \ d \ in \ m}$')
plt.ylabel(r'$\mathrm{Intensität \ I/A}$')
plt.savefig('doppelspalt.pdf')

print((d_dopp[:51]))
print()
