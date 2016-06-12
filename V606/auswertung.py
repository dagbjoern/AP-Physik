import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const



#def gaus1(x,mu,sigma):
#    return((1/(sigma*np.sqrt(2*np.pi))*np.exp(-0.5*((x-mu)/sigma)**2)))

def gaus2(x,a,mu,sigma):
  	return(a*np.exp(-0.5*((x-mu)/sigma)**2))

#Gausplot
v , U_A =np.genfromtxt('a).txt', unpack=True)
v=1000*v

n = len(v)                          #the number of data
mean = np.sum(v*U_A)/n                   #note this correction
sigma_test = np.sqrt(np.sum(U_A*(v-mean)**2)/n)
print(mean)
print(sigma_test)
params, covariance = curve_fit(gaus2,v,U_A,p0=[1.14,35300,300])
x=np.linspace(30000,40000,10000)
plt.figure(1)
print(params)
plt.plot(v,U_A,'rx',label=r'$Messwerte$')
#plt.plot(x,gaus2(x,1.14,35300,290),'k-',label=r'$Messwerte$')
plt.plot(x,gaus2(x,*params),'b-',label=r'$Ausgleichsfunktion$')
plt.legend(loc='best')
plt.xlabel(r'$Frequenz \ \nu /Hz$')
plt.ylabel(r'$frac{U_A}{U_E} $')
plt.savefig('a).pdf')
#
def fit(x):
    return(gaus(x,*params)-(1/np.sqrt(2))*np.max(U_A))

    print(scipy.optimize.fsolve(fit(x,*params),0))

Uvor_ND , Rvor_ND, Unach_ND , Rnach_ND=np.genfromtxt('b)ND.txt',unpack=True)
#
Uvor_DY , Rvor_DY, Unach_DY , Rnach_DY=np.genfromtxt('b)DY.txt',unpack=True)
#
Uvor_GD , Rvor_GD, Unach_GD , Rnach_GD=np.genfromtxt('b)GD.txt',unpack=True)

Rvor_ND=Rvor_ND*5e-3
Rvor_DY=Rvor_DY*5e-3
Rvor_GD=Rvor_GD*5e-3

Rnach_ND=Rnach_ND*5e-3
Rnach_DY=Rnach_DY*5e-3
Rnach_GD=Rnach_GD*5e-3


Delta_U_ND=Uvor_ND-Unach_ND
Delta_U_DY=Uvor_DY-Unach_DY
Delta_U_GD=Uvor_GD-Unach_GD



Delta_R_ND=Rvor_ND-Rnach_ND
Delta_R_DY=Rvor_DY-Rnach_DY
Delta_R_GD=Rvor_GD-Rnach_GD
print('\n \nRvorND',Rvor_ND,'\n RnachND', Rnach_ND,'\n Delta_R_ND',Delta_R_ND,'\n Delta_U_ND',Delta_U_ND )
print('\n \nRvorGD',Rvor_GD,'\n RnachGD', Rnach_GD,'\n Delta_R_GD',Delta_R_GD,'\n Delta_U_GD',Delta_U_GD)
print('\n \nRvorDY',Rvor_DY,'\n RnachDY', Rnach_DY,'\n Delta_R_DY',Delta_R_DY,'\n Delta_U_DY',Delta_U_DY)

#Dichten der Stoffel

#R3
R3=998
#daten der Messspule
n=250
F=86.6e-6
l=125e-3
R=0.7


def Qreal(M,L,rho):
    return(M/(L*rho))


#eignschafeten von den proben
#ND203
rhow_ND = 7240
M_ND=0.0095# masse der probe
l_ND=0.16 #länge der probe
Q_ND=Qreal(M_ND,l_ND,rhow_ND)

#GD203
rhow_GD = 7400
M_GD=0.01408
l_GD=0.16
Q_GD=Qreal(M_GD,l_GD,rhow_GD)
#DJ203
rhow_DY = 7800
M_DY=0.0185
l_DY=0.16
Q_DY=Qreal(M_DY,l_DY,rhow_DY)


print('Qreal',Q_ND,Q_GD,Q_DY)

def Fall1(U_br,Q):
    w=2*np.pi*35300
    X=(U_br/100)*(4*l/(w*const.mu_0*(n**2)*Q))*np.sqrt((R**2)+(w**2)*(const.mu_0*(n**2)*F/l)**2)
    print('\n Test mit Vereinfacherung \n',(4*F*U_br)/(Q*100))
    return X


def Fall2(Delta_R,Q):
    return(2*Delta_R*F/(R3*Q))

def Fall3(N,J,gj):
    mu_B=(const.e*const.hbar)/(2*const.m_e)
    return(((const.mu_0*(mu_B**2)*(gj**2)*N*J)*(J+1))/(3*const.k*293))


J=np.array([9/2,7/2,15/2])
gj=np.array([8/11,2,4/3])




#momente pro Volumeneinheit
Mol=np.array([0.33684,0.3625,0.372998])
rho_alle=np.array([7240,7400,7800])
N=2*const.N_A*rho_alle/Mol


theo=Fall3(N,J,gj)
print('X für theorie',Fall3(N,J,gj))
print('\n X für ND Fall 1',Fall1(Unach_ND,Q_ND),
'\n X für ND Fall 2',Fall2(Delta_R_ND,Q_ND),
'\n X für ND Fall 1 Mittelwert',np.mean(Fall1(Unach_ND,Q_ND)),np.std(Fall1(Unach_ND,Q_ND)),
'\n X für ND Fall 2 Mittelwert',np.mean(Fall2(Delta_R_ND,Q_ND)),np.std(Fall2(Delta_R_ND,Q_ND)),
'\n relative Abweichung fall 1',(theo[0]-np.mean(Fall1(Unach_ND,Q_ND)))/theo[0],
'\n relative Abweichung fall 2',(theo[0]-np.mean(Fall2(Delta_R_ND,Q_ND)))/theo[0],)



print('\nX für GD Fall 1',Fall1(Unach_GD,Q_GD),
'\n X für GD Fall 2',Fall2(Delta_R_GD,Q_GD),
'\n X für GD Fall 1 Mittelwert',np.mean(Fall1(Unach_GD,Q_GD)),np.std(Fall1(Unach_GD,Q_GD)),
'\n X für ND Fall 2 Mittelwert',np.mean(Fall2(Delta_R_GD,Q_GD)),np.std(Fall2(Delta_R_GD,Q_GD)),
'\n relative Abweichung fall 1',(theo[1]-np.mean(Fall1(Unach_GD,Q_GD)))/theo[2],
'\n relative Abweichung fall 2',(theo[1]-np.mean(Fall2(Delta_R_GD,Q_GD)))/theo[2],)

print('\n X für DY Fall 1',Fall1(Unach_DY,Q_DY),
'\n X für DY Fall 2',Fall2(Delta_R_DY,Q_DY),
'\n X für DY Fall 1 Mittelwert',np.mean(Fall1(Unach_DY,Q_DY)),np.std(Fall1(Unach_DY,Q_DY)),
'\n X für DY Fall 2 Mittelwert',np.mean(Fall2(Delta_R_DY,Q_DY)),np.std(Fall2(Delta_R_DY,Q_DY)),
'\n relative Abweichung fall 1',(theo[2]-np.mean(Fall1(Unach_DY,Q_DY)))/theo[2],
'\n relative Abweichung fall 2',(theo[2]-np.mean(Fall2(Delta_R_DY,Q_DY)))/theo[2],)
