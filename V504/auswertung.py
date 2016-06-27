import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const




# #Kennlinen
U, I1 , I2 , I3 , I4 , I5 =np.genfromtxt('messung1.txt',unpack=True)

U1=6.5
I1_leucht=2.5
U2=6
I2_leucht=2.4
U3=5.5
I3_leucht=2.3
U4=5.2
I4_leucht=2.2
U5=5
I5_leucht=2.1

I1=I1*1e-6
I2=I2*1e-6
I3=I3*1e-6
I4=I4*1e-6
I5=I5*1e-6
print(U)



plt.figure(1)
plt.plot(U,I1*1e6,'rx',label=r'$\mathrm{Messwerte \  für \ I_1}$')
plt.plot(U,I2*1e6,'gx',label=r'$\mathrm{Messwerte \  für \ I_2}$')
plt.plot(U,I3*1e6,'bx',label=r'$\mathrm{Messwerte \  für \ I_3}$')
plt.plot(U,I4*1e6,'cx',label=r'$\mathrm{Messwerte \  für \ I_4}$')
plt.plot(U,I5*1e6,'yx',label=r'$\mathrm{Messwerte \  für \ I_5}$')

#plt.plot([100,100],[0,100],'b--',label=r'$\mathrm{Theorie}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Anodenspannung \ U/V}$')
plt.ylabel(r'$\mathrm{Anodenstrom \ I/mA}$')
plt.savefig('plot1.pdf')


#Fit für 3/2
def V_23(V,b,a):
    #c=(4/9)*const.epsilon_0*np.sqrt(2*const.e/const.m_e)*b
    return(a*(V**b))

U_fit=U[:np.size(U)-8]
print(U_fit)
I1_fit=I1[:np.size(U)-8]
test=3/2
params_1, covariance_1 =curve_fit(V_23,U_fit,I1_fit)
errors_1 = np.sqrt(np.diag(covariance_1))
print(params_1)
print(errors_1)
x=np.linspace(0,190)
plt.figure(2)
plt.plot(U,I1*1e6,'cx',label=r'$\mathrm{Messwerte \ außerhalb \ des \ Raumladungsgebietes} $')
plt.plot(U_fit,I1_fit*1e6,'rx',label=r'$\mathrm{Messwerte \ innerhalb \ des\ Raumladungsgebietes }$')
plt.plot(x,V_23(x,*params_1)*1e6,'b-',label=r'$\mathrm{Ausgleichsfunktion}$')
plt.plot(U_fit,I1_fit*1e6,'rx')

#plt.plot(x,V_23(x,3/2,400),'b-',label=r'$Ausgleichsfunktion$')

#plt.plot([100,100],[0,100],'b--',label=r'$\mathrm{Theorie}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Anodenspannung \ U/V}$')
plt.ylabel(r'$\mathrm{Anodenstrom \ I/mA}$')
plt.savefig('plotfitt.pdf')




U_anlauf ,I_anlauf =np.genfromtxt('messung2.txt',unpack=True)

I_anlauf=I_anlauf*1e-9

r_i=1e6

U_korr=U_anlauf-I_anlauf*r_i

print('u_korr',U_korr)

m , b , r ,p ,std =stats.linregress(U_korr,np.log(I_anlauf))

def fehler_b(std_m,x):
    return(std_m*np.sqrt(np.mean(x**2)))
M=unp.uarray(m,std)
B=unp.uarray(b,fehler_b(std,U_korr))

x=np.linspace(0,1)
plt.figure(3)
plt.plot(U_korr,I_anlauf,'cx',label=r'$\mathrm{Messwerte}$')
plt.plot(x,np.exp(m*x+b),'b-',label=r'$\mathrm{Ausgleichsfunktion}$')
plt.yscale('log')
#plt.plot(x,V_23(x,3/2,400),'b-',label=r'$Ausgleichsfunktion$')
#plt.plot([100,100],[0,100],'b--',label=r'$\mathrm{Theorie}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Gegenspannung \ U_k/V}$')
plt.ylabel(r'$\mathrm{ln(I_A)}$')
plt.savefig('plotanlauf.pdf')


T=-const.e/(M*const.k)
print('Steigung m',M)

print(' b',B)
print(T)
#
# U1=6.5
# I1_leucht=2.5

f=0.32
o=5.7e-12
n=0.28
n_wl=unp.uarray(0.95,0.05)

T_w1=((U1*I1_leucht-n_wl)/(f*n*o))**(1/4)
T_w2=((U2*I2_leucht-n_wl)/(f*n*o))**(1/4)
T_w3=((U3*I3_leucht-n_wl)/(f*n*o))**(1/4)
T_w4=((U4*I4_leucht-n_wl)/(f*n*o))**(1/4)
T_w5=((U5*I5_leucht-n_wl)/(f*n*o))**(1/4)



print('T_w1',T_w1 )
print('T_w2',T_w2 )
print('T_w3',T_w3 )
print('T_w4',T_w4 )
print('T_w5',T_w5 )





#N_zu

I_s1=3.110
I_s2=2.440
I_s3=1.294
I_s4=0.721
I_s5=0.250


Austritt_1=-unp.log((I_s1/f)*(const.h**3)/(const.e*const.m_e*(const.k**2)*(T_w1**2)))*const.k*T_w1/const.e
Austritt_2=-unp.log((I_s2/f)*(const.h**3)/(const.e*const.m_e*(const.k**2)*(T_w2**2)))*const.k*T_w2/const.e
Austritt_3=-unp.log((I_s3/f)*(const.h**3)/(const.e*const.m_e*(const.k**2)*(T_w3**2)))*const.k*T_w3/const.e
Austritt_4=-unp.log((I_s4/f)*(const.h**3)/(const.e*const.m_e*(const.k**2)*(T_w4**2)))*const.k*T_w4/const.e
Austritt_5=-unp.log((I_s5/f)*(const.h**3)/(const.e*const.m_e*(const.k**2)*(T_w5**2)))*const.k*T_w5/const.e




#Austritt_Leistung=-unp.log((I_s1/f)*(const.h**3)/(const.e*const.m_e*(const.k**2)*(T_w1**2)))*const.k*T_w1/const.e

print('Austritt_1',Austritt_1)
print('Austritt_2',Austritt_2)
print('Austritt_3',Austritt_3)
print('Austritt_4',Austritt_4)
print('Austritt_5',Austritt_5)


print('Mittelwert normal',np.mean([Austritt_1,Austritt_2,Austritt_3,Austritt_4,Austritt_5]),'+-',np.std([noms(Austritt_1),noms(Austritt_2),noms(Austritt_3),noms(Austritt_4),noms(Austritt_5)]))
#print('Mittelwert np',np.mean([Austritt_2,Austritt_1]),'+-',np.std([noms(Austritt_2),noms(Austritt_1)]))
