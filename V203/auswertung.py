import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit


Tf , Tg , p1 =np.genfromtxt('1messung.txt',unpack=True)
T2 , p2 =np.genfromtxt('2messung.txt',unpack=True)
R=8.3144598

p1=p1*1e-3
Tg=Tg+273.15
T2=T2+273.15
p2=8.14+p2



m, b, r, p, stdsm=stats.linregress(Tg,np.log(p1))

def F(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d

params, covariance =curve_fit(F,T2,p2)


a = params[0]
b = params[1]
c = params[2]
d = params[3]


print('a =', params[0])
print('b =', params[1])
print('c =', params[2])
print('d =', params[3])




x1=np.linspace(20+273.15,100+273.15)
plt.figure(1)
plt.plot(Tg, np.log(p1),'rx')
plt.plot(x1,m*x1+b,'b-')
plt.savefig('plot1.pdf')


x2=np.linspace(20+273.15,200+273.15)
plt.figure(2)
plt.plot(T2, p2,'rx')
plt.plot(x2,F(x2,*params),'b-')
plt.savefig('plot2.pdf')

def p_t(T):
    return a*T**3+b*T**2+c*T+d

def p_t_abl(T):
    return 3*a*T**2+2*b*T+c

def L_max(T):
    return T*((R*T/(2*p_t(T)))+np.sqrt((R*T/(2*p_t(T)))**2-(0.9/p_t(T))))*p_t_abl(T)


def L_min(T):
    return T*((R*T/(2*p_t(T)))-np.sqrt((R*T/(2*p_t(T)))**2-(0.9/p_t(T))))*p_t_abl(T)



plt.figure(3)
plt.plot(T2,L_min(T2),'rx')
plt.plot(x2,L_min(x2))
plt.savefig('plot3.pdf')

plt.figure(4)
plt.plot(T2,L_max(T2),'rx')
plt.plot(x2,L_max(x2))
plt.savefig('plot4.pdf')


np.savetxt('1tabelle.txt',np.column_stack((Tg,p1)),fmt='%r',delimiter=' & ')
np.savetxt('2tabelle.txt',np.column_stack((T2,p2)),fmt='%r',delimiter=' & ')
