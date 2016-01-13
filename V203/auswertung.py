import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit


Tf , Tg , p1 =np.genfromtxt('1messung.txt',unpack=True)
T2 , p2 =np.genfromtxt('2messung.txt',unpack=True)


p1=p1*1e-3
Tg=Tg
T2=T2
p2=8.14+p2



m, b, r, p, stdsm=stats.linregress(Tg,np.log(p1))

def F(x,a,b,c,d):
    return a*x**3+b*x**2+c*x+d

params, covariance =curve_fit(F,T2,p2)

print('a =', params[0])
print('b =', params[1])
print('c =', params[2])
print('d =', params[3])




x1=np.linspace(20,100)
plt.figure(1)
plt.plot(Tg, np.log(p1),'rx')
plt.plot(x1,m*x1+b,'b-')
plt.savefig('plot1.pdf')


x2=np.linspace(20,200)
plt.figure(2)
plt.plot(T2, p2,'rx')
plt.plot(x2,F(x2,*params),'b-')
plt.savefig('plot2.pdf')

np.savetxt('1tabelle.txt',np.column_stack((Tg,p1)),fmt='%r',delimiter=' & ')
np.savetxt('2tabelle.txt',np.column_stack((T2,p2)),fmt='%r',delimiter=' & ')
