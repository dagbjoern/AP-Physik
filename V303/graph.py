import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

x, y = np.genfromtxt('daten.txt', unpack=True)

x1, y1 ,y2=np.genfromtxt('1).txt', unpack=True)

plt.figure(1)
plt.xlim(0,110)
plt.plot(x,y, 'rx')
#def f(r , a):
#    return a/r**2
#params, pcov=curve_fit(f,x,y)
#sex=np.linspace(0,110)
#sey=np.linspace(0,50)
#plt.plot(sex, f(sex,*params), 'b-')
plt.savefig('plot.pdf')


plt.figure(2)
plt.plot((np.pi/180)*x1, y1, 'xr')
plt.plot((np.pi/180)*x1, y2, '+b')
xlin=np.linspace(0,2*np.pi)
plt.plot(xlin,(2/np.pi)*2*np.cos(xlin),'--k')
plt.savefig('plot1.pdf')
