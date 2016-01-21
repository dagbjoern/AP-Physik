import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.signal as sc

def D(T):
    Tb=[]
    Tb=np.array(Tb)
    x=0
    while x < (len(T)-1):
        x=x+4
        Tb=np.append(Tb,T[x])
    return Tb

def hoch(t,T,hoch1):
  thoch=[]
  thoch=np.array(thoch)
  i=0
  y=0
  while y < len(hoch1):
    while T[i] != hoch1[y]:
          i=i+1
    thoch=np.append(thoch,t[i])
    y=y+1
  return thoch

def tief(t,T,tief1):
    ttief=[]
    ttief=np.array(ttief)
    i=len(T)-1
    y=len(tief1)-1
    while y >= 0:
        while T[i] != tief1[y]:
            i=i-1
            #print('2schleife')
        ttief=np.append(ttief,t[i])
        y=y-1
        #print('1schleife')
    return np.sort(ttief, axis=None)

def k(p,c,Thnah,Ttnah,thnah,ttnah,Thfern,Ttfern,thfern,ttfern):
 k=(p*c*(0.03)**2)/(2*Phasenver(thnah,ttnah,thfern,ttfern)*np.log(Ampver(Thnah,Ttnah)/Ampver(Thfern,Ttfern)))
 return k

def Ampver(Th,Tt):
     return np.average((Th-Tt)/2)

def Phasenver(thnah,ttnah,thfern,ttfern):
    dt1=np.average(thfern-thnah)
    dt2=np.average(ttfern-ttnah)
    dt=np.average([dt1,dt2])
    return dt

c_mes=385.
c_alu=830.
c_stahl=400.

p_mes=8520.
p_alu=2800.
p_stahl=8000.

t1,T11,T12,T13,T14,T15,T16,T17,T18=np.genfromtxt('messwerte1.txt',unpack=True)
t2,T21,T22,T23,T24,T25,T26,T27,T28=np.genfromtxt('messwerte2.txt',unpack=True)
t3,T31,T32,T33,T34,T35,T36,T37,T38=np.genfromtxt('messwerte3.txt',unpack=True)

A2Th , A2Tt= np.genfromtxt('extrema1.txt',unpack=True)
A5Th , A5Tt=np.genfromtxt('extremaAlu1.txt',unpack=True)
A6Th , A6Tt=np.genfromtxt('extremaAlu2.txt',unpack=True)
A7Th , A7Tt ,A8Th, A8Tt =np.genfromtxt('extremaEdel7u8.txt',unpack=True)
A1th=np.array([76.00,143.00,221.50,300.50,377.50,457.50,536.50,616.00,696.00,779.50])
A1Th=np.array([32.04,35.85,38.83,41.11,42.93,44.50,45.80,46.86,47.80,48.83])
A1tt=np.array([87.00,170.50,252.00,333.50,413.50,493.50,573.50,653.00,734.50,814.00])
A1Tt=np.array([32.03,35.60,38.36,40.43,42.14,43.61,44.78,45.76,46.66,47.66])



A2th=hoch(t2,T22,A2Th)
A2tt=tief(t2,T22,A2Tt)
A5th=hoch(t2,T25,A5Th)
A5tt=tief(t2,T25,A5Tt)
A6th=hoch(t2,T26,A6Th)
A6tt=tief(t2,T26,A6Tt)
A7th=hoch(t3,T37,A7Th)
A7tt=tief(t3,T37,A7Tt)
A8th=hoch(t3,T38,A8Th)
A8tt=tief(t3,T38,A8Tt)

#def k(p,c,Thnah,Ttnah,thnah,ttnah,Thfern,Ttf..
#def Phasenver(thnah,ttnach,thfern,ttfern)
#def Ampver(Th,Tt):
print('messing=',k(p_mes,c_mes,A2Th,A2Tt,A2th,A2tt,A1Th,A1Tt,A1th,A1tt))
print('alu=',k(p_alu,c_alu,A6Th,A6Tt,A6th,A6tt,A5Th,A5Tt,A5th,A5tt))
print('stahl=',k(p_stahl,c_stahl,A7Th,A7Tt,A7th,A7tt,A8Th,A8Tt,A8th,A8tt))

#print(Phasenver(A2th,A2tt,A1th,A1tt))
# dt1=np.average(thfern-thnah)
# dt2=np.average(ttfern-ttnah)
print(A2Th)
print(A2Tt)
print(Ampver(A2Th,A2Tt))

plt.figure(1)#T1T4
plt.plot(t1,T11+273.15,'b-',label=r'$Temperatur\ am\ Punkt\ T1 $')
plt.plot(t1,T14+273.15,'r-',label=r'$Temperatur\ am\ Punkt\ T4 $')
plt.xlabel(r'$Zeit\ \frac{t}{s}$')
plt.ylabel(r'$Temperatur\ \frac{T}{K} $')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plotT1T4.pdf')

plt.figure(2)#T5T8
plt.plot(t1,T15+273.15,'b-',label=r'$Temperatur\ am\ Punkt\ T5 $')
plt.plot(t1,T18+273.15,'r-',label=r'$Temperatur\ am\ Punkt\ T8 $')
plt.xlabel(r'$Zeit\ \frac{t}{s}$')
plt.ylabel(r'$Temperatur\ \frac{T}{K} $')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plotT5T8.pdf')

plt.figure(3)#T7-T8
plt.plot(t1,T17-T18,'b-',label=r'$Temperaturdifferents\ T7-T8 $')
plt.xlabel(r'$Zeit\ \frac{t}{s}$')
plt.ylabel(r'$Temperatur\ \frac{T}{K} $')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plotT7-T8.pdf')





plt.figure(4)#T2-T1
plt.plot(t1,T12-T11,'r-',label=r'$Temperaturdifferents\ T2-T1 $')
plt.xlabel(r'$Zeit\ \frac{t}{s}$')
plt.ylabel(r'$Temperatur\ \frac{T}{K} $')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plotT2-T1.pdf')

plt.figure(5)#T1T2
plt.plot(D(t2),D(T21+273.15),'b-',label=r'$Temperatur\ am\ Punkt\ T1 $')
plt.plot(t2,T22+273.15,'r-',label=r'$Temperatur\ am\ Punkt\ T2 $')
plt.xlabel(r'$Zeit\ \frac{t}{s}$')
plt.ylabel(r'$Temperatur\ \frac{T}{K} $')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plotT1T2.pdf')

plt.figure(6)#T7T8
plt.plot(t3,T37+273.15,'b-',label=r'$Temperatur\ am\ Punkt\ T7 $')
plt.plot(t3,T38+273.15,'r-',label=r'$Temperatur\ am\ Punkt\ T8 $')
plt.xlabel(r'$Zeit\ \frac{t}{s}$')
plt.ylabel(r'$Temperatur\ \frac{T}{K} $')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plotT7T8.pdf')

#ausgeben der Temperatur nach 700 s
print('Temperatur nach 700s für T1:',T11[700/0.2])
print('Temperatur nach 700s für T5:',T15[700/0.2])
print('Temperatur nach 700s für T4:',T14[700/0.2])
print('Temperatur nach 700s für T8:',T18[700/0.2])






#T2max=argrelextrema(T22, np.greater)

#T1min=argrelextrema(T22, np.less)
#T2min=argrelextrema(T22, np.less)
