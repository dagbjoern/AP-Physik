import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit

r31, r41  ,c21 = np.genfromtxt('werte1.txt',unpack=True)
r311, r411 ,r211= np.genfromtxt('werte11.txt',unpack=True)
r33, r43  ,c23 = np.genfromtxt('werte3.txt',unpack=True)
r28,r38,r48,  c28    = np.genfromtxt('werte8.txt',unpack=True)
r314,r414,r214= np.genfromtxt('werte14.txt',unpack=True)
r217,r317,r417,l217= np.genfromtxt('werte17.txt',unpack=True)
r217m,r317m,r417m,c417m= np.genfromtxt('werte17Max.txt',unpack=True)
f, u = np.genfromtxt('werteTT.txt',unpack=True)


def linfehler(r3,r4):
    rx=unp.uarray((r3/r4),(r3/r4)*0.005)
    return rx


R3_41=linfehler(r31,r41)
R31=unp.uarray(r31,r31*0.03)
R41=unp.uarray(r41,r41*0.03)
C21=unp.uarray(c21,c21*0.002)
R3_411=linfehler(r311,r411)
R311=unp.uarray(r311,r311*0.03)
R411=unp.uarray(r411,r411*0.03)
R211=unp.uarray(r211,r211*0.03)
R3_43=linfehler(r33,r43)
R33=unp.uarray(r33,r33*0.03)
R43=unp.uarray(r43,r43*0.03)
C23=unp.uarray(c23,c23*0.002)
R28=unp.uarray(r28,r28*0.03)
R3_48=linfehler(r38,r48)
R38=unp.uarray(r38,r38*0.03)
R48=unp.uarray(r48,r48*0.03)
C28=unp.uarray(c28,c28*0.002)
R3_414=linfehler(r314,r414)
R314=unp.uarray(r314,r314*0.03)
R414=unp.uarray(r414,r414*0.03)
R214=unp.uarray(r214,r214*0.03)
R217=unp.uarray(r217,r217*0.03)
R3_417=linfehler(r317,r417)
print(R3_417)
R317=unp.uarray(r317,r317*0.03)
R417=unp.uarray(r417,r417*0.03)
L217=unp.uarray(l217,l217*0.002)
R217m=unp.uarray(r217m,r217m*0.03)
R317m=unp.uarray(r317m,r317m*0.03)
R417m=unp.uarray(r417m,r417m*0.03)
C417m=unp.uarray(c417m,c417m*0.002)


def wheat(R2,R3_4):
    Rx=R2*(R3_4)
    return Rx


def kapaziR(R2,R3_4):
    Rx=R2*(R3_4)
    return Rx

def kapaziC(C2,R3_4):
    Cx=C2*(1/R3_4)
    return Cx

def indukR(R2,R3_4):
    Rx=R2*(R3_4)
    return Rx

def indukL(L2,R3_4):
    Lx=L2*(R3_4)
    return Lx

def maxR(R2,R3,R4):
    Rx=(R2*R3)/R4
    return Rx

def maxL(R2,R3,C4):
    Lx=R2*R3*C4
    return Lx

#a)
rx14=wheat(R214,R3_414)
Rx14=rx14.mean()
print('14=',Rx14)
rx11=wheat(R211,R3_411)
Rx11=rx11.mean()
print('11=',Rx11)

#b)
cx3=kapaziC(C23,R3_43)
Cx3=cx3.mean()
print('3=',Cx3)


cx1=kapaziC(C21,R3_41)
Cx1=cx1.mean()
print('1=',Cx1)

rx8=kapaziR(R28,R3_48)
cx8=kapaziC(C28,R3_48)

Rx8=rx8.mean()
Cx8=cx8.mean()
print('C8=',Cx8)
print('R8=',Rx8)
#c)
rx17=indukR(R217,R3_417)
Rx17=rx17.mean()
print('17R',Rx17)

lx17=indukL(L217,R3_417)
Lx17=lx17.mean()
print('17L',Lx17)

#d)
rx17m=maxR(R217m,R317m,R417m)
Rx17m=rx17m.mean()
print('Rx17m',Rx17m)
lx17m=maxL(R317m,R417m,C417m)
Lx17m=lx17m.mean()
print('Lx17m',Lx17m)

RTT=unp.uarray(1000,1000*0.03)
CTT=Cx3
print(Cx3)
w_0=1/(RTT*CTT)
print('W_o',w_0)
print('v_o',w_0/(2*np.pi))


i=0
while u[i] != u.min():
    i=i+1

v_0=f[i]

print(2*np.pi*v_0)



def fk(O):

    return np.sqrt(((O**2-1)**2)/((1-O**2)**2 + 16*O**2))

def TT_B(w,w_0):
    O=w/w_0
    return np.sqrt(noms(((O**2-1)**2)/((1-O**2)**2 + 16*O**2)))

u_s=1.65
x=np.linspace(20*np.pi*2,6000*np.pi*2,2000)

plt.figure(1)
plt.plot(noms(x/(w_0)),TT_B(x,w_0),'b-',label=r'$Theoriekurve$')
plt.plot(noms(f/v_0),u/u_s,'rx',label=r'$Messwerte$')
plt.xlabel(r'$\Omega$')
plt.ylabel(r'$U_{Br}/U_S$')
plt.legend(loc='best')
plt.grid(True)
plt.savefig('plot2.pdf')

U_2=0.005/fk(2)
print('U_2',U_2)
print('k',U_2/u_s)
#np.savetxt('tabelle14.txt',np.column_stack((r314,r414,R3_414,R214,rx14)),fmt='%r',delimiter=' & ')
#np.savetxt('tabelle11.txt',np.column_stack((r311,r411,R3_411,R211,rx11)),fmt='%r',delimiter=' & ')
#np.savetxt('tabelle3.txt',np.column_stack((r33,r43,R3_43,C23,cx3)),fmt='%r',delimiter=' & ')
np.savetxt('tabelle1.txt',np.column_stack((r31,r41,R3_41,C21,cx1)),fmt='%r',delimiter=' & ')
#np.savetxt('tabelle8.txt',np.column_stack((r38,r48,R3_48,R28,C28,rx8,cx8)),fmt='%r',delimiter=' & ')
#np.savetxt('tabelle17.txt',np.column_stack((r317,r417,R3_417,R217,L217,rx17,lx17)),fmt='%r',delimiter=' & ')
#np.savetxt('tabelle17m.txt',np.column_stack((R317m,R417m,R217m,C417m,rx17m ,lx17m)),fmt='%r',delimiter=' & ')

print(lx17)
