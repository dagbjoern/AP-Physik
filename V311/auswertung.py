import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const

#formeln

#Dicke der Platte
def d(I,B,U,A_h):
    d=unp.uarray(np.mean((A_h*I*B)/(U)),np.std((A_h*I*B)/(U)))
    return d

#Wiederstand
def URI(I,U):
    return U/I


#ladungsträger pro Volumen
def n(I,d,B,U):
    n=(B*I)/(U*const.e*d)
    n=np.mean(n)
    return n


#Ladungsträger pro Atom
def z(n,Vm):
    return n*Vm/const.N_A
#mittlere Flugzeit
def t(n,R,Q,L):
    t=(2*const.m_e*L)/(R*(const.e**2)*n*Q)
    return t

#mittlere driftgeschwindigkeit für j =1
def v(n):
    v=-1/(const.e*n*(10**(-6)))
    return v

#fermie Energie
def Ef(n):
    return((const.h**2)/(2*const.m_e))*((3*n)/(8*np.pi))**(2/3)

#mittlere freie Wellenlänge
def l(t,n):
 return t*((2*Ef(n))/const.m_e)**(1/2)

#Beweglichkeit
def mu(t):
    return -(const.e*t)/(2*const.m_e)

#Totalgeschwindigkeit
def v_total(n):
 return (2*Ef(n)/const.m_e)**(1/2)

#hysterese:
I_hys, B_hys=np.genfromtxt('hysterese.txt',unpack=True)
B=np.max(B_hys)
print(np.max(B_hys))

plt.figure(1)
plt.plot(I_hys,B_hys,'rx',label=r'$Messwerte$')
plt.legend(loc='best')
plt.xlabel(r'$Strom I/A$')
plt.ylabel(r'$Flussdichte B/T $')
plt.savefig('hysterese.pdf')


#zink
I_z, U_z =np.genfromtxt('messungRzink.txt',unpack=True)
I_hall_z, U_hall_z=np.genfromtxt('hallzink.txt',unpack=True)
b_z=0.024#breite
l_z=0.034#länge
print('länge',l_z)
A_hz=6.4*10**(-11)  #hallkonstante
Vm_z=9.16*10**(-6) #molares volumen
U_z=U_z*10**(-3)
U_hall_z=U_hall_z*10**(-3)
R_z=URI(I_z,U_z)
print('\nWerte Für Zink','R_Z=',R_z)
R_z=unp.uarray(np.mean(R_z),np.std(R_z))
print('Wiederstand Mittelwert',R_z)


d_z=d(I_hall_z,B,U_hall_z,A_hz)
print('Dicke d',d_z)
Q_z=d_z*b_z

#def n(I,d,B,U):
n_z=n(I_hall_z,d_z,B,U_hall_z)
print('n von zink',n_z)

#def z(n,Vm):
print('ladungsträger pro Atom',z(n_z,Vm_z))

#def t(n,R,Q,L):
t_z=t(n_z,R_z,Q_z,l_z)
print('tau von zink',t_z)

#def v(n):
v_z=v(n_z)
print('mittlere Driftgeschwindigkeit',v_z)

#def Ef(n):
Ef_z=Ef(n_z)
print('fermie Energie',Ef_z)

#def v_total(n):
v_total_z=v_total(n_z)
print('v Total',v_total_z)

#def l(t,n):
l_z=l(t_z,n_z)
print('mittlere Wellenlänge',l_z)

#def mu(t):
mu_z=mu(t_z)
print('Beweglichkeit mu',mu_z)

print('Querschnitt',Q_z)

#wolfram
I_w, U_w =np.genfromtxt('messungRwolfram.txt',unpack=True)
I_hall_w, U_hall_w=np.genfromtxt('hallwolfram.txt',unpack=True)
b_w=0.026#breite
l_w=0.04#länge
print('länge',l_w)
Vm_w=9.47*10**(-6)#molares Volumen Wolfram
A_hw=11.8*10**(-11)
U_w=U_w*10**(-3)
U_hall_w=U_hall_w*10**(-3)
R_w=URI(I_w,U_w)
print('\nWerte für wolfram','R_w=',R_w)
R_w=unp.uarray(np.mean(R_w),np.std(R_w))
print('gemittelter wiederstand ',R_w)

d_w=d(I_hall_w,B,U_hall_w,A_hw)
print('Dicke d',d_w)
Q_w=d_w*b_w

#def n(I,d,B,U):
n_w=n(I_hall_w,d_w,B,U_hall_w)
print('n von wolfram',n_w)

#def z(n,Vm):
print('ladungsträger pro Atom',z(n_w,Vm_w))

#def t(n,R,Q,L):
t_w=t(n_w,R_w,Q_w,l_w)
print('tau von wolfram',t_w)

#def v(n):
v_w=v(n_w)
print('mittlere Driftgeschwindigkeit',v_w)

#def Ef(n):
Ef_w=Ef(n_w)
print('fermie Energie',Ef_w)

#def v_total(n):
v_total_w=v_total(n_w)
print('v Total',v_total_w)

#def l(t,n):
l_w=l(t_w,n_w)
print('mittlere Wellenlänge',l_w)

#def mu(t):
mu_w=mu(t_w)
print('Beweglichkeit mu',mu_w)

print('Querschnitt',Q_w)
