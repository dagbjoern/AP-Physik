import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit
import scipy.constants as const


def fehler_b(std_m,x):
    return(std_m*np.sqrt(np.mean(x**2)))


# #Kennlinen
Ud1, Ud2 , Ud3 , Ud4 , Ud5 =np.genfromtxt('V501mess.txt',unpack=True)
I1, I2 , I3 , I4 , I5 =np.genfromtxt('V502mess.txt',unpack=True)

käst=np.array([0,1,2,3,4,5,6,7,8])


I4=I4[:8]
D_sonst=käst*2.54e-2*(1/4)

D_U1=D_sonst
D_I1=unp.uarray(D_sonst,0.002)

D_U2=D_sonst
D_I2=unp.uarray(D_sonst,0.002)

D_U3=D_sonst
D_I3=unp.uarray(D_sonst,0.002)

D_U4=käst*2.54e-2*(1/8)
D_I4=unp.uarray(käst[:8]*2.54e-2*(1/4),0.002)

D_U5=D_U4
D_I5=unp.uarray(käst*2.54e-2*(1/8),0.002)

print(2.54e-2)
print(D_I4)

m1 , b1 , r1 ,p1 ,std1 =stats.linregress(Ud1,D_U1)
m2 , b2 , r2 ,p2 ,std2 =stats.linregress(Ud2,D_U2)
m3 , b3 , r3 ,p3 ,std3 =stats.linregress(Ud3,D_U3)
m4 , b4 , r4 ,p4 ,std4 =stats.linregress(Ud4,D_U4)
m5 , b5 , r5 ,p5 ,std5 =stats.linregress(Ud5,D_U5)


M1=unp.uarray(m1,std1)
B1=unp.uarray(b1,fehler_b(std1,Ud1))
M2=unp.uarray(m2,std2)
B2=unp.uarray(b2,fehler_b(std2,Ud2))
M3=unp.uarray(m3,std3)
B3=unp.uarray(b3,fehler_b(std3,Ud3))
M4=unp.uarray(m4,std4)
B4=unp.uarray(b4,fehler_b(std4,Ud4))
M5=unp.uarray(m5,std5)
B5=unp.uarray(b5,fehler_b(std5,Ud5))





print('M1',M1)
print('M2',M2)
print('M3',M3)
print('M4',M4)
print('M5',M5)


Ub1=180
Ub2=250
Ub3=300
Ub4=400
Ub5=500
#print(D_U1)

x_U=np.linspace(-40,20)
plt.figure(1)
plt.plot(Ud1,D_U1,'rx',label=r'$\mathrm{Messwerte \  für \ U_b=180\mathrm{V}}$')
plt.errorbar(Ud1,D_U1,xerr=0,yerr=0.002, fmt='rx')
plt.plot(Ud3,D_U3,'yx',label=r'$\mathrm{Messwerte \  für \ U_b=300\mathrm{V}}$')
plt.errorbar(Ud3,D_U3,xerr=0,yerr=0.002, fmt='yx')

plt.plot(Ud5,D_U5,'bx',label=r'$\mathrm{Messwerte \  für \ U_b=500\mathrm{V}}$')
plt.errorbar(Ud5,D_U5,xerr=0,yerr=0.002, fmt='bx')


plt.plot(x_U,m1*x_U+b1,'r-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=180V  }$')

plt.plot(x_U,m3*x_U+b3,'y-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=300V  }$')

plt.plot(x_U,m5*x_U+b5,'b-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=500V  }$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Ablenkspannung \ U_d/V}$')
plt.ylabel(r'$\mathrm{Verschiebung \ D/m}$')
plt.savefig('plotV5011.pdf')

plt.figure(2)

plt.plot(Ud2,D_U2,'cx',label=r'$\mathrm{Messwerte \  für \ U_b=250\mathrm{V}}$')
plt.errorbar(Ud2,D_U2,xerr=0,yerr=0.002, fmt='cx')

plt.plot(Ud4,D_U4,'gx',label=r'$\mathrm{Messwerte \  für \ U_b=400\mathrm{V}}$')
plt.errorbar(Ud4,D_U4,xerr=0,yerr=0.002, fmt='gx')

plt.plot(x_U,m2*x_U+b2,'c-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=250V  }$')
plt.plot(x_U,m4*x_U+b4,'g-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=400V  }$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Ablenkspannung \ U_d/V}$')
plt.ylabel(r'$\mathrm{Verschiebung \ D/m}$')
plt.savefig('plotV5012.pdf')

np.savetxt('1tabelle.txt',np.column_stack((D_U1,Ud1,Ud2,Ud3)),fmt='%r',delimiter=' & ')


Mges=unp.uarray([m1,m2,m3,m4,m5],[stds(M1),stds(M2),stds(M3),stds(M4),stds(M5)])
Ubges=np.array([Ub1,Ub2,Ub3,Ub4,Ub5])


print(noms(Mges))
ma , ba , ra ,pa ,stda =stats.linregress(1/Ubges,noms(Mges))

x_ub=np.linspace(0,0.006)
plt.figure(3)
plt.plot(1/Ubges,noms(Mges),'rx',label=r'$\mathrm{Messwerte}$')
plt.errorbar(1/Ubges,noms(Mges),xerr=0,yerr=stds(Mges), fmt='rx')
plt.plot(x_ub,ma*x_ub+ba,'c-',label=r'$\mathrm{Ausgleichsfunkion}$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{U_b^{-1}/V^{-1}}$')
plt.ylabel(r'$D\,U_d^{-1}\mathrm{/m V^{-1} }$')
plt.savefig('plotV501a.pdf')


Ma=unp.uarray(ma,stda)
Ba=unp.uarray(ba,fehler_b(stda,1/Ubges))


print('ma',Ma)
print('ba',Ba)

L=143e-3
p=19e-3
d=3.8e-3
a=(p*L)/(2*d)

D_max=unp.uarray(2.54e-2*(3/8),2e-3)
U_bmax=300
U_max=300*D_max/a
print('a ',a)
print('abweich',(a-ma)/a)
print('U_max sinus',U_max)




#############-B-Feld   #########



U_b11=250
U_b22=300
U_b33=350
U_b44=400
U_b55=500


U_bges2=np.array([U_b11,U_b22,U_b33,U_b44,U_b55])


def f(D):
    L=143e-3
    return(D/(L**2 +  D**2))



def B(I):
    n=20
    r=0.282
    return(const.mu_0*8*(125**(-0.5))*n*I/r)



m11 , b11 , r11 ,p11 ,std11 =stats.linregress(B(I1),noms(f(D_I1)))
m22 , b22 , r22 ,p22 ,std22 =stats.linregress(B(I2),noms(f(D_I2)))
m33 , b33 , r33 ,p33 ,std33 =stats.linregress(B(I3),noms(f(D_I3)))
m44 , b44 , r44 ,p44 ,std44 =stats.linregress(B(I4),noms(f(D_I4)))
m55 , b55 , r55 ,p55 ,std55 =stats.linregress(B(I5),noms(f(D_I5)))




M11=unp.uarray(m11,std11)
B11=unp.uarray(b11,fehler_b(std11,B(I1)))
M22=unp.uarray(m22,std22)
B22=unp.uarray(b22,fehler_b(std22,B(I2)))
M33=unp.uarray(m33,std33)
B33=unp.uarray(b33,fehler_b(std33,B(I3)))
M44=unp.uarray(m44,std44)
B44=unp.uarray(b44,fehler_b(std44,B(I4)))
M55=unp.uarray(m55,std55)
B55=unp.uarray(b55,fehler_b(std55,B(I5)))


Mges2=unp.uarray([m11,m22,m33,m44,m55],[stds(M11),stds(M22),stds(M33),stds(M44),stds(M55)])
Bges2=unp.uarray([b11,b22,b33,b44,b55],[stds(B11),stds(B22),stds(B33),stds(B44),stds(B55)])




print('M11',M11)
print('M22',M22)
print('M33',M33)
print('M44',M44)
print('M55',M55)






x_I=np.linspace(0,0.0002)

plt.figure(4)
plt.plot(B(I1),noms(f(D_I1)),'rx',label=r'$\mathrm{Messwerte \  für \ U_b=250\mathrm{V}}$')
plt.errorbar(B(I1),noms(f(D_I1)),xerr=0,yerr=stds(f(D_I1)), fmt='rx')

plt.plot(B(I3),noms(f(D_I3)),'bx',label=r'$\mathrm{Messwerte \  für \ U_b=350\mathrm{V}}$')
plt.errorbar(B(I3),noms(f(D_I3)),xerr=0,yerr=stds(f(D_I3)), fmt='bx')

plt.plot(B(I5),noms(f(D_I5)),'mx',label=r'$\mathrm{Messwerte \  für \ U_b=500\mathrm{V}}$')
plt.errorbar(B(I5),noms(f(D_I5)),xerr=0,yerr=stds(f(D_I5)), fmt='mx')

plt.plot(x_I,m11*x_I+b11,'r-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U-B=250V  }$')
plt.plot(x_I,m33*x_I+b33,'b-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_B=350V  }$')
plt.plot(x_I,m55*x_I+b55,'m-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_B=500V  }$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Feldstärke \  B/T}$')
plt.ylabel(r'$\mathrm{ D/(L^2+D^2)\ /m}$')
plt.savefig('plotV5021.pdf')



plt.figure(5)
plt.plot(B(I2),noms(f(D_I2)),'cx',label=r'$\mathrm{Messwerte \  für \ U_B=300\mathrm{V}}$')
plt.errorbar(B(I2),noms(f(D_I2)),xerr=0,yerr=stds(f(D_I2)), fmt='cx')

plt.plot(B(I4),noms(f(D_I4)),'gx',label=r'$\mathrm{Messwerte \  für \ U_B=400\mathrm{V}}$')
plt.errorbar(B(I4),noms(f(D_I4)),xerr=0,yerr=stds(f(D_I4)), fmt='gx')

plt.plot(x_I,m22*x_I+b22,'c-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=300V  }$')
plt.plot(x_I,m44*x_I+b44,'g-',label=r'$\mathrm{Ausgleichsfunkion \ für\ U_b=400V  }$')
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Flussdichte \  B/T}$')
plt.ylabel(r'$\mathrm{ D/(L^2+D^2)\ /m}$')
plt.savefig('plotV5022.pdf')

def e0_m0(U_B,m):
    return(8*U_B*(m**2))

print('e0/m0',e0_m0(U_bges2,Mges2))
print('e_0/m_0 mittel', np.mean(e0_m0(U_bges2,Mges2)))


np.savetxt('tabelleB.txt',np.column_stack((U_bges2,Mges2,Bges2,e0_m0(U_bges2,Mges2)  )),fmt='%r',delimiter=' & ')


I_erde=0.16
B_hor=B(I_erde)

print('B_hor',B_hor)
phi=75*2*np.pi/360

B_erde=B_hor/np.cos(phi)

print('B_erde',B_erde)





#
#
# #Fit für 3/2
# def V_23(V,b,a):
#     #c=(4/9)*const.epsilon_0*np.sqrt(2*const.e/const.m_e)*b
#     return(a*(V**b))
#
# U_fit=U[:np.size(U)-8]
# print(U_fit)
# I1_fit=I1[:np.size(U)-8]
# test=3/2
# params_1, covariance_1 =curve_fit(V_23,U_fit,I1_fit)
# errors_1 = np.sqrt(np.diag(covariance_1))
#
#
#
#
#
#
#
