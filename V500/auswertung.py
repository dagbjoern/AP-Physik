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
I_all , U_O , U_G , U_T , U_B ,U_V  =np.genfromtxt('V500messa.txt',unpack=True)
U_b , I_b =np.genfromtxt('V500messb.txt',unpack=True)






I_all=I_all*10**(-12)

I_b=I_b*10**(-12)

I_allT=I_all[0:9]

U_T=U_T[0:9]


mB , bB , rB ,pB ,stdB =stats.linregress(U_B,np.sqrt(I_all))
mG , bG , rG ,pG ,stdG =stats.linregress(U_G,np.sqrt(I_all))
mO , bO , rO ,pO ,stdO =stats.linregress(U_O,np.sqrt(I_all))
mT , bT , rT ,pT ,stdT =stats.linregress(U_T,np.sqrt(I_allT))
mV , bV , rV ,pV ,stdV =stats.linregress(U_V,np.sqrt(I_all))



MB=unp.uarray(mB,stdB)
BB=unp.uarray(bB,fehler_b(stdB,U_B))
MG=unp.uarray(mG,stdG)
BG=unp.uarray(bG,fehler_b(stdG,U_G))
MO=unp.uarray(mO,stdO)
BO=unp.uarray(bO,fehler_b(stdO,U_O))
MT=unp.uarray(mT,stdT)
BT=unp.uarray(bT,fehler_b(stdT,U_T)) # könnte wegen der länge scjeiteren
MV=unp.uarray(mV,stdV)
BV=unp.uarray(bV,fehler_b(stdV,U_V))

M_ges=unp.uarray([noms(MO),noms(MG),noms(MT),noms(MB),noms(MV)],
[stds(MO),stds(MG),stds(MT),stds(MB),stds(MV)])

B_ges=unp.uarray([noms(BO),noms(BG),noms(BT),noms(BB),noms(BV)],
[stds(BO),stds(BG),stds(BT),stds(BB),stds(BV)])




def Ug(M,B):
    return(-B/M)

U_gO=Ug(MO,BO)
U_gG=Ug(MG,BG)
U_gT=Ug(MT,BT)
U_gB=Ug(MB,BB)
U_gV=Ug(MV,BV)


U_gGes=unp.uarray([noms(U_gO),noms(U_gG),noms(U_gT),noms(U_gB),noms(U_gV)],
[stds(U_gO),stds(U_gG),stds(U_gT),stds(U_gB),stds(U_gV)])


x_U=np.linspace(0.9,1)

plt.figure(1)
plt.plot(U_B,np.sqrt(I_all)*1e6,'bx',label=r'$\mathrm{Messwerte \  für \ Blau}$')
plt.plot(x_U,(mB*x_U+bB)*1e6,'b-',label=r'$\mathrm{Ausgleichsfunkion \ für\ Blau  }$')
plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\mathrm{Spannung \ U/V}$')
plt.ylabel(r'$\mathrm{\sqrt{I} \cdot 10^{-6} }$')
plt.savefig('plotV500aBlau.pdf')
#
#


x_U=np.linspace(0.25,0.5)

plt.figure(2)
plt.plot(U_O,np.sqrt(I_all)*1e6,'yx',label=r'$\mathrm{Messwerte \  für \ Gelb}$')
plt.plot(x_U,(mO*x_U+bO)*1e6,'y-',label=r'$\mathrm{Ausgleichsfunkion \ für\ Gelb  }$')
plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\mathrm{Spannung \ U/V}$')
plt.ylabel(r'$\mathrm{\sqrt{I} \cdot 10^{-6} }$')
plt.savefig('plotV500aGelb.pdf')


x_U=np.linspace(0.4,0.6)

plt.figure(3)
plt.plot(U_G,np.sqrt(I_all)*1e6,'gx',label=r'$\mathrm{Messwerte \  für \ Grün}$')
plt.plot(x_U,(mG*x_U+bG)*1e6,'g-',label=r'$\mathrm{Ausgleichsfunkion \ für\ Grün  }$')
plt.grid()
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Spannung \ U/V}$')
plt.ylabel(r'$\mathrm{\sqrt{I} \cdot 10^{-6} }$')
plt.savefig('plotV500aGrun.pdf')


x_U=np.linspace(0.9,1.1)

plt.figure(4)
plt.plot(U_V,np.sqrt(I_all)*1e6,'mx',label=r'$\mathrm{Messwerte \  für \ Violett}$')
plt.plot(x_U,(mV*x_U+bV)*1e6,'m-',label=r'$\mathrm{Ausgleichsfunkion \ für\ Violett  }$')
plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\mathrm{Spannung \ U/V}$')
plt.ylabel(r'$\mathrm{\sqrt{I} \cdot 10^{-6} }$')
plt.savefig('plotV500aviolett.pdf')


x_U=np.linspace(0,0.9)

plt.figure(5)
plt.plot(U_T,np.sqrt(I_allT)*1e6,'cx',label=r'$\mathrm{Messwerte \  für \ Blau-Grün}$')
plt.plot(x_U,(mT*x_U+bT)*1e6,'c-',label=r'$\mathrm{Ausgleichsfunkion \ für\ Blau-Grün  }$')
plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\mathrm{Spannung \ U/V}$')
plt.ylabel(r'$\mathrm{\sqrt{I} \cdot 10^{-6} }$')
plt.savefig('plotV500agrunblau.pdf')











lamO=578e-9
lamG=546e-9
lamT=492e-9
lamB=435e-9
lamV=405e-9

lam_ges=np.array([lamO,lamG,lamT,lamB,lamV])




v_ges=const.c/lam_ges


#np.savetxt('tabelleU_g.txt',np.column_stack((lam_ges,v_ges,M_ges,B_ges,U_gGes)),fmt='%r',delimiter=' & ')






m1 , b1 , r1 ,p1 ,std1 =stats.linregress(v_ges,noms(U_gGes))

x_v=np.linspace(0,9e14)

plt.figure(6)
#plt.plot(v_ges,noms(U_gGes),'rx',label=r'$\mathrm{Messwerte }$')
plt.errorbar(v_ges*1e-12,noms(U_gGes),yerr=stds(U_gGes), fmt='rx',label=r'$\mathrm{Messwerte }$')
plt.plot(x_v*1e-12,m1*x_v+b1,'--',label=r'$\mathrm{Ausgleichsfunkion}$')
plt.grid()
plt.legend(loc='best')
plt.xlabel(r'$\mathrm{Frequenz \ v/THz}$')
plt.ylabel(r'$\mathrm{Grenzspannung\ U_g/V}$')

plt.savefig('plotV500b.pdf')


M1=unp.uarray(m1,std1)
B1=unp.uarray(b1,fehler_b(std1,v_ges))


print('M1',M1)
print('h/e0',const.h/const.e)
print('Ak in eV',B1)



plt.figure(7)
plt.plot(U_b,I_b*1e12,'rx',label=r'$\mathrm{Messwerte}$')
plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\mathrm{Grenzspannung\ U_g/V}$')
plt.ylabel(r'$\mathrm{Photostrom\ I/V}$')
plt.savefig('plotV500c.pdf')


print(U_b.size)
U_b_red=U_b[30-13:30+12]
I_b_red=I_b[30-13:30+12]
print(U_b_red)
#um null herrum
plt.figure(8)
plt.plot(U_b_red,I_b_red*1e12,'rx',label=r'$\mathrm{Messwerte}$')
plt.legend(loc='best')
plt.grid()
plt.xlabel(r'$\mathrm{Grenzspannung\ U_g/V}$')
plt.ylabel(r'$\mathrm{Photostrom\ I/pA}$')
plt.savefig('plotV500c0.pdf')

#--------------------------
#
# I_b   , U_b
#
#
#
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
