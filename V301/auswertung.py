
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

x1 , y1 = np.genfromtxt('b).txt',unpack=True)
x2 , y2 = np.genfromtxt('d)recht.txt',unpack=True)
x3 , y3 = np.genfromtxt('d)sin.txt',unpack=True)
x4 , y4 = np.genfromtxt('c).txt',unpack=True)



def FehlerB(x,sa):
    b=sa*(sum(x**2)/11)**(1/2)
    return b

a1 ,b1 ,r1 ,p1 ,sa1 =stats.linregress(x1, y1)
a2 ,b2 ,r2 ,p2, sa2 =stats.linregress(x2, y2)
a3 ,b3 ,r3 ,p3, sa3 =stats.linregress(x3, y3)
a4 ,b4 ,r4 ,p4, sa4 =stats.linregress(x4, y4)
print('b)',a1, b1,r1,p1 , sa1 )
print('d)s',a2, b2,r2,p2,sa2)
print('d)r',a3, b3,r3,p3,sa3)
print('c)',a4, b4,r4,p4,sa4)


sb1=FehlerB(x1,sa1)
sb2=FehlerB(x2,sa2)
sb3=FehlerB(x3,sa3)
sb4=FehlerB(x4,sa4)


print('Steigung=',a1,'y-Ach=', b1 ,'Fehler der Steigung=', sa1 ,'Fehler des y-Achsenabschnitts= ', sb1 )
print('Steigung=',a2,'y-Ach=', b2 ,'Fehler der Steigung=', sa2 ,'Fehler des y-Achsenabschnitts= ', sb2 )
print('Steigung=',a3,'y-Ach=', b3 ,'Fehler der Steigung=', sa3 ,'Fehler des y-Achsenabschnitts= ', sb3 )
print('Steigung=',a4,'y-Ach=', b4 ,'Fehler der Steigung=', sa4 ,'Fehler des y-Achsenabschnitts= ', sb4 )

p=(x1*y1)
ra=(y1/x1)
deltaP=((y1*0.03)**2 +(x1*0.045)**2)**(1/2)
deltaRa=((0.045/x1)**2 +(y1*0.03/(x1**2))**2)**(1/2)




print (p, deltaP, ra, deltaRa)
plt.figure(1)
plt.errorbar(y1/x1,x1*y1,yerr=deltaP, xerr=deltaRa,fmt='rx')
plt.plot(y1/x1,x1*y1,'xr')
ra1=np.linspace(0,50)
plt.plot(ra1,(b1**2/(ra1-a1)**2)*ra1)
plt.xlabel(  r'$Lastwiederstand\ R_a / \Omega$')
plt.ylabel(r'$Leistung\ P/W$')
plt.savefig('b)leistung.jpg')






plt.figure(2)
plt.plot(y2/x2,x2*y2,'xr')
ra2=np.linspace(0,250)
plt.plot(ra2,(b2**2/(ra2-a2)**2)*ra2)
plt.savefig('d)rleistung.pdf')




plt.figure(3)
plt.plot(y3/x3,x3*y3,'xr')
ra3=np.linspace(0.1,5000)
plt.plot(ra3,(b3**2/(ra3-a3)**2)*ra3)
plt.savefig('d)sleistung.pdf')

plt.figure(4)
x1lin=np.linspace(0,0.15)
plt.plot(x1lin,a1*x1lin+b1,'-b')
plt.errorbar(x1,y1,yerr=0.45, xerr=0.03, fmt='xr')
plt.plot(x1,y1,'kx')
plt.savefig('b).jpg')

plt.figure(5)
x2lin=np.linspace(0,0.008)
plt.plot(x2lin,a2*x2lin+b2,'-b')
plt.errorbar(x2,y2,yerr=0.015, xerr=0.0003, fmt='rx')
plt.plot(x2,y2,'kx')
plt.savefig('d)recht.jpg')

plt.figure(6)
x3lin=np.linspace(0,0.0012)
plt.plot(x3lin,a3*x3lin+b3,'-b')
plt.errorbar(x3,y3,yerr=0.015, xerr=0.0003, fmt='rx')
plt.plot(x3,y3,'kx')
plt.savefig('d)sin.jpg')

plt.figure(7)
x4lin=np.linspace(0,0.2)
plt.plot(x4lin,a4*x4lin+b4,'-b')
plt.errorbar(x4,y4,yerr=0.15, xerr=0.03, fmt='rx')
plt.plot(x4,y4,'kx')
plt.savefig('c).jpg')


deltaP=((y1*0.03)**2 +(x1*0.045)**2)**(1/2)
