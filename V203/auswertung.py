import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit


Tf , Tg , p1 =np.genefromtxt('1messung.txt',unpack=True)
T2 , p2 =np.genefromtxt('2messung.txt',unpack=True)


p1=p1*1e-3
Tg=Tg+273.15
T2=T2+273.15
p2=8.14+p2





np.savetxt('1tabelle.txt',np.column_stack((Tg,p1)),fmt='%r',delimiter=' & ')
np.savetxt('2tabelle.txt',np.column_stack((T,p2)),fmt='%r',delimiter=' & ')
