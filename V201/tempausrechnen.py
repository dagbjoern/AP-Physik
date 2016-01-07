import numpy as np
import matplotlib.pyplot as plt
import uncertainties.unumpy as unp
from uncertainties.unumpy import (nominal_values as noms,
                                  std_devs as stds)
from scipy import stats
from scipy.optimize import curve_fit

cw=4.18
m_wk=379.83-139.77
m_ww=466.88-139.77
T_wk=21.49
T_ww=85.25
T_m=54.67




while 1==1:
 U_th=float(input('spannung?'))
 T=25.157*U_th-0.19*(U_th**2)
 print('temperatur=',T)


#cg=(cw*m_ww*(T_ww-T_m)-cw*m_wk*(T_m-T_wk))/(T_m-T_wk)
#print(cg)
