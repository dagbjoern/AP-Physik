
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

x1 , y1 = np.genfromtxt('b).txt', unpack=True)
x2 , y2 = np.genfromtxt('d)recht.txt',unpack=True)
x3 , y3 = np.genfromtxt('d)sin.txt', unpack=True)
x4 , y4 = np.genfromtxt('c).txt', unpack=True)
x1mit=np.sum(x1)/11
x2mit=np.sum(x2)/11
x3mit=np.sum(x3)/11
x4mit=np.sum(x4)/11

y1mit=np.sum(y1)/11
y2mit=np.sum(y2)/11
y3mit=np.sum(y3)/11
y4mit=np.sum(y4)/11

stdx1=np.std(x1)
stdx2=np.std(x2)
stdx3=np.std(x3)
stdx4=np.std(x4)

stdy1=np.std(y1)
stdy2=np.std(y2)
stdy3=np.std(y3)
stdy4=np.std(y4)


x1sumq=np.sum(x1**2)
print(x1sumq)
print('b) mittx=',x1mit, 'mitty=',y1mit, 'stdx=',stdx1, 'stdy=',stdy1)
print('d)recht mittx=',x2mit, 'mitty=',y2mit, 'stdx=',stdx2, 'stdy=',stdy2)
print('d)sin mittx=',x3mit, 'mitty=',y3mit, 'stdx=',stdx3, 'stdy=',stdy3)
print('c) mittx=',x4mit, 'mitty=',y4mit, 'stdx=',stdx4, 'stdy=',stdy4)


#np.sum(np.sum)
