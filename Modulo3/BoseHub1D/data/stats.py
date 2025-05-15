import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#FIT CORRELAZIONI NEL FSDMRG

#Load data from text file
dataMI = np.loadtxt("./corr.txt")
dataSF = np.loadtxt("./corrCrit.txt")

xMI = dataMI[:,0]
yMI = dataMI[:,1]
xSF = dataSF[:,0]
ySF = dataSF[:,1]

#Fitting data
def power(x, c, a):
    return c*x**a

def exponential(x, c, t):
    return c*np.exp(-t*x)

poptsMI, covMI = curve_fit(exponential, xMI, yMI)
poptsSF, covSF = curve_fit(power, xSF[5:50], ySF[5:50])

print('exponential opt parameters ( c* exp(-t*x)) rispettivamente c e t:')
print(poptsMI[0], '+/-', np.sqrt(covMI[0,0]))
print(poptsMI[1], '+/-', np.sqrt(covMI[1,1]))
print('power opt parameters ( c* x**a rispettivamente c e a:')
print(poptsSF[0], '+/-', np.sqrt(covSF[0,0]))
print(poptsSF[1], '+/-', np.sqrt(covSF[1,1]))

#PLOTTINGS FITS
plt.errorbar(xMI, yMI, fmt='.', label = 'MI (μ,J) = (0.5, 0.05)', linestyle='None', capsize=3, markersize = 4)
x = np.linspace(xMI[0], xMI[-1], 100)
plt.plot(x, exponential(x, *poptsMI), color='red', label='Exponential fit')
plt.xlabel('|i-j|')
plt.ylabel('C(|i-j|)')
plt.title('Fit correlazioni')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()

plt.errorbar(xSF, ySF, fmt='.', label = 'SF (μ,J) = (0.5, 0.5)', linestyle='None', capsize=3, markersize = 4)
x = np.linspace(xSF[0], xSF[-1], 100)
plt.plot(x, power(x, *poptsSF), color='red', label='Power law fit')
plt.xlabel('|i-j|')
plt.ylabel('C(|i-j|)')
plt.title('Fit correlazioni')
plt.grid(True)
plt.legend()
plt.xlim(0,60)
plt.show()
