import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#Load data from text file
data10 = np.loadtxt("autoCr10.txt")
data20 = np.loadtxt("autoCr20.txt")
data40 = np.loadtxt("autoCr40.txt")
data60 = np.loadtxt("autoCr60.txt")
data80 = np.loadtxt("autoCr80.txt")
data100 = np.loadtxt("autoCr100.txt")

#Name columns
k10 = data10[:, 0]
c10 = data10[:,1]
k20 = data20[:, 0]
c20 = data20[:,1]
k40 = data40[:, 0]
c40 = data40[:,1]
k60 = data60[:, 0]
c60 = data60[:,1]
k80 = data80[:, 0]
c80 = data80[:,1]
k100 = data100[:, 0]
c100 = data100[:,1]

#FIT ESPONENZIALE PER LA FUNZIONE DI AUTOCORRELAZIONE  

def exp_func(x, A, B, tau):
        return A*np.exp(-x/tau) + B

#p0 = ()
popt10, pcov10 = curve_fit(exp_func, k10[:30], c10[:30]/c10[0], maxfev = 3000)
popt20, pcov20 = curve_fit(exp_func, k20[:30], c20[:30]/c20[0], maxfev = 3000)
popt40, pcov40 = curve_fit(exp_func, k40[:30], c40[:30]/c40[0], maxfev = 3000)
popt60, pcov60 = curve_fit(exp_func, k60[:30], c60[:30]/c60[0], maxfev = 3000)
popt80, pcov80 = curve_fit(exp_func, k80[:30], c80[:30]/c80[0], maxfev = 3000)
popt100, pcov100 = curve_fit(exp_func, k100[:30], c100[:30]/c100[0], maxfev = 3000)

#PRINT FIT PARAMETERS
print('tau10 = ', popt10[2], '+/-', np.sqrt(pcov10[2,2]))
print('A10 = ', popt10[0], '+/-', np.sqrt(pcov10[0,0]))
print('B10 = ', popt10[1], '+/-', np.sqrt(pcov10[1,1]))


#PLOTTING AUTOCORRELATION FUNCIONT C(k) FOR EACH L
plt.errorbar(k10, c10/c10[0], fmt='.', linestyle = '--', capsize = 3, markersize = 4, label = 'L=10')
plt.errorbar(k20, c20/c20[0],fmt='.', linestyle = '--', capsize = 3, markersize = 4, label = 'L=20')
plt.errorbar(k40, c40/c40[0], fmt='.', linestyle = '--', capsize = 3, markersize = 4, label = 'L=40')
plt.errorbar(k60, c60/c60[0], fmt='.', linestyle = '--', capsize = 3, markersize = 4, label = 'L=60')
plt.errorbar(k80, c80/c80[0], fmt='.', linestyle = '--', capsize = 3, markersize = 4, label = 'L=80')
plt.errorbar(k100, c100/c100[0], fmt='.', linestyle = '--', capsize = 3, markersize = 4, label = 'L=100')
plt.xlim(0, 50)
plt.xlabel('k')
plt.ylabel('C(k)')
plt.title('Funzione di autocorrelazione per m')
plt.grid(True)
plt.legend()
plt.show()

A,B,tau = popt10
#PLOTTING EXPONENTIAL FIT FOR L=10
plt.errorbar(k10, c10/c10[0], fmt='v', linestyle = 'None', capsize = 3, markersize = 4, label = 'L=10')
plt.errorbar(k10, exp_func(k10, *popt10), linestyle='--', color = 'darkorange', capsize=3, markersize = 4,  label = 'exp fit: Aexp(-t/τ) + B')
#plt.errorbar(k20, c20/c20[0],fmt='v', linestyle = 'None', capsize = 3, markersize = 4, label = 'L=20')
#plt.errorbar(k20, exp_func(k20, *popt20), linestyle='--', color = 'green', capsize=3, markersize = 4,  label = 'exp.fit for L=20')
#plt.errorbar(k40, c40/c40[0], fmt='v', linestyle = 'None', capsize = 3, markersize = 4, label = 'L=40')
#plt.errorbar(k40, exp_func(k40, *popt40), linestyle='--', color = 'red', capsize=3, markersize = 4,  label = 'exp.fit for L=40')
#plt.errorbar(k60, c60/c60[0], fmt='v', linestyle = 'None', capsize = 3, markersize = 4, label = 'L=60')
#plt.errorbar(x, exp_func(x, *popt60), linestyle='--', color = 'purple', capsize=3, markersize = 4,  label = 'exp.fit for L=60')
#plt.errorbar(k60, c80/c80[0], fmt='v', linestyle = 'None', capsize = 3, markersize = 4, label = 'L=80')
#plt.errorbar(x, exp_func(x, *popt80), linestyle='--', color = 'brown', capsize=3, markersize = 4,  label = 'exp.fit for L=80')
#plt.errorbar(k100, c100/c100[0], fmt='v', linestyle = 'None', capsize = 3, markersize = 4, label = 'L=100')
#plt.errorbar(x, exp_func(x, *popt100), linestyle='--', color = 'pink', capsize=3, markersize = 4,  label = 'exp.fit for L=100')
plt.xlim(0, 20)
plt.xlabel('k')
plt.ylabel('C(k)')
plt.title('Funzione di autocorrelazione per m')
plt.grid(True)
plt.legend()
plt.show()

#PLOTTING TAU AS FUNCTION OF SIZE L
x = [10,20,40,60,80,100]
y = [popt10[2], popt20[2], popt40[2], popt60[2], popt80[2], popt100[2]]
yerr = [np.sqrt(pcov10[2,2]), np.sqrt(pcov20[2,2]), np.sqrt(pcov40[2,2]), np.sqrt(pcov60[2,2]), np.sqrt(pcov80[2,2]), np.sqrt(pcov100[2,2])]
plt.errorbar(x, y, fmt='None', yerr = yerr, linestyle = 'None', capsize = 3, markersize = 4, label = 'L=100')
plt.xlabel('L')
plt.ylabel('tau_exp')
plt.title('Tau as function of size L')
plt.grid(True)
plt.legend()
plt.show()