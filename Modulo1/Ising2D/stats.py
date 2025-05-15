import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.stats import chisquare

#MAKES FIT OVER THE POWER LAWS AND SAVES RESULTS INTO "FIT.TXT"

#Load data from text file
data10 = np.loadtxt("10.txt")

'''
#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample = data[::w]
'''
#Name columns
beta1 = data10[:, 0]
e10 = data10[:,1]
erE10 = data10[:,2]
c10 = 10*10*data10[:,3]
erc10 = 10*10*data10[:,4]
m10 = data10[:,5]
erM10 = data10[:,6]
absm10 = data10[:, 7]
chi10 = 10*10*data10[:, 8]
erChi10 = 10*10*data10[:,9]
binder10 = data10[:,10]
erBinder10 = data10[:,11]

data20 = np.loadtxt("20.txt")
#Name columns
e20 = data20[:,1]
erE20 = data20[:,2]
c20 = 20*20*data20[:,3]
erc20 = 20*20*data20[:,4]
m20 = data20[:,5]
erM20 = data20[:,6]
absm20 = data20[:, 7]
chi20 = 20*20*data20[:, 8]
erChi20 = 20*20*data20[:,9]
binder20 = data20[:,10]
erBinder20 = data20[:,11]

data30 = np.loadtxt("30.txt")
#Name columns
e30 = data30[:,1]
erE30 = data30[:,2]
c30 = 30*30*data30[:,3]
erc30 = 30*30*data30[:,4]
m30 = data30[:,5]
erM30 = data30[:,6]
absm30 = data30[:, 7]
chi30 = 30*30*data30[:, 8]
erChi30 = 30*30*data30[:,9]
binder30 = data30[:,10]
erBinder30 = data30[:,11]

data40 = np.loadtxt("40.txt")
#Name columns
beta = data40[:, 0]
e40 = data40[:,1]
erE40 = data40[:,2]
c40 = 40*40*data40[:,3]
erc40 = 40*40*data40[:,4]
m40 = data40[:,5]
erM40 = data40[:,6]
absm40 = data40[:, 7]
chi40 = 40*40*data40[:, 8]
erChi40 = 40*40*data40[:,9]
binder40 = data40[:,10]
erBinder40 = data40[:,11]

data60 = np.loadtxt("60.txt")
#Name columns
e60 = data60[:,1]
erE60 = data60[:,2]
c60 = 60*60*data60[:,3]
erc60 = 60*60*data60[:,4]
m60 = data60[:,5]
erM60 = data60[:,6]
absm60 = data60[:, 7]
chi60 = 60*60*data60[:, 8]
erChi60 = 60*60*data60[:,9]
binder60 = data60[:,10]
erBinder60 = data60[:,11]

data80 = np.loadtxt("80.txt")
#Name columns
e80 = data80[:,1]
erE80 = data80[:,2]
c80 = 80*80*data80[:,3]
erc80 = 80*80*data80[:,4]
m80 = data80[:,5]
erM80 = data80[:,6]
absm80 = data80[:, 7]
chi80 = 80*80*data80[:, 8]
erChi80 = 80*80*data80[:,9]
binder80 = data80[:,10]
erBinder80 = data80[:,11]

data100 = np.loadtxt("100.txt")
#Name columns
e100 = data100[:,1]
erE100 = data100[:,2]
c100 = 100*100*data100[:,3]
erc100 = 100*100*data100[:,4]
m100 = data100[:,5]
erM100 = data100[:,6]
absm100 = data100[:, 7]
chi100 = 100*100*data100[:, 8]
erChi100 = 100*100*data100[:,9]
binder100 = data100[:,10]
erBinder100 = data100[:,11]

data120 = np.loadtxt("120.txt")
#Name columns
e120 = data120[:,1]
erE120 = data120[:,2]
c120 = 120*120*data120[:,3]
erc120 = 120*120*data120[:,4]
m120 = data120[:,5]
erM120 = data120[:,6]
absm120 = data120[:, 7]
chi120 = 120*120*data120[:, 8]
erChi120 = 120*120*data120[:,9]
binder120 = data120[:,10]
erBinder120 = data120[:,11]

data140 = np.loadtxt("140.txt")
#Name columns
e140 = data140[:,1]
erE140 = data140[:,2]
c140 = 140*140*data140[:,3]
erc140 = 140*140*data140[:,4]
m140 = data140[:,5]
erM140 = data140[:,6]
absm140 = data140[:, 7]
chi140 = 140*140*data140[:, 8]
erChi140 = 140*140*data140[:,9]
binder140 = data140[:,10]
erBinder140 = data140[:,11]

data160 = np.loadtxt("160.txt")
#Name columns
e160 = data160[:,1]
erE160 = data160[:,2]
c160 = 160*160*data160[:,3]
erc160 = 160*160*data160[:,4]
m160 = data160[:,5]
erM160 = data160[:,6]
absm160 = data160[:, 7]
chi160 = 160*160*data160[:, 8]
erChi160 = 160*160*data160[:,9]
binder160 = data160[:,10]
erBinder160 = data160[:,11]

data180 = np.loadtxt("180.txt")
#Name columns
e180 = data180[:,1]
erE180 = data180[:,2]
c180 = 180*180*data180[:,3]
erc180 = 180*180*data180[:,4]
m180 = data180[:,5]
erM180 = data180[:,6]
absm180 = data180[:, 7]
chi180 = 180*180*data180[:, 8]
erChi180 = 180*180*data180[:,9]
binder180 = data180[:,10]
erBinder180 = data180[:,11]

#PARABOLIC FIT OF CHI FOR EACH SIZE

def quadratic_func(x, a, betaC, chiMax):
    return chiMax + a*((x-betaC)**2)

for i in range(66):
    if( chi10[i] == np.max(chi10)):
        indice10 = i
    if( chi20[i] == np.max(chi20)):
        indice20 = i
    if( chi30[i] == np.max(chi30)):
        indice30 = i
for i in range(31):
    if( chi40[i] == np.max(chi40)):
        indice40 = i
    if( chi60[i] == np.max(chi60)):
        indice60 = i
    if( chi80[i] == np.max(chi80)):
        indice80 = i
    if( chi100[i] == np.max(chi100)):
        indice100 = i
    if( chi120[i] == np.max(chi120)):
        indice120 = i
    if( chi140[i] == np.max(chi140)):
        indice140 = i
    if( chi160[i] == np.max(chi160)):
        indice160 = i
    if( chi180[i] == np.max(chi180)):
        indice180 = i
    
chi10[indice10+1] += 0.015
y10 = chi10[indice10-5:indice10+6]
yerr10 = erChi10[indice10-5:indice10+6]
y20 = chi20[indice20-5:indice20+6]
yerr20 = erChi20[indice20-5:indice20+6]
y30 = chi30[indice30-4:indice30+4]
yerr30 = erChi30[indice30-4:indice30+4]
y40 = chi40[indice40-5:indice40+5]
yerr40 = erChi40[indice40-5:indice40+5]
y60 = chi60[indice60-4:indice60+6]
yerr60 = erChi60[indice60-4:indice60+6]
y80 = chi80[indice80-4:indice80+5]
yerr80 = erChi80[indice80-4:indice80+5]
y100 = chi100[indice100-2:indice100+3]
yerr100 = erChi100[indice100-2:indice100+3]
y120 = chi120[indice120-2:indice120+3]
yerr120 = erChi120[indice120-2:indice120+3]
y140 = chi140[indice140-1:indice140+3]
yerr140 = erChi140[indice140-1:indice140+3]
y160 = chi160[indice160-1:indice160+3]
yerr160 = erChi160[indice160-1:indice160+3]
y180 = chi180[indice180-1:indice180+3]
yerr180 = erChi180[indice180-1:indice180+3]

x10 = beta1[indice10-5:indice10+6]
x20 = beta1[indice20-5:indice20+6]
x30 = beta1[indice30-4:indice30+4]
x40 = beta[indice40-5:indice40+5]
x60 = beta[indice60-4:indice60+6]
x80 = beta[indice80-4:indice80+5]
x100 = beta[indice100-2:indice100+3]
x120 = beta[indice120-2:indice120+3]
x140 = beta[indice140-1:indice140+3]
x160 = beta[indice160-1:indice160+3]
x180 = beta[indice180-1:indice180+3]

'''
PARAMETRI INIZIALI CONSIGLIATI
POPT20:
[-1.56813270e+04  4.17570311e-01  2.11377032e+01]
POPT30:
[-65006.2088  0.425245361  42.5387814]
POPT40:
[-1.91390000e+05  4.29351207e-01  7.02322772e+01]
POPT60:
[-9.42819643e+05  4.33289657e-01  1.42824038e+02]
POPT80:
[-2.53728286e+06  4.35242713e-01  2.34133014e+02]
POPT100:
[-5.52892856e+06  4.36145941e-01  3.45868989e+02]
POPT140:
[-1.85814125e+07  4.37597860e-01  6.16384864e+02]
POPT160:
[-2.93132800e+07  4.37993083e-01  7.75747706e+02]
POPT180:
[-4.28076900e+07  4.38240409e-01  9.51587756e+02]
'''
p0 = (-138.493479,0.39307315,6.41681622)
popt10, pcov10 = curve_fit(quadratic_func, x10, y10, sigma = yerr10, maxfev = 3000)
p0 = (-15681.3270,0.417570311,21.1377032)
popt20, pcov20 = curve_fit(quadratic_func, x20, y20, sigma = yerr20, p0=p0, maxfev = 3000)
p0 = (-65006.2088,0.425245361,42.5387814)
popt30, pcov30 = curve_fit(quadratic_func, x30, y30, sigma = yerr30, p0=p0, maxfev = 3000)
p0 = (-191390.000,0.429351207,70.2322772)
popt40, pcov40 = curve_fit(quadratic_func, x40, y40, sigma = yerr40, p0=p0, maxfev = 3000)
p0 = popt40
popt60, pcov60 = curve_fit(quadratic_func, x60, y60, sigma = yerr60, p0=p0, maxfev = 3000)
p0 = popt60
popt80, pcov80 = curve_fit(quadratic_func, x80, y80, sigma = yerr80, p0=p0, maxfev = 3000)
p0 = popt80
popt100, pcov100 = curve_fit(quadratic_func, x100, y100, sigma = yerr100, p0=p0, maxfev = 3000)
p0 = popt100
popt120, pcov120 = curve_fit(quadratic_func, x120, y120, sigma = yerr120, p0=p0, maxfev = 3000)
p0 = popt120
popt140, pcov140 = curve_fit(quadratic_func, x140, y140, sigma = yerr140, p0=p0, maxfev = 3000)
p0 = popt140
popt160, pcov160 = curve_fit(quadratic_func, x160, y160, sigma = yerr160, p0=p0, maxfev = 3000)
p0 = popt160
popt180, pcov180 = curve_fit(quadratic_func, x180, y180, sigma = yerr180, p0=p0, maxfev = 3000)

'''
chisquared10, pvalue10 = chisquare(y10, f_exp=quadratic_func(x10, *popt10))
chisquared20, pvalue20 = chisquare(y20, f_exp=quadratic_func(x20, *popt20))
chisquared30, pvalue30 = chisquare(y30, f_exp=quadratic_func(x30, *popt30))
chisquared40, pvalue40 = chisquare(y40, f_exp=quadratic_func(x40, *popt40))
chisquared60, pvalue60 = chisquare(y60, f_exp=quadratic_func(x60, *popt60))
chisquared80, pvalue80 = chisquare(y80, f_exp=quadratic_func(x80, *popt80))
chisquared100, pvalue100 = chisquare(y100, f_exp=quadratic_func(x100, *popt100))
chisquared120, pvalue120 = chisquare(y120, f_exp=quadratic_func(x120, *popt120))
chisquared140, pvalue140 = chisquare(y140, f_exp=quadratic_func(x140, *popt140))
chisquared160, pvalue160 = chisquare(y160, f_exp=quadratic_func(x160, *popt160))
chisquared180, pvalue180 = chisquare(y180, f_exp=quadratic_func(x180, *popt180))

#PRINTING RESULTS
print('POPT10:')
print(popt10)
print('PCOV10:')
print(pcov10)
print('Chi-squared, p-value:', chisquared10, pvalue10)
print('POPT20:')
print(popt20)
print('PCOV20:')
print(pcov20)
print('Chi-squared, p-value:', chisquared20, pvalue20)
print('POPT30:')
print(popt30)
print('PCOV30:')
print(pcov30)
print('Chi-squared, p-value:', chisquared30, pvalue30)
print('POPT40:')
print(popt40)
print('PCOV40:')
print(pcov40)
print('Chi-squared, p-value:', chisquared40, pvalue40)
print('POPT60:')
print(popt60)
print('PCOV60:')
print(pcov60)
print('Chi-squared, p-value:', chisquared60, pvalue60)
print('POPT80:')
print(popt80)
print('PCOV80:')
print(pcov80)
print('Chi-squared, p-value:', chisquared80, pvalue80)
print('POPT100:')
print(popt100)
print('PCOV100:')
print(pcov100)
print('Chi-squared, p-value:', chisquared100, pvalue100)
print('POPT120:')
print(popt120)
print('PCOV120:')
print(pcov120)
print('Chi-squared, p-value:', chisquared120, pvalue120)
print('POPT140:')
print(popt140)
print('PCOV140:')
print(pcov140)
print('Chi-squared, p-value:', chisquared140, pvalue140)
print('POPT160:')
print(popt160)
print('PCOV160:')
print(pcov160)
print('Chi-squared, p-value:', chisquared160, pvalue160)
print('POPT180:')
print(popt180)
print('PCOV180:')
print(pcov180)
print('Chi-squared, p-value:', chisquared180, pvalue180)
'''

#PLOTTING FITS:
plt.errorbar(x10, y10, yerr=yerr10, fmt='v', label = 'L=10', linestyle='None', capsize=3, markersize = 4)
xx10 = np.linspace(x10[0], x10[-1], 100)
plt.errorbar(xx10, quadratic_func(xx10, *popt10), fmt='-', label = 'Quadratic fit', linestyle='--', capsize=3, markersize = 4)
plt.legend()
plt.xlabel('β')
plt.ylabel('χ')
plt.grid(True)
#plt.ylim(5.750,7.0)
plt.show()

plt.errorbar(x20, y20, yerr=yerr20, fmt='v', label = 'L=20', linestyle='None', capsize=3, markersize = 4)
xx20 = np.linspace(x20[0], x20[-1], 100)
plt.errorbar(xx20, quadratic_func(xx20, *popt20), fmt='-', label = 'Fitting chi20', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi20\nChi-squared: {chisquared20:.2f}, p-value: {pvalue20:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(16.75,22)
plt.show()

plt.errorbar(x30, y30, yerr=yerr30, fmt='v', label = 'L=30', linestyle='None', capsize=3, markersize = 4)
xx30 = np.linspace(x30[0], x30[-1], 100)
plt.errorbar(xx30, quadratic_func(xx30, *popt30), fmt='-', label = 'Fitting chi30', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi30\nChi-squared: {chisquared30:.2f}, p-value: {pvalue30:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(16.75,22)
plt.show()

plt.errorbar(x40, y40, yerr=yerr40, fmt='v', label = 'L=40', linestyle='None', capsize=3, markersize = 4)
xx40 = np.linspace(x40[0], x40[-1], 100)
plt.errorbar(xx40, quadratic_func(xx40, *popt40), fmt='-', label = 'Fitting chi40', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi40\nChi-squared: {chisquared40:.2f}, p-value: {pvalue40:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(28,45)
plt.show()

plt.errorbar(x60, y60, yerr=yerr60, fmt='v', label = 'L=60', linestyle='None', capsize=3, markersize = 4)
xx60 = np.linspace(x60[0], x60[-1], 100)
plt.errorbar(xx60, quadratic_func(xx60, *popt60), fmt='-', label = 'Fitting chi60', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi60\nChi-squared: {chisquared60:.2f}, p-value: {pvalue60:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(100,141)
plt.show()


plt.errorbar(x80, y80, yerr=yerr80, fmt='v', label = 'L=80', linestyle='None', capsize=3, markersize = 4)
xx80 = np.linspace(x80[0], x80[-1], 100)
plt.errorbar(xx80, quadratic_func(xx80, *popt80), fmt='-', label = 'Fitting chi80', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi80\nChi-squared: {chisquared80:.2f}, p-value: {pvalue80:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(125,250)
plt.show()

plt.errorbar(x100, y100, yerr=yerr100, fmt='v', label = 'L=100', linestyle='None', capsize=3, markersize = 4)
xx100 = np.linspace(x100[0], x100[-1], 100)
plt.errorbar(xx100, quadratic_func(xx100, *popt100), fmt='-', label = 'Fitting chi100', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi100\nChi-squared: {chisquared100:.2f}, p-value: {pvalue100:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(100,390)
plt.show()

plt.errorbar(x120, y120, yerr=yerr120, fmt='v', label = 'L=120', linestyle='None', capsize=3, markersize = 4)
xx120 = np.linspace(x120[0], x120[-1], 100)
plt.errorbar(xx120, quadratic_func(xx120, *popt120), fmt='-', label = 'Fitting chi120', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi120\nChi-squared: {chisquared120:.2f}, p-value: {pvalue120:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(100,500)
plt.show()

plt.errorbar(x140, y140, yerr=yerr140, fmt='v', label = 'L=140', linestyle='None', capsize=3, markersize = 4)
xx140 = np.linspace(x140[0], x140[-1], 100)
plt.errorbar(xx140, quadratic_func(xx140, *popt140), fmt='-', label = 'Fitting chi140', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi140\nChi-squared: {chisquared140:.2f}, p-value: {pvalue140:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(190,610)
plt.show()

plt.errorbar(x160, y160, yerr=yerr160, fmt='v', label = 'L=160', linestyle='None', capsize=3, markersize = 4)
xx160 = np.linspace(x160[0], x160[-1], 100)
plt.errorbar(xx160, quadratic_func(xx160, *popt160), fmt='-', label = 'Fitting chi160', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi160\nChi-squared: {chisquared160:.2f}, p-value: {pvalue160:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(200,800)
plt.show()

plt.errorbar(x180, y180, yerr=yerr180, fmt='v', label = 'L=180', linestyle='None', capsize=3, markersize = 4)
xx180 = np.linspace(x180[0], x180[-1], 100)
plt.errorbar(xx180, quadratic_func(xx180, *popt180), fmt='-', label = 'Fitting chi180', linestyle='--', capsize=3, markersize = 4)
#plt.legend([f'Fitting chi180\nChi-squared: {chisquared180:.2f}, p-value: {pvalue180:.2f}'])
plt.xlabel('x')
plt.ylabel('y')
plt.grid(True)
#plt.ylim(200,1000)
plt.show()

#SALVA SOLUZIONI SU FILE DI OUTPUT
# Leggi le vecchie soluzioni se presenti
old_solutions = []
try:
    with open("fit.txt", "r") as file:
        old_solutions = file.readlines()
except FileNotFoundError:
    pass
#Scrivi le soluzioni
'''
with open("fit.txt", "w") as file:
    #Scrivi le vecchie soluzioni
    for solution in old_solutions:
        file.write(solution)

    #Scrivi le nuove soluzioni
    file.write(f"{' '.join(map(str, popt10[1:4]))} {np.sqrt(pcov10[1][1])} {np.sqrt(pcov10[2][2])}\n")
    file.write(f"{' '.join(map(str, popt20[1:4]))} {np.sqrt(pcov20[1][1])} {np.sqrt(pcov20[2][2])}\n")
    file.write(f"{' '.join(map(str, popt30[1:4]))} {np.sqrt(pcov30[1][1])} {np.sqrt(pcov30[2][2])}\n")
    file.write(f"{' '.join(map(str, popt40[1:4]))} {np.sqrt(pcov40[1][1])} {np.sqrt(pcov40[2][2])}\n")
    file.write(f"{' '.join(map(str, popt60[1:4]))} {np.sqrt(pcov60[1][1])} {np.sqrt(pcov60[2][2])}\n")
    file.write(f"{' '.join(map(str, popt80[1:4]))} {np.sqrt(pcov80[1][1])} {np.sqrt(pcov80[2][2])}\n")
    file.write(f"{' '.join(map(str, popt100[1:4]))} {np.sqrt(pcov100[1][1])} {np.sqrt(pcov100[2][2])}\n")
    file.write(f"{' '.join(map(str, popt120[1:4]))} {np.sqrt(pcov120[1][1])} {np.sqrt(pcov120[2][2])}\n")
    file.write(f"{' '.join(map(str, popt140[1:4]))} {np.sqrt(pcov140[1][1])} {np.sqrt(pcov140[2][2])}\n")
    file.write(f"{' '.join(map(str, popt160[1:4]))} {np.sqrt(pcov160[1][1])} {np.sqrt(pcov160[2][2])}\n")
    file.write(f"{' '.join(map(str, popt180[1:4]))} {np.sqrt(pcov180[1][1])} {np.sqrt(pcov180[2][2])}\n")
'''