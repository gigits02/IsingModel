import numpy as np
import matplotlib.pyplot as plt


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

#PLOTTING ENERGIES
#plt.errorbar(beta1, e10, yerr=erE10/10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, e20, yerr=erE20/10, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, e30, yerr=erE30/10, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e40, yerr=erE40/10, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e60, yerr=erE60/10, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e80, yerr=erE80/10, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e100, yerr=erE100/10, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e120, yerr=erE120/10, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e140, yerr=erE140/10, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e160, yerr=erE160/10, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, e180, yerr=erE180/10, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('β')
plt.ylabel('ε')
plt.title('Densità di energia')
plt.grid(True)
plt.legend()
plt.xlim(0.41,0.47)
plt.show()

#PLOTTING SPECIFIC HEAT
plt.errorbar(beta1, c10, yerr=erc10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, c20, yerr=erc20, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, c30, yerr=erc30, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c40, yerr=erc40, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c60, yerr=erc60, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c80, yerr=erc80, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c100, yerr=erc100, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c120, yerr=erc120, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c140, yerr=erc140, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c160, yerr=erc160, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, c180, yerr=erc180, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('β')
plt.ylabel('c')
plt.title('Calore specifico')
plt.grid(True)
plt.legend()
#plt.xlim(0.41,0.47)
plt.show()

#PLOTTING HISTOGRAM ON MAGNETIZATION (sample of L=10, beta = 0.340 < betaC)
rawdata = np.loadtxt("./L10/beta0340.txt")
m = rawdata[:, 1]
plt.hist(m, bins=1000, edgecolor='orange')
plt.xlabel('m')
plt.ylabel('Occurence')
plt.title('m (β < βC )')
plt.show()
#PLOTTING HISTOGRAM ON MAGNETIZATION (sample of L=10, beta = 0.470 > betaC)
rawdata = np.loadtxt("./L10/beta0470.txt")
m = rawdata[:, 1]
plt.hist(m, bins=1000, edgecolor='orange')
plt.xlabel('m')
plt.ylabel('Occurence')
plt.title('m (β > βC)')
plt.show()

#PLOTTING AVERAGE MAGNETIZATION
#plt.errorbar(beta1, m10, yerr=erM10, fmt='.', label = 'L=10', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta1, m20, yerr=erM20, fmt='.', label = 'L=20', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta1, m30, yerr=erM30, fmt='.', label = 'L=30', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m40, yerr=erM40, fmt='.', label = 'L=40', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m60, yerr=erM60, fmt='.', label = 'L=60', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m80, yerr=erM80, fmt='.', label = 'L=80', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m100, yerr=erM100, fmt='.', label = 'L=100', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m120, yerr=erM120, fmt='.', label = 'L=120', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m140, yerr=erM140, fmt='.', label = 'L=140', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m160, yerr=erM160, fmt='.', label = 'L=160', linestyle='None', capsize=3, markersize = 4)
plt.errorbar(beta, m180, yerr=erM180, fmt='.', label = 'L=180', linestyle='None', capsize=3, markersize = 4)
plt.xlabel('β')
plt.ylabel('<m>')
plt.title('Magnetizzazione')
plt.grid(True)
#plt.legend()
plt.xlim(0.41,0.47)
plt.show()

#PLOTTING ABS(M)
plt.errorbar(beta1, absm10, yerr=erM10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, absm20, yerr=erM20, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, absm30, yerr=erM30, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm40, yerr=erM40, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm60, yerr=erM60, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm80, yerr=erM80, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm100, yerr=erM100, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm120, yerr=erM120, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm140, yerr=erM140, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm160, yerr=erM160, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, absm180, yerr=erM180, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('β')
plt.ylabel('|m|')
plt.title('Valore assoluto della magnetizzazione')
plt.grid(True)
plt.legend()
plt.xlim(0.41,0.47)
plt.show()

#PLOTTING ORDER PARAM.
plt.errorbar(beta1, chi10, yerr=erChi10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, chi20, yerr=erChi20, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, chi30, yerr=erChi30, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi40, yerr=erChi40, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi60, yerr=erChi60, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi80, yerr=erChi80, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi100, yerr=erChi100, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi120, yerr=erChi120, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi140, yerr=erChi140, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi160, yerr=erChi160, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, chi180, yerr=erChi180, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('β')
plt.ylabel('χ')
plt.title('Suscettività')
plt.grid(True)
plt.legend(loc = 'upper left')
#plt.xlim(0.4,0.5)
plt.show()

#PLOTTING BINDER CUMULANT
plt.errorbar(beta1, binder10, yerr=erBinder10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, binder20, yerr=erBinder20, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta1, binder30, yerr=erBinder30, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder40, yerr=erBinder40, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder60, yerr=erBinder60, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder80, yerr=erBinder80, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder100, yerr=erBinder100, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder120, yerr=erBinder120, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder140, yerr=erBinder140, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder160, yerr=erBinder160, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(beta, binder180, yerr=erBinder180, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('β')
plt.ylabel('U4')
plt.title('Cumulante di Binder')
plt.grid(True)
plt.legend()
plt.xlim(0.41,0.47)
plt.show()

#FSS
#(for pasting: βγν)
b = 1./8.
nu = 1.
betac = 0.5*np.log(1+np.sqrt(2))
gamma = 7./4.

t10 = (beta1-betac)*(10**(1/nu))
t20 = (beta1-betac)*(20**(1/nu))
t30 = (beta1-betac)*(30**(1/nu))
t40 = (beta-betac)*(40**(1/nu))
t60 = (beta-betac)*(60**(1/nu))
t80 = (beta-betac)*(80**(1/nu))
t100 = (beta-betac)*(100**(1/nu))
t120 = (beta-betac)*(120**(1/nu))
t140 = (beta-betac)*(140**(1/nu))
t160 = (beta-betac)*(160**(1/nu))
t180 = (beta-betac)*(180**(1/nu))

f10 = absm10*(10**(b/nu))
f20 = absm20*(20**(b/nu))
f30 = absm30*(30**(b/nu))
f40 = absm40*(40**(b/nu))
f60 = absm60*(60**(b/nu))   
f80 = absm80*(80**(b/nu))
f100 = absm100*(100**(b/nu))
f120 = absm120*(120**(b/nu))
f140 = absm140*(140**(b/nu))
f160 = absm160*(160**(b/nu))
f180 = absm180*(180**(b/nu))

g10 = chi10/(10**(gamma/nu))
g20 = chi20/(20**(gamma/nu))
g30 = chi30/(30**(gamma/nu))
g40 = chi40/(40**(gamma/nu))
g60 = chi60/(60**(gamma/nu))
g80 = chi80/(80**(gamma/nu))
g100 = chi100/(100**(gamma/nu))
g120 = chi120/(120**(gamma/nu))
g140 = chi140/(140**(gamma/nu))
g160 = chi160/(160**(gamma/nu))
g180 = chi180/(180**(gamma/nu))

#FSS FOR ABS(M)
plt.errorbar(t10, f10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t20, f20, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t30, f30, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t40, f40, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t60, f60, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t80, f80, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t100, f100, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t120, f120, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t140, f140, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t160, f160, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t180, f180, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('(β - βC)*L^(1/ν)')
plt.ylabel('|m|*L^(β/ν)')
plt.title('Finite Size Scaling per |m|')
plt.grid(True)
plt.legend()
#plt.xlim(-5000,5000)
plt.show()

#FSS FOR CHI
plt.errorbar(t10, g10, fmt='.', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t20, g20, fmt='.', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t30, g30, fmt='.', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t40, g40, fmt='.', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t60, g60, fmt='.', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t80, g80, fmt='.', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t100, g100, fmt='.', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t120, g120, fmt='.', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t140, g140, fmt='.', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t160, g160, fmt='.', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t180, g180, fmt='.', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('(β - βC)*L^(1/ν)')
plt.ylabel('χ/L^(γ/ν)')
plt.title('Finite Size Scaling per la suscettività')
plt.grid(True)
plt.legend()
#plt.xlim(-5000,5000)
plt.show()

#FSS FOR BINDER CUMULANT
plt.errorbar(t10, binder10, yerr=erBinder10, fmt='-', label = 'L=10', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t20, binder20, yerr=erBinder20, fmt='-', label = 'L=20', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t30, binder30, yerr=erBinder30, fmt='-', label = 'L=30', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t40, binder40, yerr=erBinder40, fmt='-', label = 'L=40', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t60, binder60, yerr=erBinder60, fmt='-', label = 'L=60', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t80, binder80, yerr=erBinder80, fmt='-', label = 'L=80', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t100, binder100, yerr=erBinder100, fmt='-', label = 'L=100', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t120, binder120, yerr=erBinder120, fmt='-', label = 'L=120', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t140, binder140, yerr=erBinder140, fmt='-', label = 'L=140', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t160, binder160, yerr=erBinder160, fmt='-', label = 'L=160', linestyle='--', capsize=3, markersize = 4)
plt.errorbar(t180, binder180, yerr=erBinder180, fmt='-', label = 'L=180', linestyle='--', capsize=3, markersize = 4)
plt.xlabel('(β - βC)*L^(1/ν)')
plt.ylabel('U4')
plt.title('Finite Size Scaling per il cumulante di Binder')
plt.grid(True)
plt.legend()
#plt.xlim(-5000,5000)
plt.show()

