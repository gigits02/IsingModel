import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

#LEVEL SPACE STATISTICS AND ENERGY SCALING

#Load data from text file
energy010 = np.loadtxt("./N10/eh0.txt")
energy10 = np.loadtxt("./N10/eopen.txt")

data3 = np.loadtxt("./N3/0solutions.txt")
data4 = np.loadtxt("./N4/0solutions.txt")
data5 = np.loadtxt("./N5/0solutions.txt")
data6 = np.loadtxt("./N6/0solutions.txt")
data7 = np.loadtxt("./N7/0solutions.txt")
data8 = np.loadtxt("./N8/0solutions.txt")
data9 = np.loadtxt("./N9/0solutions.txt")
data10 = np.loadtxt("./N10/0solutions.txt")

'''
#Sample data with given frequency
w = 1  # (Modifica questo valore a tuo piacimento)
sample = data[::w]
'''

#LSS

#Name columns
g = data3[:,4]

#CONSIDERO LE DIFFERENZE DI ENERGIE TRA STATI SUCCESSIVI
#CIASCUN ARRAY DIVISO PER LA MEDIA DELLO SPACING

e010 = energy010[9,1:]
e10 = energy10[19,1:]
lss010 = np.diff(e010)
lss010 /= np.mean(lss010)
lss10 = np.diff(e10)
lss10 /= np.mean(lss10)
print('vettore lss NON integrabile:')
#STIMA DELLA BINSIZE PER GLI ISTOGRAMMI
std_dev0 = np.std(lss010)
num_observations0 = len(lss010)
bin_width0 = (3.5 * std_dev0) / (num_observations0 ** (1/3))

std_dev = np.std(lss10)
num_observations = len(lss10)
bin_width = (3.5 * std_dev) / (num_observations ** (1/3))

# FITTING POISSONIAN
# Calcola gli istogrammi per lss010
counts0, bin_edges0 = np.histogram(lss010, bins = int(len(lss010) / bin_width0), density=True)
bin_centers0 = (bin_edges0[:-1] + bin_edges0[1:]) / 2


def decay(x,c,lamb):
    return c*np.exp(-lamb*x)

params, cov = curve_fit(decay, bin_centers0, counts0)
print('C*exp(-lambda x) (rispettivamente parametri c e lambda):')
print(params[0], '+/-', np.sqrt(cov[0,0]))
print(params[1], '+/-', np.sqrt(cov[1,1]))

#FITTING WIGNER-DYSON
# Calcola gli istogrammi per lss10
counts, bin_edges = np.histogram(lss10, bins = int(len(lss10) / bin_width), density=True)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
def wigner_dysonGOE(s, a, b):
    return (a * np.pi * s / 2) * np.exp(-b * np.pi * s**2 / 4)

def wigner_dyson2GUE(s,a,b):
    return (a * 32 * s**2 )/(np.pi**2) * np.exp(-b * 4 * s**2 / np.pi)

def wigner_dyson2DEF(s,a,b,c):
    return (a * s**c ) * np.exp(-b * s**2 )

params_wd, cov_wd = curve_fit(wigner_dyson2DEF, bin_centers, counts)
#print('(A*pi*x/2) * exp(- A * pi * (x^2)/4) parametro A:')
#print(params_wd[0], '+/-', np.sqrt(cov_wd[0,0]))

print('(A* s**c) * exp(- b x^2) parametro A:')
print(params_wd[0], '+/-', np.sqrt(cov_wd[0,0]))
print('(A* s**c) * exp(- b x^2) parametro B:')
print(params_wd[1], '+/-', np.sqrt(cov_wd[1,1]))
print('(A* s**c) * exp(- b x^2) parametro C:')
print(params_wd[2], '+/-', np.sqrt(cov_wd[2,2]))

#PLOTTING HISTOGRAMS STATS
plt.hist(lss010, bins = int(len(lss010) / bin_width0), density=True, label='(h,g) = (0, 0.5)')
#plt.plot(bin_centers0, decay(bin_centers0, *params), color='red', label='Fit Poisson')
x = np.linspace(0.09, 3.0, 100)
#plt.plot(x, decay(x, *params), color='red', label='Fit Poisson')
plt.xlabel('Level Spacing')
plt.ylabel('Probability Density')
plt.title('Level Spacing Statistics per N=10, sistema integrabile')
plt.xlim(0.0, 3.0)
plt.ylim(0.0, 4.5)
plt.legend()
plt.show()

plt.hist(lss10, bins = int(len(lss10) / bin_width), density=True, label='h,g = (1, 1)')
#plt.plot(bin_centers, wigner_dyson2DEF(bin_centers, *params_wd), color='red', label='Fit Wigner-Dyson')
plt.xlabel('Level Spacing')
plt.ylabel('Probability Density')
plt.title('Level Spacing Statistics per N=10, sistema non integrabile')
#plt.xlim(0.0, 200.0)
plt.xlim(0.0, 7.0)
plt.legend()
plt.show()

#CHECK NORMALIZZAZIONE ISTOGRAMMI
counts1, bin_edges1 = np.histogram(lss010, bins=20, density=True)
bin_widths1 = np.diff(bin_edges1)
area1 = np.sum(counts1 * bin_widths1)
print(f"Area totale sotto l'istogramma: {area1}")
counts2, bin_edges2 = np.histogram(lss10, bins=20, density=True)
bin_widths2 = np.diff(bin_edges2)
area2 = np.sum(counts2 * bin_widths2)
print(f"Area totale sotto l'istogramma: {area2}")



#ENERGY SCALING

#POWER LAW FIT OF GAP FOR EACH SIZE N
def power_func(N, a, c):
    return c*N**(-a)  

N = [3,4,5,6,7,8,9,10]
# a g = 1 non torna bene il fit (torna meglio quando g = 0.75)
deltaE = [data3[19,2], data4[19,2], data5[19,2], data6[19,2], data7[19,2], data8[19,2], data9[19,2], data10[19,2]]
#deltaE = [data3[14,2], data4[14,2], data5[14,2], data6[14,2], data7[14,2], data8[14,2], data9[14,2], data10[14,2]]

p0 = (1.0, 0.0)
popt, pcov = curve_fit(power_func, N, deltaE, p0, maxfev = 3000)
z, a = popt
print('Parametri del fit ottimizzati (rispettivamente z, c):')
print(z, '+/-', np.sqrt(pcov[0,0]), a, '+/-', np.sqrt(pcov[1,1]))
#PLOTTING FIT
plt.errorbar(N, deltaE, fmt='v', linestyle='None', capsize=3, markersize = 4,  label = 'Energy scaling')
N = np.linspace(N[0],N[-1],100)
plt.errorbar(N, power_func(N, *popt), fmt='', linestyle='--', capsize=3, markersize = 4, label = f'Fitting ΔE: y = {a:.2f}N^(-{z:.2f})')
plt.xlabel('N')
plt.ylabel('ΔE')
plt.title('Energy scaling (h, g = (0 , 1.0))')
plt.grid(True)
plt.legend()
#plt.xlim(0.4,0.5)
plt.show()
