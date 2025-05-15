import numpy as np
import matplotlib.pyplot as plt

#MAKES FIT OVER THE POWER LAWS AND SAVES RESULTS INTO "FIT.TXT"

#Load data from text file
data1 = np.loadtxt("./L30/beta0340.txt")
data2 = np.loadtxt("./L30/beta0470.txt")
#Sample data with given frequency
w = 100  # (Modifica questo valore a tuo piacimento)
sample1 = data1[::w]
sample2 = data2[::w]


plt.scatter(range(len(sample1[:,1])), sample1[:,1], s = 1, color = 'darkorange')
plt.xlabel('Iterazione')
plt.ylabel('m')
plt.title('Storico Montecarlo (β = 0.34, L = 30)')
plt.show()

plt.scatter(range(len(sample2[:,1])),sample2[:,1], s = 1, color = 'darkorange')
plt.xlabel('Iterazione')
plt.ylabel('m')
plt.title('Storico Montecarlo (β = 0.47, L = 30)')
plt.show()
