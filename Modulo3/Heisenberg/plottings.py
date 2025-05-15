import numpy as np
import matplotlib.pyplot as plt

#Load data from text file
dati = np.loadtxt("dati1.txt")

#Name columns
syslenght = dati[:, 0]
En = dati[:, 1]
Ebond = dati[:, 2]
Edens = dati[:, 3]
truncErr = dati[:, 4]
ent = dati[:, 5]

#Plots
plt.plot(syslenght, ent, marker='.', linestyle='-')
plt.xlabel('N')
plt.ylabel('S')
plt.title('Bipartite entanglement')
plt.grid(True)
plt.show()

plt.plot(syslenght, Edens, marker='.', linestyle='-' )
plt.xlabel('N')
plt.ylabel('Edens')
plt.title('Energy per site')
plt.grid(True)
plt.show()