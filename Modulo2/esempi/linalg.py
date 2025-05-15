import numpy as np

# Definizione della dimensione della matrice
ROWS = 2
COLS = 2

# Creazione della matrice complessa iniziale M
M = np.array([[1 + 0j, 0 + 0j],
              [0 + 0j, -1 + 0j],])

print("Matrice M iniziale:")
print(M)

# Decomposizione agli autovalori e autovettori
eigval, eigvec = np.linalg.eig(M)

print("\nAutovalori della matrice:")
print(eigval)
print("\nAutovettori della matrice:")
print(eigvec)

# Creazione della matrice degli autovettori (reshape)
V = eigvec.reshape(ROWS, COLS)

print("\nMatrice V degli autovettori:")
print(V)


# Costruzione della matrice diagonale degli autovalori
D = np.diag(eigval)

print("\nMatrice D degli autovalori:")
print(D)

# Calcolo dell'approssimazione di M
M_approx = V @ D @ np.conj(V.T)

# Confronto tra M e M_approx
print("\nDifferenza tra M e VDV':")
print(M - M_approx)
