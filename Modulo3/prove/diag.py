import numpy as np

# Definisci la matrice da diagonalizzare
matrice = np.array([[1, 2, 3],
                    [4, 5, 6],
                    [7, 8, 9]])

# Calcola gli autovalori e gli autovettori
autovalori, autovettori = np.linalg.eig(matrice)

print("Autovalori:")
print(autovalori)

print("\nAutovettori:")
print(autovettori)