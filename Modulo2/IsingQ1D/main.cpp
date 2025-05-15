#include <iostream>
#include <cmath>
#include <fstream>
#include<math.h>
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<time.h>
#include<math.h>
#include<stdint.h>
#include<complex>
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Eigenvalues"
#include "mylib.h"

using namespace std;
using namespace Eigen;


/*
	PROGRAMMA PER LA DIAGONALIZZAZIONE ESATTA DEL MODELLO DI ISING-1D-QUANTISTICO;
	DIAGONALIZZAZIONE OTTIMIZZATA TRAMITE LE ROUTINE DI EIGEN, PER MATRICI SPARSE;
	VENGONO USATE P.B.C;

	H=-J*\sum_{i} \sigma_{i}^{z}\sigma_{i+1}^{z} -h*\sum_{i} \sigma_{i}^{z} -g*\sum_{i}\sigma_{i}^{x}; 

	SONO PRESENTI IL TERMINE TRASVERSO E QUELLO LONGITUDINALE, IN PARTICOLARE J=1, IL PARAMETRO RILEVANTE è g==g/J;

*/


int main()
{

    SparseMatrix<int> base(rows,N);

	// COSTRUZIONE DELLA BASE COMPUTAZIONALE;	
	costruisciBase(base);
	double g, h;
	
	// DICHIARAZIONE DEI PARAMETRI HAMILTONIANI g E h;
	
    g = 0.05;
	h =	1.00;

	// CREAZIONE DEI FILE DI OUTPUT PER IL GROUND_STATE E PER LE ENERGIE;
	ofstream file1("N" + to_string(N) + "/open.txt");
	ofstream file2("N" + to_string(N) + "/eopen.txt");
	
	// APERTURA DEI FILE SU CUI SCRIVERE;
	if (!file1.is_open() || !file2.is_open()) {
        cerr << "Impossibile aprire uno o entrambi i file di output." << endl;
        return 1; // Esci con un codice di errore
    }

	// INIZIO DEL CICLO CHE RESTITUISCE I VARI GROUND_STATE ED ENERGIE AL VARIARE DEI PARAMETRI HAMILTONIANI;
	while( g <= 1.55 )
	{
		SparseMatrix<double> H(rows,rows);
		hamiltonianaOPEN(H, base, g, h);

		//VIENE CREATA L'HAMILTONIANA DEL SISTEMA E VIENE CONTROLLATA LA SUA HERMITIANITA'

		// CREAZIONE DEGLI AUTOVALORI E DEGLI AUTOVETTORI;
		SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(H);

		VectorXd eigenvaluesH = eigensolver.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
        MatrixXd eigenvectorsH = eigensolver.eigenvectors(); // MARICE DEGLI AUTOVETTORI;
        
        // ORGANIZZAZIONE DEI DATI;

        VectorXd sortedEigenvaluesH(eigenvaluesH.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
        MatrixXd sortedEigenvectorsH(eigenvectorsH.rows(), eigenvaluesH.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;

		// ORDINA LO SPETTRO ENERGETICO E CONTEMPORANEAMENTE RIORDINA LA BASE;
        ORDINE_ASC(eigenvaluesH,eigenvectorsH,sortedEigenvaluesH,sortedEigenvectorsH);


		// SALVA SULLA PRIMA COLONNA IL VALORE DI g/h E SULLE RIGHE I VARI GROUND STATE;
		
		file1 << g << " ";
		for(int i = 0; i < rows; i++)
		{
			file1 << sortedEigenvectorsH(i,0) << " ";
		}
		file1 << endl;
		
		// SALVA LO SPETTRO ENERGETICO IN MODO ORDINATO

		file2 << g << " ";
		for(int i = 0; i < rows; i++)
		{
			file2 << sortedEigenvaluesH(i) << " ";
		}
		file2 << endl;

		g += 0.05;
	}


	file1.close();
	file2.close();	

	return 0;

}