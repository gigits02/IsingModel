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

double GAP(SparseMatrix<double>& H, VectorXd& gs,double g,double hI,SparseMatrix<int>& base)
{
    hamiltoniana(H, base, g, hI);
    SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(H);

    VectorXd eigenvaluesH = eigensolver.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
    MatrixXd eigenvectorsH = eigensolver.eigenvectors(); // MARICE DEGLI AUTOVETTORI;
    
    // ORGANIZZAZIONE DEI DATI;

    VectorXd sortedEigenvaluesH(eigenvaluesH.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
    MatrixXd sortedEigenvectorsH(eigenvectorsH.rows(), eigenvaluesH.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;

    // ORDINA LO SPETTRO ENERGETICO E CONTEMPORANEAMENTE RIORDINA LA BASE;
    ORDINE_ASC(eigenvaluesH,eigenvectorsH,sortedEigenvaluesH,sortedEigenvectorsH);
    gs = eigenvectorsH.col(0);
    return(sortedEigenvaluesH(1)-sortedEigenvaluesH(0));

}


void GROUND_STATE(SparseMatrix<double>& H,VectorXd& gs,double g,double hI,SparseMatrix<int>& base)
{
    hamiltoniana(H, base, g, hI);
            
    SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(H);

    VectorXd eigenvaluesH = eigensolver.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
    MatrixXd eigenvectorsH = eigensolver.eigenvectors(); // MARICE DEGLI AUTOVETTORI;
    
    // ORGANIZZAZIONE DEI DATI;

    VectorXd sortedEigenvaluesH(eigenvaluesH.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
    MatrixXd sortedEigenvectorsH(eigenvectorsH.rows(), eigenvaluesH.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;

    // ORDINA LO SPETTRO ENERGETICO E CONTEMPORANEAMENTE RIORDINA LA BASE;
    ORDINE_ASC(eigenvaluesH,eigenvectorsH,sortedEigenvaluesH,sortedEigenvectorsH);

    gs=eigenvectorsH.col(0);
    
}

int main()
{
    SparseMatrix<int> base(rows, N);
	SparseMatrix<double> H0(rows, rows);
    VectorXd gs(rows);
	costruisciBase(base);

    double g = 0.5, h = -0.020;
    
    double gap = GAP(H0,gs,g,0.0,base);
    double K = Mz(gs,base)*N/gap;
    
    // CREAZIONE DEL FILE DI OUTPUT CONTENENTE VALORI DI h*K E Mz
    ofstream solutions("N" + to_string(N) + "/primaSp" + to_string(N) + ".txt");
	if (!solutions.is_open())
    {
        cerr << "Impossibile aprire il file di output." << endl;
        return 1; // Esci con un codice di errore
    }
    
    while(h <= 0.021)
    {
	    SparseMatrix<double> H(rows, rows);
        GROUND_STATE(H,gs,g,h,base);
        solutions << -h*K << " " << Mz(gs,base) << endl;
        h += 0.001;
    }
    
    solutions.close();
    cout << "Ho compilato il file di output primSp.txt" << endl;
    return 0;
}