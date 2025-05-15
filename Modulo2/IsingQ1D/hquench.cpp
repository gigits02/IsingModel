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
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Eigenvalues"
#include<complex>
#include "mylib.h"

using namespace std;
using namespace Eigen;

using namespace std;
int main()
{
    SparseMatrix<int> base(rows,N);

	// COSTRUZIONE DELLA BASE COMPUTAZIONALE;	
	costruisciBase(base);
	double g, h;
	// DICHIARAZIONE DEI PARAMETRI HAMILTONIANI g E h;
    //g = static_cast<double>(1.0-1./N);
    g = 1.0;
    h = -pow(N, -15./8.);

    SparseMatrix<double> H(rows,rows);
    hamiltoniana(H, base, g, h);
    VectorXd gs(rows);
    //PARTIAMO DAL GROUND STATE CON I VALORI DEI PARAMETRI INIZIALI INSERITI
    SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(H);
    VectorXd eigenvaluesH = eigensolver.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
    MatrixXd eigenvectorsH = eigensolver.eigenvectors(); // MARICE DEGLI AUTOVETTORI;
    VectorXd sortedEigenvaluesH(eigenvaluesH.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
    MatrixXd sortedEigenvectorsH(eigenvectorsH.rows(), eigenvaluesH.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;
    // ORDINA LO SPETTRO ENERGETICO E CONTEMPORANEAMENTE RIORDINA LA BASE;
    ORDINE_ASC(eigenvaluesH,eigenvectorsH,sortedEigenvaluesH,sortedEigenvectorsH);
    
    for(int i = 0; i < rows; i++)
        gs[i] = sortedEigenvectorsH(i,0);


	//DINAMICA (quench su h = 0.0)
    double quench = 0.0;
	double t = 0.0;
	double dt = 0.01;   
    VectorXcd z(rows);
    z = gs;
    SparseMatrix<double> Hquench(rows,rows);
    hamiltoniana(Hquench, base, g, quench);
    ofstream solutions("dinamica/hdinamica" + to_string(N) + ".txt");
    if (!solutions.is_open())
    {
        cerr << "Impossibile aprire il file di output." << endl;
        return 1; // Esci con un codice di errore
    }
	
	for (int i = 0; i < 30000; ++i)
	{
        cout << "iterazione " << i << "° ..." << endl;
	    solutions << t << " " << Mz_complex(z, base) << endl;
		rungeKutta(z, Hquench, t, dt);
		t += dt;
    }    

    solutions.close();
    
    cout << "Le soluzioni sono state salvate correttamente nel file di output nel seguente ordine:" << endl;
    cout << "t   mz(t)  " << endl;

    /*
    //CALCOLA LE DERIVATE PARZIALI PER TROVARE SUSCETTIVITÀ LONG E TRASV
    ofstream susc("zdinamica/suscTrasv" + to_string(N) + ".txt");

    if (!susc.is_open())
    {
        cerr << "Impossibile aprire il file di output." << endl;
        return 1; // Esci con un codice di errore
    }
    for(int i = 1; i < dim-1; i++)
    {    
        VectorXd GSsuc(rows);
        VectorXd GSprec(rows);
        double partial = 0.0;
        GSsuc=GS.row(i+1);
        GSprec=GS.row(i-1);

        for(int j = 0; j < rows; j++)
        {
            
            // Calcola il valore di y in t + deltat
            double Mx_plus_delta = Mx(GSsuc, base);
    
            // Calcola il valore di mx in g - delta_g
            double Mx_minus_delta = Mx(GSprec, base);
    
            // Calcola la derivata parziale utilizzando la formula delle differenze finite
            partial = (Mx_plus_delta - Mx_minus_delta) / (2 * 0.05);

        }

        susc << 0.05*i << " " << partial << endl;
    }

    cout << "Le soluzioni per suscettività sono state salvate correttamente secondo l'ordine:" << endl;
    cout << "g  suscTr" << endl;

    susc.close();
    */
    return 0;
}