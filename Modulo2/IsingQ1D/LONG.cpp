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

void GROUND_STATE(SparseMatrix<double>&  H,VectorXd& gs,double g,double hI,SparseMatrix<int>& base)
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



//CALCOLA I VALORI DELLA SUSCETTIVITÀ (LONGITUDINALE) AL VARIARE DI g (per h piccolo) E AL VARIARE DI h (per g = 0.40)

int main()
{

	
    SparseMatrix<int> base(rows,N);
    VectorXd gs(rows);

	costruisciBase(base);
    double g=0.05,hI=-0.02,hF=0.02;
    const int DIM = 50;
    double dh= static_cast<double>((hF-hI)/(DIM+2));
    double mzDef[DIM+2];
    double mz[DIM+2];
    double chiZDef[DIM];
    double chiZ[DIM];

    // CREAZIONE DEL FILE DI OUTPUT
    /*ofstream soluzionig("dinamica/suscLong" + to_string(N) + ".txt");
	if (!soluzionig.is_open())
    {
        cerr << "Impossibile aprire uno o entrambi i file di output." << endl;
        return 1; // Esci con un codice di errore
    }*/
    ofstream soluzioniG("dinamica/suscLongDEF" + to_string(N) + ".txt");
	if (!soluzioniG.is_open())
    {
        cerr << "Impossibile aprire uno o entrambi i file di output." << endl;
        return 1; // Esci con un codice di errore
    }



    for(int i=0;i<40;i++)
    {
        // CHECK 2
        cout<<i<<endl;


        //soluzionig<<g<<" ";
        soluzioniG<<g<<" ";
        hI=-0.02;

        // COMPILAZIONE DEL VETTORE DELLE MAGNETIZZAZIONI AL VARIARE DI hI
        for(int j=0;j<DIM+2;j++)
        {
        	SparseMatrix<double> H(rows, rows);
            GROUND_STATE(H,gs,g,hI,base);
            mzDef[j]=MzDef(gs,base);
            mz[j]=Mz(gs,base);
            hI += dh;
        }

        // CALCOLO CON UNO SCHEMA ALLE DIFFERENZE FINITE DELLA SUSCETTIVITÀ
        for(int k=0;k<DIM;k++)
        {   
            chiZ[k]=(mz[k+2]-mz[k])/(2*dh); 
            chiZDef[k]=(mzDef[k+2]-mzDef[k])/(2*dh); 
            //soluzionig<<chiZ[k]<<" ";
            soluzioniG<<chiZDef[k]<<" ";
        }
        //soluzionig << endl;
        soluzioniG << endl;
        g += 0.05;
    }
    //soluzionig.close();
    soluzioniG.close();
    return 0;   
}