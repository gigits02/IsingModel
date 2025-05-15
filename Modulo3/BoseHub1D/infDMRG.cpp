#include <iostream>
#include <cmath>
#include <fstream>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <stdint.h>
#include <iomanip>
#include <random>
#include <complex>
#include <algorithm>
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;


/*
    INFINITE-SYSTEM - DMRG METHODS FOR BOSE-HUBBARD MODEL IN 1D;

    H=-J*\sum_{i}[b_{i}^{dagger}b_{i+1} + h.c]+U/2*\sum_{i}[n_{i}(n_{i}-1)-\mu*\sum_{i}n_{i}];

*/

const int Nmax = 4; // HALF FILLING;
const int ITER = 20; // MAX NUMBER OF ITERATIONS;
int H_SPACE_DIM = 5; // HIGHEST DIMENTION OF HILBERT SPACE USABLE;

// HAMILTONIAN PARAMETERS SCALED BY U WHICH WAS SET TO 1;
const double J = 0.05; 
const double mu = 0.5;

/*

    DICHIARAZIONE E DEFINIZIONE DELLE FUNZIONI CHE SVILUPPANO PRODOTTI E RIDUZIONI MATRICIALI;

*/

// CALCOLA IL PRODOTTO DI KRONECKER PER MATRICI SPARSE;
SparseMatrix<double> kron(const SparseMatrix<double>& A, const SparseMatrix<double>& B)
{
    // Calcolo delle dimensioni della matrice risultante
    const int rows_result = A.rows() * B.rows();
    const int cols_result = A.cols() * B.cols();

    // Costruzione della matrice risultante come matrice sparsa
    SparseMatrix<double> result(rows_result, cols_result);

    // Iterazione sui coefficienti non nulli della matrice A
    for (int kA = 0; kA < A.outerSize(); ++kA) {
        for (SparseMatrix<double>::InnerIterator itA(A, kA); itA; ++itA) {
            // Iterazione sui coefficienti non nulli della matrice B
            for (int kB = 0; kB < B.outerSize(); ++kB) {
                for (SparseMatrix<double>::InnerIterator itB(B, kB); itB; ++itB) {
                    // Calcolo dell'indice della riga e della colonna nel risultato
                    int i_result = itA.row() * B.rows() + itB.row();
                    int j_result = itA.col() * B.cols() + itB.col();
                    // Inserimento del prodotto del coefficiente di A e di B nella posizione corrispondente nella matrice risultante
                    result.insert(i_result, j_result) = itA.value() * itB.value();
                }
            }
        }
    }

    // Finalizzazione della costruzione della matrice risultante
    result.makeCompressed();

    return result;
}

// DEFINIZIONE DEGLI OPERATORI DI SALITA E DISCESA BOSONICI E DUNQUE L'OPERATORE NUMERO A DATO NUMERO DI OCCUPAZIONE MASSIMO (ON SITE) n;
void LadderOps( SparseMatrix<double>& A, SparseMatrix<double>& B, SparseMatrix<double>& C,SparseMatrix<double>& I)
{
    for(int i=0;i<Nmax+1;i++)
    {
        if(i != Nmax)
        {
            A.insert(i+1,i)=(sqrt(i+1));
            B.insert(i,i+1)=(sqrt(i+1));
        }
        
            C.insert(i,i)=i;
            I.insert(i,i)=1;
    }

    A.makeCompressed();
    B.makeCompressed();
    C.makeCompressed();
    I.makeCompressed();
}

// FUNZIONE CHE RIORDINA GLI AUTOVETTORI E AUTOVALORI IN ORDINE CRESCENTE IN ENERGIA;
void ORDINE_ASC(VectorXd& eigenvalues, MatrixXd& eigenvectors,VectorXd& sortedEigenvalues, MatrixXd& sortedEigenvectors)
{
    // Creazione di un vettore di indici per tenere traccia dell'ordine originale degli autovalori
    vector<int> indices(eigenvalues.size());
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        indices[i] = i;
    }

    // Ordinamento degli indici in base agli autovalori corrispondenti
    sort(indices.begin(), indices.end(), [&](int i, int j) { return eigenvalues[i] < eigenvalues[j]; });

    
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        sortedEigenvalues(i) = eigenvalues[indices[i]];
        sortedEigenvectors.col(i) = eigenvectors.col(indices[i]);
    }

}

void ORDINE_DISC(VectorXd& eigenvalues, MatrixXd& eigenvectors,VectorXd& sortedEigenvalues, MatrixXd& sortedEigenvectors)
{
    // Creazione di un vettore di indici per tenere traccia dell'ordine originale degli autovalori
    vector<int> indices(eigenvalues.size());
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        indices[i] = i;
    }

    // Ordinamento degli indici in base agli autovalori corrispondenti
    sort(indices.begin(), indices.end(), [&](int i, int j) { return eigenvalues[i] > eigenvalues[j]; });

    
    for (int i = 0; i < eigenvalues.size(); ++i)
    {
        sortedEigenvalues(i) = eigenvalues[indices[i]];
        sortedEigenvectors.col(i) = eigenvectors.col(indices[i]);
    }

}

// FUNZIONE CHE IMPLEMENTA IL RESHAPE E DA LA MATRICE DENSITA';
void costruisciRo(VectorXd GS, MatrixXd& Ro)
{
    //cout<<"sono dentro la funzione che crea ro"<<endl;
    int size = static_cast<int>(sqrt(GS.rows()));
    //cout << "SIZE DI PSI:" << size << endl;
    MatrixXd PsiT(size,size);  // IL METODO RESTITUISCE IL TRASPOSTO DI QUELLO CHE VORREMMO;
    MatrixXd Psi(size,size);
    PsiT = Map<MatrixXd>(GS.data(),size,size);
    Psi = PsiT.transpose();
    Ro = Psi*PsiT;
    //cout << "RIGHE DI RO:" << Ro.rows() << endl;
    //cout << "COLONNE DI RO:" << Ro.cols() << endl;
}

// IMPLEMENTA L'ALGORITMO DEL NRG;
void DMRG(SparseMatrix<double>& Block_b, SparseMatrix<double>& Block_b_dagger, SparseMatrix<double>& Block_n, SparseMatrix<double>& Block_I, SparseMatrix<double>& Block_H, SparseMatrix<double>& b, SparseMatrix<double>& b_dagger, SparseMatrix<double>& n, SparseMatrix<double>& I)
{   
    // APERTURA FILE;
    ofstream solutions("data/infMI2.txt");
    SparseMatrix<double>Hsuper;
    SparseMatrix<double>HsuperT;
    int LENGHT=0;
    double Energy=0.0;
    double LastEnergy=0.0;
    double En_Density=0.0;
    double En_Bond=0.0;
    for(int i=1;i<ITER+1;i++)
    {
        LENGHT= 2*i+2;
        //cout<<"SONO NEL CICLO FOR"<<endl;
        // COSTRUZIONE DEL BLOCCO SINISTRO i.e A--o;
        Block_H = kron(Block_H,I)-J*( kron(Block_b,b_dagger) + kron(Block_b_dagger,b) )+0.5*( kron( (Block_n*(Block_n-Block_I)), I )+kron( Block_I, (n*(n-I)) ) )-mu*(kron(Block_n, I)+kron(Block_I, n));
        Block_b = kron(Block_I,b);
        Block_b_dagger = kron(Block_I,b_dagger);
        Block_n = kron(Block_I,n);
        Block_I = kron(Block_I,I);
        // COSTRUZIONE HAMILTONIANO SUPERBLOCCO i.e. A--o--o--B;
        Hsuper= kron(Block_H,Block_I)+kron(Block_I,Block_H)-J*(kron(Block_b_dagger,Block_b)+kron(Block_b,Block_b_dagger))+0.5*(kron((Block_n*(Block_n-Block_I)),Block_I)+kron(Block_I,(Block_n*(Block_n-Block_I))))-mu*(kron(Block_n,Block_I)+kron(Block_I,Block_n));
        // CHECK HERMITIANITA';
        HsuperT = Hsuper.transpose();
        Hsuper = 0.5*(Hsuper+HsuperT);
        LastEnergy = Energy;
        //cout<<"HO CREATO HSUPER"<<endl;
        //cout<<"Numero colonne hsuper:"<<Hsuper.cols()<<endl;
        //cout<<"Numero righe hsuper:"<<Hsuper.rows()<<endl;
        // DIAGONALIZZAZIONE DI H;
        SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(Hsuper); //DEFINIZIONE DELL'OGETTO PER LA DIAGONALIZZAZIONE;

       
        if (eigensolver.info() != Success)
        {
            cout<< "Errore durante la diagonalizzazione!" <<endl; //CONTROLLO SU POSSIBILI ERRORI;
            exit(1);
        }

        VectorXd eigenvaluesH = eigensolver.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
        MatrixXd eigenvectorsH = eigensolver.eigenvectors(); // METRICE DEGLI AUTOVETTORI;
        
        // ORGANIZZAZIONE DEI DATI;

        VectorXd sortedEigenvaluesH(eigenvaluesH.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
        MatrixXd sortedEigenvectorsH(eigenvectorsH.rows(), eigenvaluesH.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;
        ORDINE_ASC(eigenvaluesH,eigenvectorsH,sortedEigenvaluesH,sortedEigenvectorsH);
        //cout<<"HO ORDINATO I DATI"<<endl;
        //cout<<"Numero colonne autovettori:"<<sortedEigenvectorsH.cols()<<endl;
        //cout<<"Numero righe autovettori:"<<sortedEigenvectorsH.rows()<<endl;


        // UPDATE DELLE OSSERVABILI DEL SISTEMA;
        Energy=sortedEigenvaluesH(0);
        En_Density = Energy/LENGHT;
        En_Bond=(Energy-LastEnergy)/2;

        // TRCCIA PARZIALE TRAMITE RESHAPING;
        VectorXd GS = sortedEigenvectorsH.col(0);
        //cout<<"GS:"<<GS.rows()<<endl;
        MatrixXd Ro(static_cast<int>(sqrt(sortedEigenvaluesH.size())),static_cast<int>(sqrt(sortedEigenvaluesH.size())));
        costruisciRo(GS,Ro);

        // DIAGONALIZZO RO ED ESEGUO IL TRONCAMENTO DELLO SPAZIO DI HILBERT;
        SelfAdjointEigenSolver<SparseMatrix<double>> eigensolverRo(Ro); //DEFINIZIONE DELL'OGETTO PER LA DIAGONALIZZAZIONE;
        VectorXd eigenvaluesRO = eigensolverRo.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
        MatrixXd eigenvectorsRO = eigensolverRo.eigenvectors(); // METRICE DEGLI AUTOVETTORI;
        //cout << "HO DIAGONALIZZATO RO" << endl;
        // ORGANIZZAZIONE DEI DATI (Parte 2);

        VectorXd sortedEigenvaluesRO(eigenvaluesRO.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
        MatrixXd sortedEigenvectorsRO(eigenvectorsRO.rows(), eigenvaluesRO.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;
        ORDINE_DISC(eigenvaluesRO,eigenvectorsRO,sortedEigenvaluesRO,sortedEigenvectorsRO);
        //cout << "HO ORDINATO I DATI PER RO" << endl;

        // RIDUZIONE DELLO SPAZIO DI HILBERT;
        int n_keep= min(static_cast<int>(sortedEigenvaluesRO.size()),H_SPACE_DIM); // SELEZIONA IL MINIMO TRA IL NUMERO DI AUTOVALORI == DIAG RO E LA DIMENSIONE INIZIALE H_SPACE_DIM;
        //cout << "IL MINIMO È :" << n_keep << endl;
        MatrixXd O = sortedEigenvectorsRO.block(0,0,sortedEigenvectorsRO.rows(),n_keep); // SCREMA LA MATRICE DI AUTOVETTORI SELEZIONANDO I PRIMI n_keep (N.B.: HO UNA MATRICE RETTANGOLARE);
        MatrixXd OT = O.transpose();
        //cout << "MATRICE O DI DIMENSIONI:" << O.rows() << "X" << O.cols() << endl;
        //cout << "MATRICE OTRASPOSTA DI DIMENSIONI" << OT.rows() << "X" << OT.cols() << endl;

        // CLACOLO DELL'ERRORE DI TRONCAMENTO ED ENTROPIA DI VON NEUMANN
        double Tronc_error = 1 - sortedEigenvaluesRO.head(n_keep).sum();
        double ent=0;
        for(int i=0;i<n_keep;i++)
        {
            ent += -sortedEigenvaluesRO(i)*log2(sortedEigenvaluesRO(i));
        }

        solutions<<LENGHT<<" "<<Energy<<" "<<En_Bond<<" "<<En_Density<<" "<<Tronc_error<<" "<<ent<<endl;

        // RIDUZIONE DELLE MATRICI;
        Block_H = (OT*Block_H*O).sparseView();
        Block_b = (OT*Block_b*O).sparseView();
        Block_b_dagger = (OT*Block_b_dagger*O).sparseView();
        Block_n = (OT*Block_n*O).sparseView();
        Block_I = (OT*Block_I*O).sparseView();

    }
    solutions.close();
}

int main()
{
    cout << "Simulazione per J = " << J << ", mu = " << mu << endl;
    cout <<"Sto salvando nel file di output ... " << endl;
    cout << "LENGHT" << " " << "Energy" << " " << "En_Bond" << " " << "En_Density" << " " << "TruncErr" << " " << "Ent" << endl;
    //INIZIALIZZAZIONE OPERATORI ON-SITE;
    SparseMatrix<double> b(Nmax+1,Nmax+1);
    SparseMatrix<double> b_dagger(Nmax+1,Nmax+1);
    SparseMatrix<double> n(Nmax+1,Nmax+1);
    SparseMatrix<double> I(Nmax+1,Nmax+1);
    SparseMatrix<double> Block_H(Nmax+1,Nmax+1);
    Block_H.setZero();
    
    LadderOps(b_dagger,b,n,I);

    //INIZIALIZZAZIONE BLOCCHI ASSUMENDO A--o--o--B CON o--B SPECULARE AD A--o;
    SparseMatrix<double> Block_b=b;
    SparseMatrix<double> Block_b_dagger=b_dagger;
    SparseMatrix<double> Block_n=n;
    SparseMatrix<double> Block_I=I;

    DMRG(Block_b,Block_b_dagger,Block_n,Block_I,Block_H,b,b_dagger,n,I);
    return 0;
}