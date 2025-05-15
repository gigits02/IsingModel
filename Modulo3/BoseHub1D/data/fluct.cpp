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
#include <vector>
#include <sstream>
#include "../eigen/Eigen/Sparse"
#include "../eigen/Eigen/Dense"
#include "../eigen/Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;


/*
    COMPUTING "ORDER PARAMETER" (FLUCTUATION OF TOTAL NUMBER OF BOSONS) FOR DIFFERENT J.
    AND ITS SCALING WITH M (NUMBER OF TOTAL SITES CONSIDERED IN THE SUM)
    FOR BOSE-HUBBARD MODEL IN 1D.

    H=-J*\sum_{i}[b_{i}^{dagger}b_{i+1} + h.c]+U/2*\sum_{i}[n_{i}(n_{i}-1)-\mu*\sum_{i}n_{i}];

*/


const int ITER = 20; // MAX NUMBER OF ITERATIONS;
int H_SPACE_DIM = 6; // HIGHEST DIMENTION OF HILBERT SPACE USABLE;


/*

    CREAZIONE DI STRUTTURE DASTI E VARIABILI DI AMBIENTE PER L'ESECUZIONE DELLO SWEEPING;

*/

struct BlockStructure 
{
    SparseMatrix<double> b;
    SparseMatrix<double> b_dagger;
    SparseMatrix<double> n;
    SparseMatrix<double> I;
    SparseMatrix<double> H;
};

vector<BlockStructure> history_left;
vector<BlockStructure> history_right;

BlockStructure initializeBlockStructure(SparseMatrix<double>& b, SparseMatrix<double>& b_dagger, SparseMatrix<double>& n, SparseMatrix<double>& I, SparseMatrix<double>& H)
{
    BlockStructure STRUCTURE;
    STRUCTURE.b = b;
    STRUCTURE.b_dagger = b_dagger;
    STRUCTURE.n = n;
    STRUCTURE.I = I;
    STRUCTURE.H = H;
    return STRUCTURE;
}

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
void LadderOps( SparseMatrix<double>& A, SparseMatrix<double>& B, SparseMatrix<double>& C,SparseMatrix<double>& I, int Nmax)
{
    for(int i=0;i < Nmax+1;i++)
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
    int size = static_cast<int>(sqrt(GS.rows()));
    MatrixXd PsiT(size,size);  // IL METODO RESTITUISCE IL TRASPOSTO DI QUELLO CHE VORREMMO;
    MatrixXd Psi(size,size);
    PsiT = Map<MatrixXd>(GS.data(),size,size);
    Psi = PsiT.transpose();
    Ro = Psi*PsiT;
}

// IMPLEMENTA L'ALGORITMO DEL DMRG E SALVA VALORI DI MU CRITICI;
void DMRG(BlockStructure& STRUCTURE, SparseMatrix<double>& b, SparseMatrix<double>& b_dagger, SparseMatrix<double>& n, SparseMatrix<double>& I, double mu, double J, int M)
{  

    // APERTURA FILE;
    ifstream solutions("fluct" + to_string(ITER - M) + ".txt");
    
    // Leggi le vecchie soluzioni se presenti
    vector<string> oldSolutions;
    string line;
    if (solutions.is_open())
    {
        while(getline(solutions, line))
        {
            oldSolutions.push_back(line);
        }
        solutions.close();
    }

    // Scrivi le soluzioni su file
    ofstream outputFile("fluct" + to_string(ITER-M) + ".txt");
    if(outputFile.is_open())
    {
        // Scrivi le vecchie soluzioni
        for (const auto& solution : oldSolutions)
        {
            outputFile << solution << '\n';
        }
    }

    SparseMatrix<double>Hsuper;
    SparseMatrix<double>HsuperT;
    double Energy = 0.0;
    double N_M = 0.0;
    SparseMatrix<double> Nmedio;
    for(int i=1;i<ITER+1;i++)
    {    
        cout << i << "°: effettuando il DMRG ... " << endl;
        // CREAZIONE OPERATORI PER CALCOLO DELLE FLUTTUAZIONI DI n (PARAMETRO D'ORDINE)
        //SARÒ INTERESSATO AL NUMERO DI BOSONI DEL SITO J E VOGLIO SOMMARE FINO A M
        if( i == M)    
            Nmedio = kron(STRUCTURE.I, n);
        
        if( i > M && i <= ITER)
            Nmedio = kron(Nmedio, I) + kron(STRUCTURE.I, n);
        
        if( i == ITER)
            Nmedio = kron(Nmedio, kron(I, STRUCTURE.I));

        // COSTRUZIONE DEL BLOCCO SINISTRO i.e A--o;
        STRUCTURE.H = kron(STRUCTURE.H,I)-J*( kron(STRUCTURE.b,b_dagger) + kron(STRUCTURE.b_dagger,b) )+0.5*( kron( (STRUCTURE.n*(STRUCTURE.n-STRUCTURE.I)), I )+kron( STRUCTURE.I, (n*(n-I)) ) )-mu*(kron(STRUCTURE.n, I)+kron(STRUCTURE.I, n));
        STRUCTURE.b = kron(STRUCTURE.I,b);
        STRUCTURE.b_dagger = kron(STRUCTURE.I,b_dagger);
        STRUCTURE.n = kron(STRUCTURE.I,n);
        STRUCTURE.I = kron(STRUCTURE.I,I);
        // COSTRUZIONE HAMILTONIANO SUPERBLOCCO i.e. A--o--o--B;
        Hsuper = kron(STRUCTURE.H,STRUCTURE.I)+kron(STRUCTURE.I,STRUCTURE.H)-J*(kron(STRUCTURE.b_dagger,STRUCTURE.b)+kron(STRUCTURE.b,STRUCTURE.b_dagger))+0.5*(kron((STRUCTURE.n*(STRUCTURE.n-STRUCTURE.I)),STRUCTURE.I)+kron(STRUCTURE.I,(STRUCTURE.n*(STRUCTURE.n-STRUCTURE.I))))-mu*(kron(STRUCTURE.n,STRUCTURE.I)+kron(STRUCTURE.I,STRUCTURE.n));
        // CHECK HERMITIANITA';
        HsuperT = Hsuper.transpose();
        Hsuper = 0.5*(Hsuper+HsuperT);
        // DIAGONALIZZAZIONE DI H;
        SelfAdjointEigenSolver<SparseMatrix<double>> eigensolver(Hsuper); //DEFINIZIONE DELL'OGETTO PER LA DIAGONALIZZAZIONE;

       
        if (eigensolver.info() != Success)
        {
            cout << "Errore durante la diagonalizzazione!" <<endl; //CONTROLLO SU POSSIBILI ERRORI;
            exit(1);
        }

        VectorXd eigenvaluesH = eigensolver.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
        MatrixXd eigenvectorsH = eigensolver.eigenvectors(); // METRICE DEGLI AUTOVETTORI;
        
        // ORGANIZZAZIONE DEI DATI;

        VectorXd sortedEigenvaluesH(eigenvaluesH.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
        MatrixXd sortedEigenvectorsH(eigenvectorsH.rows(), eigenvaluesH.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;
        ORDINE_ASC(eigenvaluesH,eigenvectorsH,sortedEigenvaluesH,sortedEigenvectorsH);


        // UPDATE DELLE OSSERVABILI DEL SISTEMA;
        Energy = sortedEigenvaluesH(0);

        // TRACCIA PARZIALE TRAMITE RESHAPING E CALCOLO ULTERIORI OSSERVABILI
        VectorXd GS = sortedEigenvectorsH.col(0);
        MatrixXd Ro(static_cast<int>(sqrt(sortedEigenvaluesH.size())),static_cast<int>(sqrt(sortedEigenvaluesH.size())));
        costruisciRo(GS,Ro);
        //n_medio = (Ro*NmedioA).trace();

        // DIAGONALIZZO RO ED ESEGUO IL TRONCAMENTO DELLO SPAZIO DI HILBERT;
        SelfAdjointEigenSolver<SparseMatrix<double>> eigensolverRo(Ro); //DEFINIZIONE DELL'OGETTO PER LA DIAGONALIZZAZIONE;
        VectorXd eigenvaluesRO = eigensolverRo.eigenvalues(); //MATRICE DEGLI AUTOVALORI;
        MatrixXd eigenvectorsRO = eigensolverRo.eigenvectors(); // METRICE DEGLI AUTOVETTORI;

        // ORGANIZZAZIONE DEI DATI (Parte 2);
        VectorXd sortedEigenvaluesRO(eigenvaluesRO.size()); // CREAZIONE DEL VETTORE DI AUTOVALORI ORDINATI;
        MatrixXd sortedEigenvectorsRO(eigenvectorsRO.rows(), eigenvaluesRO.size()); // CREAZIONE DELLA MATRICE DI AUTOVETTORI ORDINATI;
        ORDINE_DISC(eigenvaluesRO,eigenvectorsRO,sortedEigenvaluesRO,sortedEigenvectorsRO);

        // RIDUZIONE DELLO SPAZIO DI HILBERT;
        int n_keep= min(static_cast<int>(sortedEigenvaluesRO.size()),H_SPACE_DIM); // SELEZIONA IL MINIMO TRA IL NUMERO DI AUTOVALORI == DIAG RO E LA DIMENSIONE INIZIALE H_SPACE_DIM;
        MatrixXd O = sortedEigenvectorsRO.block(0,0,sortedEigenvectorsRO.rows(),n_keep); // SCREMA LA MATRICE DI AUTOVETTORI SELEZIONANDO I PRIMI n_keep (N.B.: HO UNA MATRICE RETTANGOLARE);
        MatrixXd OT = O.transpose();

        if( i == ITER)
        {
            N_M = GS.dot(Nmedio*GS);
            double Nsq_medio = GS.dot(Nmedio*Nmedio*GS);
            double deltan = Nsq_medio - (N_M*N_M);
            outputFile << J << " " << deltan << endl;
        }

        // MEMORIZZAZIONE STRUTTURA (SERVIRÀ PER LO SWEEPING)
        history_right.push_back({STRUCTURE});

        // RIDUZIONE DELLE MATRICI;
        STRUCTURE.H = (OT*STRUCTURE.H*O).sparseView();
        STRUCTURE.b = (OT*STRUCTURE.b*O).sparseView();
        STRUCTURE.b_dagger = (OT*STRUCTURE.b_dagger*O).sparseView();
        STRUCTURE.n = (OT*STRUCTURE.n*O).sparseView();
        STRUCTURE.I = (OT*STRUCTURE.I*O).sparseView();
        
        //MI SERVE TRONCARE ANCHE QUESTA OSSERVABILE PER PASSARE AL SITO SUCCESSIVO
        if( i >= M && i < ITER)
            Nmedio = (OT*Nmedio*O).sparseView();
    }
    outputFile.close();
}

void filterMonotonic(const string& filename)
{
    ifstream infile(filename);
    if (!infile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    vector<pair<double, double>> data;
    string line;
    bool startFound = false;
    double lastY = -1;

    // Read the file line by line
    while (getline(infile, line)) {
        istringstream iss(line);
        double x, y;
        if (!(iss >> x >> y)) { break; }

        if (startFound) {
            if (y >= lastY) {
                data.emplace_back(x, y);
                lastY = y;
            }
        } else {
            data.emplace_back(x, y);
            if (abs(x-0.03) <= 0.01) {
                startFound = true;
                lastY = y;
            }
        }
    }

    infile.close();

    // Overwrite the original file with filtered data
    ofstream outfile(filename);
    if (!outfile.is_open()) {
        cerr << "Error opening file: " << filename << endl;
        return;
    }

    for (const auto& dp : data) {
        outfile << dp.first << " " << dp.second << endl;
    }

    outfile.close();
}
int main()
{

    int Nmax = 4; // FILLING;

    // HAMILTONIAN PARAMETERS (U WAS SET TO 1);
    double J = 0.2; 
    double mu = 0.5; // CHEMICAL POTENTIAL;

    //INIZIALIZZAZIONE OPERATORI ON-SITE;
    SparseMatrix<double> b(Nmax+1,Nmax+1);
    SparseMatrix<double> b_dagger(Nmax+1,Nmax+1);
    SparseMatrix<double> n(Nmax+1,Nmax+1);
    SparseMatrix<double> I(Nmax+1,Nmax+1);
    SparseMatrix<double> H(Nmax+1,Nmax+1);
    H.setZero();
    LadderOps(b_dagger,b,n,I, Nmax);

    // INIZIALIZZAZIONE BLOCCHI ASSUMENTO A--a--b--B CON b--B SPECULARE AD A--a;
    BlockStructure STRUCTURE = initializeBlockStructure(b,b_dagger,n,I,H);

    //NUMERO DI SITI CONSIDERATI NELLA SOMMA PER IL CALCOLO DELLE FLUTTUAZIONI DI BOSONI
    
    int M = ITER-1;
    while( M >= 0)
    {
        J = 0.0;
        
        while(J < 0.41)
        {
            cout << "Simulazione per Nmax = " << Nmax <<", J = " << J << ", mu = " << mu << ", M = " << M << endl;
            DMRG(STRUCTURE, b, b_dagger, n, I, mu, J, M);     
            J += 0.02;
        }
                 
        cout <<"Ho salvato nel file di output ... " << endl;
        cout << "J  deltan" << endl;
        //string filename = "fluct" + to_string(ITER-M) + ".txt";
        //filterMonotonic(filename);
        M = M - 3;
    }
    return 0;

}
