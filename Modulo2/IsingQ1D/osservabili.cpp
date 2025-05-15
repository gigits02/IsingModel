#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <string>
#include <complex>
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Eigenvalues"
#include "mylib.h"

using namespace std;
using namespace Eigen;

void leggiFile(const string& nomeFile, VectorXd& primaColonna, vector<vector<double>>& dati) {
    ifstream file(nomeFile);
    if (!file.is_open()) {
        cerr << "Impossibile aprire il file " << nomeFile << endl;
        exit(EXIT_FAILURE);
    }

    primaColonna.resize(0); // Svuota il vettore della prima colonna
    dati.clear();

    string riga;
    while(getline(file, riga))
    {
        stringstream ss(riga);
        double valore;
        ss >> valore; // Leggi il valore della prima colonna
        primaColonna.conservativeResize(primaColonna.size() + 1);
        primaColonna(primaColonna.size() - 1) = valore;

        vector<double> rigaDati;
        while (ss >> valore)
        {
            rigaDati.push_back(valore);
        }
        dati.push_back(rigaDati);
    }

    file.close();
}

void costruisciMatrici(const string& fileGS, const string& fileE, MatrixXd& GS, MatrixXd& E, VectorXd& first_col) {
    vector<vector<double>> datiGS, datiE;

    leggiFile(fileGS, first_col, datiGS);
    leggiFile(fileE, first_col, datiE);

    // Inserisci i dati in una matrice
    int numRighe = datiGS.size();
    int numColonne = datiGS.empty() ? 0 : datiGS[0].size(); // Numero di colonne dei dati
    GS.resize(numRighe, numColonne);
    E.resize(numRighe, numColonne);
    for (int i = 0; i < numRighe; ++i)
    {
        for (int j = 0; j < numColonne; ++j)
        {
            GS(i, j) = datiGS[i][j];
            E(i, j) = datiE[i][j];
        }
    }
}


int main() 
{
    SparseMatrix<int> base(rows,N);
    costruisciBase(base);

    MatrixXd GROUND_STATE;
    MatrixXd E;
    VectorXd g;
    costruisciMatrici("N" + to_string(N) + "/h0.txt", "N" + to_string(N) + "/eh0.txt", GROUND_STATE, E,g);

    ofstream solutions("N" + to_string(N) + "/0solutions.txt");

    if (!solutions.is_open()) {
        cerr << "Impossibile aprire il file di output." << endl;
        return 1; // Esci con un codice di errore
    }
    
   
    
    for(int i=0 ; i<g.size();i++ )
    {
        VectorXd Psi = GROUND_STATE.row(i);
        double U4 = M4(Psi,base)/(M2(Psi,base)*M2(Psi,base));
        solutions << Mx(Psi,base) << " " << MzDef(Psi,base) << " " << E(i,1)-E(i,0) << " " << E(i,2)-E(i,0) << " " << g(i) <<" "<<U4<<endl;
    }


    cout << "Le soluzioni sono state correttamente salvate secondo l'ordine seguente:" << endl;
    cout << "mx mz (E1-E0) (E2-E1) g U4" << endl;

    solutions.close();

    return 0;
}