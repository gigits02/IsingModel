#include <iostream>
#include <armadillo>
#include <complex>

using namespace std;
using namespace arma;

#define ROWS 3
#define COLS 3

void stampaMatrice( complex<double> M[][COLS] )
{							
	for( int i = 0; i < ROWS; i++ )
	{	
		for( int j = 0; j < COLS; j++ )
		{	
			cout << M[i][j] << " ";
		}
		cout << endl;
	}	

}

void stampaVettore( complex<double> vet[ROWS])
{
    for(int i = 0; i < ROWS; i++)
        cout << vet[i] << " ";  
    
    cout << endl;
}

int main()
{

    complex<double> M[ROWS][COLS] = {0};
    stampaMatrice(M);
    cx_mat A( &M[0][0], ROWS, COLS);
    //stampaMatrice(M);

    cx_vec eigval; // eigenvalues may well be complex
    cx_mat eigvec;

    // eigenvalue decomposition for general dense matrices
    
    eig_gen(eigval, eigvec, A);

    cout << "Autovalori della matrice: " << endl;
    cout << eigval << endl;
    cout << "Matrice degli autovettori: " << endl;
    cout << eigvec << endl;
    
    cout << endl << endl << endl;
    for(int i = 0; i < ROWS; i++)
    {
        for(int j = 0; j < COLS; j++)
        {
               M[i][j] = eigvec(i,j);
        }
    }
    cout << "Matrice degli autovettori M: " << endl;
    stampaMatrice(M);


    return 0;

}
