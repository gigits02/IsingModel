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
#include <complex>
#include <random>
#include <chrono>
#include "mylib.h"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Sparse"
#include "eigen/Eigen/Eigenvalues"

using namespace std;
using namespace Eigen;

//DEF VARIABILI DI AMBIENTE

//COSTRUZIONI DI OPERATORI
//COSTRUZIONE DELLA BASE
uniform_real_distribution<> dis(0.0, 1.0);
void costruisciBase(SparseMatrix<int>& base)					
{    
    for (int i = 0; i < rows; i++)
	{
        int value = i;
        for (int j = N-1; j >= 0; j--)
		{
			if(value%2 == 1)
			{
				base.insert(i,j)=1;
			}
    
            value /= 2;
        }
    }
	base.makeCompressed();
}

void hamiltoniana(SparseMatrix<double>& H, SparseMatrix<int>& base, double g, double h)
{
	for(int i = 0; i < rows; i++)
	{	 
		//Costruiamo pezzo di H diagonale, salvo in riga l'indice relativo all'elemento di base considerato
		if( abs(h - 0.00) < 0.0000001 ) //se h è nullo non gli faccio chiamare la funzione del termine longitudinale
		{
			cout << "ATTENZIONE" << endl;
			//inserisco += anche se entro una sola volta per ogni i e j perché sono termini tutti diagonali e H è inzialmente azzerata.
			//qualora volessi aggiungere qualche termine non diagonale il += torna utile.
			H.insert(i,i) = -accoppiamento(base, i);
		}
		else
		{
			H.insert(i,i) = -accoppiamento(base, i) - h*longitudinale(base, i);
		}
	}
	if( abs(g - 0.00) > 0.00001 )
	{
		trasverso(base,H,g);
	}
	SparseMatrix<double> HT = H.transpose();
	H= 0.5*(H+HT);
	H.makeCompressed();

}

void hamiltonianaOPEN(SparseMatrix<double>& H, SparseMatrix<int>& base, double g, double h)
{
	for(int i = 0; i < rows; i++)
	{	 
		//Costruiamo H non integrabile, aggiungendo un campo longitudinale che varia randomicamente (noise)
		{
			H.insert(i,i) = -accoppiamentoOPEN(base, i);
		}
		
	}
	if( abs(g - 0.00) > 0.00001 )
	{
		trasverso(base,H,g);
	}
	
	insertNoise(base, H);

	SparseMatrix<double> HT = H.transpose();
	H= 0.5*(H+HT);
	H.makeCompressed();
}

void insertNoise(SparseMatrix<int>& base, SparseMatrix<double>& H)
{
	//Inizializzo generatore numeri casuali attraverso algoritmo Mersenne Twister
	// Usa l'orologio di sistema per generare un seed casuale
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    mt19937 generator(seed);
	uniform_real_distribution<> dis(-1, 1);

	for(int i = 0; i < N; i++)
	{
		double c = dis(generator);
		for(int j = 0; j < rows; j++)
		{
			if(base.coeff(j,i) == 0)
			{
				H.coeffRef(j,j) += c;
				cout << "sono nell'IF e c vale: " << c << endl;
			}
			else
			{
				H.coeffRef(j,j) += -c;
				cout << "sono nell'ELSE e c vale: " << c << endl;
			}
		}
	}
}

double longitudinale(SparseMatrix<int>& base, int i)
{
	double somma = 0;
	for(int k = 0; k < N; k++)
	{
		if(base.coeff(i,k) == 0)
		{
			somma += 1;
		}
		else
		{
			somma -= 1;
		}
	}
		
	return somma;
}

double accoppiamento(SparseMatrix<int>& base, int i)
{
	double somma = 0;
	
	for(int k = 0; k < N; k++)
	{		
		if(base.coeff(i,k) != base.coeff(i,((k+1)%N)))
		{
			somma -= 1;
		}
		else
		{
			somma += 1;
		}	
	}
	
	return somma;
}

double accoppiamentoOPEN(SparseMatrix<int>& base, int i)
{
	double somma = 0;
	
	for(int k = 0; k < N-1; k++)
	{		
		if(base.coeff(i,k) != base.coeff(i,((k+1))))
		{
			somma -= 1;
		}
		else
		{
			somma += 1;
		}	
	}
	
	return somma;
}

void trasverso(SparseMatrix<int>& base, SparseMatrix<double>& H, double g)
{
	
	int indice = 0;
	for(int j = 0; j < rows; j++)
		for(int k = 0; k < N; k++)
		{
			if(base.coeff(j,k) == 0)
			{
				indice = j + static_cast<int>((pow(2, (N-1-k))));
			}
			else
			{
				indice = j - static_cast<int>((pow(2, (N-1-k))));
			}
			
			H.insert(indice,j) = -g;
		}	
		
}

double mtrasversa(VectorXd& v, SparseMatrix<int>& base)
{
	double mx = 0.0;
	int indice = 0;

	for( int j = 0; j < rows; j++)
	{
		if(base.coeff(j,0) == 0)
		{
			indice = j + static_cast<int>((pow(2, (N-1))));
		}
		else
		{
			indice = j - static_cast<int>((pow(2, (N-1))));
		}
		
		mx += v[indice]*v[j];
	}
	return mx;
}

double mlongitudinale(VectorXd& v, SparseMatrix<int>& base)
{
	double mz = 0.0;
	for( int i = 0; i < rows; i++)
	{
		if(base.coeff(i,0) == 0)
		{
			mz += v[i]*v[i];
		}
		else
		{
			mz -= v[i]*v[i];
		}
	}

	return mz;
}

double Mx(VectorXd& v, SparseMatrix<int>& base)
{
	double mx = 0.0;
	int indice = 0;
	for(int i = 0; i < rows; i++)
		for( int j = 0; j < N; j++)
		{
			if(base.coeff(i,j) == 0)
			{
				indice = i + static_cast<int>((pow(2, (N-1-j))));
			}
			else
			{
				indice = i - static_cast<int>((pow(2, (N-1-j))));
			}
			
			mx += v[indice]*v[i];
		}
	return mx/N;
}

double Mz(VectorXd& v, SparseMatrix<int>& base)
{
	double mz = 0.0;
	for( int i = 0; i < rows; i++)
	{
		for(int j = 0; j < N; j++)
		{
			if(base.coeff(i,j) == 0)
			{
				mz += v[i]*v[i];
			}
			else
			{
				mz -= v[i]*v[i];
			}
		}
	}
	return mz/N;
}

double Mz_complex(VectorXcd& v, SparseMatrix<int>& base)
{
	double mz = 0.0;
	for( int i = 0; i < rows; i++)
	{
		for(int j = 0; j < N; j++)
		{
			if(base.coeff(i,j) == 0)
			{
				mz += norm(v[i]);
			}
			else
			{
				mz -= norm(v[i]);
			}
		}
	}
	return mz/N;
}

double MzDef(VectorXd& v, SparseMatrix<int>& base)
{
	double mzdef = 0.0;
	double sum = 0.0;
	for( int i = 0; i < rows; i++)
	{
		sum = 0.0;
		for(int j = 0; j < N; j++)
		{	
			if(base.coeff(i,j) == 0)
			{
				sum += 1;
			}
			else
			{
				sum -= 1;
			}
		}
		mzdef += abs(sum)*v[i]*v[i];
	}
	return mzdef/N;
}

double MzDef_complex(VectorXcd& v, SparseMatrix<int>& base)
{
	double mzdef = 0.0;
	double sum = 0.0;
	for( int i = 0; i < rows; i++)
	{
		sum = 0.0;
		for(int j = 0; j < N; j++)
		{	
			if(base.coeff(i,j) == 0)
			{
				sum += 1;
			}
			else
			{
				sum -= 1;
			}
		}
		mzdef += abs(sum)*norm(v[i]);
	}
	return mzdef/N;
}

double M4(VectorXd& v, SparseMatrix<int>& base)
{
	double m4 = 0.0;
	double sum = 0.0;
	int sign = 1;
	int sign1 = 1;
	int sign2 = 1;
	int sign3 = 1;
	int sign4 = 1;

	for(int riga =0 ; riga<rows ;riga++)
	{
		sign = 1;
		
		for(int i=0 ; i<N ;i++)
		{
			if(base.coeff(riga,i)==0)
			{
				sign1 = sign;
			}
			else
			{
				sign1 = -sign;
			}

			for(int j=0 ; j<N ;j++)
			{
				if(base.coeff(riga,j)==0)
				{
					sign2 = sign1;
				}
				else
				{
					sign2 = -sign1;
				}

				for(int k=0 ; k<N; k++)
				{
					if(base.coeff(riga,k)==0)
					{
						sign3 = sign2;
					}
					else
					{
						sign3 = -sign2;
					}

					for(int l=0 ; l<N ;l++)
					{
						if(base.coeff(riga,l)==0)
						{
							sign4 = sign3;
						}
						else
						{
							sign4 = -sign3;
						}

						sum += sign4*v[riga]*v[riga];
					}
				}
			}
		}
	}

	m4= sum/(N*N*N*N);
	return m4;

}

double M2(VectorXd& v, SparseMatrix<int>& base)
{
	double m2 = 0.0;
	double sum = 0.0;
	int sign = 1;
	int sign1 = 1;
	int sign2 = 1;

	for(int riga =0 ; riga<rows ;riga++)
	{
		sign = 1;
		
		for(int i=0 ; i<N ;i++)
		{
			if(base.coeff(riga,i)==0)
			{
				sign1 = sign;
			}
			else
			{
				sign1 = -sign;
			}

			for(int j=0 ; j<N ;j++)
			{
				if(base.coeff(riga,j)==0)
				{
					sign2 = sign1;
				}
				else
				{
					sign2 = -sign1;
				}

				sum += sign2*v[i]*v[i];
			}
		}
	}

	m2= sum/(N*N);
	return m2;

}
//GROUND STATE E AFFINI 
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

//DINAMICA
void rungeKutta(VectorXcd& y, SparseMatrix<double>& H, double t, double dt)
{
	VectorXcd k1(rows), k2(rows), k3(rows), k4(rows);
	VectorXcd tmp(rows);

	k1= H*y;
	tmp = y + k1*(dt/2.0);
	k2= H*tmp;
	tmp= y + k2*(dt/2.0);
	k3= H*tmp;
	tmp= y + k3*(dt/2.0);
	k4= H*tmp;
	
	y = y -I*(k1 + 2.0 * k2 + 2.0 * k3 + k4) * dt/6.0;
	
	//normalizzazione del vettore
	y.normalize();
}