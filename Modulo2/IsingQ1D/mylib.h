#ifndef MYLIB_H
#define MYLIB_H

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

using namespace std;
using namespace Eigen;

//DEF VARIABILI DI AMBIENTE
const int N = 10;
const int rows = 1024;
const double hbar = 1.0;
const complex<double> I = {0.0,1.0};

//COSRTUZIONI DI OPERATORI
void costruisciBase(SparseMatrix<int>& base);
void hamiltoniana(SparseMatrix<double>& H, SparseMatrix<int>& base, double g, double h);
void hamiltonianaOPEN(SparseMatrix<double>& H, SparseMatrix<int>& base, double g, double h);
void insertNoise(SparseMatrix<int>& base, SparseMatrix<double>& H);
double longitudinale(SparseMatrix<int>& base, int i);
double accoppiamento(SparseMatrix<int>& base, int i);
double accoppiamentoOPEN(SparseMatrix<int>& base, int i);
void trasverso(SparseMatrix<int>& base, SparseMatrix<double>& H, double g);
double mtrasversa(VectorXd& v, SparseMatrix<int>& base);
double mlongitudinale(VectorXd& v, SparseMatrix<int>& base);
double Mx(VectorXd& v, SparseMatrix<int>& base);
double Mz(VectorXd& v, SparseMatrix<int>& base);
double Mz_complex(VectorXcd& v, SparseMatrix<int>& base);
double MzDef(VectorXd& v, SparseMatrix<int>& base);
double MzDef_complex(VectorXcd& v, SparseMatrix<int>& base);
double M4(VectorXd& v, SparseMatrix<int>& base);
double M2(VectorXd& v, SparseMatrix<int>& base);



//GROUND STATE E AFFINI
void ORDINE_ASC(VectorXd& eigenvalues, MatrixXd& eigenvectors,VectorXd& sortedEigenvalues, MatrixXd& sortedEigenvectors);

//DINAMICA
void rungeKutta(VectorXcd& y, SparseMatrix<double>& H, double t, double dt);

#endif