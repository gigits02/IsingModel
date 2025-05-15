#include <iostream>
#include <cmath>
#include <algorithm>
#include "eigen/Eigen/Sparse"

using namespace std;
using namespace Eigen;

MatrixXd kron(const MatrixXd& A, const MatrixXd& B) {
    int rowsA = A.rows();
    int colsA = A.cols();
    int rowsB = B.rows();
    int colsB = B.cols();
    MatrixXd result(rowsA * rowsB, colsA * colsB);
    for (int i = 0; i < rowsA; ++i) {
        for (int j = 0; j < colsA; ++j) {
            result.block(i * rowsB, j * colsB, rowsB, colsB) = A(i, j) * B;
        }
    }
    return result;
}

int main()
{
    
    MatrixXd BlockH(2, 2);
    BlockH.setZero();
    Matrix2d I = Matrix2d::Identity();
    MatrixXd risultato = kron(I,BlockH);
    cout << "Prodotto di kronecker tra BlockH e Identità:" << endl;
    cout << risultato << endl;


    return 0;

}