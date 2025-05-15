#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>
#include <algorithm>

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
    // Definizione delle costanti
    double Delta = 0.0; // Anisotropia ZZ
    int m = 10;         // Numero di stati mantenuti m
    int NIter = 10;     // Numero di iterazioni. La dimensione finale del reticolo è 2*NIter + 2

    // Inizializzazione degli operatori locali
    Matrix2d I = Matrix2d::Identity();
    Matrix2d Sz, Sp, Sm;
    Sz << 0.5, 0,
          0, -0.5;
    Sp << 0, 0,
          1, 0;
    Sm << 0, 1,
          0, 0;

    // Blocchi iniziali (si assume simmetria di riflessione)
    Matrix2d BlockSz = Sz;
    Matrix2d BlockSp = Sp;
    Matrix2d BlockSm = Sm;
    Matrix2d BlockI = I;
    MatrixXd BlockH(2, 2);
    BlockH.setZero();
    double Energy = 0.0;

    // Numerical Renormalization Group (NRG)
    for (int l = 1; l <= NIter; ++l) {
        int SystSize = pow(2, l);

        // Matrice Hamiltoniana del superblocco
        MatrixXd H_super = kron(BlockH, BlockI) + kron(BlockI, BlockH)
                         - Delta * kron(BlockSz, BlockSz)
                         + 0.5 * (kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp));
        H_super = 0.5 * (H_super + H_super.transpose()); // Assicura che H sia simmetrica

        double LastEnergy = Energy;

        // Diagonalizzazione dell'Hamiltoniana
        SelfAdjointEigenSolver<MatrixXd> es(H_super);
        double Energy = es.eigenvalues()[0];
        double Ener2 = Energy / SystSize;
        double EnergyPerBond = (Energy - LastEnergy) / (SystSize / 2);

        // Costruzione dell'operatore di troncamento
        int NKeep = min(static_cast<int>(es.eigenvalues().size()), m);
        MatrixXd Omatr = es.eigenvectors().leftCols(NKeep);

        cout << SystSize << "\t" << Energy << "\t" << EnergyPerBond << "\t" << Ener2 << endl;

        BlockSz = kron(BlockI, BlockSz);
        BlockSp = kron(BlockI, BlockSp);
        BlockSm = kron(BlockI, BlockSm);
        BlockI = kron(BlockI, BlockI);

        // Trasforma gli operatori del blocco nella base troncata
        BlockH = Omatr.transpose() * H_super * Omatr;
        BlockSz = Omatr.transpose() * BlockSz * Omatr;
        BlockSp = Omatr.transpose() * BlockSp * Omatr;
        BlockSm = Omatr.transpose() * BlockSm * Omatr;
        BlockI = Omatr.transpose() * BlockI * Omatr;
    }

    // Density Matrix Renormalization Group (DMRG) - Algoritmo infinito
    for (int l = 1; l <= NIter; ++l) {
        int SystSize = 2 * l + 2;

        // Ottieni gli operatori dimensionali 2m per il blocco + sito
        BlockH = kron(BlockH, I) - Delta * kron(BlockSz, Sz)
               + 0.5 * (kron(BlockSp, Sm) + kron(BlockSm, Sp));
        BlockSz = kron(BlockI, Sz);
        BlockSp = kron(BlockI, Sp);
        BlockSm = kron(BlockI, Sm);
        BlockI = kron(BlockI, I);

        // Matrice Hamiltoniana del superblocco
        MatrixXd H_super = kron(BlockH, BlockI) + kron(BlockI, BlockH)
                         - Delta * kron(BlockSz, BlockSz)
                         + 0.5 * (kron(BlockSp, BlockSm) + kron(BlockSm, BlockSp));
        H_super = 0.5 * (H_super + H_super.transpose()); // Assicura che H sia simmetrica

        // Diagonalizzazione dell'Hamiltoniana
        double LastEnergy = Energy;
        SelfAdjointEigenSolver<MatrixXd> es(H_super);
        double Energy = es.eigenvalues()[0];
        double EnergyPerBond = (Energy - LastEnergy) / 2;
        double Ener2 = Energy / SystSize;

        // Forma la matrice densità ridotta
        MatrixXd Psi = es.eigenvectors().col(0);
        MatrixXd PsiMatrix = Map<MatrixXd>(Psi.data(), es.eigenvalues().size(), 1);
        MatrixXd Rho = PsiMatrix * PsiMatrix.transpose();

        // Diagonalizza la matrice densità
        SelfAdjointEigenSolver<MatrixXd> es_rho(Rho);
        VectorXd D = es_rho.eigenvalues();
        MatrixXd V = es_rho.eigenvectors();
        D /= D.sum(); // Normalizza i valori

        // Costruzione dell'operatore di troncamento
        int NKeep = min(static_cast<int>(D.size()), m);
        MatrixXd Omatr = V.leftCols(NKeep);
        double TruncationError = 1 - D.head(NKeep).sum();

        cout << SystSize << "\t" << Energy << "\t" << EnergyPerBond << "\t" << Ener2 << "\t" << TruncationError << endl;

        // Trasforma gli operatori del blocco nella base troncata
        BlockH = Omatr.transpose() * BlockH * Omatr;
        BlockSz = Omatr.transpose() * BlockSz * Omatr;
        BlockSp = Omatr.transpose() * BlockSp * Omatr;
        BlockSm = Omatr.transpose() * BlockSm * Omatr;
        BlockI = Omatr.transpose() * BlockI * Omatr;
    }

    return 0;
}