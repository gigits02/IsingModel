#include <iostream>
#include <fstream>
#include <vector>

using namespace std;

// Funzione per calcolare l'autocorrelazione per un dato ritardo k
double autocorrelation_k(const vector<double>& data, int k, double mean) {
    int n = data.size();
    double acf_sum = 0;
    for (int i = 0; i < n - k; ++i) {
        acf_sum += (data[i] - mean)*(data[i + k] - mean);
    }
    acf_sum /= (n - k);
    return acf_sum;
}

int main() {
    int t = 100000; // Tempo di termalizzazione (numero di dati da scartare)
    vector<double> seconda_colonna;

    // Carica i dati dal file di testo
    ifstream file("../L10/beta0340.txt");
    double value;
    for (int i = 0; i < t; ++i) {
        file.ignore(); // Ignora il delimitatore (presumendo che sia uno spazio)
        file >> value; // Leggi il valore della seconda colonna
    }
    while (file >> value) {
        file.ignore(); // Ignora il delimitatore (presumendo che sia uno spazio)
        file >> value; // Leggi il valore della seconda colonna
        seconda_colonna.push_back(abs(value));
    }
    file.close();

    // Calcola l'autocorrelazione per vari valori di k
    int max_lag = 10000; // Modifica questo valore in base alle tue esigenze
    vector<double> autocorrs;
    double mean = 0; // Media della seconda colonna
    for (double val : seconda_colonna) {
        mean += abs(val);
    }
    mean /= seconda_colonna.size();

    for (int k = 0; k < max_lag; ++k) {
        double autocorr = autocorrelation_k(seconda_colonna, k, mean);
        autocorrs.push_back(autocorr);
    }

    // Salva i risultati in un file di output
    ofstream outfile("auto10.txt");
    if (outfile.is_open()) {
        for (int k = 0; k < max_lag; ++k) {
            outfile << k << " " << autocorrs[k] << "\n";
        }
        outfile.close();
        cout << "Risultati salvati correttamente in auto10.txt" << endl;
    } else {
        cerr << "Impossibile aprire il file di output!" << endl;
        return 1;
    }

    return 0;
}