#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;

int main() {
    for (int i = 410; i < 472; i += 2) {
        ifstream inputFile("L180/beta0" + to_string(i) + ".txt");
        if (!inputFile.is_open()) {
            cerr << "Error: Unable to open file L10/beta0" + to_string(i) + ".txt\n";
            continue;
        }

        // Termalizzazione time  
        int t = 100000;

        // Nome colonne
        vector<double> En;
        vector<double> m;
        double value;
        while (inputFile >> value) {
            En.push_back(value);
            inputFile >> value;
            m.push_back(value);
        }
        inputFile.close();

        // Definisci binsize e binnumber
        int dim = En.size();
        int binsize = 2000;
        int binnumber = dim / binsize;

        //  BINNING
        vector<double> mediaE(binnumber);
        vector<double> quadratoE(binnumber);
        vector<double> mediaM(binnumber);
        vector<double> moduloM(binnumber);
        vector<double> quadratoM(binnumber);
        vector<double> quartaM(binnumber);
        vector<double> M_escluso(binnumber), E_escluso(binnumber), c_escluso(binnumber), sigmaM_escluso(binnumber), chi_escluso(binnumber), binder_escluso(binnumber);


        for (int j = 0; j < binnumber; j++) {
            mediaE[j] = 0;
            quadratoE[j] = 0;
            mediaM[j] = 0;
            moduloM[j] = 0;
            quadratoM[j] = 0;
            quartaM[j] = 0;

            for (int k = 0; k < binsize; ++k) {
                mediaE[j] += En[j * binsize + k];
                quadratoE[j] += En[j * binsize + k] * En[j * binsize + k];
                mediaM[j] += m[j * binsize + k];
                moduloM[j] += abs(m[j * binsize + k]);
                quadratoM[j] += m[j * binsize + k] * m[j * binsize + k];
                quartaM[j] += m[j * binsize + k] * m[j * binsize + k] * m[j * binsize + k] * m[j * binsize + k];
            }

            mediaE[j] /= binsize;
            quadratoE[j] /= binsize;
            mediaM[j] /= binsize;
            moduloM[j] /= binsize;
            quadratoM[j] /= binsize;
            quartaM[j] /= binsize;
            
        }

        // Calcola la media per energia e magnetizzazione
        double mediaE_media = 0;
        double quadratoE_media = 0;
        double mediaM_media = 0;
        double moduloM_media = 0;
        double quadratoM_media = 0;
        double quartaM_media = 0;

        for (int j = 0; j < binnumber; ++j) {
            mediaE_media += mediaE[j];
            quadratoE_media += quadratoE[j];
            mediaM_media += mediaM[j];
            moduloM_media += moduloM[j];
            quadratoM_media += quadratoM[j];
            quartaM_media += quartaM[j];
        }
        mediaE_media /= binnumber;
        quadratoE_media /= binnumber;
        mediaM_media /= binnumber;
        moduloM_media /= binnumber;
        quadratoM_media /= binnumber;
        quartaM_media /= binnumber;


        // JACKKNIFE
        double c = 0.0, ei = 0.0, mi = 0.0, chi = 0.0, binder = 0.0;
        double erc = 0.0, devstE = 0.0, devstM = 0.0, erChi = 0.0, erBinder = 0.0;
        
        for (int j = 0; j < binnumber; j++) {
            // Calcola le medie escludendo il bin j
            double mediaE_escluso = (mediaE_media * binnumber - mediaE[j]) / (binnumber - 1);
            double quadratoE_escluso = (quadratoE_media * binnumber - quadratoE[j]) / (binnumber - 1);
            double mediaM_escluso = (mediaM_media * binnumber - mediaM[j]) / (binnumber - 1);
            double moduloM_escluso = (moduloM_media * binnumber - moduloM[j]) / (binnumber - 1);
            double quadratoM_escluso = (quadratoM_media * binnumber - quadratoM[j]) / (binnumber - 1);
            double quartaM_escluso = (quartaM_media * binnumber - quartaM[j]) / (binnumber - 1);

            // Aggiorna le funzioni d'interesse
            c += quadratoE_escluso - mediaE_escluso*mediaE_escluso;
            ei += mediaE_escluso;
            mi += mediaM_escluso;
            chi += quadratoM_escluso - moduloM_escluso*moduloM_escluso;
            binder += quartaM_escluso/(quadratoM_escluso*quadratoM_escluso);

            //Tiene traccia degli esclusi
            c_escluso[j]=quadratoE_escluso - mediaE_escluso*mediaE_escluso;
            E_escluso[j]=mediaE_escluso;
            M_escluso[j]=mediaM_escluso;
            chi_escluso[j]=quadratoM_escluso - moduloM_escluso*moduloM_escluso;
            binder_escluso[j]=quartaM_escluso/(quadratoM_escluso*quadratoM_escluso);    
        }

        // Calcola le varianze e deviazioni standard finali
        c /= (binnumber-1);
        ei /= (binnumber-1);
        mi /= (binnumber-1);
        chi /= (binnumber-1);
        binder /= (binnumber-1);

        for(int i=0; i<binnumber; i++)
        {
            erc += (c_escluso[i]-c)*(c_escluso[i]-c);
            devstE += (E_escluso[i]-ei)*(E_escluso[i]-ei);
            devstM += (M_escluso[i]-mi)*(M_escluso[i]-mi);
            erChi += (chi_escluso[i]-chi)*(chi_escluso[i]-chi);
            erBinder += (binder_escluso[i]-binder)*(binder_escluso[i]-binder);
        }

        erc = erc *(binnumber/(binnumber-1));
        erc = sqrt(erc);
        devstE = devstE *(binnumber/(binnumber-1));
        devstE = sqrt(devstE);
        devstM = devstM *(binnumber/(binnumber-1));
        devstM = sqrt(devstM);
        erChi = erChi *(binnumber/(binnumber-1));
        erChi = sqrt(erChi);
        erBinder = erBinder *(binnumber/(binnumber-1));
        erBinder = sqrt(erBinder);

        // Leggi le vecchie soluzioni se presenti
        ifstream oldSolutionsFile("180.txt");
        vector<string> oldSolutions;
        string line;
        if (oldSolutionsFile.is_open()) {
            while (getline(oldSolutionsFile, line)) {
                oldSolutions.push_back(line);
            }
            oldSolutionsFile.close();
        }

        // Scrivi le soluzioni su file
        ofstream outputFile("180.txt");
        if (outputFile.is_open()) {
            // Scrivi le vecchie soluzioni
            for (const auto& solution : oldSolutions) {
                outputFile << solution << '\n';
            }

            // Scrivi le nuove soluzioni
            outputFile << i / 1000.0 << ' ' << mediaE_media << ' ' << devstE << ' ' << c << ' ' << erc << ' ' << mediaM_media << ' ' << devstM << ' ' << moduloM_media << ' ' << chi << ' ' << erChi << ' ' << binder << ' ' << erBinder << '\n';
            outputFile.close();
        } else {
            cerr << "Error: Unable to write to file 30.txt\n";
        }
    }

    return 0;
}