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
#include <vector>

using namespace std;

int main()
{
    // Apri il file di dati in input
    ifstream inputFile("dati.txt"); // Sostituisci "dati.txt" con il nome del tuo file

    if(!inputFile.is_open())
    {
        cout << "Impossibile aprire il file." << endl;
        return 1;
    }

    vector<double> colonna1;
    vector<double> colonna2;
    double valore1, valore2;
    
    // Leggi i dati dal file e inseriscili nelle colonne
    while(inputFile >> valore1 >> valore2)
    {
        colonna1.push_back(valore1);
        colonna2.push_back(valore2);
    }

    inputFile.close();

    // Calcola la media e la deviazione standard per ciascuna colonna
    double mediaColonna1 = 0.0;
    double mediaColonna2 = 0.0;
    double deviazioneStandardColonna1 = 0.0;
    double deviazioneStandardColonna2 = 0.0;

    for(double dato : colonna1)
    {
        mediaColonna1 += dato;
    }
    
    for(double dato : colonna2)
    {
        mediaColonna2 += dato;
    }

    mediaColonna1 /= colonna1.size();
    mediaColonna2 /= colonna2.size();

    for(double dato : colonna1)
    {
        double scarto = dato - mediaColonna1;
        deviazioneStandardColonna1 += scarto * scarto;
    }

    for(double dato : colonna2)
    {
        double scarto = dato - mediaColonna2;
        deviazioneStandardColonna2 += scarto * scarto;
    }

    deviazioneStandardColonna1 = sqrt(deviazioneStandardColonna1 / colonna1.size());
    deviazioneStandardColonna2 = sqrt(deviazioneStandardColonna2 / colonna2.size());

    // Stampare i risultati per entrambe le colonne
    cout << "Colonna 1 - Media: " << mediaColonna1 << ", Deviazione standard: " << deviazioneStandardColonna1 << endl;
    cout << "Colonna 2 - Media: " << mediaColonna2 << ", Deviazione standard: " << deviazioneStandardColonna2 << endl;

    return 0;
}





