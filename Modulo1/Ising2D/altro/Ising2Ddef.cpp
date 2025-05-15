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
#include <iomanip>
#include <random>

using namespace std;


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//-------------------------------
// *Really* minimal PCG32 code / (c) 2014 M.E. O'Neill / pcg-random.org
// Licensed under Apache License 2.0 (NO WARRANTY, etc. see website)

typedef struct { uint64_t state;  uint64_t inc; } pcg32_random_t;

uint32_t pcg32_random_r(pcg32_random_t* rng)
    {
    uint64_t oldstate = rng->state;
    // Advance internal state
    rng->state = oldstate * 6364136223846793005ULL + (rng->inc|1);
    // Calculate output function (XSH RR), uses old state for max ILP
    uint32_t xorshifted = (uint32_t) ( ((oldstate >> 18u) ^ oldstate) >> 27u );
    uint32_t rot = (uint32_t) ( oldstate >> 59u );
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
    }

void pcg32_srandom_r(pcg32_random_t* rng, uint64_t initstate, uint64_t initseq)
    {
    rng->state = 0U;
    rng->inc = (initseq << 1u) | 1u;
    pcg32_random_r(rng);
    rng->state += initstate;
    pcg32_random_r(rng);
    }
//-----------------------------------


//----------------my wrapper for pcg32

// random number internal state
pcg32_random_t pcg32_random_state;

// initialization
void myrand_init(unsigned long int initstate, unsigned long int initseq)
  {
  pcg32_srandom_r(&pcg32_random_state, (uint64_t) initstate, (uint64_t) initseq);
  }

// number in [0,1)
double myrand(void)
  {
  return (double) pcg32_random_r(&pcg32_random_state)/(pow(2.0, 32.0));
  }
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


const int dim = 4;

// Funzione per generare un numero casuale INTERO tra 0 e dim
int generateRandomInt()
{
    // Inizializzazione del generatore di numeri casuali con un seme basato sull'orologio di sistema
    static std::random_device rd;
    static std::mt19937 gen(rd()); // Utilizza il seme ottenuto da random_device

    // Definizione del range desiderato per i numeri casuali (0, 1)
    static std::uniform_int_distribution<> dis(0,dim-1);

    // Generazione di un numero casuale nell'intervallo (0, 1)
    return dis(gen);
}


// Funzione per generare un numero casuale REALE tra 0 e 1
double generateRandomNumber()
{
    // Inizializzazione del generatore di numeri casuali con un seme basato sull'orologio di sistema
    static std::random_device rd;
    static std::mt19937 gen(rd()); // Utilizza il seme ottenuto da random_device

    // Definizione del range desiderato per i numeri casuali (0, 1)
    static std::uniform_real_distribution<> dis(0.0, 1.0);

    // Generazione di un numero casuale nell'intervallo (0, 1)
    return dis(gen);
}

bool accetto(long int lattice[dim][dim], int rx, int ry,double sigma[],int lenght,double beta)
{   
    int pv;
    if(rx - 1 < 0 )
        rx = dim;
    if( ry - 1 < 0)
        ry = dim;
    
    pv = lattice[rx][(ry+1)%dim] + lattice[(rx+1)%dim][ry] + lattice[rx][(ry-1)%dim] + lattice[(rx-1)%dim][ry];
    int deltaE = 2*lattice[rx][ry]*pv;
    
    
    if(-(double)(beta*deltaE) >= 0)
        return true;
    
    
    else
    {
        if(abs(pv)==2)
        {
            //cout << "SONO NEL CASO pv = 2" << endl;
            if(generateRandomNumber() /*myrand()*/ < sigma[0])
            {
                //cout << "L'IF È SODDISFATTO" << endl;
                return true;
            }
        }
        if(abs(pv)==4)
        {   
            //cout << "SONO NEL CASO pv = 4" << endl;
            if(generateRandomNumber() /*myrand()*/ < sigma[1])
            {    
                //cout << "L'IF È SODDISFATTO" << endl;
                return true;
            }
        }
    }
        

    return false;
}


void compila(long int lattice[dim][dim])
{
    for(int i = 0 ;i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            lattice[i][j] = 1;
        }
    }
}

void stampa(long int lattice[dim][dim])
{
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            cout << lattice[i][j] << " ";
        }
        cout << endl;
    }
}

double energia(long int lattice[dim][dim])
{
    long double E = 0;
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
                E += -(lattice[i][j]*lattice[(i+1)%dim][j] + lattice[i][j]*lattice[i][(j+1)%dim]);
        }
    }
    return (double)(E/(dim*dim));
}

double magnetizzazione(long int lattice[dim][dim])
{
    long double M = 0;
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
                M += lattice[i][j];
        }
    }
    M = (double)(M/(dim*dim));
    return M;
}


void Metropolis(long int lattice[dim][dim],double sigma[], int lenght,double beta)
{

    long int rx, ry;
    for(int i = 0; i < dim; i++)
        for(int j = 0; j < dim; j++)
        {
            rx = generateRandomInt();
            ry = generateRandomInt();
             /*rx=(int)(dim* (double)generateRandomNumber());
             ry=(int)(dim* (double)generateRandomNumber());*/
            
            if( accetto(lattice, rx, ry, sigma, lenght,beta) == true)
                lattice[rx][ry] = -lattice[rx][ry];
            
        }
        
}

int main()
{
    const unsigned long int seed1=(unsigned long int) time(NULL);
    const unsigned long int seed2=seed1+127;
    myrand_init(seed1,seed2);
    long int lattice[dim][dim];
    compila(lattice);
    srand((unsigned)time(NULL));
    
    double beta;
    cout << "Inserisci valore di beta per il quale effettuare la simulazione:" << endl;
    cin >> beta; 
    
    char nfile[256];
    cout << "Inserisci nome del file di output: " << endl;
    cin >> nfile;
    
    ofstream fout;
   
    fout.open( nfile );
    
    //fout<<"t,E,M"<<endl;

    const int lenght=2;
    double sigma[lenght];
    sigma[0]=exp(-4*beta);
    sigma[1]=exp(-8*beta);
    
    for(int k = 0; k < 1000000; k++)
    {
        Metropolis(lattice,sigma, lenght,beta);
        /*cout << "Reticolo aggiornata: " << endl;
        stampa(lattice);
        cout << endl << endl; 
        
        cout << energia(lattice) << " " << magnetizzazione(lattice) << endl;
        */
        fout << energia(lattice) << " " << magnetizzazione(lattice) << endl;
    }

    fout.close();

    return 0;
}
