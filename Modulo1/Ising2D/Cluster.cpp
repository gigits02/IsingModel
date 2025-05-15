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


const int dim =30;


void ClusterCreation(int lattice[dim][dim],double ProbAcc,int ClusterMatrix[dim][dim],int a,int b)
{
    if(lattice[a][b]==lattice[(a+1)%dim][b] && ClusterMatrix[(a+1)%dim][b]==0)
    {
        if(myrand()<ProbAcc)
        {
            ClusterMatrix[(a+1)%dim][b]=1;
            ClusterCreation(lattice,ProbAcc,ClusterMatrix,(a+1)%dim,b);
        }
    }
    if(lattice[a][b]==lattice[a][(b+1)%dim] && ClusterMatrix[a][(b+1)%dim]==0)
    {
        if(myrand()<ProbAcc)
        {
             ClusterMatrix[a][(b+1)%dim]=1;
            ClusterCreation(lattice,ProbAcc,ClusterMatrix,a,(b+1)%dim);
        }
    }

    if(a==0 || b==0)
    {
        if(a==0 & b!=0)
        {
            if(lattice[a][b]==lattice[(dim-1)%dim][b] && ClusterMatrix[(dim-1)%dim][b]==0)
                {
                     if(myrand()<ProbAcc)
                        {
                            ClusterMatrix[(dim-1)%dim][b]=1;
                            ClusterCreation(lattice,ProbAcc,ClusterMatrix,(dim-1)%dim,b);
                        }
                }
            if(lattice[a][b]==lattice[a][(b-1)%dim] && ClusterMatrix[a][(b-1)%dim]==0)
                {
                    if(myrand()<ProbAcc)
                        {
                            ClusterMatrix[a][(b-1)%dim]=1;
                            ClusterCreation(lattice,ProbAcc,ClusterMatrix,a,(b-1)%dim);
                        }
                }
        }
        if(b==0 & a!=0)
        {
            if(lattice[a][b]==lattice[(a-1)%dim][b] && ClusterMatrix[(a-1)%dim][b]==0)
                {
                    if(myrand()<ProbAcc)
                        {
                            ClusterMatrix[(a-1)%dim][b]=1;
                            ClusterCreation(lattice,ProbAcc,ClusterMatrix,(a-1)%dim,b);
                        }
                }
            if(lattice[a][b]==lattice[a][(dim-1)%dim] && ClusterMatrix[a][(dim-1)%dim]==0)
                {
                    if(myrand()<ProbAcc)
                        {
                            ClusterMatrix[a][(dim-1)%dim]=1;
                            ClusterCreation(lattice,ProbAcc,ClusterMatrix,a,(dim-1)%dim);
                        }
                }
        }
        if(a == 0 & b == 0)
        {
            if(lattice[a][b] == lattice[(dim-1)%dim][b] && ClusterMatrix[(dim-1)%dim][b] == 0)
                {
                    if(myrand() < ProbAcc)
                        {
                            ClusterMatrix[(dim-1)%dim][b] = 1;
                            ClusterCreation(lattice, ProbAcc, ClusterMatrix, (dim-1)%dim, b);
                        }
                }
            if(lattice[a][b] == lattice[a][(dim-1)%dim] && ClusterMatrix[a][(dim-1)%dim] == 0)
                {
                    if(myrand() < ProbAcc)
                        {
                            ClusterMatrix[a][(dim-1)%dim] = 1;
                            ClusterCreation(lattice, ProbAcc, ClusterMatrix, a, (dim-1)%dim);
                        }
                }

        }
    }
    else
    {

        if(lattice[a][b] == lattice[(a-1)%dim][b] && ClusterMatrix[(a-1)%dim][b] == 0)
        {
            if(myrand()<ProbAcc)
            {
                ClusterMatrix[(a-1)%dim][b]=1;
                ClusterCreation(lattice,ProbAcc,ClusterMatrix,(a-1)%dim,b);
            }
        }
        if(lattice[a][b] == lattice[a][(b-1)%dim] && ClusterMatrix[a][(b-1)%dim] == 0)
        {
            if(myrand() < ProbAcc)
            {
                ClusterMatrix[a][(b-1)%dim] = 1;
                ClusterCreation(lattice, ProbAcc, ClusterMatrix, a, (b-1)%dim);
            }
        }
    }
}

void ClusterFlipping(int lattice[dim][dim],int ClusterMatrix[dim][dim])
{
    for(int i=0;i<dim;i++)
    {
        for(int j=0;j<dim;j++)
        {
            if(ClusterMatrix[i][j]==1)
            {
                lattice[i][j]=-lattice[i][j];
            }
        }
    }
}

void compila(int lattice[dim][dim])
{
    for(int i = 0 ;i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
            lattice[i][j] = 1;
        }
    }
}

void azzera(int cluster[dim][dim])
{
    for(int i = 0 ;i < dim; i++)
        for(int j = 0; j < dim; j++)
            cluster[i][j] = 0;
}

void stampa(int lattice[dim][dim])
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

double energia(int lattice[dim][dim])
{
    double E = 0;
    int tmp=0;
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
                tmp += -(lattice[i][j]*lattice[(i+1)%dim][j] + lattice[i][j]*lattice[i][(j+1)%dim]);
        }
    }
    E=static_cast<double>(tmp);
    return E/(dim*dim);
}

double magnetizzazione(int lattice[dim][dim])
{
    double M = 0;
    int tmp=0;
    for(int i = 0; i < dim; i++)
    {
        for(int j = 0; j < dim; j++)
        {
                tmp += lattice[i][j];
        }
    }
    M =static_cast<double>(tmp);
    return M/(dim*dim);
}

int main()
{
    myrand_init(3,334555);
    srand((unsigned)time(NULL));

    cout<<dim<<endl;
    int lattice[dim][dim];
    compila(lattice);

/*
    cout << "Reticolo iniziale: " << endl;
    stampa(lattice);
    cout << endl; 
*/
    double beta;
    //cout << "Inserire il valore di beta: " << endl;
    //cin>>beta;

    char nfile[256];
    //cout << "Inserire il nome del file (ex. beta1.txt): " << endl;
    //cin >> nfile;
    
    for (double i = 0.340; i <= 0.472; i+=0.002)
    {
        beta = i;
        std::string nfile = "beta0" + std::to_string(int(i*1000)) + ".txt";
    
        int ClusterMatrix[dim][dim];
        azzera(ClusterMatrix);

        double ProbAcc=(1-exp(-2*beta));

        int a=0,b=0;

        ofstream fout;
        fout.open( nfile );
        for(int k=0;k<1000000;k++)
        {
            a=static_cast<int>(dim*myrand());
            b=static_cast<int>(dim*myrand());

            azzera(ClusterMatrix);

            ClusterCreation(lattice, ProbAcc,ClusterMatrix,a,b);
            ClusterFlipping(lattice,ClusterMatrix);

            /*cout<<"reticolo risultante (dopo il flip): "<<endl;
            stampa(lattice);
            cout<<endl;*/

            fout << energia(lattice) << " " << magnetizzazione(lattice) << endl;        

        }
        
    }

    return 0;

}