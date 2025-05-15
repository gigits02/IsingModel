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
using namespace std;

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

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////




void gengauss(double delta, long int Npassi,double x0)
{
    double xn,xt,E,r,t;
    char nfile[256]="dati.txt";
    ofstream fout;
    fout.open( nfile );
    fout<<"# t  x"<<endl;
    fout << 0 << "  " << x0 << endl;
    xn = x0;
    for(int i = 0; i < Npassi; i++)
    {
        xt=xn+delta*(1-2*myrand()); 
        E=-0.5*(xt*xt-xn*xn);
        if(E >= 0)
        {
            xn = xt;
            fout << i << "  " << xn <<endl;
        }
        else
        {
            if(1 < exp(E))
            {
                t = 1;
            }
            else
            {
                t = exp(E);
            }
            r = myrand();
            if(r < t)
            {
                xn = xt;
                fout<< i << "  " << xn <<endl;
            }
            else
            {
                fout<< i << "  "<< xn << endl;
            }
        }
    }
fout.close();
}

int main()
{
    double x0;
    long int Npassi;
    double delta;

    cout << "inserisci il punto di inizio"<<endl;
    cin >> x0;
    cout << "inserisci il numero di iterazioni"<<endl;
    cin >> Npassi;
    cout<<"inserisci lo step"<<endl;
    cin >> delta;

    gengauss(delta,Npassi,x0);


    return 0;
}