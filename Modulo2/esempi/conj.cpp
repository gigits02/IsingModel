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
#include<armadillo>
#include<complex>

using namespace std;
using namespace arma;

int main()
{
    complex<double> array[3] = { {1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0} };
     
    cout << "array complesso: " << endl;
    for(int i = 0; i < 3; i++)
        cout << array[i] << " ";
    cout << endl;

    cout << "array complesso coniugato: " << endl; 
    for(int i = 0; i < 3; i++)
        cout << conj(array[i]) << " ";
    cout << endl;
    
    cout << "conj(v[i])*v[i]:" << endl;
    for(int i = 0; i < 3; i++)
        cout << conj(array[i])*array[i] << " ";
    cout << endl;

    //abbiamo appurato che norm(v[i]) è il modulo quadro dell'elemento i dell'array complesso.
    cout << "norm(v[i]):" << endl;
    for(int i = 0; i < 3; i ++)
        cout << norm(array[i]) << " ";
    cout << endl;

    return 0;
}