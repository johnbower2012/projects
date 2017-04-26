#ifndef STATESET_H
#define STATESET_H	

#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>
#include "armadillo"
#include "time.h"
#include "memory.h"

class stateset{
		int cutoff, dimension, spin,
			info;
	public:
		double hbar, omega;
		int** state;
		int states;

		stateset(int,int,double,double);
		~stateset();

		void print();
};
class numbers{
	public:
		int value;

		numbers(){value=0;}

		int factorial(int);
		int choose(int,int);
};

void sp_energies(const stateset&,int,double,arma::mat&);		
void twobody(const stateset&,double*&);
void twobody(const stateset&, matrix4D<double>& V);


#endif
