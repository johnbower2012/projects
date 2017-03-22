#ifndef HARTREEFOCK_H
#define HARTREEFOCK_H

#include<iostream>
#include<iomanip>
#include<fstream>
#include<cmath>
#include<string>
#include<sstream>
#include "armadillo"
#include "time.h"
#include "memory.h"

#define TOLERANCE 1e-8
#define maxITERATIONS 30

void sp_energies(arma::mat&,int,arma::mat&);
void compute_densityMatrix(arma::mat&,int,arma::mat&);
void solve_iterations(arma::mat&,matrix4D<double>&,int,arma::vec&,arma::mat&,arma::vec&);

#endif
