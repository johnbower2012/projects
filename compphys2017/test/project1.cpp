#include<iostream>
#include<iomanip>
#include<cmath>
#include<stdlib.h>
#include<cstdlib>
#include "armadillo"
#include "Coulomb_Functions.hpp"
#include "hartreefock.hpp"
#include "stateset.hpp"
#include "memory.h"

int main(int argc, char* argv[]){
	int part, shells, s=1, states;

	double hbar,omega,hbaromega;

	if(argc<5){
		std::cout << "Bad usage. Enter also 'part# shells hbar omega' on same line." << std::endl;
		exit(1);
	}
	else{
		part = atoi(argv[1]);
		shells = atoi(argv[2]);
		hbar = atof(argv[3]);
		omega = atof(argv[4]);
		hbaromega = hbar*omega;
	}

	stateset test(shells,s,hbar,omega);	
	std::cout << std::endl;
	test.print();

	states = test.states;

	matrix4D<double> V(states,states,states,states);
	arma::mat	H0 = arma::zeros<arma::mat>(states,states),
				H = H0;
	arma::vec	E = arma::zeros<arma::vec>(states);

	sp_energies(test,states,hbaromega,H0);
	twobody(test,V);

	solve_iterations(H0, V, states, part, H, E);

	std::cout << E << std::endl;

	return 0;
}
