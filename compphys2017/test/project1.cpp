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

	double hbar,omega,hbaromega,
			E0_sp=0.0,E0_hf1=0.0,E0_hf2=0.0,
			V_sp=0.0, V_hf=0.0;

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
				densityMatrix = H0;
	arma::vec	E = arma::zeros<arma::vec>(states);

	sp_energies(test,states,hbaromega,H0);
	twobody(test,V);

	solve_iterations(H0, V, states, part, densityMatrix, E);

	for(int i=0;i<part;i++){
		V_sp=0.0;
		for(int j=0;j<part;j++){
			V_sp += V.memory[i][j][i][j];
		}
		E0_sp += H0(i,i) + 0.5*V_sp;
		E0_hf1 += E(i);
	}
	for(int i=0;i<states;i++){
		E0_hf2 += densityMatrix(i,i)*H0(i,i);
		for(int j=0;j<states;j++){
			for(int k=0;k<states;k++){
				for(int l=0;l<states;l++){
					V_hf += densityMatrix(i,k)*densityMatrix(j,l)*V.memory[i][j][k][l];
				}
			}
		}
	}
	E0_hf1 -= 0.5*V_hf;
	E0_hf2 += 0.5*V_hf;
	

	std::cout << "E0_sp:    " << E0_sp << std::endl;
	std::cout << "E0_hf1:    " << E0_hf1 << std::endl;
	std::cout << "E0_hf2:    " << E0_hf2 << std::endl << std::endl;

	return 0;
}
