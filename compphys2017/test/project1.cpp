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

/* declaring proper integer values
	part = particles
	shells = shells
	s = spin
	states --> to be calculated from shells */
	int part, shells, s=1, states;

/* declaring proper double values
	hbar = physical constant
	omega = HO frequency
	E0_etc to be used in final result */
	double hbar=1.0, omega, hbaromega,
			E0_sp=0.0, E0_hf1=0.0, E0_hf2=0.0,
			V_sp=0.0, V_hf=0.0;

/* check for proper inputs from command line */
	if(argc<4){
		std::cout << "Bad usage. Enter also 'part# shells omega' on same line." << std::endl;
		exit(1);
	}
	else{
		part = atoi(argv[1]);
		shells = atoi(argv[2]);
		omega = atof(argv[3]);
		hbaromega = hbar*omega;
	}

/* declare proper object 'stateset'
	contains a list of sp states 0,1,2... for spin s particles in a given number of shells
	print list to screen */
	stateset HarmonicOsc(shells,s,hbar,omega);	
	std::cout << std::endl;
	HarmonicOsc.print();

/* define number of sp states in stateset object
	define two-body V and hamiltonian/density matrix */
	states = HarmonicOsc.states;
	matrix4D<double> V(states,states,states,states);
	arma::mat	H0 = arma::zeros<arma::mat>(states,states),
					densityMatrix = H0;
	arma::vec	E = arma::zeros<arma::vec>(states);

/* calculate sp energies */
	std::cout << "Creating H0...\n";
	sp_energies(HarmonicOsc,states,hbaromega,H0);

/* calculate two-body potential */
	std::cout << "Creating V....";
	std::cout << std::endl;
	twobody(HarmonicOsc,V);

/* conduct HF solver
	prints each E vector iteration to screen */
	std::cout << "Beginning iterations...\n" << std::endl;	
	solve_iterations(H0, V, states, part, densityMatrix, E);

/* calculate final filled-shell energies for both HO and HF bases */ 
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
	
/* print energies to screen */
	std::cout << std::endl;
	std::cout << "E0_sp:     " << E0_sp << std::endl;
	std::cout << "E0_hf1:    " << E0_hf1 << std::endl;
	std::cout << "E0_hf2:    " << E0_hf2 << std::endl << std::endl;

	return 0;
}
