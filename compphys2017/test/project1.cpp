#include<iostream>
#include<iomanip>
#include<cmath>
#include<stdlib.h>
#include<cstdlib>
#include "Coulomb_Functions.hpp"
#include "hartreefock.hpp"

class stateset{
		int cutoff, dimension, spin,
			info;
		double hbar, omega;
	public:
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
stateset::stateset(int max_dist_energy, int s, double hbar_input, double omega_input){
	hbar = hbar_input;
	omega = omega_input;

	dimension = 2;
	info = 2 + dimension;
	cutoff = max_dist_energy;
	spin = s;

	numbers calc;
	states=0;
	for(int i=0;i<cutoff+1;i++){
		states += calc.choose(i+dimension-1,i);
	}
	states *= (spin+1);

	int count=0, energy, n, m, ms;

	state = new int*[states];
	for(int i=0;i<states;i++){
		state[i] = new int[info];
	}
	for(int i=0;i<cutoff+1;i++){
		energy = i;
		for(int j=0;j<energy+1;j+=2){
			n = j;
			m = energy - n;
			n /= 2;
			for(int k=-spin;k<spin+1;k+=2){
				ms = k;
				state[count][0] = n;
				state[count][1] = m;
				state[count][2] = spin;
				state[count][3] = ms;
				count++;
				if(m!=0){
					state[count][0] = n;
					state[count][1] = -m;
					state[count][2] = spin;
					state[count][3] = ms;
					count++;
				}
			}
		}
	}
}
stateset::~stateset(){
	for(int i=0;i<cutoff+1;i++){
		delete[] state[i];
	}
	delete[] state;
}
void stateset::print(){
	double TE;
	std::cout << "E cut off:	" << cutoff << std::endl;
	std::cout << "State count: 	" << states << std::endl;
	std::cout << std::endl;
	std::cout << std::setw(8) << "TE";
	std::cout << std::setw(8) << "n" << std::setw(8) << "m";
	std::cout << std::setw(8) << "spin" << std::setw(8) << "ms" << std::endl;
	std::cout << std::endl;
	for(int i=0;i<states;i++){
		TE = 2*state[i][0] + abs(state[i][1]) + 1;
		TE *= hbar*omega;
		std::cout << std::setw(8) << TE;
		for(int j=0;j<info;j++){
			std::cout << std::setw(8) << state[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
int numbers::factorial(int n){
	int x=1;
	if(n<0){
		x=0;
	}
	else if(n==0){
		x=1;
	}
	else{
		for(int i=1;i<n;i++){
			x *= i;
		}
	}

	value = x;
	return x;
}
int numbers::choose(int a, int b){
	int x=1;
	if(b>a){
		x=0;
	}
	else if(a<0||b<0){
		x=0;
	}
	else{
		int c = a-b;
		if(b<c){
			c=b;
		}
		for(int i=0;i<c;i++){
			x *= (a-i);
			x /= (i+1);
		}
	}
	value = x;
	return x;
}
	
int main(int argc, char* argv[]){
	int 	n,s,states,
		ni, mi, nj, mj, nk, mk, nl, ml,
		msi, msj, msk, msl,
		Minit, Mfinal,
		inum, jnum, knum, lnum,
		states3, states2,
		place, number=0;
	double hbar,omega,hbaromega;
	if(argc<5){
		std::cout << "Bad usage. Enter also 'n s hbar omega' on same line." << std::endl;
		exit(1);
	}
	else{
		n = atoi(argv[1]);
		s = atoi(argv[2]);
		hbar = atof(argv[3]);
		omega = atof(argv[4]);
	}

	hbaromega = hbar*omega;
	stateset test(n+1,s,hbar,omega);

	states = test.states;
	states2 = states*states;
	states3 = states2*states;

	int size = pow(states,4.0);
	double* V = new double[size];
	for(int i=0;i<size;i++){
		V[i] = 0.0;
	}

	for(int i=0;i<states;i++){
		inum = i*states3;
		ni = test.state[i][0];
		mi = test.state[i][1];
		msi = test.state[i][3];
		for(int j=0;j<states;j++){
			jnum = j*states2;
			nj = test.state[j][0];
			mj = test.state[j][1];
			msj = test.state[j][3];
			Minit = mi + mj;
			for(int k=0;k<states;k++){
				knum = k*states;
				nk = test.state[k][0];
				mk = test.state[k][1];
				msk = test.state[k][3];
				if(msi==msk){
					for(int l=0;l<states;l++){
						lnum = l;
						nl = test.state[l][0];
						ml = test.state[l][1];
						msl = test.state[l][3];
						Mfinal = mk + ml;
						place = inum + jnum + knum + lnum;
						
						if(msj==msl){
							if(Minit==Mfinal){
								V[place] = Coulomb_HO(hbaromega,ni,mi,nj,mj,nl,ml,nk,mk);
								//std::cout << number << std::setw(10) << i << std::setw(10) << j << std::setw(10) << k << std::setw(10) << l << std::setw(10) << V[place] << std::endl;
								number++;
							}
						}
					}
				}
			}
		}
	}

	std::cout << std::endl;
	test.print();
	std::cout << std::endl << number << std::endl;

	delete[] V;

	return 0;
}