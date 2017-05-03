#include "stateset.hpp"
#include "Coulomb_Functions.hpp"
#include "memory.h"

stateset::stateset(int shells, int s, double hbar_input, double omega_input){
	hbar = hbar_input;
	omega = omega_input;

	dimension = 2;
	info = 2 + dimension;
	cutoff = shells - 1;
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
/*
	for(int k=-spin;k<spin+1;k+=2){
		for(int i=-cutoff;i<cutoff+1;i++){
			m = i;
			energy = cutoff - abs(m);
			energy /= 2;
			for(int j=0;j<energy+1;j++){
				n = j;
				ms = k;
				state[count][0] = n;
				state[count][1] = m;
				state[count][2] = spin;
				state[count][3] = ms;
				count++;
			}
		}
	}
*/
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
	std::cout << "E cut off:	" << cutoff+1 << std::endl;
	std::cout << "State count: 	" << states << std::endl;
	std::cout << std::endl;
	std::cout << std::setw(8) << " ";
	std::cout << std::setw(8) << "TE";
	std::cout << std::setw(8) << "n" << std::setw(8) << "m";
	std::cout << std::setw(8) << "spin" << std::setw(8) << "ms" << std::endl;
	std::cout << std::endl;
	for(int i=0;i<states;i++){
		TE = 2*state[i][0] + abs(state[i][1]) + 1;
		TE *= hbar*omega;
		std::cout << std::setw(8) << i << std::setw(8) << TE;
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

void sp_energies(const stateset &object, int size, double hbaromega, arma::mat& H0){
	int i,j;
	for(i=0;i<size;i++){
		for(j=i+1;j<size;j++){
			H0(i,j) = H0(j,i) = 0.0;
		}
		H0(i,i) = hbaromega*(2*object.state[i][0] + abs(object.state[i][1]) + 1.0);
	}
}
void twobody(const stateset &object, double*& V){
	
	int states,	states2, states3,
		ni, mi, nj, mj, nk, mk, nl, ml,
		msi, msj, msk, msl,
		Minit, Mfinal, MSinit, MSfinal,
		inum, jnum, knum, lnum,
		place/*,number=0*/;

	double hbaromega = object.hbar*object.omega,
			dir, exch;

	states = object.states;
	states2 = states*states;
	states3 = states2*states;

	int size = pow(states,4.0);
	for(int i=0;i<size;i++){
		V[i] = 0.0;
	}

	for(int i=0;i<states;i++){
		inum = i*states3;
		ni = object.state[i][0];
		mi = object.state[i][1];
		msi = object.state[i][3];
		for(int j=0;j<states;j++){
			jnum = j*states2;
			nj = object.state[j][0];
			mj = object.state[j][1];
			msj = object.state[j][3];
			Minit = mi + mj;
			MSinit = msi + msj;
			for(int k=0;k<states;k++){
				knum = k*states;
				nk = object.state[k][0];
				mk = object.state[k][1];
				msk = object.state[k][3];
				for(int l=0;l<states;l++){
					lnum = l;
					nl = object.state[l][0];
					ml = object.state[l][1];
					msl = object.state[l][3];
					Mfinal = mk + ml;
					MSfinal = msk + msl;
					if(Minit==Mfinal){
						if(MSinit==MSfinal){
							place = inum + jnum + knum + lnum;
							dir = exch = 0.0;
							if(msi==msk){
								dir = Coulomb_HO(hbaromega,ni,mi,nj,mj,nl,ml,nk,mk);
							}
							if(msi==msl){
								exch = Coulomb_HO(hbaromega,ni,mi,nj,mj,nk,mk,nl,ml);
							}
							V[place] = dir - exch;
							//std::cout << number << std::setw(10) << i << std::setw(10) << j << std::setw(10) << k << std::setw(10) << l << std::setw(10) << V[place] << std::endl;
							//number++;
						}
					}
				}
			}
		}
	}
}
void twobody(const stateset &object, matrix4D<double>& V){
	
	int states,
		ni, mi, nj, mj, nk, mk, nl, ml,
		msi, msj, msk, msl,
		Minit, Mfinal, MSinit, MSfinal/*,
		number=0*/;

	double hbaromega = object.hbar*object.omega,
			dir, exch;

	states = object.states;

	for(int i=0;i<states;i++){
		for(int j=0;j<states;j++){
			for(int k=0;k<states;k++){
				for(int l=0;l<states;l++){
					V.memory[i][j][k][l] = 0.0;
				}
			}
		}
	}

	for(int i=0;i<states;i++){
		ni = object.state[i][0];
		mi = object.state[i][1];
		msi = object.state[i][3];
		for(int j=0;j<states;j++){
			nj = object.state[j][0];
			mj = object.state[j][1];
			msj = object.state[j][3];
			Minit = mi + mj;
			MSinit = msi + msj;
			for(int k=0;k<states;k++){
				nk = object.state[k][0];
				mk = object.state[k][1];
				msk = object.state[k][3];
				for(int l=0;l<states;l++){
					nl = object.state[l][0];
					ml = object.state[l][1];
					msl = object.state[l][3];
					Mfinal = mk + ml;
					MSfinal = msk + msl;
					if(Minit==Mfinal){
						if(MSinit==MSfinal){
							dir = 0.0;
							exch = 0.0;
							if(msi==msk){
								dir = Coulomb_HO(hbaromega,ni,mi,nj,mj,nk,mk,nl,ml);
							}
							if(msi==msl){
								exch = Coulomb_HO(hbaromega,ni,mi,nj,mj,nl,ml,nk,mk);
							}
							V.memory[i][j][k][l] = dir - exch;
/*
							if(V.memory[i][j][k][l]!=0.0){
								number++;
								std::cout << number << " " << i << j << k << l << " "  << mi << mj << mk << ml << " "  << msi << msj << msk << msl << " "  << Minit << Mfinal << " " << MSinit << MSfinal << std::setw(10) << dir << std::setw(10) << exch << std::setw(10) << V.memory[i][j][k][l] << std::endl;


							}
*/
						}
					}
				}
			}
		}
	}
}
void twobodywrite(const stateset &object, matrix4D<double>& V){
	
	int states,
		ni, mi, nj, mj, nk, mk, nl, ml,
		msi, msj, msk, msl,
		Minit, Mfinal, MSinit, MSfinal,
		number=0, shell;

/* create file name then open */
	std::string outfile;
	std::stringstream shel,omeg; 
	shell = object.cutoff+1;
	shel << shell;
	outfile = "twobodys" + shel.str();
	omeg << object.omega;
	outfile += "o" + omeg.str() + ".dat";
	std::ofstream ofile;

	ofile.open(outfile);

	double hbaromega = object.hbar*object.omega,
			dir, exch;

	states = object.states;

	for(int i=0;i<states;i++){
		for(int j=0;j<states;j++){
			for(int k=0;k<states;k++){
				for(int l=0;l<states;l++){
					V.memory[i][j][k][l] = 0.0;
				}
			}
		}
	}

	for(int i=0;i<states;i++){
		ni = object.state[i][0];
		mi = object.state[i][1];
		msi = object.state[i][3];
		for(int j=0;j<states;j++){
			nj = object.state[j][0];
			mj = object.state[j][1];
			msj = object.state[j][3];
			Minit = mi + mj;
			MSinit = msi + msj;
			for(int k=0;k<states;k++){
				nk = object.state[k][0];
				mk = object.state[k][1];
				msk = object.state[k][3];
				for(int l=0;l<states;l++){
					nl = object.state[l][0];
					ml = object.state[l][1];
					msl = object.state[l][3];
					Mfinal = mk + ml;
					MSfinal = msk + msl;
					if(Minit==Mfinal){
						if(MSinit==MSfinal){
							dir = 0.0;
							exch = 0.0;
							if(msi==msk){
								dir = Coulomb_HO(hbaromega,ni,mi,nj,mj,nk,mk,nl,ml);
							}
							if(msi==msl){
								exch = Coulomb_HO(hbaromega,ni,mi,nj,mj,nl,ml,nk,mk);
							}
							V.memory[i][j][k][l] = dir - exch;
							if(V.memory[i][j][k][l]!=0.0){
								number++;
							}

						}
					}
				}
			}
		}
	}
	ofile << number << std::endl;
	for(int i=0;i<states;i++){
		for(int j=0;j<states;j++){
			for(int k=0;k<states;k++){
				for(int l=0;l<states;l++){
					if(V.memory[i][j][k][l]!=0.0){
						ofile << i << std::setw(15) << j << std::setw(15) << k << std::setw(15) << l << std::setw(15) << V.memory[i][j][k][l] << std::endl;
					}
				}
			}
		}
	}
	ofile.close();
}
void twobodyread(const stateset &object, matrix4D<double>& V){
	int states;
	states = object.states;
	for(int i=0;i<states;i++){
		for(int j=0;j<states;j++){
			for(int k=0;k<states;k++){
				for(int l=0;l<states;l++){
					V.memory[i][j][k][l] = 0.0;
				}
			}
		}
	}
	int shell, length, i, j, k, l;

/* create file name then open */
	std::string infile;
	std::stringstream shel,omeg; 
	shell = object.cutoff+1;
	shel << shell;
	infile = "twobodys" + shel.str();
	omeg << object.omega;
	infile += "o" + omeg.str() + ".dat";
	std::ifstream ifile;

	ifile.open(infile);

	ifile >> length;
	for(int num=0;num<length;num++){
		ifile >> i;
		ifile >> j;
		ifile >> k;
		ifile >> l;
		ifile >> V.memory[i][j][k][l];
	}
	

	ifile.close();
}
