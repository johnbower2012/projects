#include "hartreefock.hpp"

void sp_energies(arma::mat& input, int size, arma::mat& H0){
	int i;
	for(i=0;i<size;i++){
		H0(i,i) = 10*(2*input(i,0) + abs(input(i,1) + 1.0);
	}
}
void compute_densityMatrix(arma::mat& C, int size, arma::mat& densityMatrix){
	int i, j, k;
	double rho;	
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			rho = 0;
			for(k=0;k<part;k++){
				rho += C(i,k)*C(j,k);
			}
			densityMatrix(i,j) = rho;
		}
	}
}
void solve_iterations(arma::mat& H0, matrix4D<double>& V, int size, arma::vec& mask, arma::mat& H, arma::vec& E){	
	int i, j, k, l,
		I, J, K, L,
		iterations = 0;
	double energy, time;

	arma::mat 	C = arma::eye<arma::mat>(size,size),
				densityMatrix = C;
	arma::vec	prevE = arma::zeros<arma::vec>(size),
				diff = prevE;
	
	clock_t start, finish;
	start = clock();	

	while(iterations < maxITERATIONS){
		H = arma::zeros<arma::mat>(size,size);
		diff = arma::zeros<arma::vec>(size);
		compute_densityMatrix(C,size,densityMatrix);
		for(i=0;i<size;i++){
			for(j=i;j<size;j++){
				energy = 0;
				for(k=0;k<size;k++){
					for(l=0;l<size;l++){
						energy += densityMatrix(k,l)*V.memory[i][k][j][l];
					}
				}
				H(i,j) = H(j,i) = energy + H0(i,j);
			}
		}
		arma::eig_sym(E,C,H);
		diff = prevE - E;
		prevE = E;
		if(fabs(diff.max()) < TOLERANCE){
			break;
		}
		iterations++;
	}

	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	std::cout << std::endl;
	std::cout << "Convergence in " << time << " seconds after " << iterations;
	std::cout << " iterations with " << fabs(diff.max()) << " < " << TOLERANCE << std::endl;
	std::cout << std::endl;
	
	E = arma::sort(E);
}
