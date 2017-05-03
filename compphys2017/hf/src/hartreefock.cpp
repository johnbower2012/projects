#include "hartreefock.hpp"

#define TOLERANCE 1e-8
#define maxITERATIONS 500

void compute_densityMatrix(arma::mat& C, int size, int particles, arma::mat& densityMatrix){
	int i, j, k;
	double rho;
	for(i=0;i<size;i++){
		for(j=0;j<size;j++){
			rho = 0;
			for(k=0;k<particles;k++){
				rho += C(i,k)*C(j,k);
			}
			densityMatrix(i,j) = rho;
		}
	}
}
void solve_iterations(arma::mat& H0, matrix4D<double>& V, int size, int particles, arma::mat& densityMatrix, arma::vec& E){	
	int i, j, k, l,
		iterations = 0;
	double energy, time, avg;

	arma::mat 	C = arma::eye<arma::mat>(size,size),
				H = C;
				densityMatrix = H;
	arma::vec	prevE = arma::zeros<arma::vec>(size),
				diff = prevE;
	
	clock_t start, finish;
	start = clock();	

	while(iterations < maxITERATIONS){
		H = arma::zeros<arma::mat>(size,size);
		diff = arma::zeros<arma::vec>(size);
		compute_densityMatrix(C,size,particles,densityMatrix);

		for(i=0;i<size;i++){
			for(j=0;j<size;j++){
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
		avg=0.0;
/* does not function properly as a method:
		for(i=0;i<size;i++){
			avg += abs(diff(i));
		}
		avg /= (double) size;
*/
		avg = fabs(diff.max());
		std::cout << std::endl << "   Iteration  " << iterations << " with diff= " << avg << std::endl << std::setw(5) << "[";
		for(i=0;i<size;i++){
			if(i%10==0){
				std::cout << std::endl << "     ";
			}
			std::cout << std::setw(10) << E(i);
		}
		std::cout << std::endl << std::setw(5) << "]" << std::endl;
		if(avg < TOLERANCE){
			break;
		}
		iterations++;
	}

	finish = clock();
	time = (finish - start)/(double) CLOCKS_PER_SEC;
	E = arma::sort(E);

	if(iterations<maxITERATIONS){
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "Convergence in " << time << " seconds after " << iterations;
		std::cout << " iterations with " << avg << " < " << TOLERANCE << std::endl;
		std::cout << std::endl;
	}
	else{
		std::cout << std::endl;
		std::cout << std::endl;
		std::cout << "Did not converge in " << time << " seconds after " << iterations;
		std::cout << " iterations with " << avg << " > " << TOLERANCE << std::endl;
		std::cout << std::endl;
	}	
}
