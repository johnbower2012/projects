#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include "time.h"

double func_prob(double alpha, double rho){
	double value = alpha*rho;
	return value*value*exp(-2*value);
}
double func_energy(double alpha, double rho){
	double value = 1.0/rho;
	return -value + alpha*(value-alpha/2.0);
}

int main(int argc, char* argv[]){

	double	alpha, random, step, variance, E, ratio,
			rho, rho_new, prob, prob_new, energy, energy_new,
			e=0.0, esq=0.0,
			time;

	int 	mcsteps, accsteps=0;
	unsigned seed_move, seed_accept;
	clock_t start, finish;


	if(argc<4){
		std::cout << "\tBad usage. Enter also 'alpha step mcsteps' on same line." << std::endl;
		std::cout << "\talpha==0.0 will loop over alpha[0.6,1.5]" << std::endl;
		exit(1);
	}
	else{
		alpha = atof(argv[1]);
		step = atof(argv[2]);
		mcsteps = atoi(argv[3]);
	}

	seed_move = std::chrono::system_clock::now().time_since_epoch().count();
	seed_accept = std::chrono::system_clock::now().time_since_epoch().count();

	std::default_random_engine generator_move (seed_move);
	std::default_random_engine generator_accept (seed_accept);

	std::uniform_real_distribution<double> dist_move(-0.5,0.5);
	std::uniform_real_distribution<double> dist_accept(0.0,1.0);

	rho = 10.0*(dist_move(generator_move)+0.5);

	prob = func_prob(alpha, rho);
	energy = func_energy(alpha, rho);

	if(alpha==0.0){
		for(int i=0;i<10;i++){
			accsteps = 0.0;
			esq = e = 0.0;
			alpha = 0.6 + (double) i/10.0;

			start = clock();
			for(int j=0;j<mcsteps;j++){
				random = dist_move(generator_move);
				rho_new = fabs(rho + random*step);
				prob_new = func_prob(alpha,rho_new);
				ratio = prob_new/prob;
				if(ratio<1.0){
					random = dist_accept(generator_accept);
					if(random > ratio){
						rho_new = rho;
						prob_new = prob;
						energy_new = energy;
					}
					else{
						rho = rho_new;
						prob = prob_new;	
						energy = energy_new = func_energy(alpha, rho);
						accsteps++;
					}
				}
				else{
					rho = rho_new;
					prob = prob_new;
					energy = energy_new = func_energy(alpha, rho);
					accsteps++;
				}
				e += energy;
				esq += energy*energy;
			}
			finish = clock();
			time = (finish - start)/(double) CLOCKS_PER_SEC;

			e /= (double) mcsteps;
			esq /= (double) mcsteps;
			variance = esq - e*e;

			std::cout << std::endl;
			std::cout << std::setw(15) << "accepted steps" << std::setw(15) << "alpha" << std::setw(15) << "<E>";
			std::cout << std::setw(15) << "var" << std::setw(15) << "var/sqrt(N)" << std::setw(15) << "time" << std::endl;
			std::cout << std::setw(15) << (double) accsteps/(double) mcsteps << std::setw(15) << alpha << std::setw(15) << e << std::setw(15);
			std::cout << variance << std::setw(15) << variance/sqrt((double) mcsteps) << std::setw(15) << time << "s" << std::endl;
			std::cout << std::endl;
		}
	}
	else{
		esq = e = 0.0;
		start = clock();
		for(int j=0;j<mcsteps;j++){
			random = dist_move(generator_move);
			rho_new = fabs(rho + random*step);
			prob_new = func_prob(alpha,rho_new);
			ratio = prob_new/prob;
			if(ratio<1.0){
		 		random = dist_accept(generator_accept);
		 		if(random > ratio){
					rho_new = rho;
					prob_new = prob;
					energy_new = energy;
			 	}
				else{
					rho = rho_new;
					prob = prob_new;	
					energy = energy_new = func_energy(alpha, rho);
					accsteps++;
				}
			}
			else{
				rho = rho_new;
				prob = prob_new;
				energy = energy_new = func_energy(alpha, rho);
				accsteps++;
			}

			e += energy;
			esq += energy*energy;

		}
		finish = clock();
		time = (finish - start)/(double) CLOCKS_PER_SEC;

		e /= (double) mcsteps;
		esq /= (double) mcsteps;
		variance = esq - e*e;

		std::cout << std::endl;
		std::cout << std::setw(15) << "accepted steps" << std::setw(15) << "alpha" << std::setw(15) << "<E>";
		std::cout << std::setw(15) << "var" << std::setw(15) << "var/sqrt(N)" << std::setw(15) << "time" << std::endl;
		std::cout << std::setw(15) << accsteps << std::setw(15) << alpha << std::setw(15) << e << std::setw(15);
		std::cout << variance << std::setw(15) << variance/sqrt((double) mcsteps) << std::setw(15) << time << "s" << std::endl;
		std::cout << std::endl;
	}

	return 0;
}	
