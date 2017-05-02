#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<string>
#include "time.h"

#define STEP 2.0

double func_wave(double x[6]){
	double beta = x[5], alpha = x[4], omega = x[3],
			r1 = x[0], r2 = x[1], r12 = x[2];
	return exp(-alpha*omega*(r1*r1 + r2*r2)/2.0)*exp(r12/(1.0+beta*r12));
}
double func_prob(double x[6]){
	double beta = x[5], alpha = x[4], omega = x[3],
			r1 = x[0], r2 = x[1], r12 = x[2];
	return exp(-alpha*omega*(r1*r1 + r2*r2))*exp(2.0*r12/(1.0+beta*r12));
}
double distance(double x, double y){
	return sqrt(x*x + y*y);
}
double func_energy(double x[6]){
	double beta = x[5], alpha = x[4], omega = x[3],
			r1 = x[0], r2 = x[1], r12 = x[2],
			alphasq = alpha*alpha,
			rsq = r1*r1 + r2*r2,
			onebetar12 = 1.0/(1.0 + beta*r12),
			obr12sq = onebetar12*onebetar12,
			omegasq = omega*omega,
			a = 1.0;
	return 2.0*alphasq*omegasq*rsq - 4.0*alpha*omega - obr12sq*(2.0*alpha*a*r12*omega + 2.0*a*(a*obr12sq + 1.0/r12 - 2.0*beta*onebetar12));
}

int main(int argc, char* argv[]){

	double	r[6], r_new[6], x[4], x_new[4], X, Y, r12, 
			alpha, beta, omega, random, step, variance, E, ratio,
			rho, rho_new, prob, prob_new, energy, energy_new,
			e=0.0, esq=0.0,
			time;

	std::string choice;
	int 	mcsteps, accsteps=0, particles = 2, dim = 2;
	unsigned seed_move, seed_accept;
	clock_t start, finish;


	if(argc<6){
		std::cout << "\tBad usage. Enter also 'alpha beta omega step mcsteps' on same line." << std::endl;
		std::cout << "\talpha==0.0 or beta==0.0 will loop over alpha[0.6,1.5] or beta[0.6,1.5]" << std::endl;
		exit(1);
	}
	else{
		alpha = atof(argv[1]);
		beta = atof(argv[2]);
		omega = atof(argv[3]);
		step = atof(argv[4]);
		mcsteps = atoi(argv[5]);
	}

	seed_move = std::chrono::system_clock::now().time_since_epoch().count();
	seed_accept = std::chrono::system_clock::now().time_since_epoch().count();

	std::default_random_engine generator_move (seed_move);
	std::default_random_engine generator_accept (seed_accept);

	std::uniform_real_distribution<double> dist_move(-0.5,0.5);
	std::uniform_real_distribution<double> dist_accept(0.0,1.0);

	if(alpha==0.0){	
		for(int loop=0;loop<10;loop++){
			alpha = 0.6 + (double) loop/10.0;
			for(int i=0;i<particles;i++){
				for(int j=0;j<dim;j++){
					x[dim*i+j] = x_new[dim*i+j] = STEP*dist_move(generator_move);
				}
				X = x[dim*i]; Y = x[dim*i+1];
				r[i] = r_new[i] = distance(X,Y);
			}
			r[3] = r_new[3] = omega;
			r[4] = r_new[4] = alpha;
			r[5] = r_new[5] = beta;
			prob = func_prob(r);
			energy = func_energy(r);

			accsteps = 0;
			esq = e = 0.0;
			start = clock();
			for(int i=0;i<mcsteps;i++){
				for(int j=0;j<particles;j++){
					for(int k=0;k<dim;k++){
						x_new[dim*j+k] = step*dist_move(generator_move);
					}
				
					X = x_new[dim*j]; Y = x_new[dim*j+1];
					r_new[j] = distance(X,Y);

					prob_new = func_prob(r_new);
					ratio = prob_new/prob;

					if(ratio<1.0){
			 			random = dist_accept(generator_accept);
			 			if(random > ratio){
							for(int k=0;k<dim;k++){
								x_new[dim*j+k] = x[dim*j+k];
							}
							r_new[j] = r[j];
							prob_new = prob;
							energy_new = energy;
					 	}
						else{
							for(int k=0;k<dim;k++){
								x_new[dim*j+k] = x[dim*j+k];
							}
							r[j] = r_new[j];
							prob = prob_new;	
							energy = energy_new = func_energy(r);
							accsteps++;
						}
					}
					else{
						for(int k=0;k<dim;k++){
							x_new[dim*j+k] = x[dim*j+k];
						}
						r[j] = r_new[j];
						prob = prob_new;	
						energy = energy_new = func_energy(r);
						accsteps++;
					}
	
					e += energy;
					esq += energy*energy;
				}
			}
			finish = clock();
			time = (finish - start)/(double) CLOCKS_PER_SEC;

			e /= (double) (mcsteps*particles);
			esq /= (double) (mcsteps*particles);
			variance = esq - e*e;

			std::cout << std::endl;
			std::cout << std::setw(15) << "accepted steps" << std::setw(15) << "alpha" << std::setw(15) << "<E>";
			std::cout << std::setw(15) << "var" << std::setw(15) << "var/sqrt(N)" << std::setw(15) << "time" << std::endl;
			std::cout << std::setw(15) << (double) accsteps/(double) (mcsteps*particles) << std::setw(15) << alpha << std::setw(15) << e << std::setw(15);
			std::cout << variance << std::setw(15) << variance/sqrt((double) (mcsteps*particles)) << std::setw(15) << time << "s" << std::endl;
			std::cout << std::endl;
		}
	}
	else{				
		for(int i=0;i<particles;i++){
			for(int j=0;j<dim;j++){
				x[dim*i+j] = x_new[dim*i+j] = STEP*dist_move(generator_move);
			}
			X = x[dim*i]; Y = x[dim*i+1];
			r[i] = r_new[i] = distance(X,Y);
		}
		for(int i=0;i<particles-1;i++){
			for(int j=i+1;j<particles;j++){
				X = x[dim*i] - x[dim*j];
				Y = x[dim*i+1] - x[dim*j+1];
				r12 = distance(X,Y);
			}
		}

		for(int i=0;i<particles-1;i++){
			for(int j=i+1;j<particles;j++){
				X = x[dim*i] - x[dim*j];
				Y = x[dim*i+1] - x[dim*j+1];
				r12 = distance(X,Y);
			}
		}
		r[2] = r12;
		r[3] = r_new[3] = omega;
		r[4] = r_new[4] = alpha;
		r[5] = r_new[5] = beta;

		prob = func_prob(r);
		energy = func_energy(r);

		accsteps = 0;
		esq = e = 0.0;
		start = clock();
		for(int i=0;i<mcsteps;i++){
			for(int j=0;j<particles;j++){
				for(int k=0;k<dim;k++){
					x_new[dim*j+k] = step*dist_move(generator_move);
				}
			
				X = x_new[dim*j]; Y = x_new[dim*j+1];
				r_new[j] = distance(X,Y);

				for(int k=0;k<particles-1;k++){
					for(int l=l+1;l<particles;l++){
						X = x[dim*k] - x[dim*l];
						Y = x[dim*k+1] - x[dim*l+1];
						r12 = distance(X,Y);
					}
				}
				r_new[2] = r12;

				prob_new = func_prob(r_new);
				ratio = prob_new/prob;

				if(ratio<1.0){
		 			random = dist_accept(generator_accept);
		 			if(random > ratio){
						for(int k=0;k<dim;k++){
							x_new[dim*j+k] = x[dim*j+k];
						}
						r_new[j] = r[j];
						r_new[2] = r[2];
						prob_new = prob;
						energy_new = energy;
				 	}
					else{
						for(int k=0;k<dim;k++){
							x_new[dim*j+k] = x[dim*j+k];
						}
						r[j] = r_new[j];
						r[2] = r_new[2];
						prob = prob_new;	
						energy = energy_new = func_energy(r);
						accsteps++;
					}
				}
				else{
					for(int k=0;k<dim;k++){
						x_new[dim*j+k] = x[dim*j+k];
					}
					r[j] = r_new[j];
					r[2] = r_new[2];
					prob = prob_new;	
					energy = energy_new = func_energy(r);
					accsteps++;
				}

				e += energy;
				esq += energy*energy;
			}
		}
		finish = clock();
		time = (finish - start)/(double) CLOCKS_PER_SEC;

		e /= (double) (mcsteps*particles);
		esq /= (double) (mcsteps*particles);
		variance = esq - e*e;

		std::cout << std::endl;
		std::cout << std::setw(15) << "accepted steps" << std::setw(15) << "alpha" << std::setw(15) << "<E>";
		std::cout << std::setw(15) << "var" << std::setw(15) << "var/sqrt(N)" << std::setw(15) << "time" << std::endl;
		std::cout << std::setw(15) << (double) accsteps/(double) (mcsteps*particles) << std::setw(15) << alpha << std::setw(15) << e << std::setw(15);
		std::cout << variance << std::setw(15) << variance/sqrt((double) (mcsteps*particles)) << std::setw(15) << time << "s" << std::endl;
		std::cout << std::endl;
	}

	return 0;
}	
