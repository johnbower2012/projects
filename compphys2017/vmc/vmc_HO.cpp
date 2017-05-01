#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>
#include<random>
#include<chrono>
#include<string>
#include "time.h"

#define STEP 2.0

double func_wave(double x[4]){
	double alpha = x[3], omega = x[2],
			r1 = x[0], r2 = x[1];
	return exp(-alpha*omega*(r1*r1 + r2*r2)/2.0);
}
double func_prob(double x[4]){
	double alpha = x[3], omega = x[2],
			r1 = x[0], r2 = x[1];
	return exp(-alpha*omega*(r1*r1 + r2*r2));
}
double distance(double x, double y){
	return sqrt(x*x + y*y);
}
double func_energy(double x[4]){
	double alpha = x[3], omega = x[2],
			r1 = x[0], r2 = x[1],
			oneminusalphasq = 1.0 - alpha*alpha,
			rsq = r1*r1 + r2*r2;
	return 0.5*omega*omega*rsq*oneminusalphasq + 2.0*alpha*omega;
}
/*
double func_energy_rep(double x[4], double r12){
	double alpha = x[3], omega = x[2],
			r1 = x[0], r2 = x[1],
			oneminusalphasq = 1.0 - alpha*alpha,
			rsq = r1*r1 + r2*r2;
	return 0.5*omega*omega*rsq*oneminusalphasq + 2.0*alpha*omega + 1.0/r12;
}
if(repulsion==true){
	for(int i=0;i<particles-1;i++){
		for(int j=i+1;j<particles;j++){
			X = x[dim*i] - x[dim*j];
			Y = x[dim*i+1] - x[dim*j+1];
			r12 = distance(X,Y);
		}
	}
	energy = energy_new = func_energy_rep(r,r12);
}
else{
	energy = energy_new = func_energy(r);
}
*/

int main(int argc, char* argv[]){

	double	r[4], r_new[4], x[4], x_new[4], X, Y, r12, 
			alpha, omega, random, step, variance, E, ratio,
			rho, rho_new, prob, prob_new, energy, energy_new,
			e=0.0, esq=0.0,
			time;

	std::string choice;
	int 	mcsteps, accsteps=0, particles = 2, dim = 2;
	unsigned seed_move, seed_accept;
	clock_t start, finish;


	if(argc<5){
		std::cout << "\tBad usage. Enter also 'alpha omega step mcsteps' on same line." << std::endl;
		std::cout << "\talpha==0.0 will loop over alpha[0.6,1.5]" << std::endl;
		exit(1);
	}
	else{
		alpha = atof(argv[1]);
		omega = atof(argv[2]);
		step = atof(argv[3]);
		mcsteps = atoi(argv[4]);
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
			r[2] = r_new[2] = omega;
			r[3] = r_new[3] = alpha;
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

		r[2] = r_new[2] = omega;
		r[3] = r_new[3] = alpha;

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

	return 0;
}	
