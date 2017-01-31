#include<iostream>
#include<iomanip>
#include<stdlib.h>
#include<cmath>
#include<fstream>
#include "state_function.h"



harmonic_oscillator_func::harmonic_oscillator_func(){
	mass = 0;
	omega = 0;
	degree = 0;
	hermite_coeff = NULL;
	prefactors = NULL;
}
harmonic_oscillator_func::harmonic_oscillator_func(int n){
	mass = 1;
	omega = 1;
	degree = n;
	hermite_polynomial temp(degree+1);
	hermite_coeff = new int*[degree+1];
	for(int i=0;i<degree+1;i++){
		hermite_coeff[i] = new int[degree+1];
		for(int j=0;j<degree+1;j++){
			hermite_coeff[i][j] = temp.coeff[i][j];
		}
	}

	prefactors = new double[degree+1];
	prefactors[0] = pow(mass*omega/pi/hbar,0.25);
	for(int i=1;i<degree+1;i++){
		prefactors[i] = prefactors[i-1]/sqrt(2.0*(double)i);
	}
}
harmonic_oscillator_func::harmonic_oscillator_func(int n, double m, double o){
	mass = m;
	omega = o;
	degree = n;
	hermite_polynomial temp(degree+1);
	hermite_coeff = new int*[degree+1];
	for(int i=0;i<degree+1;i++){
		hermite_coeff[i] = new int[degree+1];
		for(int j=0;j<degree+1;j++){
			hermite_coeff[i][j] = temp.coeff[i][j];
		}
	}

	prefactors = new double[degree+1];
	prefactors[0] = pow(mass*omega/pi/hbar,0.25);
	for(int i=1;i<degree+1;i++){
		prefactors[i] = prefactors[i-1]/sqrt(2.0*(double)i);
	}
}
harmonic_oscillator_func::~harmonic_oscillator_func(){
	mass = 0;
	omega = 0;
	for(int i=0;i<degree+1;i++){
		delete[] hermite_coeff[i];
	}
	delete[] hermite_coeff;
	delete[] prefactors;
	degree = 0;
	hermite_coeff = NULL;
	prefactors = NULL;
}
double harmonic_oscillator_func::value(int n, double x){
	double value=0;
	double X = sqrt(mass*omega/hbar)*x;
	for(int i=0;i<n+1;i++){
		value += hermite_coeff[n][i]*pow(X,i);
	}
	value *= prefactors[n];
	value *= exp(-X*X/2.0);
	return value;
} 
hermite_polynomial::hermite_polynomial(int n){
	H00 = 1; H10 = 0; H11 = 2;
	degree = n + 1;
	coeff = new int*[degree];
	for(int i=0;i<degree;i++){
		coeff[i] = new int[degree];
		for(int j=0;j<degree;j++){
			coeff[i][j] = 0.0;
		}
	}

	if(n==0){
		coeff[0][0] = H00;
	}
	else if(n>0){
		coeff[0][0] = H00;
		coeff[1][0] = H10;
		coeff[1][1] = H11;
		for(int i=2;i<degree;i++){
			coeff[i][0] = -2*(i-1)*coeff[i-2][0];
			for(int j=0;j<i+1;j++){
				coeff[i][j] = 2*coeff[i-1][j-1] - 2*(i-1)*coeff[i-2][j];
			}
		}
	}
}
hermite_polynomial::~hermite_polynomial(){
	for(int i=0;i<degree;i++){
		delete[] coeff[i];
	}
	delete[] coeff;
}
double hermite_polynomial::value(int n, double x){
	double value = 0;
	int parity = n%2;
	for(int i=parity;i<n+1;i+=2){
		value += coeff[n][i]*pow(x,i);
	}
	return value;
}
void hermite_polynomial::print(){
	std::cout << std::endl;
	for(int i=0;i<degree;i++){
		for(int j=0;j<degree;j++){
			std::cout << std::setw(5) << coeff[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}
factorial::factorial(int a){
	int x = 1;
	n=a;

	if(a<0){
		x=0;
	}
	else if(a==0){
		x=1;
	}
	else{
		for(int i=1;i<a+1;i++){
			x*=i;
		}
	}
	value=x;
}
int factorial::result(){
	return value;
}