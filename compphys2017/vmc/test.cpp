#include<iostream>
#include<iomanip>
#include<cmath>
#include<cstdlib>

double func1(double x[2]){
	return x[0]*x[0]*x[0]*x[1];
}
double func2(double x[2]){
	return 3.0*x[1]*x[0]*x[0];
}
double func3(double x[2]){
	return 6.0*x[1]*x[0];
}
double der_first(double func(double[2]), double x[2], double step){
	double xp[2], xm[2];
	xp[1] = xm[1] = x[1];
	xp[0] = x[0] + step;
	xm[0] = x[0] - step;	
	return (func(xp) - func(xm))/step/2.0;
}
double der_second(double func(double[2]), double x[2], double step){
	double xp[2], xm[2];
	xp[1] = xm[1] = x[1];
	xp[0] = x[0] + step;
	xm[0] = x[0] - step;
	return (func(xp) - 2.0*func(x)+func(xm))/step/step;
}/*
double func_energy1(double x[2], double h){
	return -0.5*der_second(func1,x,h) - func1(x)/x;
}*/

int main(int argc, char* argv[]){
	double x[2], step;

	if(argc<4){
		std::cout << "Bad usage. Enter also 'x alpha h' on same line." << std::endl;
		exit(1);
	}
	else{
		x[0] = atof(argv[1]);
		x[1] = atof(argv[2]);
		step = atof(argv[3]);
	}	
	double xp[2], xm[2];
	xp[1] = xm[1] = x[1];
	xp[0] = x[0] + step;
	xm[0] = x[0] - step;

	std::cout << func1(x) << std::endl;
	std::cout << func2(x) << std::endl;
	std::cout << func3(x) << std::endl;
	std::cout << std::endl;
	std::cout << der_first(func1,x,step) << std::endl;
	std::cout << der_second(func1,x,step) << std::endl;

	return 0;
}
