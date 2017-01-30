#include<iostream>
#include<iomanip>
#include<cmath>
#include<fstream>
#include "state_function.cpp"

int main(int argc,char* argv[]){
	int nx,ny,m;
	double x,y;
	if(argc<6){
		std::cout << "Bad usage. Enter also 'nx ny x y m' on same line." << std::endl;
		exit(1);
	}
	else{
		nx = atoi(argv[1]);
		ny = atoi(argv[2]);
		x = atof(argv[3]);
		y = atof(argv[4]);
		m = atoi(argv[5]);
	}

	double xder, yder;
	double xspacing = x/(double) m;
	double yspacing = y/(double) m;
	double integral = 0.0;
	double** function = new double*[2*m+1];
	for(int i=0;i<2*m+1;i++){
		function[i] = new double[2*m+1];
		for(int j=0;j<2*m+1;j++){
			function[i][j] = 0.0;
		}
	}
	double* xfunction = new double[2*m+1];
	double* yfunction = new double[2*m+1];

	harmonic_oscillator_func HOFx(nx);
	harmonic_oscillator_func HOFy(ny);
	for(int i=-m;i<m+1;i++){
		xfunction[i+m] = HOFx.value(nx,i*xspacing);
		for(int j=-m;j<m+1;j++){
			yfunction[j+m] = HOFy.value(ny,j*yspacing);
			function[i+m][j+m] = xfunction[i+m]*yfunction[j+m];
		}
	}
	for(int i=1;i<2*m;i++){
		x = (i-m)*xspacing;
		xder = (xfunction[i-1] + xfunction[i+1] - 2.0*xfunction[i])/xspacing/xspacing;
		for(int j=1;j<2*m;j++){	
			y = (j-m)*yspacing;
			yder = (yfunction[j-1] + yfunction[j+1] - 2.0*yfunction[j])/yspacing/yspacing;
			integral += function[i][j]*(0.5*(x*x+y*y)*function[i][j] - 0.5*(yfunction[j]*xder + yder*xfunction[i]));
		}
	}
	integral *= xspacing*yspacing;
	std::cout << integral << std::endl;

	for(int i=0;i<2*m+1;i++){
		delete[] function[i];
	}
	delete[] function;
	delete[] xfunction;
	delete[] yfunction;
	
	return 0;
}	
