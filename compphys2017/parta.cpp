#include<iostream>
#include<iomanip>
#include<stdlib.h>

int main(int argc, char* argv[]){
	int energy, n, nx, ny, spin, ms,
		states, info=5, degeneracy, count=0;

	if(argc<2){
		std::cout << "Bad usage. Enter also 'energy_cut_off==n' on same line." << std::endl;
		exit(1);
	}
	else{
		n = atoi(argv[1]);
		states = n*(n+1);
	}

	int** state = new int*[states];
	for(int i=0;i<states;i++){
		state[i] = new int[info];
	}
	for(int i=0;i<n;i++){
		energy = i+1;
		degeneracy = 2*(energy);
		for(int j=0;j<degeneracy;j++){
			state[count][0] = energy;
			state[count][1] = energy - j/2 - 1;
			state[count][2] = j/2;
			state[count][3] = 1;
			state[count][4] = (count%2)*2 - 1;
			count++;
		}
	}

	std::cout << "E cut off:	" << n << std::endl;
	std::cout << "State count: 	" << states << std::endl;
	std::cout << std::endl;
	std::cout << std::setw(5) << "TE" << std::setw(5) << "nx" << std::setw(5) << "ny" << std::setw(5) << "spin" << std::setw(5) << "ms" << std::endl;
	std::cout << std::endl;

	for(int i=0;i<states;i++){
		for(int j=0;j<info;j++){
			std::cout << std::setw(5) << state[i][j];
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;

	for(int i=0;i<states;i++){
		delete[] state[i];
	}
	delete[] state;

	return 0;
}
