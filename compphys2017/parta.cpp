#include<iostream>
#include<iomanip>
#include<stdlib.h>

class stateset{
		int cutoff, dimension, spin,
			info;
		double hbar, omega;
	public:
		int** state;
		int states;

		stateset(int,int,int);
		~stateset();

		void print();
};
class numbers{
	public:
		int value;

		numbers(){value=0;}

		int factorial(int);
		int choose(int,int);
};
stateset::stateset(int n, int d, int s){
	hbar = 1.0;
	omega = 1.0;

	info = 2+d;
	cutoff = n;
	dimension = d;
	spin = s;

	numbers calc;
	states=0;
	for(int i=0;i<cutoff+1;i++){
		states += calc.choose(i+d-1,i);
	}
	states *= (spin+1);

	int count=0, ecount, dcount, penergy, position,
		energy, degeneracy, check;
	int* energies;

	state = new int*[states];
	for(int i=0;i<states;i++){
		state[i] = new int[info];
	}
	for(int i=0;i<n+1;i++){
		energy = i;
		degeneracy=calc.choose(energy+dimension-1,energy);
		energies = new int[energy+1];
		for(int j=0;j<energy;j++){
			energies[j]=j;
		}
		energies[energy]=-1;
	//Read in energy distribution from energies
		for(int deg=0;deg<degeneracy;deg++){
			dcount=0; ecount=0; position=0;
			while(dcount<dimension){
				penergy=0;
				while(energies[ecount]==position){
					position++;
					ecount++;
					penergy++;
				}
				state[count][dcount]=penergy;
				position++;
				dcount++;
			}
			state[count][dimension]=spin;
			state[count][dimension+1]=-spin;
			count++;

			for(int j=0;j<spin;j++){
				for(int k=0;k<dimension;k++){
					state[count][k]=state[count-1][k];
				}
				state[count][dimension]=spin;
				state[count][dimension+1]=state[count-1][dimension+1]+2;
				count++;
			}
			check=0;
			if(i>0){
				if(energies[i-1]<i+dimension-2){
					energies[i-1]++;
				}
				else{
					while(energies[i-1-check]==i+dimension-2-check){
						check++;
						if(check==i-1){break;}
					}
					energies[i-1-check]++;
					for(int j=0;j<check;j++){
						energies[i-check+j] = energies[i-check+j-1] + 1;
					}
				}
			}
		}
	}
}
stateset::~stateset(){
	for(int i=0;i<states;i++){
		delete[] state[i];
	}
	delete[] state;
	cutoff=0;
	info=0;
	states=0;
}
void stateset::print(){
	double TE;
	std::cout << "E cut off:	" << cutoff << std::endl;
	std::cout << "State count: 	" << states << std::endl;
	std::cout << std::endl;
	std::cout << std::setw(8) << "TE";
	for(int i=0;i<dimension;i++){
		std::cout << std::setw(7) << "n" << i+1;
	}
	std::cout << std::setw(8) << "spin" << std::setw(8) << "ms" << std::endl;
	std::cout << std::endl;
	for(int i=0;i<states;i++){
		TE = 0.0;
		for(int j=0;j<dimension;j++){
			TE += state[i][j];
		}
		TE += dimension/2.0;
		TE *= hbar*omega;
		std::cout << std::setw(8) << TE;
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
	
int main(int argc, char* argv[]){
	int n,d,s;
	if(argc<4){
		std::cout << "Bad usage. Enter also 'n d s' on same line." << std::endl;
		exit(1);
	}
	else{
		n = atoi(argv[1]);
		d = atoi(argv[2]);
		s = atoi(argv[3]);
	}

	stateset test(n,d,s);
	test.print();

	return 0;
}
