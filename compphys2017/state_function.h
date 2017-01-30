class hermite_polynomial{
		int	H00, H10, H11;
	public:
		int** coeff;
		int degree;

		hermite_polynomial(int);
		~hermite_polynomial();
	
		double value(int n, double x);	
		void print();
};
class factorial{
	public:
		int n, value;

		factorial(){n=1;value=1;};
		factorial(int a);

		int result();
};
class harmonic_oscillator_func{
		double mass,
				omega,
				pi = acos(-1),
				hbar=1;
	public:
		int degree;
		int** hermite_coeff;
		double* prefactors;

		harmonic_oscillator_func();
		harmonic_oscillator_func(int);
		harmonic_oscillator_func(int,double,double);
		~harmonic_oscillator_func();

		double value(int, double);	
};
struct particle{
	public:
		double	mass,
				charge;
		int		spin;
};
class state_function{
		int pnum,
			dnum;
	public:
		particle* part;
		harmonic_oscillator_func* HOF;
};	
