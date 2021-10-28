#include <iostream>
#include <fstream>
#include <cmath>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

#include <ode_solver.h>

// This solver is based on a book 'Numerical Recipies in C'
// Sample regular problem: Van ner Pol oscillator

#define MU 0.4
#define OM 2

class Van_der_Pol_RK4: public RK4_stepper{
public:
	double mu;
	double om;

	Van_der_Pol_RK4(): RK4_stepper(2), om(OM), mu(MU) {};
	
	void derivatives(const double x, const double *y, double *dy) const{
		dy[0] = y[1];
		dy[1] = -om*y[0] + mu*(1 - y[0]*y[0])*y[1];
	}	
};

class Van_der_Pol_DP: public RKDP_stepper{
public:
	double mu;
	double om;

	Van_der_Pol_DP(): RKDP_stepper(2), om(OM), mu(MU) {};
	
	void derivatives(const double x, const double *y, double *dy) const{
		dy[0] = y[1];
		dy[1] = -om*om*y[0] + mu*(1 - y[0]*y[0])*y[1];
	}	
};

void save_grid(const char *file_name, const int nfunc, const int npoints, const double *xarr, const double *yarr){
	std::ofstream out;
	out.open(file_name);

	out << "x";
	for(int i = 0; i < nfunc; i++) out << "," << "y" << i;
	out << std::endl;

	const double *ycurr = yarr;
	for(int i = 0; i <= npoints; i++){
		out << xarr[i];
		for(int i = 0; i < nfunc; i++) out << "," << ycurr[i];
		out << std::endl;
		ycurr += nfunc;
	}

	out.close();
	return;
}

void integrate_VdP(int argc, char const *argv[]){
	int npoints = 5000;
	double xst = 0.0;
	double xend = 50.0;

	double *xarr = new double[npoints+1];
	double h = (xend - xst)/(1.0*npoints), xtmp = 0;
	for(int i = 0; i <= npoints; i++, xtmp += h) xarr[i] = xst + xtmp;
	double *yarr = new double[2*(npoints+1)];

	yarr[0] = 0.0;
 	yarr[1] = 0.1;

	Van_der_Pol_DP vdp;
	vdp.integrate_grid(npoints, xarr, yarr);

	save_grid("Data/solution.dat", 2, npoints, xarr, yarr);

	delete [] xarr;
	delete [] yarr;

	return;
}

// Testing accuracy

class Accuracy_DP: public RKDP_stepper{
public:
	double om;

	Accuracy_DP(): RKDP_stepper(2), om(OM) {};
	
	void derivatives(const double x, const double *y, double *dy) const{
		dy[0] = y[1];
		dy[1] = -om*om*y[0];
	}

	void solution(const double x, const double *y0, double *y) const{
		y[0] = y0[0]*std::cos(om*x) + y0[1]/om*std::sin(om*x);
		y[1] = y0[1]*std::cos(om*x) - om*y0[0]*std::sin(om*x);
	}
};

void test_accuracy(int argc, char const *argv[]){
	int nfunc = 2;
	int npoints = 5000;
	double xst = 0.0;
	double xend = 50.0;
	double h = (xend - xst)/(1.0*npoints);

	std::ofstream out_st, out_ext;
	out_st.open("Data/accuracy_st.dat"); //Accuracy of the stepper
	out_ext.open("Data/accuracy_int.dat"); //Accuracy of the extrapolator

	out_st << "x";
	out_ext << "x";
	for(int i = 0; i < nfunc; i++){
		out_st << "," << "y" << i << " solution error";
		out_ext << "," << "y" << i << " interpolation error";
	}
	out_st << std::endl;
	out_ext << std::endl;

	//Initial conditions
	double *yic = new double[nfunc];
	yic[0] = 1.0;
 	yic[1] = 1.0;

	Accuracy_DP acc_test;
	acc_test.set_initial_conditions(xst, yic);
	acc_test.set_tolerances(1.0e-10, 1.0e-10);

	double xcurr = h;
	double *ycurr = new double[nfunc];
	while(acc_test.x <= xend){
		acc_test.step();

		acc_test.solution(acc_test.x, yic, ycurr);
		out_st << acc_test.x;
		for(int i = 0; i < nfunc; i++) out_st << "," << std::fabs(ycurr[i] - acc_test.y[i]);
		out_st << std::endl;
		
		while(xcurr < acc_test.x){
			acc_test.solution(xcurr, yic, ycurr);
			out_ext << acc_test.x;
			for(int i = 0; i < nfunc; i++) out_ext << "," << std::fabs(ycurr[i] - acc_test.dense_out(i, xcurr));
			out_ext << std::endl;

			xcurr += h;
		}
	}
	
	delete [] ycurr;
	delete [] yic;
}

// This routine runs ordinary 4th order Runge-Kutta solver

int main(int argc, char const *argv[]){
	struct stat sttus = {0};
	if (stat("Data", &sttus) == -1) {
		mkdir("Data", 0700);
	}

	integrate_VdP(argc, argv);
	test_accuracy(argc, argv);
	return 0;
}