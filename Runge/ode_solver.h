#include <cmath>

//===============================================================================
//
//	These methods were mostly copied from Numerical Recipies 3rd edition
//
//===============================================================================

#ifndef ODE_SOLVER_H
#define ODE_SOLVER_H

class ODE_parameters{
public:
	// Number of functions
	const int nfunc;
	// RK stepsize
	double h;
	// Runtime values of x, the functions and derivatives
	double x;
	double *y;
	double *dy;

	ODE_parameters(int nfunci): nfunc(nfunci){
		y = new double[nfunc];
		dy = new double[nfunc];
	}

	~ODE_parameters(){
		delete [] y;
		delete [] dy;
	}

	virtual void derivatives(const double xd, const double *yd, double *dyd) const = 0;
	
	virtual void set_initial_conditions(const double xi, const double *yi){
		x = xi;
		for(int i = 0; i < nfunc; i++) y[i] = yi[i];
		derivatives(xi, yi, dy);
	}
};

class RK4_stepper: public ODE_parameters{
public:
	// Temporary arrays
	double *dym;
	double *dyt;
	double *yt;
	// Next step array
	double *ynext;

public:
	RK4_stepper(int nfunc): ODE_parameters(nfunc){
		h = 1.0e-4;//Default stepsize
		dym = new double[nfunc];
		dyt = new double[nfunc];
		yt = new double[nfunc];
		ynext = new double[nfunc];
	}

	~RK4_stepper(){
		delete [] dym;
		delete [] dyt;
		delete [] yt;
		delete [] ynext;
	}

	void rk4_step(void){
		double hh = h*0.5;
		double h6 = h/6.0;
		double xh = x + hh;

		for(int i = 0; i < nfunc; i++) yt[i] = y[i] + hh*dy[i];
		derivatives(xh, yt, dyt); //Second step
		for(int i = 0; i < nfunc; i++) yt[i] = y[i] + hh*dyt[i];
		derivatives(xh, yt, dym); //Third step
		for(int i = 0; i < nfunc; i++){
			yt[i] = y[i] + h*dym[i];
			dym[i] += dyt[i];}
		derivatives(x + h, yt, dyt);//Fourth step
		for(int i = 0; i < nfunc; i++) //Adding up
			ynext[i] = y[i] + h6*(dy[i] + dyt[i] + 2.0*dym[i]);
		return;
	}

	void step(void){
		rk4_step();
		x += h;
		for(int i = 0; i < nfunc; i++) y[i] = ynext[i];
		derivatives(x, y, dy);
	}

	// Integrates equations and saves function set on a grid xarr
	// xarr[0] and yarr[0] ... yarr[nfunc-1] serve as initial conditions
	// For every step xarr[i] ... xarr[i+1] makes substeps substeps
	void integrate_grid(const int npoints, const int substeps, const double *xarr, double *yarr){
		double *ycurr = yarr;
		set_initial_conditions(xarr[0], yarr);
		
		for(int i = 0; i < npoints; i++){
			ycurr += nfunc; 
			h = (xarr[i+1] - xarr[i])/(1.0*substeps);
			for(int j = 0; j < substeps; j++) step();
			for(int j = 0; j < nfunc; j++) ycurr[j] = y[j];
		}
	}
};

class RKDP_stepper: public ODE_parameters{
public:
	// Temporary arrays for RK steps
	double *ytemp;
	double *k2;
	double *k3;
	double *k4;
	double *k5;
	double *k6;
	// Temporary arrays for interpolation
	double *rcont1;
	double *rcont2;
	double *rcont3;
	double *rcont4;
	double *rcont5;
	// Next step function and derivative array
	double *ynext;
	double *dynext;
	// Error parameters
	double *yerr; // Error array
	double atol;
	double rtol;

	// Error controller
	struct Controller {
		double hnext, errold;
		bool reject;
		Controller() : reject(false), errold(1.0e-4) {};
		bool success(const double err, double &h);
	};
	Controller con;

	// Does a single Runge-Kutta Dormand-Price step
	void rkdp_step(void);
	// Evaluates scaled Runge-Kutta Dormand-Price error
	double scaled_error(void);

	// Dense interpolation variables and functions
	bool dense;
	double xold;
	double hdid;
	void prepare_dense(void);

public:
	RKDP_stepper(int nfunc): ODE_parameters(nfunc),
							 atol(1.0e-6), rtol(1.0e-6), dense(true)// Some default tolerance parameters
	{
		h = 1.0e-4;
		ytemp = new double[nfunc];
		k2 = new double[nfunc];
		k3 = new double[nfunc];
		k4 = new double[nfunc];
		k5 = new double[nfunc];
		k6 = new double[nfunc];

		rcont1 = new double[nfunc];
		rcont2 = new double[nfunc];
		rcont3 = new double[nfunc];
		rcont4 = new double[nfunc];
		rcont5 = new double[nfunc];

		ynext = new double[nfunc];
		dynext = new double[nfunc];
		yerr = new double[nfunc];
	}

	~RKDP_stepper(){
		delete [] ytemp;
		delete [] k2;
		delete [] k3;
		delete [] k4;
		delete [] k5;
		delete [] k6;

		delete [] rcont1;
		delete [] rcont2;
		delete [] rcont3;
		delete [] rcont4;
		delete [] rcont5;

		delete [] ynext;
		delete [] dynext;
		delete [] yerr;
	}

	// Setting tolerances
	void set_tolerances(const double atoli, const double rtoli){atol = atoli, rtol = rtoli;}

	void step();
	void fixed_step(const double hi);

	// Interpolation routines within the stepsize
	double dense_out(const int i, const double x);

	// Integrates equations and saves function set on a grid xarr
	// xarr[0] and yarr[0] ... yarr[nfunc-1] serve as initial conditions
	void integrate_grid(const int npoints, const double *xarr, double *yarr);
};

#endif