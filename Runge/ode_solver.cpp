#include <ode_solver.h>
#include <iostream>

#define EPS 1.0e-12

#ifndef FMAX
#define FMAX(a, b) (a > b ? a : b)
#endif

#ifndef FMIN
#define FMIN(a, b) (a < b ? a : b)
#endif

#ifndef SQR
#define SQR(a) a*a
#endif

//Dormand-Prince routines

void RKDP_stepper::step(){
	while(1){
		rkdp_step(); // Take a step.
		if(con.success(scaled_error(), h)) break;
		//Step rejected. Try again with reduced h set by controller.
		if(std::fabs(h) <= std::fabs(x)*EPS) throw("stepsize underflow in RKDP_stepper");
	}
	if(dense) prepare_dense(); // Step succeeded. Compute coeﬃcients for dense

	for(int i = 0; i < nfunc; i++){
		y[i] = ynext[i];
		dy[i] = dynext[i];
	}

	// Used for dense output
	xold = x;
	x += (hdid = h);
	
	h = con.hnext;
}

void RKDP_stepper::fixed_step(const double hi){
	h = hi;
	rkdp_step(); // Take a step.
	if(dense) prepare_dense(); // Step succeeded. Compute coeﬃcients for dense

	for(int i = 0; i < nfunc; i++){
		y[i] = ynext[i];
		dy[i] = dynext[i];
	}
	// Used for dense output
	xold = x;
	x += (hdid = h);
}

void RKDP_stepper::rkdp_step(){
	static const double c2 = 0.2, c3 = 0.3, c4 = 0.8, c5 = 8.0/9.0, a21 = 0.2, a31 = 3.0/40.0;
	static const double a32 = 9.0/40.0, a41 = 44.0/45.0, a42 = -56.0/15.0, a43 = 32.0/9.0, a51 = 19372.0/6561.0;
	static const double a52 = -25360.0/2187.0, a53 = 64448.0/6561.0, a54 = -212.0/729.0, a61 = 9017.0/3168.0;
	static const double a62 = -355.0/33.0, a63 = 46732.0/5247.0, a64 = 49.0/176.0, a65 = -5103.0/18656.0;
	static const double a71 = 35.0/384.0, a73 = 500.0/1113.0, a74 = 125.0/192.0, a75 = -2187.0/6784.0;
	static const double a76 = 11.0/84.0, e1 = 71.0/57600.0, e3 = -71.0/16695.0, e4 = 71.0/1920.0;
	static const double e5 = -17253.0/339200.0, e6 = 22.0/525.0, e7 = -1.0/40.0;

	// First step.
	for(int i = 0; i < nfunc; i++)
		ytemp[i] = y[i] + h*a21*dy[i];
	derivatives(x + c2*h, ytemp, k2);
	// Second step.
	for(int i=0; i < nfunc; i++)
		ytemp[i] = y[i] + h*(a31*dy[i] + a32*k2[i]);
	derivatives(x + c3*h, ytemp, k3);
	// Third step.
	for(int i = 0; i < nfunc; i++)
		ytemp[i] = y[i] + h*(a41*dy[i] + a42*k2[i] + a43*k3[i]);
	derivatives(x + c4*h, ytemp, k4);
	// Fourth step.
	for(int i = 0; i < nfunc; i++)
		ytemp[i] = y[i] + h*(a51*dy[i] + a52*k2[i] + a53*k3[i] + a54*k4[i]);
	derivatives(x + c5*h, ytemp, k5);
	// Fifth step.
	for(int i = 0; i < nfunc; i++)
		ytemp[i] = y[i] + h*(a61*dy[i] + a62*k2[i] + a63*k3[i] + a64*k4[i] + a65*k5[i]);
	double xph = x + h;
	derivatives(xph, ytemp, k6);
	// Sixth step.
	// Accumulate increments with proper weights.
	for(int i = 0; i < nfunc; i++)
		ynext[i] = y[i] + h*(a71*dy[i] + a73*k3[i] + a74*k4[i] + a75*k5[i] + a76*k6[i]);
	derivatives(xph, ynext, dynext);
	// Will also be ﬁrst evaluation for next step.
	// Estimate error as diﬀerence between fourth- and ﬁfth-order methods.
	for(int i = 0; i < nfunc; i++)
		yerr[i] = h*(e1*dy[i] + e3*k3[i] + e4*k4[i] + e5*k5[i] + e6*k6[i] + e7*dynext[i]);
}

void RKDP_stepper::prepare_dense(void){
	static const double d1 = -12715105075.0/11282082432.0;
	static const double d3 = 87487479700.0/32700410799.0;
	static const double d4 = -10690763975.0/1880347072.0;
	static const double d5 = 701980252875.0/199316789632.0;
	static const double d6 = -1453857185.0/822651844.0;
	static const double d7 = 69997945.0/29380423.0;
	double ydiff, bspl;
	
	for(int i=0; i < nfunc; i++){
		rcont1[i] = y[i];
		ydiff = ynext[i] - y[i];
		rcont2[i] = ydiff;
		bspl = h*dy[i] - ydiff;
		rcont3[i] = bspl;
		rcont4[i] = ydiff - h*dynext[i] - bspl;
		rcont5[i] = h*(d1*dy[i] + d3*k3[i] + d4*k4[i] + d5*k5[i] + d6*k6[i] + d7*dynext[i]);
	}
}

double RKDP_stepper::dense_out(const int i, const double xi){
	double s = (xi - xold)/hdid;
	if(s < 0 || s > 1.0){
		std::cout << "Trying to interpolate solution outside the step" << std::endl;
		throw("Trying to interpolate solution outside the step");
	}
	double s1 = 1.0 - s;
	return rcont1[i] + s*(rcont2[i] + s1*(rcont3[i] + s*(rcont4[i] + s1*rcont5[i])));
}

double RKDP_stepper::scaled_error(void){
	double err = 0.0, sk;
	for(int i = 0; i < nfunc; i++){
		sk = atol + rtol*FMAX(std::fabs(y[i]), std::fabs(ynext[i]));
		err += SQR(yerr[i]/sk);
	}
	return std::sqrt(err/(1.0*nfunc));
}

// Returns true if err \leq 1, false otherwise. If step was successful, sets hnext to the estimated
// optimal stepsize for the next step. If the step failed, reduces hc appropriately for another try.
bool RKDP_stepper::Controller::success(const double err, double &hc){
	static const double beta = 0.0, alpha = 0.2 - beta*0.75, safe = 0.9;
	static const double minscale = 0.2, maxscale = 10.0;

	// Set beta to a nonzero value for PI control. beta D 0:04–0.08 is a good default.
	double scale;
	if(err <= 1.0){ // Step succeeded. Compute hnext.
		if(err == 0.0) scale = maxscale;
		else{// PI control if beta \neq 0.
			scale = safe*std::pow(err, -alpha)*std::pow(errold, beta);
			if(scale < minscale) scale = minscale; // Ensure minscale \leq hnext = hc \leq maxscale.
			if(scale > maxscale) scale = maxscale;
		}
	
		if(reject) hnext = hc*FMIN(scale, 1.0); //Don’t let step increase if last one was rejected.
		else hnext = hc*scale;
	
		// Bookkeeping for next call.
		errold = FMAX(err, 1.0e-4);
		reject = false;
		return true;
	} 
	else{// Truncation error too large, reduce stepsize.
		scale = FMAX(safe * std::pow(err, -alpha), minscale);
		hc *= scale;
		reject = true;
		return false;
	}
}

// Integrates equations and saves function set on a grid xarr
// xarr[0] and yarr[0] ... yarr[nfunc-1] serve as initial conditions
void RKDP_stepper::integrate_grid(const int npoints, const double *xarr, double *yarr){
	double *ycurr = yarr;
	set_initial_conditions(xarr[0], yarr);

	int st = 1;
	while(x + h < xarr[npoints]){
		step();
		// std::cout << "Made a step to x = " << x << std::endl;
		while(x > xarr[st]){
			// std::cout << "Saving step at x = " << xarr[st] << std::endl;
			ycurr += nfunc;
			for(int i = 0; i < nfunc; i++) ycurr[i] = dense_out(i, xarr[st]);
			st++;
		}
	}

	fixed_step(xarr[npoints] - x);
	// std::cout << "Made a step to x = " << x << std::endl;
	for(; st <= npoints; st++){
		// std::cout << "Saving step at x = " << xarr[st] << std::endl;
		ycurr += nfunc;
		for(int i = 0; i < nfunc; i++) ycurr[i] = dense_out(i, xarr[st]);
	}
}