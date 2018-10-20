#include<stdio.h>
#include <iostream>
#include "mesolver.h"


void RK4(H, rho, dt, *cOps, *coeff) {
  .......
}

void lindbladME(H, rho, *cOps, *coeff) {
  .......
}



//inclue the Armadillo c++ library
//define lindbladME(t,tf.xf):
cx_mat lindbladME(double t, cx_mat H, cx_mat rho,*cOps, *coeff)
{
	const   complex<double> i(0.0,1.0); 
	cx_mat mx;
	mx = -i/h_bar*[H,rho]+(*coeff)*[(*cOps)*rho*(*cOps.dag())-1/2*{(*cOps.dag())*(*cOps),rho}];//this is pseudocode,,,,,
	return(mx);
}

//define the Runge-Kutta 4-th order method:
cx_mat RK4(double t, double tf, cx_mat xf, double dt, cx_mat H1, cx_mat H2)
{
	cx_mat lilndbladME(double t,cx_mat H, cx_mat rho,*cOps, *coeff);
  lilndbladME(t,H,rho,*cOps, *coeff)
	cx_mat k1 = lilndbladME(t,H,rho,*cOps, *coeff) * hs;
	cx_mat rho1 = rho + k1 * 0.5;
	cx_mat k2 = lilndbladME(t + hs/2, H, rho1, *cOps, *coeff) * hs;
	cx_mat rho2 = rho + k2 * 0.5;
	cx_mat k3 = fx(t + hs/2, H, rho2, *cOps, *coeff) * hs;
	cx_mat rho3 = rho + k3;
	cx_mat k4 = fx(t + hs, H, rho3, *cOps, *coeff) * hs;
	cx_mat xx = rho + (k1 + 2*(k2 + k3) + k4)/6;

	return(xx); 
}
