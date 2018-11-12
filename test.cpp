// #include <iostream>
#include "mesolver.h"

using namespace std;
using namespace arma;

float pulse1(float t, vector<float> params) {
	return(params[0]*sin(PI*t/params[1]));
}

float pulse2(float t, vector<float> params) {
	return(params[0]*sin(2.0*PI*t/params[1]));
}

int main() {

	mesolver me;

	cx_mat rho0, H0, a, ad, rho1, sigmax, sigmaz;
	vector<cx_mat> cOps, tOps, Ht;
	float dt, tmax;
	vector<float> coeff;
	vector<vector<float>> params;
	bool isSparse = 0;

	a = zeros<cx_mat>(2,2);
	ad = zeros<cx_mat>(2,2);
	rho0 = zeros<cx_mat>(2,2);
	rho1 = zeros<cx_mat>(2,2);
	a(0, 1) = 1; ad(1, 0) = 1;
	rho0(0, 0) = 1; rho1(1, 1) = 1;
	sigmax = a + ad;
	sigmaz = rho0 - rho1;
	H0 = sigmaz;
	cOps = {a};
	tOps = {rho0, rho1};
	Ht = {sigmax};
	dt = 0.1; tmax = 10;
	coeff = {2.3};
	params = {{0.2, tmax}};
	vector<float(*)(float, vector<float>)> f = {pulse1};

	cout << H0 << endl;

	fmat data;

	data = me.evolveState(rho0, H0, cOps, tOps, Ht, params, f, dt, tmax, coeff, isSparse);

	cout << data << endl;

	return 0;
}
