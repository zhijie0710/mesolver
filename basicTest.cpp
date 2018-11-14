#include <fstream>
#include "mesolver.h"

using namespace std;
using namespace arma;

float pulse1(float t, vector<float> params) {
	return(params[0]*sin(PI*t/params[1]));
}

int main() {

	mesolver me;
	ofstream fout;

	fout.open("output.dat");

	cx_mat rho0, H0, a, ad, rho1, sigmax, sigmaz, someState;
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
	
	dt = 0.001; tmax = 15;
	coeff = {2.0};
	params = {{0.2, tmax}};
	vector<float(*)(float, vector<float>)> f = {pulse1};

	fmat data;
	someState = ones<cx_mat>(2,2);
	someState *= 0.5;

	vector<float> tvec;
	for(int i = 0; i < 10000; i++) {
		tvec.push_back(dt*i);
	}

	data = me.evolveState(rho0, H0, cOps, tOps, Ht, params, f, dt, tvec, coeff, isSparse);

	fout << data;

	fout.close();

	return 0;
}
