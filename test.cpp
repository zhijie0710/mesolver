#include <fstream>
#include "timeEvolution.h"

using namespace std;
using namespace arma;

float pulse1(float t, vector<float> params) {
	return(params[0]*sin(PI*t/(2*params[1])));
}

int main() {

	timeEvolution te;

	cx_mat rho0, H0, a, ad, rho1, sigmax, sigmaz;
	cx_vec ket0(2), ket1(2);
	vector<cx_mat> cOps, tOps, Ht;
	float dt, tmax;
	vector<float> coeff;
	vector<vector<float>> params;

	a = zeros<cx_mat>(2,2);
	ad = zeros<cx_mat>(2,2);

	rho0 = zeros<cx_mat>(2,2);
	rho1 = zeros<cx_mat>(2,2);
	ket0 = {1, 0};
	ket1 = {0, 1};

	a(0, 1) = 1; ad(1, 0) = 1;
	rho0(0, 0) = 1; rho1(1, 1) = 1;
	sigmax = a + ad;
	sigmaz = rho0 - rho1;
	
	H0 = sigmaz;
	cOps = {a};
	tOps = {rho0, rho1};
	Ht = {sigmax};
	
	dt = 0.01; tmax = 15;
	coeff = {0.1};
	params = {{10.0, tmax}};
	vector<float(*)(float, vector<float>)> f = {pulse1};

	vector<float> tvec;
	for(int i = 0; i < tmax/dt; i++) {
		tvec.push_back(dt*i);
	}

	fmat data;
	data = te.qtEvolveState(ket0, H0, cOps, tOps, Ht, params, f, tvec, coeff, 1);
	// data = te.meEvolveState(rho0, H0, cOps, tOps, Ht, params, f, tvec, coeff, 1);

	ofstream fout;
	fout.open("data/qtbasicTest.dat");
	// fout.open("data/basicTest.dat");
	fout << data;
	fout.close();

	return 0;
}
