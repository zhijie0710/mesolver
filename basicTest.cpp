#include <fstream>
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
	ofstream fout;

	fout.open("output.dat");

	cx_mat I, r00, r01, r10, r11, r20, r21, target, currentState,
		   II, III, s0, s1, ss0, ss1, ss2, a, aa, a1, a2,
		   H0, HX;
	vector<cx_mat> cOps, tOps, Ht;

	II = eye<cx_mat>(2,2); III = eye<cx_mat>(3,3);
	s0 = zeros<cx_mat>(2,2); s0(0,0) = 1;
	s1 = zeros<cx_mat>(2,2); s1(1,1) = 1;
	ss0 = zeros<cx_mat>(3,3); ss0(0,0) = 1;
	ss1 = zeros<cx_mat>(3,3); ss1(1,1) = 1;
	ss2 = zeros<cx_mat>(3,3); ss2(2,2) = 1;
	a = zeros<cx_mat>(2,2); a(0,1) = 1;
	aa = zeros<cx_mat>(3,3); aa(0,1) = 1; aa(1,2) = sqrt(2);

	a1 = kron(aa, II); a2 = kron(III, a);
	r00 = kron(ss0, s0); r01 = kron(ss0, s1); r10 = kron(ss1, s0);
	r11 = kron(ss1, s1); r20 = kron(ss2, s0); r21 = kron(ss2, s1);

	H0 = -0.1*a1.t()*a1.t()*a1*a1;
	HX = a1.t()*a2.t() + a1*a2;

	cOps = {a1, a2};
	tOps = {r00, r11};

	float dt = 0.001;
	vector<float> tvec;
	for(int i = 0; i <= 10000; i++) {
		tvec.push_back(dt*i);
	}

	fmat data;
	data = me.evolveState(r00, H0, cOps, tOps, Ht, params, f, dt, tvec, coeff, isSparse);

	/*

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
	// for(int i = 0; i < tvec.size(); i++)
	// 	cout << tvec[i] << " ";
	// cout << endl;

	data = me.evolveState(rho0, H0, cOps, tOps, Ht, params, f, dt, tvec, coeff, isSparse);
	// cout << data << endl;

	fout << data;
*/
	fout.close();

	return 0;
}