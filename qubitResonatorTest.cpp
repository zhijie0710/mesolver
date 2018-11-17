#include <fstream>
#include "mesolver.h"

using namespace std;
using namespace arma;

float pulse(float t, vector<float> params) {
	return(params[0]*sin(PI*t/params[1]));
}

int main() {

	mesolver me;
	ofstream fout;
	bool isSparse = 0;

	fout.open("qubitResonator_correct_output.dat");

	cx_mat I, r00, r01, r10, r11, r20, r21, target, CurrentState, II, III, s0, s1, ss0, ss1, ss2, a, aa, a1, a2, H0, HX, HY;
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

	H0 = -0.2*PI*a1.t()*a1.t()*a1*a1;
	HX = a1.t()*a2.t() + a1*a2;
	HY = IM*(a1.t()*a2.t() - a1*a2);

	cOps = {a1, a2};
	tOps = {r00, r11};
	Ht = {HX, HY};

	float dt = 0.001;
	float tmax = 20;
	vector<float> tvec;
	for(int i = 0; i <= 10*tmax/dt; i++)
		tvec.push_back(dt*i);
	vector<vector<float>> params = {{0.04*PI, tmax},{0.0, tmax}};
	vector<float(*)(float, vector<float>)> f = {pulse, pulse};
	vector<float> coeff = {0.0};


	fmat data;
	data = me.evolveState(r00, H0, cOps, tOps, Ht, params, f, dt, tvec, coeff, isSparse);

	fout << data;
	fout.close();

	return 0;
}