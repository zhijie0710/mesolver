#include <fstream>
#include "time.h"
#include "mesolver.h"

using namespace std;
using namespace arma;

void printRunTime(time_t t0, time_t t1) {
    int t_run, d, h, m, s; 
    t_run = difftime(t1, t0);
    d = floor(t_run/(24*3600));
    h = floor(t_run/3600);
    m = floor(t_run/60); m = m%60;
    s = t_run%60;
    cout << "TIME TO RUN:\n"
         << d << "d " 
         << h << "h " 
         << m << "m " 
         << s << "s " << endl;
}

cx_mat tensor(vector<cx_mat> matrix_list) {
    cx_mat output = kron(matrix_list[0], matrix_list[1]);
    for(int i = 2; i < matrix_list.size(); i++)
		output = kron(output, matrix_list[i]);
    return output;
}

fmat integrate(int N, vector<float> h, vector<float> Jx, vector<float> Jy, vector<float> Jz, cx_mat psi0, vector<float> tlist, vector<float> gamma) {

	cx_mat si(2,2), sx(2,2), sy(2,2), sz(2,2);
	si = eye<cx_mat>(2,2);
	sx(0,1) = 1; sx(1,0) = 1;
	sy(0,1) = -IM; sy(1,0) = IM;
	sz(0,0) = 1; sz(1,1) = -1;

	vector<cx_mat> sx_list, sy_list, sz_list;

	for(int n = 0; n < N; n++) {
		vector<cx_mat> op_list;
		for(int m = 0; m < N; m++)
			op_list.push_back(si);
		
		op_list[n] = sx;
		sx_list.push_back(tensor(op_list));

		op_list[n] = sy;
		sy_list.push_back(tensor(op_list));

		op_list[n] = sz;
		sz_list.push_back(tensor(op_list));
	}

	cx_mat H = 0*sx_list[0];

	for(int n = 0; n < N; n++)
		H += -0.5*h[n]*sz_list[n];

	for(int n = 0; n < N-1; n++) {
		H += -0.5*Jx[n]*sx_list[n]*sx_list[n + 1];
		H += -0.5*Jy[n]*sy_list[n]*sy_list[n + 1];
		H += -0.5*Jz[n]*sz_list[n]*sz_list[n + 1];
	}

	vector<cx_mat> c_op_list;

	for(int n = 0; n < N; n++) {
		if(gamma[n] > 0.0)
			c_op_list.push_back(sz_list[n]);
	}

	mesolver me;
	fmat result;
	result = me.evolveState(psi0, H, c_op_list, sz_list, {}, {}, {}, tlist, gamma, 1);
	return result;
}

int main() {

	fstream fout;
	time_t t0, t1;

	fout.open("output.dat");

	int N = 10;

	vector<float> h(N), Jx(N), Jy(N), Jz(N), gamma(N);
	for(int i = 0; i < N; i++) {
		h[i] = 1.0*2*PI;
		Jx[i] = 0.1*2*PI;
		Jy[i] = 0.1*2*PI;
		Jz[i] = 0.1*2*PI;
		gamma[i] = 0.01;
	}

	cx_mat s0(2,2), s1(2,2);
	s0(0,0) = 1; s1(1,1) = 1;

	vector<cx_mat> rho_list;
	rho_list.push_back(s1);
	for(int n = 0; n < N-1; n++)
		rho_list.push_back(s0);
	cx_mat rho0 = tensor(rho_list);

	vector<float> tlist;
	float dt = 50.0/200;
	for(int i = 0; i <= 200; i++)
		tlist.push_back(dt*i);

	fmat sz_expt;
	time(&t0);
	sz_expt = integrate(N, h, Jx, Jy, Jz, rho0, tlist, gamma);
	time(&t1);

	fout << sz_expt;
	cout << sz_expt;

	printRunTime(t0, t1);

	fout.close();

	return 0;
}