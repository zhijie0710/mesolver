#include <iostream>
#include <vector>
#include <armadillo>

using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

#ifndef MESOLVER
#define MESOLVER

class mesolver {
	public:
		mesolver();
		~mesolver();
		fmat evolveState(cx_mat, cx_mat, vector<cx_mat>, vector<cx_mat>, vector<cx_mat>, vector<vector<float>>, vector<float(*)(float, vector<float>)>, float, float, vector<float>, bool);

	private:
		cx_mat lindbladME(cx_mat, cx_mat, vector<cx_mat>, vector<float>);
		cx_mat RK4(float, cx_mat, cx_mat, vector<cx_mat>, vector<float>);
		// inline void RK4(H, rho, dt, *cOps, *coeff);
		// inline void lindvladME(H, rho, *cOps, *coeff);
		// inline void RK4(cx_mat&, cx_mat);
		// inline void lindbladME(cx_mat&, cx_mat);

};

#endif // MESOLVER

// inline void mesolver::RK4(H, rho, dt, *cOps, *coeff) {
//     .......
// }

// inline void mesolver::lindbladME(H, rho, *cOps, *coeff) {
//     .......
// }

// inline void mesolver::RK4(cx_mat& rho, cx_mat H) {
// 	lindbladME(rho, H);
// } 

// inline void mesolver::lindbladME(cx_mat& rho, cx_mat H) {
// 	rho = 2.0*IM*PI*(rho*H - H*rho);
// }