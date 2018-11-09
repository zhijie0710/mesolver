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
		fmat evolveState(cx_mat, cx_mat, vector<cx_mat>, vector<cx_mat>, vector<cx_mat>, vector<vector<float>>, float(*)(float, vector<float>), float, float, vector<float>, bool);

	private:
		// inline void RK4(H, rho, dt, *cOps, *coeff);
		// inline void lindvladME(H, rho, *cOps, *coeff);
		inline void RK4(cx_mat&, cx_mat);
		inline void lindbladME(cx_mat&, cx_mat);
<<<<<<< HEAD
//anychksjfoa
=======
	// David edits
>>>>>>> ea3f27dfede201e5799adeb3784d932a10f70602

};

#endif // MESOLVER


// inline void mesolver::RK4(H, rho, dt, *cOps, *coeff) {
//     .......
// }

// inline void mesolver::lindbladME(H, rho, *cOps, *coeff) {
//     .......
// }

inline void mesolver::RK4(cx_mat& rho, cx_mat H) {
	lindbladME(rho, H); something
} 

inline void mesolver::lindbladME(cx_mat& rho, cx_mat H) {
	rho = 2.0*IM*PI*(rho*H - H*rho);
}
/*
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
*/
