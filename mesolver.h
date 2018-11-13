#include <iostream>
#include <vector>
#include <armadillo>
#include <Eigen>

using namespace std;
using namespace arma;
using namespace Eigen;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

#ifndef MESOLVER
#define MESOLVER

class mesolver {
	public:
		mesolver();
		~mesolver();
		fmat evolveState(cx_mat, cx_mat, vector<cx_mat>, vector<cx_mat>, vector<cx_mat>, vector<vector<float>>, vector<float(*)(float, vector<float>)>, float, vector<float>, vector<float>, bool);
		// MatrixXf evolveState(MatrixXcd, MatrixXcd, vector<MatrixXcd>, vector<MatrixXcd>, vector<MatrixXcd>, vector<vector<float>>, vector<float(*)(float, vector<float>)>, float, float, vector<float>, bool);

	private:
		cx_mat lindbladME(cx_mat, cx_mat, vector<cx_mat>, vector<float>);
		cx_mat RK4(float, cx_mat, cx_mat, vector<cx_mat>, vector<float>);
		// MatrixXcd lindbladME(MatrixXcd, MatrixXcd, vector<MatrixXcd>, vector<float>);
		// MatrixXcd RK4(float, MatrixXcd, MatrixXcd, vector<MatrixXcd>, vector<float>);

};

#endif // MESOLVER