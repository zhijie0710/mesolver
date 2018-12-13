#include <iostream>
#include <vector>
#include <random>
#include <armadillo>

using namespace std;
using namespace arma;

#define PI 3.14159265358979323846
#define IM std::complex<double> (0.0, 1.0)

#ifndef TIME_EVOLUTION
#define TIME_EVOLUTION

class timeEvolution {
	public:
		timeEvolution();
		~timeEvolution();
		fmat meEvolveState(cx_mat, cx_mat, vector<cx_mat>, vector<cx_mat>, vector<cx_mat>, vector<vector<float>>, vector<float(*)(float, vector<float>)>, vector<float>, vector<float>, bool);
		fmat qtEvolveState(cx_vec, cx_mat, vector<cx_mat>, vector<cx_mat>, vector<cx_mat>, vector<vector<float>>, vector<float(*)(float, vector<float>)>, vector<float>, vector<float>, bool);

	private:
		inline void lindbladMEDense(cx_mat&, cx_mat&, vector<cx_mat>&, vector<float>&, cx_mat&);
		inline void meRK4Dense(float, cx_mat&, cx_mat&, vector<cx_mat>&, vector<float>&);
		inline void lindbladMESparse(sp_cx_mat&, sp_cx_mat&, vector<sp_cx_mat>&, vector<float>&, sp_cx_mat&);
		inline void meRK4Sparse(float, sp_cx_mat&, sp_cx_mat&, vector<sp_cx_mat>&, vector<float>&);
		inline void qtRK4Dense(float hs, cx_mat&, vector<cx_mat>&, cx_vec&);
		inline void HeffDense(cx_mat&,  vector<cx_mat>&, cx_vec&, cx_vec&);
		void makePListDense(vector<cx_mat>&, cx_vec&, vector<float>&);
		void jumpDense(vector<cx_mat>&, cx_vec&);
		inline void qtRK4Sparse(float hs, sp_cx_mat&, vector<sp_cx_mat>&, cx_vec&);
		inline void HeffSparse(sp_cx_mat&,  vector<sp_cx_mat>&, cx_vec&, cx_vec&);
		void makePListSparse(vector<sp_cx_mat>&, cx_vec&, vector<float>&);
		void jumpSparse(vector<sp_cx_mat>&, cx_vec&);
};

#endif // TIME_EVOLUTION

inline void timeEvolution::lindbladMEDense(cx_mat& H, cx_mat& rho, vector<cx_mat>& cOps, vector<float>& coeff, cx_mat& newrho)
{
	newrho = IM*(rho * H - H * rho);
	for (int i = 0; i < coeff.size() ; i++)
	{
		newrho += coeff[i]*(cOps[i]*rho*cOps[i].t() - 0.5*(cOps[i].t()*cOps[i]*rho + rho*cOps[i].t()*cOps[i]));
	}
}

inline void timeEvolution::meRK4Dense(float dt, cx_mat& H, cx_mat& rho, vector<cx_mat>& cOps, vector<float>& coeff)
{
	cx_mat k1, k2, k3, k4, rho1, rho2, rho3;
	lindbladMEDense(H, rho, cOps, coeff, k1);
	rho1 = rho + k1 * dt * 0.5;
	lindbladMEDense(H, rho1, cOps, coeff, k2);
	rho2 = rho + k2 * dt * 0.5;
	lindbladMEDense(H, rho2, cOps, coeff, k3);
	rho3 = rho + dt * k3;
	lindbladMEDense(H, rho3, cOps, coeff, k4);

	rho = rho + dt*(k1 + 2*(k2 + k3) + k4)/6.0;
}

inline void timeEvolution::lindbladMESparse(sp_cx_mat& H, sp_cx_mat& rho, vector<sp_cx_mat>& cOps, vector<float>& coeff, sp_cx_mat& newrho)
{
	newrho = IM*(rho * H - H * rho);
	for (int i = 0; i < coeff.size() ; i++)
	{
		newrho += coeff[i]*(cOps[i]*rho*cOps[i].t() - 0.5*(cOps[i].t()*cOps[i]*rho + rho*cOps[i].t()*cOps[i]));
	}
}

inline void timeEvolution::meRK4Sparse(float dt, sp_cx_mat& H, sp_cx_mat& rho, vector<sp_cx_mat>& cOps, vector<float>& coeff)
{
	sp_cx_mat k1, k2, k3, k4, rho1, rho2, rho3;
	lindbladMESparse(H, rho, cOps, coeff, k1);
	rho1 = rho + k1*dt*0.5;
	lindbladMESparse(H, rho1, cOps, coeff, k2);
	rho2 = rho + k2*dt*0.5;
	lindbladMESparse(H, rho2, cOps, coeff, k3);
	rho3 = rho + dt*k3;
	lindbladMESparse(H, rho3, cOps, coeff, k4);

	rho = rho + dt*(k1 + 2*(k2 + k3) + k4)/6.0;
}

inline void timeEvolution::HeffDense(cx_mat& H, vector<cx_mat>& cOps, cx_vec& phi, cx_vec& newphi)
{
	cx_mat Heff = H;
	for (int i = 0; i < cOps.size(); i++)
	{
		Heff += - 0.5 * i * (cOps[i].t()*cOps[i]);
	}
	newphi = -IM*Heff*phi;
}

inline void timeEvolution::qtRK4Dense(float dt, cx_mat& H, vector<cx_mat>& cOps, cx_vec& phi)
{
	cx_vec k1, k2, k3, k4, phi1, phi2, phi3;
	HeffDense(H, cOps, phi, k1);
	phi1 = phi + k1*dt*0.5;
	HeffDense(H, cOps, phi, k2);
	phi2 = phi + k2*dt*0.5;
	HeffDense(H, cOps, phi, k3);
	phi3 = phi + k3;
	HeffDense(H, cOps, phi, k4);

	phi = phi + (k1 + 2*(k2 + k3) + k4)/6.0;
}

inline void timeEvolution::HeffSparse(sp_cx_mat& H, vector<sp_cx_mat>& cOps, cx_vec& phi, cx_vec& newphi)
{
	sp_cx_mat Heff = H;
	for (int i = 0; i < cOps.size(); i++)
	{
		Heff += - 0.5 * i * (cOps[i].t()*cOps[i]);
	}
	newphi = -IM*Heff*phi;
}

inline void timeEvolution::qtRK4Sparse(float dt, sp_cx_mat& H, vector<sp_cx_mat>& cOps, cx_vec& phi)
{
	cx_vec k1, k2, k3, k4, phi1, phi2, phi3;
	HeffSparse(H, cOps, phi, k1);
	phi1 = phi + k1*dt*0.5;
	HeffSparse(H, cOps, phi, k2);
	phi2 = phi + k2*dt*0.5;
	HeffSparse(H, cOps, phi, k3);
	phi3 = phi + k3;
	HeffSparse(H, cOps, phi, k4);

	phi = phi + (k1 + 2*(k2 + k3) + k4)/6.0;
}