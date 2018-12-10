#include "mesolver.h"
#include <cstdlib>



mesolver::mesolver() {

}

mesolver::~mesolver() {

}

fmat mesolver::evolveState(cx_mat rho0, cx_mat H0, vector<cx_mat> cOps, vector<cx_mat> tOps, vector<cx_mat> Ht, vector<vector<float>> params, vector<float(*)(float, vector<float>)> f, vector<float> tvec, vector<float> coeff, bool isSparse) {
	int dataLength, dataIndex;
	float dt = tvec[1] - tvec[0];
	fmat dataList = zeros<fmat>(tOps.size() + 1, tvec.size());
	dataIndex = 0;
	if(isSparse) {
		sp_cx_mat H(H0), currentState(rho0);
		vector<sp_cx_mat> Ht_sparse(Ht.size()), tOps_sparse(tOps.size()), cOps_sparse(cOps.size());
		for(int i = 0; i < Ht.size(); i++)
			Ht_sparse[i] = sp_cx_mat(Ht[i]);
		for(int i = 0; i < tOps.size(); i++)
			tOps_sparse[i] = sp_cx_mat(tOps[i]);
		for(int i = 0; i < cOps.size(); i++)
			cOps_sparse[i] = sp_cx_mat(cOps[i]);
		for(int i = 0; i < tvec.size(); i++) {
			for(int j = 0; j < Ht.size(); j++)
				H += f[j](tvec[i], params[j])*Ht_sparse[j];
			dataList(0, i) = tvec[i];
			for(int j = 0; j < tOps.size(); j++)
				dataList(j + 1, i) = trace(real(tOps_sparse[j]*currentState.t()));
			currentState = RK4Sparse(dt, H, currentState, cOps_sparse, coeff);
			H = H0;
		}
	}
	else {
		cx_mat currentState, H;
		currentState = rho0;
		H = H0;
		for(int i = 0; i < tvec.size(); i++) {
			for(int j = 0; j < Ht.size(); j++)
				H += f[j](tvec[i], params[j])*Ht[j];
			dataList(0, i) = tvec[i];
			for(int j = 0; j < tOps.size(); j++)
				dataList(j + 1, i) = trace(real(tOps[j]*currentState.t()));
			currentState = RK4(dt, H, currentState, cOps, coeff);
			H = H0;
		}
	}
	return dataList;
}

/*
sp_cx_mat mesolver::lindbladMESparse(sp_cx_mat H, sp_cx_mat rho, vector<sp_cx_mat> cOps, vector<float> coeff)
{
	sp_cx_mat newrho;
	newrho = IM*(rho * H - H * rho);
	for (int i = 0; i < coeff.size() ; i++)
	{
		newrho += coeff[i]*(cOps[i]*rho*cOps[i].t() - 0.5*(cOps[i].t()*cOps[i]*rho + rho*cOps[i].t()*cOps[i]));
	}
	return(newrho);
}
*/



sp_cx_mat f(sp_cx_mat H,  vector<sp_cx_mat> cOps, sp_cx_mat phi)  //phi is the current state |phi(t)>.
{
	sp_cx_mat Heff;
	Heff = H;
	for (int i = 0; i < cOps.size(); i++)
	{
		Heff += - 0.5 * i * (cOps[i].t()*cOps[i]);
	}
	return(-i*Heff*phi);
}

//pick a r ~ [0,1], if |<phi|phi>| < r, apply quantum jump, then normalize |phi>, and pick a new r





sp_cx_mat mesolver::RK4Sparse(float hs, sp_cx_mat H, vector<sp_cx_mat> cOps, sp_cx_mat phi)
{
	sp_cx_mat k1 = f( H, cOps, phi) * hs;
	sp_cx_mat phi1 = phi + k1 * 0.5;
	sp_cx_mat k2 = f( H, cOps, phi) * hs;
	sp_cx_mat phi2 = phi + k2 * 0.5;
	sp_cx_mat k3 = f( H, cOps, phi) * hs;
	sp_cx_mat phi3 = phi + k3;
	sp_cx_mat k4 = f( H, cOps, phi) * hs;
	sp_cx_mat xx = phi + (k1 + 2*(k2 + k3) + k4)/6;

	return(xx); 
}