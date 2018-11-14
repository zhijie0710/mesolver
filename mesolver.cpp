#include "mesolver.h"

mesolver::mesolver() {

}

mesolver::~mesolver() {

}

fmat mesolver::evolveState(cx_mat rho0, cx_mat H0, vector<cx_mat> cOps, vector<cx_mat> tOps, vector<cx_mat> Ht, vector<vector<float>> params, vector<float(*)(float, vector<float>)> f, float dt, vector<float> tvec, vector<float> coeff, bool isSparse) {
	cx_mat currentState, H;
	int dataLength, dataIndex;
	currentState = rho0;
	H = H0;
	fmat dataList(tOps.size() + 1, tvec.size());
	dataIndex = 0;

	for(int i = 0; i < tvec.size(); i++) {
		for(int j = 0; j < Ht.size(); j++)
			H += f[j](tvec[i], params[j])*Ht[j];
		dataList(0, i) = tvec[i];
		for(int j = 0; j < tOps.size(); j++)
			dataList(j + 1, i) = trace(abs(tOps[j]*currentState.t()));
		currentState = RK4(dt, H, currentState, cOps, coeff);
		H = H0;
	}

	return dataList;
}

cx_mat mesolver::lindbladME(cx_mat H, cx_mat rho, vector<cx_mat> cOps, vector<float> coeff)
{
	cx_mat newrho;
	newrho = IM*(rho * H - H * rho);
	for (int i = 0; i < coeff.size() ; i++)
	{
		newrho += coeff[i]*(cOps[i]*rho*cOps[i].t() - 0.5*(cOps[i].t()*cOps[i]*rho + rho*cOps[i].t()*cOps[i]));
	}
	return(newrho);
}

cx_mat mesolver::RK4(float hs, cx_mat H, cx_mat rho, vector<cx_mat> cOps, vector<float> coeff)
{
	cx_mat k1 = lindbladME(H, rho, cOps, coeff) * hs;
	cx_mat rho1 = rho + k1 * 0.5;
	cx_mat k2 = lindbladME( H, rho1, cOps, coeff) * hs;
	cx_mat rho2 = rho + k2 * 0.5;
	cx_mat k3 = lindbladME( H, rho2, cOps, coeff) * hs;
	cx_mat rho3 = rho + k3;
	cx_mat k4 = lindbladME( H, rho3, cOps, coeff) * hs;
	cx_mat xx = rho + (k1 + 2*(k2 + k3) + k4)/6;

	return(xx); 
}