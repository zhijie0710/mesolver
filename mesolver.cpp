#include "mesolver.h"

mesolver::mesolver() {

}

mesolver::~mesolver() {

}

fmat mesolver::evolveState(cx_mat rho0, cx_mat H0, vector<cx_mat> cOps, vector<cx_mat> tOps, vector<cx_mat> Ht, vector<vector<float>> params, float(*f)(float, vector<float>), float dt, float tmax, vector<float> coeff, bool isSparse) {
	cx_mat currentState, H;
	int dataLength, dataIndex;
	currentState = rho0;
	H = H0;
	float t = 0;
	dataLength = tmax/dt;
	fmat dataList(tOps.size() + 1, dataLength);
	dataIndex = 0;
	while(t <= tmax && dataIndex < dataLength) {
		for(int i = 0; i < Ht.size(); i++) {
			H += f(t, params[i])*Ht[i];
		}
		dataList(0, dataIndex) = t;
		for(int i = 0; i < tOps.size(); i++) {
			dataList(i + 1, dataIndex) = trace(abs(tOps[i]*currentState.t()));
		}
		RK4(currentState, H);
		t += dt; dataIndex++;
	}
	return dataList;
}
