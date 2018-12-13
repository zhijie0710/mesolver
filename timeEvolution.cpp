#include "timeEvolution.h"

timeEvolution::timeEvolution() {

}

timeEvolution::~timeEvolution() {

}

void timeEvolution::makePListDense(vector<cx_mat>& cOps, cx_vec& phi, vector<float>& plist)
{
	int N = cOps.size();
	for (int i = 0; i < N; i++)
		plist[i] = abs(as_scalar(phi.t()*cOps[i].t()*cOps[i]*phi));//absolute value?

    for (int i = 0; i < N-1; i++)
    	plist[i+1] = plist[i] + plist[i+1];

    for (int i = 0; i < N; i++)
        plist[i] = plist[i]/plist[N-1];
}

void timeEvolution::jumpDense(vector<cx_mat>& cOps, cx_vec& phi)
{
    vector<float> plist(cOps.size());
	makePListDense(cOps, phi, plist);

	random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<float> dist(0.0, 1.0);
    float s = dist(mt);

    int i = 0;
    while(plist[i] < s) {
    	i++;
    }
    phi = cOps[i]*phi;
}

void timeEvolution::makePListSparse(vector<sp_cx_mat>& cOps, cx_vec& phi, vector<float>& plist)
{
	int N = cOps.size();
	for (int i = 0; i < N; i++)
		plist[i] = abs(as_scalar(phi.t()*cOps[i].t()*cOps[i]*phi));

    for (int i = 0; i < N-1; i++)
    	plist[i+1] = plist[i] + plist[i+1];

    for (int i = 0; i < N; i++)
        plist[i] = plist[i]/plist[N-1];
}

void timeEvolution::jumpSparse(vector<sp_cx_mat>& cOps, cx_vec& phi)
{
    vector<float> plist(cOps.size());
	makePListSparse(cOps, phi, plist);

	random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<float> dist(0.0, 1.0);
    float s = dist(mt);

    int i = 0;
    while(plist[i] < s) {
    	i++;
    }
    phi = cOps[i]*phi;
}

fmat timeEvolution::qtEvolveState(cx_vec psi0, cx_mat H0, vector<cx_mat> cOps, vector<cx_mat> tOps, vector<cx_mat> Ht, vector<vector<float>> params, vector<float(*)(float, vector<float>)> f, vector<float> tvec, vector<float> coeff, bool isDense) {
	int dataLength, dataIndex;
	float dt, r, p;
	dt = tvec[1] - tvec[0];
    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<float> dist(1e-8, 1.0 - 1e-8);
	r = dist(mt);
	fmat dataList = zeros<fmat>(tOps.size() + 1, tvec.size());
	dataIndex = 0;
	cx_vec currentState = psi0;
	if(isDense) {
		cx_mat H;
		for(int i = 0; i < tvec.size(); i++) {
			H = H0;
			for(int j = 0; j < Ht.size(); j++)
				H += f[j](tvec[i], params[j])*Ht[j];
			dataList(0, i) = tvec[i];
			for(int j = 0; j < tOps.size(); j++)
				dataList(j + 1, i) = as_scalar(abs((currentState.t()*tOps[j]*currentState)));
			qtRK4Dense(dt, H, cOps, currentState);
			p = abs(dot(currentState, currentState));
			if(p < r) {
				jumpDense(cOps, currentState);
				currentState = normalise(currentState);
				r = dist(mt);
			}
		}
	}
	else {
		sp_cx_mat H, H0s(H0);
		vector<sp_cx_mat> Ht_sparse(Ht.size()), tOps_sparse(tOps.size()), cOps_sparse(cOps.size());
		for(int i = 0; i < Ht.size(); i++)
			Ht_sparse[i] = sp_cx_mat(Ht[i]);
		for(int i = 0; i < tOps.size(); i++)
			tOps_sparse[i] = sp_cx_mat(tOps[i]);
		for(int i = 0; i < cOps.size(); i++)
			cOps_sparse[i] = sp_cx_mat(cOps[i]);
		for(int i = 0; i < tvec.size(); i++) {
			H = H0s;
			for(int j = 0; j < Ht.size(); j++)
				H += f[j](tvec[i], params[j])*Ht_sparse[j];
			dataList(0, i) = tvec[i];
			for(int j = 0; j < tOps.size(); j++)
				dataList(j + 1, i) = as_scalar(abs(currentState.t()*tOps_sparse[j]*currentState));
			qtRK4Sparse(dt, H, cOps_sparse, currentState);
			p = abs(dot(currentState, currentState));
			if(p < r) {
				jumpSparse(cOps_sparse, currentState);
				currentState = normalise(currentState);
				r = dist(mt);
			}
		}
	}
	return dataList;
}

fmat timeEvolution::meEvolveState(cx_mat rho0, cx_mat H0, vector<cx_mat> cOps, vector<cx_mat> tOps, vector<cx_mat> Ht, vector<vector<float>> params, vector<float(*)(float, vector<float>)> f, vector<float> tvec, vector<float> coeff, bool isDense) {
	int dataLength, dataIndex;
	float dt = tvec[1] - tvec[0];
	fmat dataList = zeros<fmat>(tOps.size() + 1, tvec.size());
	dataIndex = 0;
	if(isDense) {
		cx_mat currentState, H;
		currentState = rho0;
		for(int i = 0; i < tvec.size(); i++) {
			H = H0;
			for(int j = 0; j < Ht.size(); j++)
				H += f[j](tvec[i], params[j])*Ht[j];
			dataList(0, i) = tvec[i];
			for(int j = 0; j < tOps.size(); j++)
				dataList(j + 1, i) = trace(real(tOps[j]*currentState.t()));
			meRK4Dense(dt, H, currentState, cOps, coeff);
		}
	}
	else {
		sp_cx_mat H, H0s(H0), currentState(rho0);
		vector<sp_cx_mat> Ht_sparse(Ht.size()), tOps_sparse(tOps.size()), cOps_sparse(cOps.size());
		for(int i = 0; i < Ht.size(); i++)
			Ht_sparse[i] = sp_cx_mat(Ht[i]);
		for(int i = 0; i < tOps.size(); i++)
			tOps_sparse[i] = sp_cx_mat(tOps[i]);
		for(int i = 0; i < cOps.size(); i++)
			cOps_sparse[i] = sp_cx_mat(cOps[i]);
		for(int i = 0; i < tvec.size(); i++) {
			H = H0s;
			for(int j = 0; j < Ht.size(); j++)
				H += f[j](tvec[i], params[j])*Ht_sparse[j];
			dataList(0, i) = tvec[i];
			for(int j = 0; j < tOps.size(); j++)
				dataList(j + 1, i) = trace(real(tOps_sparse[j]*currentState.t()));
			meRK4Sparse(dt, H, currentState, cOps_sparse, coeff);
		}
	}
	return dataList;
}
