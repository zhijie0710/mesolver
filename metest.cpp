#include<stdio.h>
#include <iostream>
#include <vector>
#include <sstream>
#include <iomanip>
#include <armadillo>
#include <complex.h>
#include <typeinfo>
#include <random>

using namespace std;
using namespace arma;
const   complex<double> i(0.0,1.0); 



cx_vec f(sp_cx_mat H,  vector<sp_cx_mat> cOps, cx_vec phi)  //phi is the current state |phi(t)>.
{
	sp_cx_mat Heff;
	Heff = H;
	for (int i = 0; i < cOps.size(); i++)
	{
		Heff += - 0.5 * i * (cOps[i].t()*cOps[i]);
	}
	return(-i*Heff*phi);
}

cx_vec RK4(float hs, sp_cx_mat H, vector<sp_cx_mat> cOps, cx_vec phi)
{
	//cx_vec f(sp_cx_mat H,  vector<sp_cx_mat> cOps, sp_cx_mat phi);
	cx_vec k1 = f( H, cOps, phi) * hs;
	cx_vec phi1 = phi + k1 * 0.5;
	cx_vec k2 = f( H, cOps, phi) * hs;
	cx_vec phi2 = phi + k2 * 0.5;
	cx_vec k3 = f( H, cOps, phi) * hs;
	cx_vec phi3 = phi + k3;
	cx_vec k4 = f( H, cOps, phi) * hs;
	cx_vec xx = phi + (k1 + 2*(k2 + k3) + k4)/6;

	return(xx); 
}

//pick a r ~ [0,1], if |<phi|phi>| < r, apply quantum jump, then normalize |phi>, and pick a new r

//generate a list of {pi}
vector<float> makePList(vector<sp_cx_mat> cOps, cx_vec phi)
{
	int N = cOps.size();
	vector<float> coeffp(N);
	for (int i=0; i < N ; i++)
	{
		coeffp[i]=phi.t() * cOps[i].t() * cOps[i] * phi;
	}

    for (int i=0; i < N-1 ; i++)
    {
    	coeffp[i+1]=coeffp[i] + coeffp[i+1];
    }
    for (int i=0; i < N ; i++)
    {
        coeffp[i] = coeffp[i]/coeffp[N-1];
    }
return (coeffp);
}




//define the quantum jump function
cx_vec Jump(vector<float> coeff,  vector<sp_cx_mat> cOps, cx_vec phi)
{
    vector<float> coeffp (vector<sp_cx_mat> cOps, sp_cx_mat phi);
    coeff  =  coeffp(cOps, phi);


	random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    double s = dist(mt);

    int i = 0;
    while(coeff[i] < s) {
    	i++;
    }
    phi = cOps[i]*phi;
    return(phi);
/*
for (int i = 0; i < coeff.size(); i++)
{
	if(coeff[i]>s)
	{
		phi = cOps[i] * phi;
		return(phi);
	}
}*/

}



//main 
int main()

{
	float tf=10;
	float hs = 0.01;
    int a = tf/hs;

    random_device rd;
    mt19937 mt(rd());
    uniform_real_distribution<double> dist(0.0, 1.0);
    double r = dist(mt);

    for (int i=0; i < a; i++)
    {
        phi = RK4( hs,  H, cOps, phi);
        double p=abs(phi.t() * phi);
        if(p<r)
        {
    	    phi = normalise(Jump(phi));
    	    r = dist(mt);
        }
    }
return(0);
}












