#include <vector>
#include <cmath>
#include <memory>
#include <new>
#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <complex>
#include <time.h>
//#include <mpi.h>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>

#define MKL_Complex16 std::complex<double>
#include <mkl.h>
#include <mkl_spblas.h>
#include <mkl_lapacke.h>


#include "rscg_meanfields.cpp"
#include "HBdG.cpp"

using namespace std;


class Hamiltonian{
    private:
    int _Nx;
    int _Ny;
    

    RSCG_Meanfields* _rscg_meanfields;
    HBdG* _hbdg;
    
    public:
    int N;
    double U;
    vector <complex <double> > vec_delta;
    vector <complex <double> > vec_delta_old;

    
    complex <double> get_delta(int ix,int iy){
        int i = _hbdg->xy2i(ix,iy);
        return vec_delta[i];
    };

    complex <double> get_value(int i,int j){
        return _hbdg->H(i,j);
    };

    void set_value(int i, int j,complex <double> value){
        _hbdg->H(i,j) = value;
    };

    
    Hamiltonian(int Nx,int Ny,double Ui,double T,double omegamax,double mu);

    void calc_meanfields();
    void update_hamiltonian();
    double residual();

    void set_f(std::function<vector <complex <double> >(vector <complex <double> >&)> f){
        _rscg_meanfields->set_f(f);
    };

    void set_f(std::function<vector <double>(vector <double>&)> f){
        _rscg_meanfields->set_f(f);
    };

    vector <complex <double> > matvec(vector <complex <double> > &x){
        cout << "in matvec "<< std::endl;
        auto Ax = _hbdg->matvec(x);
        return Ax;
    };
};