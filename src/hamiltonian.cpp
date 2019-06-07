#include "hamiltonian.h"



void Hamiltonian::calc_meanfields(){
    
    complex <double> v;
    for(int ix = 0; ix < _Nx; ix++){
        for(int iy = 0; iy < _Ny; iy++){
            int jx = ix;
            int jy = iy;
            int i = _hbdg->xy2i(ix,iy);
            int j = _hbdg->xy2i(jx,jy)+N;
            v = _rscg_meanfields->calc_meanfield(i,j);
            //cout << "ix,iy " << ix << "," << iy << " " << v << std::endl;
            //cout << U << std::endl;
            vec_delta[i] = U*v;
            //cout << "ix,iy " << ix << "," << iy << " " << vec_delta[i]  << std::endl;
            //
        };
        //cout << vec_delta[ix] << std::endl;
    };
    //cout << vec_delta[0] << std::endl;
};

void Hamiltonian::update_hamiltonian(){
    for(int i = 0; i < N; i++){
        complex <double>  v = vec_delta[i];
        int j = i + N;
        _hbdg->H(i,j) = v;
        _hbdg->H(j,i) = conj(v);
        vec_delta_old[i] = vec_delta[i];
    };

    _hbdg ->matvecinit();
    function<vector <complex <double> >(vector <complex <double> >&)>  f = std::bind(&HBdG::matvec,_hbdg,std::placeholders::_1);
    set_f(f);

};

double Hamiltonian::residual(){
    double res = 0.0;
    double norm = 0.0;
    for(int i = 0; i < N; i++){
        res = res + abs(vec_delta[i]-vec_delta_old[i])*abs(vec_delta[i]-vec_delta_old[i]);
        norm = norm + abs(vec_delta_old[i])*abs(vec_delta_old[i]);
    };
    res = res/norm;
    return res;
};

Hamiltonian::Hamiltonian(int Nx,int Ny,double Ui,double T,double omegamax,double mu){
    _Nx = Nx;
    _Ny = Ny;
    N = Nx*Ny;
    U = Ui;


    _hbdg = new HBdG(Nx,Ny,mu);
    _hbdg ->matvecinit();
    auto f = std::bind(&HBdG::matvec,_hbdg,std::placeholders::_1);

    _rscg_meanfields= new RSCG_Meanfields(2*Nx*Ny,T,omegamax,f);


    vector <complex <double> > vec_delta0(N,0.0);
    for(int i = 0; i < N; i++){
        vec_delta0[i] = _hbdg->H(i,i+N);
    };
    vec_delta = vec_delta0;
    vec_delta_old = vec_delta0;
};
