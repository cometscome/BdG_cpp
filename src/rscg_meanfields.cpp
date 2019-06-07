#include "rscg_meanfields.h"


RSCG_Meanfields::RSCG_Meanfields(int n,double temperature, double omegamax,std::function<vector <complex <double> >(vector <complex <double> >&)> f){
        initrscg(n,temperature,omegamax);
        _rscg = new RSCG(n,f);
        make_Matsubara();
    };

RSCG_Meanfields::RSCG_Meanfields(int n,double temperature, double omegamax,std::function<vector <double>(vector <double>&)> f){
        initrscg(n,temperature,omegamax);
        _rscg = new RSCG(n,f);
        make_Matsubara();
    };

complex <double> RSCG_Meanfields::calc_meanfield(int i,int j){
    vector <complex <double> > Gij = _rscg->calc_greenfunction_complexH(i,j,_vec_iOmega_n);
    //cout << Gij[0] << std::endl; 
    //exit;

    complex <double> csum = 0.0;


//    for(int n = 0; n < _nmax; n++){
    //for(int n = -1000; n < 1000; n++){
    for(int n = -_nmax; n < _nmax; n++){        
        //complex <double> iOmegan = get_iOmegan(n);
        //cout << iOmegan  << std::endl; 
        //complex <double> omegan(0.0,pi*T*(2*n+1));
        //csum = csum +exp(omegan*0.001)/(omegan-0.001);
        //cout << n << " " << csum*T  << std::endl; 
        //exp(iOmegan*0.00000001)*1.0/iOmegan;
        csum = csum + Gij[n+_nmax];//-1.0/iOmegan;
        //csum = csum + Gij[n];//-1.0/iOmegan;
    };
    //cout << pi << std::endl; 
    //cout << T << std::endl; 
    
    csum = csum*T;//+0.5;
    //cout << csum << std::endl; 
    //abort();

    return csum;
};

