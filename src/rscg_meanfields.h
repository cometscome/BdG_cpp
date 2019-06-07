using namespace std;
#include "rscg.cpp"

class RSCG_Meanfields{
    private:
    RSCG* _rscg;
    int _n;
    vector <complex <double> > _vec_iOmega_n; //Matsubara frequencies
    int _omegamax;//Maximum Matsubara frequency
    int _nmax;
    const double pi = std::acos(-1);
    
    public:
    double T;//Temperature

    RSCG_Meanfields(int n,double temperature, double omegamax,std::function<vector <complex <double> >(vector <complex <double> >&)> f);
    RSCG_Meanfields(int n,double temperature, double omegamax, std::function<vector <double>(vector <double>&)> f);
    complex <double> calc_meanfield(int i,int j);

    void set_T(double temperature){
        T = temperature;
    };

    double get_T(){
        return T;
    };

    void set_f(std::function<vector <complex <double> >(vector <complex <double> >&)> f){
        _rscg->set_matvec(f);
    };

    void set_f(std::function<vector <double>(vector <double>&)> f){
        _rscg->set_matvec(f);
    };


    complex <double> get_iOmegan(int n){
//        return _vec_iOmega_n[n];
        return _vec_iOmega_n[n+_nmax];
    } ;

    void set_iOmegan(int n,complex <double> iOmega){
//        _vec_iOmega_n[n] = iOmega;
        _vec_iOmega_n[n+_nmax] = iOmega;
    };

    void initrscg(int n,double temperature,double omegamax){
        _n = n;
        T = temperature;
        _omegamax = omegamax;
    };

    void make_Matsubara(){
        _nmax = (int) (_omegamax/(pi*T)-1)/2;//pi T*(2n+1) = omega -> (omega/(pi*T) -1)/2
//        vector <complex <double> > vec_iOmega_n(_nmax,0.0);
        vector <complex <double> > vec_iOmega_n(2*_nmax,0.0);
        _vec_iOmega_n = vec_iOmega_n;

        for(int n = -_nmax; n < _nmax; n++){
//        for(int n = 0; n < _nmax; n++){
            complex <double> omegan(0.0,pi*T*(2*n+1));
            set_iOmegan(n,omegan);
        };
        
    };


    



};