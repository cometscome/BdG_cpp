#include "hamiltonian.cpp"
#include <functional>
using namespace std;



int main() {
    Hamiltonian* ham;
    int Nx = 6;
    int Ny = 6;
    double T = 0.05;
    double omegamax = 60*3.1415926535897932;
    double U = -2.0;
    double mu = -0.2;
    cout << "Make Hamiltonian "<< std::endl;
    ham = new Hamiltonian(Nx,Ny,U,T,omegamax,mu);
    cout << "done "<< std::endl;
    vector<complex <double> > x(Nx*Ny*2,1);
    auto Ax = ham->matvec(x);

    int itemax=100;
    double hi = 0.0;
    double eps = 1e-10;
    

    for(int i = 0; i < itemax; i++){
        ham->calc_meanfields();        
        hi = ham->residual();
        cout << i << "-th " << "residual:  " << hi << std::endl;
        cout << "Delta at the center " << ham->get_delta(Nx/2,Ny/2) << std::endl;
        if (hi < eps && i > 0){
            return 0;
        };
        
        ham->update_hamiltonian();  
    };

    

    return 0;
};