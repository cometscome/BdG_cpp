#include "HBdG.h"

int HBdG::xy2i(int ix,int iy){
    int i = _Ny*iy + ix;
    return i;
};

vector <complex <double> > HBdG::matvec(vector <complex <double> > &x){
    
    int* row_ptr;
    int* colidx;
    complex <double>* values;
    char transa = 'N';  //do not transpose matrices in mkl routines
    vector<complex <double> > Ax(NS,0.0);
//    Ax.resize(NS);


    //H.complete_index1_data();
    row_ptr = &H.index1_data()[0];
    colidx = &H.index2_data()[0];
    values = &H.value_data()[0];

    mkl_zcsrgemv(&transa, &NS, &values[0], row_ptr, colidx, &x[0], &Ax[0]);
    return Ax;
//    mkl_zcsrgemv (&transa, &NS, &values[0], row_ptr, colidx, &h_n_i_2[0], &h_n_i_1[0]);
};

HBdG::HBdG(int Nx,int Ny,double mu):H(2*(Nx*Ny),2*(Nx*Ny)){
    
    int N = Nx*Ny;
    _Nx = Nx;
    _Ny = Ny;
    NS = 2*N;

    double delta0 = 1.0;
    complex <double> v;
    double t = 1.0;
    for(int ix = 0; ix < Nx; ix++){
        for(int iy = 0; iy < Ny; iy++){
            //Normal states
            //cout << ix << " " << iy << std::endl;

            
            //+x
            int jx = ix + 1;
            int jy = iy;
            if (jx > Nx-1){
                jx = jx - Nx;
            };
            int i = xy2i(ix,iy);
            int j = xy2i(jx,jy);

            v = -t;
            H(i,j) = v; //electron space
            H(i+N,j+N) = -v; //hole space
            //-x
            jx = ix - 1;
            jy = iy;
            if (jx < 0){
                jx = jx + Nx;
            };
            i = xy2i(ix,iy);
            j = xy2i(jx,jy);
            v = -t;

            H(i,j) = v; //electron space
            H(i+N,j+N) = -v; //hole space

            //+y
            jx = ix ;
            jy = iy+1;
            if (jy > Ny-1){
                jy = jy - Ny;
            };
            i = xy2i(ix,iy);
            j = xy2i(jx,jy);
            v = -t;

            H(i,j) = v; //electron space
            H(i+N,j+N) = -v; //hole space
            //-y
            jx = ix;
            jy = iy-1;
            if (jy < 0){
                jy = jy + Ny;
            };
            i = xy2i(ix,iy);
            j = xy2i(jx,jy);

            v = -t;
            H(i,j) = v; //electron space
            H(i+N,j+N) = -v; //hole space

            //center
            jx = ix;
            jy = iy;
            i = xy2i(ix,iy);
            j = xy2i(jx,jy);

            v = -mu;
            H(i,j) = v; //electron space
            H(i+N,j+N) = -v; //hole space


            //Delta
            jx = ix;
            jy = iy;
            i = xy2i(ix,iy);
            j = xy2i(jx,jy)+N;
            v = delta0;

            H(i,j) = v;
            H(j,i) = v;
        };
    };


};