class HBdG{
    private:
    int NS;
    int _Nx;
    int _Ny;
    
    public:
    boost::numeric::ublas::compressed_matrix< complex <double>, 
        boost::numeric::ublas::row_major, 1, boost::numeric::ublas::unbounded_array<int> > H;


    void matvecinit(){
        H.complete_index1_data();
    };

    int xy2i(int ix,int iy);
    void i2xy(int i,int &ix,int &iy);


    HBdG(int Nx,int Ny,double mu);
    vector <complex <double> > matvec(vector <complex <double> > &x);
};