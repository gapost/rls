#include <iostream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;
int main(){
    typedef Matrix <float , Dynamic, Dynamic > Mat;
	typedef Matrix< float, Dynamic, 1 > Vec;
    LLT<Eigen::Matrix2f> llt;
    int n = 20;
    Vec Time(n,1);
    Mat Phi(n,2);
    Vec Y(n);
    for (int i=0; i<n; i++){
        Time(i)=i+1;
        Phi(i,0)=1;
        Phi(i,1)=Time(i);
        Y(i)=Time(i)+5+ 2.0*rand()/RAND_MAX-1.0;
    }

    Mat A=(Phi.adjoint())*Phi;
    Mat L = A.llt().matrixL();
    llt.compute(A);
    Vec Theta = llt.solve(Phi.adjoint()*Y);

    cout << "\nThe matrix A is\n"<<A<<endl;
    cout << "\nThe Cholesky factor L is" << endl << L << endl;
    cout << "\nTo check this, let us compute L * L.transpose()" << endl;
    cout << L * L.transpose() << endl;
    cout << "This should equal the matrix A" << endl;
    cout << "\nThe solution matrix theta is:\n"<<Theta<<endl;
    cout<<endl;
    return 0;
}