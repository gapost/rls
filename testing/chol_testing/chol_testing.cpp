#include <iostream>
#include <fstream>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

typedef Matrix <float , Dynamic, Dynamic > Mat;
typedef Matrix< float, Dynamic, 1 > Vec;

void nonrecurciveUpdate( Mat Phi, Vec Y)
{
    LLT<Eigen::Matrix2f> llt;

    Mat A=(Phi.adjoint())*Phi;
    llt.compute(A);             
    Vec Theta = llt.solve(Phi.adjoint()*Y);

    ofstream ThetaFile("build/Theta.txt",  ios::out | ios::app);
    ThetaFile<<"\t"<<Theta(0,0)<<"\t"<<Theta(1,0)<<"\t";
    ThetaFile.close();

   }

void recurciveUpdate(int i, LLT<Eigen::Matrix2f> llt, Mat Phi, Vec Y)
{
    Vec v(2,1);
    v(0,0)=1;    
    v(1,0)=i+1;
    
    llt.rankUpdate( v, 1);                                      //update directly L matrix instread of re-computing L through A
    Vec Theta = llt.solve(Phi.adjoint()*Y);
    
    ofstream ThetaFile("build/Theta.txt",  ios::out | ios::app);
    ThetaFile<<"\t"<<Theta(0,0)<<"\t"<<Theta(1,0)<<"\n";
    ThetaFile.close();

}

int main(){
    
    system("rm build/Theta.txt");
    system("rm build/Signal.txt");

    LLT<Eigen::Matrix2f> llt;
    int m=50;
    Vec Y(m);
    Vec Theta;
    
    //create and fill Phi matrix
    Mat Phi(m,2);               // mx2
    ofstream SignalFile("build/Signal.txt",  ios::out | ios::app);
    for (int i=0; i<m; i++){
        Phi(i,0)=1;
        Phi(i,1)=i+1;
        Y(i)= 5 + 0.02*rand()/RAND_MAX-1.0;
    }
    SignalFile<<Y<<endl;
    SignalFile.close();

    //compute Cholesky matrix -> for the first n elements
    Mat A=((Phi.topLeftCorner(1, 2)).adjoint())*Phi.topLeftCorner(1, 2);  // 4x4
    llt.compute(A);             // 4x4

    //Update Phi and Cholesky matrix
    for(int j=0; j<m; j++){
        nonrecurciveUpdate(Phi.topLeftCorner(j, 2), Y.topLeftCorner(j, 1)); 
        recurciveUpdate(j, llt, Phi.topLeftCorner(j, 2), Y.topLeftCorner(j, 1));
    }

    return 0;
}