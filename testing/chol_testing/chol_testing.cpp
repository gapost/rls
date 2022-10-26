#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <random>

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

    Mat S=Phi.bottomLeftCorner(1,2)*(Theta);              //Recreate the Signal     
    ofstream ThetaFile("build/Theta.txt",  ios::out | ios::app);
    ThetaFile<<"\t"<<Theta(0,0)<<"\t"<<Theta(1,0)<<"\t"<<S<<"\t";
    ThetaFile.close();
   
   }

void recurciveUpdate(int i, LLT<Eigen::Matrix2f> llt, Mat Phi, Vec Y)
{
    Vec Theta = llt.solve(Phi.adjoint()*Y);

    Mat S=Phi.bottomLeftCorner(1,2)*(Theta);              //Recreate the Signal     
    ofstream ThetaFile("build/Theta.txt",  ios::out | ios::app);
    ThetaFile<<"\t"<<Theta(0,0)<<"\t"<<Theta(1,0)<<"\t"<<S<<"\n";
    ThetaFile.close();

}

int main(){
    
    system("rm build/Theta.txt");
    system("rm build/Signal.txt");

    LLT<Eigen::Matrix2f> llt;

    int   l=2500, m=1500, n=800;
    float l1=4.8, l2=6.3, l3=4;
    cout<<l1<<" "<<l2<<" "<<l3<<endl;
    cout<<l<<" "<<m<<" "<<n<<endl;
    
    Vec Y(l);
    Vec Theta;
    
    // create and fill Phi matrix
    Mat Phi(l,2);               // l x 2
    ofstream SignalFile("build/Signal.txt",  ios::out | ios::app);
    
    random_device rd;  // Will be used to obtain a seed for the random number engine
    mt19937 gen(rd()); // Standard mersenne_twister_engine seeded with rd()
    uniform_real_distribution<> dis(-0.5, 0.5);

    for (int i=0; i<n; i++){
        Phi(i,0)=1;
        Phi(i,1)=i+1;
        Y(i)= l1 + dis(gen);
    }
    for (int i=n; i<m; i++){
        Phi(i,0)=1;
        Phi(i,1)=i+1;
        Y(i)= l2 + dis(gen);

    }
    for (int i=m; i<l; i++){
        Phi(i,0)=1;
        Phi(i,1)=i+1;
        Y(i)= l3 + dis(gen);
    }
  
    SignalFile<<Y<<endl;
    SignalFile.close();
    // Compute Cholesky matrix -> for the first n elements
    int start=25;
    Mat A=((Phi.topLeftCorner(start, 2)).adjoint())*Phi.topLeftCorner(start, 2);  // 4x4
    llt.compute(A);             // 4x4
    
    // Matrix used to update llt
    Vec v(2,1);
    v(0,0)=1;    
    int j=100;

    // Update Cholesky matrix and compute Theta again
    for(int j=start; j<l; j++){
        nonrecurciveUpdate(Phi.topLeftCorner(j, 2), Y.topLeftCorner(j, 1)); 
        v(1,0)=j+1;
        llt.rankUpdate( v, 1);    //update directly L matrix instread of re-computing L through A
        recurciveUpdate(j, llt, Phi.topLeftCorner(j, 2), Y.topLeftCorner(j, 1));
    }


    return 0;
}