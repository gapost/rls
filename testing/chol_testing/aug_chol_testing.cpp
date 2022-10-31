#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <random>

using namespace std;
using namespace Eigen;

typedef Matrix <float , Dynamic, Dynamic > Mat;
typedef Matrix< float, Dynamic, 1 > Vec;

void Aug_Cholesky_Method(Mat L_Aug,Mat Phi)
{
    Mat Z=(L_Aug.bottomLeftCorner(1,2)).adjoint();
    Mat L=(L_Aug.topLeftCorner(2,2));
    
    Mat L_inv_tr=L.adjoint();
    L_inv_tr=L_inv_tr.inverse();

    Vec Theta = L_inv_tr*Z;
    Mat S=Phi.bottomLeftCorner(1,2)*(Theta);              //Recreate the Signal     

    ofstream ThetaFile("Aug_Theta.txt",  ios::out | ios::app);
    ThetaFile<<Theta(0,0)<<"\t"<<Theta(1,0)<<"\t"<<S<<"\n";
    ThetaFile.close();
}

int main(){
    
    system("rm Aug_Theta.txt");

    int   l=2500, m=1500, n=800;
    float l1=4.8, l2=5.3, l3=3.4;

    
    Vec Y(l);
    //Vec Theta;
    
    // create and fill Phi matrix
    Mat Phi(l,2);               // l x 2
    
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
  
    ofstream SignalFile("Aug_Signal.txt",  ios::out | ios::app);
    SignalFile<<Y<<endl;
    SignalFile.close();

    // Augmented Cholesky
    int start=25;
    const int win=80;

    Mat Phi_Aug(l,3);
    Phi_Aug.setZero();
    Phi_Aug.topRightCorner(start,1)=Y.topLeftCorner(start,1);
    Phi_Aug.topLeftCorner(start,2)=Phi.topLeftCorner(start,2);

    Mat A_Aug=Phi_Aug.adjoint()*Phi_Aug;
    LLT<Mat> llt(A_Aug);
    Mat L_Aug =llt.matrixL();
    Aug_Cholesky_Method(L_Aug,Phi.topLeftCorner(start, 2));

    // Matrix used to update llt
    Vec v(3,1);
    v(0,0)=1;   

    // Update Cholesky matrix and compute Theta again
    for(int j=start; j<l; j++){ 
    
        if(j>win){
            v(1,0)=j;
            v(2,0)=Y(j);
            llt.rankUpdate( v, 1);    //update directly L matrix instread of re-computing L through A
            v(1,0)=j-win;
            v(2,0)=Y(j-win);
            llt.rankUpdate( v, -1);    //downdate directly L matrix instread of re-computing L through A
            L_Aug=llt.matrixL(); 
            
            Aug_Cholesky_Method(L_Aug, Phi.block<win,2>(j-win,0));
        }
        
        else{
            v(2,0)=Y(j);
            v(1,0)=j;
            llt.rankUpdate( v, 1);    //update directly L matrix instread of re-computing L through A
            L_Aug=llt.matrixL(); 
            Aug_Cholesky_Method(L_Aug,Phi.topLeftCorner(j, 2));
        }
   }

    return 0;
}