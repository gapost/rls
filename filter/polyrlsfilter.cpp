#include <iostream>
#include <string>
#include "RLS.h"
#include "PolyRLS.h"

using namespace std;

template<class _RLS>
int do_filter(_RLS& rls)
{
    while(!cin.eof())
    {
        double data;
        cin >> data;
        if (!cin.good()) return -1;
        rls.update(data);
        cout << rls.estimatedOutput();
        cout << '\t' << rls.estimatedRate();
        cout << '\t' << rls.cost();         
        const vector<double> & w = rls.estimatedPar();
        int np = w.size();
        for(int j=0; j<np; j++) cout << '\t' << w[j];
        cout << endl;
    }
    return 0;
}

int main(int argn, char** argv)
{
    int Np = 2; // poly order
    double ff = 0.98; // forgetting factor
    int M = 0; // block width
    bool sqrt_upd = false;

    const char* usage =
        "Usage:\n"
        "  polyrlsfilter [options] \n"
        "\n"
        "Recursive least squares (RLS) fitting of a polynomial \n"
        "a_0 + a_1*t + ... to a time series using a forgetting factor or\n"
        "sliding-block algorithm.\n"
        "\n"
        "Options:\n"
        "  -nX  : number of fitting parameters, X=2,3,... (default = 2) \n"
        "  -ffX : forgetting factor, 0<X<=1, (default = 0.98) \n"
        "  -wX  : block width, X=0 means no block (default)\n"
        "  -sqrt : update square root of covariance\n"
        "  -h   : display this help message\n"
        "\n"
        "polyrlsfilter takes the input signal, y_t, t=1,2,..., from stdin\n"
        "as a stream of ascii-coded real values until eof or a non-numeric input\n"
        "is encountered.\n"
        "The output is written to stdout. Each output line contains the following tab\n"
        "separated data: estimated signal, estimated rate, sum of squared residuals,\n"
        "estimated parameters\n"
        "\n"
        "Use redirection for file input/output, e.g.:\n"
        "  polyrlsfilter < in.dat > out.dat\n";

    // parse options
    for(int i=1; i<argn; i++) {
        string opt(argv[i]);
        if (opt.find("-ff")==0) { // forgetting factor
            opt.erase(0,3);
            ff = atof(opt.c_str());
            if (ff>1. || ff<=0) {
                cout << "Invalid forgetting factor value: " << ff << endl;
                cout << usage;
                return -1;
            }
        }
        else if (opt.find("-n")==0) { // poly order
            opt.erase(0,2);
            Np = atoi(opt.c_str());
            if (Np <= 1) {
                cout << "Invalid # of fitting parameters: " << Np << endl;
                cout << usage;
                return -1;
            }
        }
        else if (opt.find("-w")==0) { // block width
            opt.erase(0,2);
            M = atoi(opt.c_str());
        }
        else if (opt.find("-sqrt")==0) { // update sqrt 
            sqrt_upd = true;
        }
        else if (opt.find("-h")==0) { // block width
            cout << usage;
            return 0;
        }
        else {
            cout << "Invalid option: " << opt << endl;
            cout << usage;
            return -1;
        }
    }

    if (M>0) {
        if (sqrt_upd) {
            RLS::PolyRLS< double, RLS::SquareRootUpdate, 
                RLS::BlockRLS<double, RLS::SquareRootUpdate> > rls(Np);
            rls.setSize(M);
            return do_filter(rls);
        } else {
            RLS::PolyRLS< double, RLS::CovarianceUpdate, 
                RLS::BlockRLS<double, RLS::CovarianceUpdate> > rls(Np);
            rls.setSize(M);
            return do_filter(rls);
        }
    }
    else {
        if (sqrt_upd) {
            RLS::PolyRLS< double, RLS::SquareRootUpdate, 
                RLS::ExpWeightedRLS<double, RLS::SquareRootUpdate> > rls(Np);
            rls.setff(ff);
            return do_filter(rls);
        } else {
             RLS::PolyRLS< double, RLS::CovarianceUpdate, 
                RLS::ExpWeightedRLS<double, RLS::CovarianceUpdate> > rls(Np);
            rls.setff(ff);
            return do_filter(rls);            
        }
    }

    return 0;
}




