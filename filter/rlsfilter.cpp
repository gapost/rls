#include <iostream>
#include <string>
#include "RLS_Estimation_Object.h"

using namespace std;
using namespace arma;

template<class _RLS>
int do_filter(_RLS& rls)
{
    unsigned int i = 0;
    while(!cin.eof())
    {
        double data;
        cin >> data;
        if (!cin.good()) return -1;
        i++;
        rls.update_par(data);
        cout << rls.getEstimatedOutput();
        const vec& w = rls.getEstimatedParameters();
        int np = w.size();
        for(int j=0; j<np; j++) cout << '\t' << w(j);
        cout << endl;
    }
    return 0;
}

int main(int argn, char** argv)
{
    int Np = 2; // poly order
    double ff = 0.98; // forgetting factor
    int M = 0; // block width

    const char* usage =
        "Usage:\n"
        "  rlstst [options] \n\n"
        "Options:\n"
        "  -nX  : fitting polynomial order, X=1,2,3,... (default = 2) \n"
        "  -ffX : forgetting factor, 0<X<=1, (default = 0.98) \n"
        "  -wX  : block width, X=0 means no block (default)\n\n"
        "rlstst takes input from stdin and writes to stdout. Use redirection for\n"
        "file input/output, e.g.:\n"
        "  rlstst < in.dat > out.dat\n";

    // parse options
    for(int i=1; i<argn; i++) {
        string opt(argv[i]);
        if (opt.find("-ff")==0) { // forgetting factor
            opt.erase(0,3);
            ff = atof(opt.c_str());
        }
        else if (opt.find("-n")==0) { // poly order
            opt.erase(0,2);
            Np = atoi(opt.c_str());
        }
        else if (opt.find("-w")==0) { // block width
            opt.erase(0,2);
            M = atoi(opt.c_str());
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
        RLS::PolyRLS<double, RLS::BlockRLS<double> > rls(Np, ff, M, 1.);
        return do_filter(rls);
    }
    else {
        RLS::PolyRLS<double, RLS::RLS_Estimator<double> > rls(Np, ff, 1.);
        return do_filter(rls);
    }

    return 0;
}




