#include "2dbat.multinest.h"

void run(int IS, int mmodal, int ceff, int nlive, double tol, double efr, int ndims, int nPar, int nClsPar, 
int maxModes, int updInt, double Ztol, char root[], int seed, int *pWrap, int fb, int resume, int outfile, 
int initMPI, double logZero, int maxiter, void (*LogLike)(double *, int *, int *, double *, TR_ringParameters *), 
void (*dumper)(int *, int *, int *, double **, double **, double **, double *, double *, double *, double *, TR_ringParameters *), 
TR_ringParameters *context)
{
    int i;
    for (i = strlen(root); i < 500; i++) root[i] = ' ';

        NESTRUN(&IS, &mmodal, &ceff, &nlive, &tol, &efr, &ndims, &nPar, &nClsPar, &maxModes, &updInt, &Ztol,
        root, &seed, pWrap, &fb, &resume, &outfile, &initMPI, &logZero, &maxiter, LogLike, dumper, context);
}


/***********************************************************************************************************************/

