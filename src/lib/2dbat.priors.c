#include "2dbat.priors.h"

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
// Inverse of complimentary error function in double precision
double dierfc(double y)
{
    double _dierfc;
    double qa,qb,qc,qd,q0,q1,q2,q3,q4,pa,pb,p0,p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,p15,p16,p17,p18;
    double p19,p20,p21,p22,x,z,w,u,s,t;
    double infinity;

    infinity=5.0E0;
    qa=9.16461398268964E-01;
    qb=2.31729200323405E-01;
    qc=4.88826640273108E-01;
    qd=1.24610454613712E-01;
    q0=4.99999303439796E-01;
    q1=1.16065025341614E-01;
    q2=1.50689047360223E-01;
    q3=2.69999308670029E-01;
    q4=-7.28846765585675E-02;

    pa=3.97886080735226000E+00;
    pb=1.20782237635245222E-01;
    p0=2.44044510593190935E-01;
    p1=4.34397492331430115E-01;
    p2=6.86265948274097816E-01;
    p3=9.56464974744799006E-01;
    p4=1.16374581931560831E+00;
    p5=1.21448730779995237E+00;
    p6=1.05375024970847138E+00;
    p7=7.13657635868730364E-01;
    p8=3.16847638520135944E-01;
    p9=1.47297938331485121E-02;
    p10=-1.05872177941595488E-01;
    p11=-7.43424357241784861E-02;
    p12=2.20995927012179067E-03;
    p13=3.46494207789099922E-02;
    p14=1.42961988697898018E-02;
    p15=-1.18598117047771104E-02;
    p16=-1.12749169332504870E-02;
    p17=3.39721910367775861E-03;
    p18=6.85649426074558612E-03;
    p19=-7.71708358954120939E-04;
    p20=-3.51287146129100025E-03;
    p21=1.05739299623423047E-04;
    p22=1.12648096188977922E-03;

    if (y==0.0)
    {
        _dierfc=infinity;
        return _dierfc;
    }
    z=y;
    if (y > 1) z=2-y;
    w=qa-log(z);
    u=sqrt(w);
    s=(qc+log(u))/w;
    t=1/(u+qb);
    x=u*(1-s*(0.5E0+s*qd))-((((q4*t+q3)*t+q2)*t+q1)*t+q0)*t;
    t=pa/(pa+x);
    u=t-0.5E0;
    s=(((((((((p22*u+p21)*u+p20)*u+p19)*u+p18)*u+p17)*u+p16)*u+p15)*u+p14)*u+p13)*u+p12;
    s=((((((((((((s*u+p11)*u+p10)*u+p9)*u+p8)*u+p7)*u+p6)*u+p5)*u+p4)*u+p3)*u+p2)*u+p1)*u+p0)*t-z*exp(x*x-pb);
    x=x+s*(1+x*s);
    if (y > 1) x=-x;
    
    _dierfc=x;
    return _dierfc;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  JeffreysPrior
double JeffreysPrior(double r, double x1, double x2)
{
    double lx1, lx2,_JeffreysPrior;
    if(r < 0.)
    {
        _JeffreysPrior=-1.0E32;
    }
    else
    {
        if(x1 == 0)
            lx1 = -1.0E20;
        else
            lx1=log(x1);

        if(x2 == 0)
            lx2 = -1.0E20;
        else
            lx2=log(x2);

        _JeffreysPrior = 1.0/(r*(lx2-lx1));
    }

    return _JeffreysPrior;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]
double gaussian_prior(double r, double mu, double sigma)
{
    double _gaussian_prior;
    double SqrtTwo;
    SqrtTwo=1.414213562E0;

    if (r <= 1.0E-16 || (1.0E0-r) <= 1.0E-16)
        _gaussian_prior = -1.0E32;
    else
        _gaussian_prior = mu + sigma*SqrtTwo*dierfc(2.E0*(1.E0-r));

    return _gaussian_prior;
}



// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  Gaussian[mean=mu,variance=sigma**2]
double gaussian_prior_paincl_skew(double r, double mu, double sigma_N, double sigma_W)
{
    double _gaussian_prior;
    double SqrtTwo;
    SqrtTwo=1.414213562E0;

    if (r <= 1.0E-16 || (1.0E0-r) <= 1.0E-16)
    {
        _gaussian_prior = -1.0E32;
    }
    else
    {
        if(mu < 0.5 && r < mu) // left skewed and smaller than mu : use narrower sigma
        {
            _gaussian_prior = mu + sigma_N*SqrtTwo*dierfc(2.E0*(1.E0-r));
            if(_gaussian_prior > mu) _gaussian_prior = mu - fabs(_gaussian_prior-mu);
        }
        else if(mu < 0.5 && r >= mu) // left skewed and larger than mu : use wider sigma
        {
            _gaussian_prior = mu + sigma_W*SqrtTwo*dierfc(2.E0*(1.E0-r));
            if(_gaussian_prior < mu) _gaussian_prior = mu + fabs(mu-_gaussian_prior);
        }
        else if(mu >= 0.5 && r < mu) // right skewed and smaller than mu : use wider sigma
        {
            _gaussian_prior = mu + sigma_W*SqrtTwo*dierfc(2.E0*(1.E0-r));
            if(_gaussian_prior > mu) _gaussian_prior = mu - fabs(_gaussian_prior-mu);
        }
        else if(mu >= 0.5 && r >= mu) // right skewed and smaller than mu : use narrower sigma
        {
            _gaussian_prior = mu + sigma_N*SqrtTwo*dierfc(2.E0*(1.E0-r));
            if(_gaussian_prior < mu) _gaussian_prior = mu + fabs(mu-_gaussian_prior);
        }
    }

    return _gaussian_prior;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  Gaussian_skew
double gaussian_prior_esigma_skew(double r, double mu, double mu_max, double sigma_N, double sigma_W)
{
    double _gaussian_prior;
    double SqrtTwo;
    SqrtTwo=1.414213562E0;

    if (r <= 1.0E-16 || (1.0E0-r) <= 1.0E-16)
    {
        _gaussian_prior = -1.0E32;
    }
    else
    {
        if(r < mu/mu_max) // left skewed and smaller than mu : use narrower sigma
        {
            _gaussian_prior = mu + (sigma_N*SqrtTwo*dierfc(2.E0*(1.E0-r)));
        }
        else if(r >= mu/mu_max) // left skewed and larger than mu : use wider sigma
        {
            _gaussian_prior = mu + (sigma_W*SqrtTwo*dierfc(2.E0*(1.E0-r)));
        }
    }
    return fabs(_gaussian_prior);
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  Gaussian_skew[mean=mu,variance=sigma**2]
double gaussian_prior_nrrho_skew(double r, double mu, double sigma_N, double sigma_W)
{
    double _gaussian_prior;
    double SqrtTwo;
    SqrtTwo=1.414213562E0;

    if (r <= 1.0E-16 || (1.0E0-r) <= 1.0E-16)
    {
        _gaussian_prior = -1.0E32;
    }
    else
    {
        if(r < mu) // left skewed and smaller than mu : use narrower sigma
        {
            _gaussian_prior = mu + sigma_N*SqrtTwo*dierfc(2.E0*(1.E0-r));
        }
        else if(r >= mu) // left skewed and larger than mu : use wider sigma
        {
            _gaussian_prior = mu + sigma_W*SqrtTwo*dierfc(2.E0*(1.E0-r));
        }
    }

    return _gaussian_prior;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  Uniform[x1:x2]
double uniform_priors(double r, double x1, double x2)
{
    double _uniform_priors;

    _uniform_priors = x1 + r*(x2 - x1);
    return _uniform_priors;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  LogUniform[x1:x2]
double loguniform_priors(double r, double x1, double x2)
{
    double _loguniform_priors;
    double SqrtTwo=1.414213562E0;
    double LogPrior;
    double lx1, lx2;

    if (r <= 0)
    {
        _loguniform_priors = -1.0E32; 
    }
    else
    {
        lx1 = log10(x1);
        lx2 = log10(x2);
        _loguniform_priors = pow(10.,  (lx1 + r*(lx2-lx1)));
    }

    return _loguniform_priors;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//! Uniform[0:1]  ->  LogNormal[mode=a,width parameter=sigma]
double logNormalPrior(double r, double a, double sigma)
{
    double _logNormalPrior;
    double SqrtTwo=1.414213562E0;
    double bracket;

    bracket = sigma*sigma + sigma*SqrtTwo*dierfc(2.0*r);
    _logNormalPrior = a*exp(bracket);

    return _logNormalPrior;
}
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


// --- End of line


