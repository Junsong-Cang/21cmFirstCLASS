#ifndef History_box_DIM
// Running & Testing outside of 21cmFAST
#include <stdio.h>
#include <math.h>
#endif
#define Use_Conde_Concentration 1

double Halo_Concentration(double m, double z, double h)
{
    /*
        Halo concentration, use Ref2 if Use_Conde_Concentration, otherwise Ref1
        Ref1: F. Ziparo, S. Gallerani, A. Ferrara, and F. Vito, Mon. Not. R. Astron. Soc. 517, 1086 (2022).
        Ref2: M. A. Sánchez-Conde and F. Prada, Mon. Not. R. Astron.Soc. 442, 2271 (2014).
    */

    double c0, c1, c2, c3, c4, c5, x, r, LgC;
    
    c0 = 37.5153;
    c1 = -1.5093;
    c2 = 1.636E-2;
    c3 = 3.66E-4;
    c4 = -2.89237E-5;
    c5 = 5.32E-7;
    x = log(m*h);
    r = c0 + c1*x + c2*pow(x, 2.0)+ c3*pow(x, 3.0)+ c4*pow(x, 4.0)+ c5*pow(x, 5.0);
    r /= 1.0+z;
    if (! Use_Conde_Concentration)
    {
        LgC = 1.071 - 0.098 * (log10(m) - 12.0);
        r = pow(10.0, LgC)/(1.0+z);
    }
    return r;
}

double Halo_F_fun(double x)
{// Function F defined in Eq.(8) of Ref1
    return log(1.0+x) - x/(1.0+x);
}

void VirialStats(double z, double m, double h, double OmM, double *Tvir, double *Rvir)
{
    /*
    Halo Virial temperature (K) and radius (m)
    */
    double OmL, OmR, OMZ, d, DeltaC, u, kpc;
    OmR = 9.1E-5; // Radiation density, assuming that neutrinos are massless 
    kpc = 3.086E19;
    OmL = 1.0 - OmM - OmR;
    OMZ = OmM*pow(1.0 + z, 3.0) / (OmM*pow(1.0 + z, 3.0) + OmL);
    d = OMZ - 1.0;
    DeltaC = 18.0 * pow(M_PI, 2.0) + 82.0*d - 39.0*pow(d, 2.0);
    u = 1.22;
    *Tvir = 1.98E4 * (u/0.6) * pow(m*h/1.0E8, 2.0/3.0) * pow(OmM * DeltaC/(OMZ * 18.0 * pow(M_PI, 2.0)), 1.0/3.0) * (1.0+z)/10.0;
    *Rvir = 0.784 * pow(m*h/1.0E8, 1.0/3.0) * pow(OmM * DeltaC/(OMZ*18.0*M_PI*M_PI), -1.0/3.0) * (10.0/(1.0+z)) /h * kpc;
}

double Escape_Velocity_SQ(double x, double m, double z, double h, double OmM, double *Tvir_out)
{
    double GravConst, msun, small, HaloC, Tvir, Rvir, vcv2, cx, vc2, result;
    GravConst = 6.6740831313131E-11;
    msun = 1.98847E30;
    small = 1.0E-20;
    HaloC = Halo_Concentration(m, z, h);
    VirialStats(z, m, h, OmM, &Tvir, &Rvir);
    vcv2 = GravConst*m*msun/Rvir;
    cx = HaloC * x;
    vc2 = vcv2 * Halo_F_fun(cx)/Halo_F_fun(HaloC);
    result = 2.0 * vc2 * (Halo_F_fun(cx) + cx/(1 + cx))/(x * Halo_F_fun(cx));
    if (x < small)
    {
        result = 2 * vcv2 * HaloC/Halo_F_fun(HaloC);
    }
    *Tvir_out = Tvir;
    return result;
}

double Halo_Baryon_Profile_kernel(double z, double m, double x, double h, double OmM)
{
    /*
    Pre-normalization halo baryon density profile
    */
    double mp, kB, ve20, ve2, Tvir, u, result;
    mp = 1.67262158e-27;
    kB = 1.38064852E-23;
    ve20 = Escape_Velocity_SQ(0.0, m, z, h, OmM, &Tvir);
    ve2 = Escape_Velocity_SQ(x, m, z, h, OmM, &Tvir);
    u = 1.22;
    result = - u * mp * (ve20 - ve2)/(2.0 * kB * Tvir);
    result = exp(result);
    if (x >= 1.0)
    {
        result = 0.0;
    }
    return result;
}

void Halo_Baryon_Profile(double z, double m, double *Rvir_out, double *RhoSQ_out, double *xax, double *rho, int nx, double h, double OmM, double OmB)
{
    int idx;
    double dx, MassNorm, RhoSQ, Rvir, Tvir, msun, rho0;
    
    MassNorm = 0.0;
    RhoSQ = 0.0; // \int dV rho^2, in SI unit, i.e., kg^2/m^3
    for (idx=0; idx<nx; idx++)
    {
        rho[idx] = Halo_Baryon_Profile_kernel(z, m, xax[idx], h, OmM);
        dx = idx==nx-1? 0.0 : xax[idx+1] - xax[idx];
        MassNorm += rho[idx] * pow(xax[idx], 2.0) * dx;
        RhoSQ += pow(rho[idx]*xax[idx], 2.0) * dx;
    }
    
    // Normalize halo profile and density^2

    VirialStats(z, m, h, OmM, &Tvir, &Rvir);
    msun = 1.98847E30;
    rho0 = m * OmB / (4 * M_PI * MassNorm * OmM * pow(Rvir, 3.0)) * msun;
    RhoSQ *= 4.0 * M_PI * pow(rho0, 2.0) * pow(Rvir, 3.0);
    for (idx=0; idx<nx; idx++)
    {
        rho[idx] = rho0 * rho[idx];
    }
    *Rvir_out = Rvir;
    *RhoSQ_out = RhoSQ;
}
