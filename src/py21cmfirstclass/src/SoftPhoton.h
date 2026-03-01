/*---- Soft Photon Heating Module ----
To test this module outside of p21f, create a link, define SOFT_PHOTON_TEST_MODE and then include the file, some functions are already defined in RadioExcess.h
*/

#define xax_len 5000
#ifdef SOFT_PHOTON_TEST_MODE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RadioExcess.h"
#define print_debug_info 0 // Don't show this if running outside of p21f
#endif

// First some useful general-purpose functions
void linspace(double xmin, double xmax, double *x, int nx)
{
	/*
	Create a linspace array
	-- inputs --
	xmin: minimum of x
	xmax: maximum of x
	x: pointer of pre-created x array
	nx: array size
	*/
	int idx;
	double dx;
	dx = (xmax - xmin) / ((double)nx - 1.0);
	for (idx = 0; idx < nx; idx++)
	{
		x[idx] = xmin + ((double)idx * dx);
	}
}

void logspace(double lgx_min, double lgx_max, double *x, int nx)
{
	/*
	Create a logspace array
	-- inputs --
	lgx_min: minimum of log10(x)
	lgx_max: maximum of log10(x)
	x: pointer of pre-created x array
	nx: array size
	*/
	int idx;
	double dlx;
	dlx = (lgx_max - lgx_min) / ((double)nx - 1.0);
	for (idx = 0; idx < nx; idx++)
	{
		x[idx] = pow(10.0, lgx_min + ((double)idx * dlx));
	}
}

double Integrate(double *x, double *fx, int nx, int method)
{
	/*
	Integrate fx over x: \int dx f(x)
	-- inputs --
	method: integration method. For equally spaced x methods 0 and 1 are pretty much the same, but trapz generally outperforms for non-equally spaced x
		0 - summation, simplest method
			\int dx f(x) = \sum_{i=0}^{n-1}f_idx_i
		1 - trapz
	*/
	double dx, result, f;
	int idx;
	result = 0.0;
	if (method == 0)
	{
		for (idx = 0; idx < nx - 1; idx++)
		{
			dx = x[idx + 1] - x[idx];
			f = fx[idx];
			result += f * dx;
		}
	}
	else if (method == 1)
	{
		for (idx = 0; idx < nx - 1; idx++)
		{
			dx = x[idx + 1] - x[idx];
			f = (fx[idx] + fx[idx + 1]) / 2.0;
			result += f * dx;
		}
	}
	return result;
}

double Lambda_BR(double z, double xe, double x, double Tk, double OmBh2, double YHe)
{
	/*
	Bremsstrahlung emissivity coefﬁcient (dimensionless) following Chluba's paper: 10.1093/mnras/stae2113
	Note: Gaunt factor in Chluba's codes is bugged!!! But equation in the paper should be correct.
	---- inputs ----
	xe: ionization fraction
	x: h*f/(kB * Tcmb)
	*/
	double kB, mec2, theta, alpha, c, h, lambda_e, gaunt, Np, r, Nb, vGHz, T4;
	c = 2.99792458E8;
	kB = 1.38064852E-23;
	h = 6.626070040818181818E-34;
	alpha = 1.0 / 137.03599976;

	// Gaunt factor following 10.1093/mnras/stad1539 (or 10.1051/0004-6361/201424643), mathematically identical to Chluba's paper
	vGHz = x * kB * 2.728 * (1.0 + z) / (h * 1.0E9);
	T4 = Tk / 1.0E4;
	gaunt = 5.9598540023339722 - sqrt(3.0) / M_PI * log(vGHz * pow(T4, -1.5));
	gaunt = log(exp(gaunt) + exp(1.0));

	mec2 = 510998.9461 * 1.602176634E-19;
	theta = Tk * kB / mec2;
	Nb = 0.25186318588882056 * OmBh2 / 0.02242 * pow(1.0 + z, 3.0); // rho_cr*OmB/mp, number density of protons assuming baryons are made entirely of Hydrogen
	Np = xe * (1.0 - YHe) * Nb;
	lambda_e = h * c / mec2;
	r = alpha * pow(lambda_e, 3.0) / (2.0 * M_PI * sqrt(6.0 * M_PI)) * pow(theta, -3.5) * Np * gaunt;
	return r;
}

double SoftPhoton_inj(double x, double xe, double z, double SFRD_II, double SFRD_III, double fR_II, double fR_III, double aR_II, double aR_III, double n_inj, double H, double YHe, double OmBh2)
{
	/*
	The injection term dn_inj/dtau
	-- inputs --
	x: h*f/(kB * Tcmb)
	xe: ionization fraction
	z: z
	SFRD_II: Pop II SFRD (comoving), in Msun/yr/Mpc^-3
	SFRD_III: Pop III SFRD (comoving), in Msun/yr/Mpc^-3
	fR_II: Pop II radio efficiency
	fR_III: Pop III radio efficiency
	aR_II: Pop II radio emission power index
	aR_III: Pop III radio emission power index
	n_inj: occupation number of injected radio photons
	H: Hubble
	YHe: Helium mass fraction, default should be 0.245
	OmBh2: Omega bh2
	*/
	double c, sigmaT, ne, nH0, nH, Tcmb, h, kB, f, sII, sIII, src, r1, r2, r, fHe;
	c = 2.99792458E8;
	sigmaT = 0.665245854E-28;
	kB = 1.38064852E-23;
	h = 6.626070040818181818E-34;
	Tcmb = 2.728 * (1.0 + z);
	fHe = YHe / (4.0 * (1.0 - YHe));
	nH0 = 0.1901567053460595 * OmBh2 / 0.02242;		 // Number density of H neclei today, in m^-3
	ne = nH0 * (1.0 + fHe) * xe * pow(1.0 + z, 3.0); // electron number density
	f = kB * Tcmb * x / h;							 // Photon frequency in Hz
	sII = fR_II * pow(f / 0.15E9, -aR_II) * SFRD_II;
	sIII = fR_III * pow(f / 0.15E9, -aR_III) * SFRD_III;
	src = sII + sIII;
	r1 = 5.507E11 * src * pow((1.0+z)/f, 3.0);
	r2 = 4.0 * n_inj * H;
	r = (r1 + r2) / (sigmaT * ne * c);
	return r;
}

double DeltaS_tilte(double x, double xe, double z, double SFRD_II, double SFRD_III, double fR_II, double fR_III, double aR_II, double aR_III, double n_inj, double H, double YHe, double OmBh2, double Tk, double LBR)
{
	/*
	DeltaS in Eq(8)
	-- inputs --
	Tk: matter temperature (K)
	LBR: Lambda_BR, can be computed by Lambda_BR function, used as input here to avoid repeated call to Lambda_BR
	*others*: see SoftPhoton_inj
	*/
	double Sinj, Tcmb, Y, ex, Smid, Sff, mec2, kB, xi, exi, f1, f2, DeltaS, r;
	mec2 = 510998.9461 * 1.602176634E-19;
	kB = 1.38064852E-23;

	// Soft Photon Injection
	Sinj = SoftPhoton_inj(x, xe, z, SFRD_II, SFRD_III, fR_II, fR_III, aR_II, aR_III, n_inj, H, YHe, OmBh2);

	// Mid-term, whatever that is
	Tcmb = 2.728 * (1.0 + z);
	ex = exp(x);
	Y = x * ex / pow(ex - 1.0, 2.0) * (x * (ex + 1.0) / (ex - 1.0) - 4.0);
	Smid = kB / mec2 * (Tk - Tcmb) * Y;

	// Sff term, Eq(9)
	// LBR = Lambda_BR(z, xe, x, Tk, OmBh2, YHe);
	xi = x * Tcmb / Tk;
	exi = exp(xi);
	f1 = LBR * (1.0 - 1.0 / exi) / pow(xi, 3.0);
	f2 = 1.0 / (exi - 1.0) - 1.0 / (ex - 1.0);
	Sff = f1 * f2;

	// DeltaS
	DeltaS = Sinj + Smid + Sff;

	r = DeltaS * pow(xi, 3.0) / (LBR * (1.0 - 1.0 / exi));
	return r;
}

double dtauff_dz(double z, double H, double xe, double LBR, double xi, double OmBh2, double YHe)
{
	//\frac{d \tau_{ff}}{dz} = - \frac{\sigma_T N_e c \Lambda_{BR} (1 - e^{-x_i})}{(1+z) H x_i^3}
	double sigmaT, ne, c, nH0, fHe, r;
	c = 2.99792458E8;
	sigmaT = 0.665245854E-28;
	nH0 = 0.1901567053460595 * OmBh2 / 0.02242; // Number density of H neclei today, in m^-3
	fHe = YHe / (4.0 * (1.0 - YHe));
	ne = nH0 * xe * (1.0 + fHe) * pow(1.0 + z, 3.0); // electron number density
	r = -sigmaT * ne * c * LBR * (1 - exp(-xi)) / ((1.0 + z) * H * pow(xi, 3.0));
	return r;
}

double Build_dTffdz_Kernel(double z, double xe, double Tk, double H, double YHe, double OmBh2, double *x, double *LBR, double *DeltaN, double *DeltaN_inj, double *EMS_inj, int nx)
{
	/*
	Compute Heating Rate dTk/dz
	-- inputs --
	z: z
	xe: ionization fraction
	H: Hubble (s^-1)
	Tk: Tk (K)
	x: pointer for hf/(kB Tcmb)
	LBR: pointer for Lambda_BR
	DeltaN: pointer for DeltaN
	DeltaN_inj: pointer of DeltaN_inj, erriely similar to ninj?
	EMS_inj: pointer of comoving emissivity output
	nx: size of pointers x, LBR, DeltaN
	*/
	double rho_cmb, Tcmb, kB, hb, c, dtdz, ne_ntot, Prefix, sigmaT, Integrand[xax_len], xi, ex, exi, dTkdz, EMS_Preix, fHe, nH0, ne;
	int idx;

	if (fabs((double)xax_len - (double)nx) > 1E-10)
	{
		fprintf(stderr, "Error: unexpected x length nx\n");
		exit(1);
	}

	kB = 1.38064852E-23;
	hb = 6.626070040818181818E-34 / (2.0 * M_PI);
	c = 2.99792458E8;
	sigmaT = 0.665245854E-28;

	Tcmb = 2.728 * (1.0 + z);
	rho_cmb = pow(M_PI, 2.0) * pow(kB * Tcmb, 4.0) / (15.0 * pow(hb * c, 3.0));
	dtdz = -1 / ((1 + z) * H);
	ne_ntot = xe / (1.0 + xe); // This is ne/(nH(1+xe+fHe+xe*fHe))
	Prefix = 2 * dtdz * ne_ntot / (3 * kB);
	Prefix = Prefix * 15 * rho_cmb * sigmaT * c * pow(Tk / Tcmb, 3.0) / pow(M_PI, 4.0); // Part not dependent on x integrand

	fHe = YHe / (4.0 * (1.0 - YHe));
	nH0 = 0.1901567053460595 * OmBh2 / 0.02242;		 // Number density of H neclei today, in m^-3
	ne = nH0 * (1.0 + fHe) * xe * pow(1.0 + z, 3.0); // electron number density
	// EMS_Preix = 30.0 * rho_cmb * sigmaT * ne * c * hb * pow(Tk*(1.0+z), 3.0)/(kB * pow(M_PI*Tcmb, 3.0) * Tcmb);
	EMS_Preix = 30.0 * rho_cmb * sigmaT * ne * c * hb * pow(Tk/(1.0+z), 3.0)/(kB * pow(M_PI*Tcmb, 3.0) * Tcmb);
	
	// Defining the integrand
	for (idx = 0; idx < nx; idx++)
	{
		ex = exp(x[idx]);
		xi = x[idx] * Tcmb / Tk;
		exi = exp(xi);
		Integrand[idx] = LBR[idx] * (1.0 - 1.0 / exi) * (DeltaN[idx] + 1.0 / (ex - 1.0) - 1.0 / (exi - 1));
		EMS_inj[idx] = EMS_Preix * LBR[idx] * (1.0 - 1.0/exi) * DeltaN_inj[idx];
		/*
		if (isnan(EMS_inj[idx]))
		{
			printf("VarX = %15E		%15E\n", DeltaN_inj[idx], xe);
		}
		*/
	}
	dTkdz = Prefix * Integrate(x, Integrand, nx, 1);
	return dTkdz;
}

double Compute_dTffdz(double *zax, double *dTdz, double *Hax, double *xe_ax, double *Tkax, double *SFRD_II_ax, double *SFRD_III_ax, double *dT_Radio_out, double fR_II, double fR_III, double aR_II, double aR_III, double OmBh2, double YHe, double redshift, double delta_z_TRadio, int nz)
{
	/*
	Compute Soft-Photon heating rate dT/dz for an array of z, many inputs must be passed externally, such that the code can be easily tested outside of 21cmFAST (e.g., with python) and it can
	merge streamelessly to 21cmFAST. Code computes dT/dz at all redshifts, which can be used by e.g. python for testing. Time used for Python to initialize input arrays seems subdominant compared
	to the min solver
	-- inputs --
	zax: z axis
	dTdz: axis for dT/dz output
	delta_z_TRadio: how instantaneous is the deposition, integrate in this z length when computing dT_Radio. To turn off this feature, set delta_z_TRadio to a large value, e.g., 1.0E10
	*/
	double xax[xax_len], LBR_ax[xax_len], DeltaN_ax[xax_len], DST_ax[xax_len], n_inj[xax_len], DeltaN_inj[xax_len], EMS_inj[xax_len];
	double z, xe, Tk, dz, dtff, xi, Tcmb, SFRD_II, SFRD_III, exp_dtff_inv, H, dninj_dtff, S_inj, lgx_min, lgx_max, f21;
	double x21, hb, kB, dT_Radio, EMS21, dT_Radio_PreFix, c;

	int xid, zid;
	FILE *OutputFile;

	// Initializing
	hb = 6.626070040818181818E-34 / (2.0 * M_PI);
	kB = 1.38064852E-23;
	c = 2.99792458E8;
	f21 = 1.42E9; // 21cm frequency in Hz
	lgx_min = -8.0;
	lgx_max = 2.0;
	x21 = f21 * 2.0 * M_PI * hb / (kB * 2.728 * (1.0 + redshift)); //log10(x21) is in [-4.78, -1.60] for redshift in [0, 1500], no need to worry about interpolation range overflow
	logspace(lgx_min, lgx_max, xax, xax_len);
	dT_Radio_PreFix = pow((1.0 + redshift)*c, 3.0) / (8.0 * M_PI * kB * f21 * f21);
	z = zax[0];
	if (zax[0] < zax[nz - 1])
	{
		fprintf(stderr, "Error: z needs to be descending, z0 = %7f, z_end = %7f\n", zax[0], zax[nz - 1]);
		exit(1);
	}
	if (zax[nz-1] < redshift)
	{
		// This can happen in test
		fprintf(stderr, "Your redshift is larger than smallest zax.\n");
		exit(1);
	}
	dTdz[0] = 0.0;
	for (xid = 0; xid < xax_len; xid++)
	{
		DeltaN_ax[xid] = 0.0;
		n_inj[xid] = 0.0;
		DeltaN_inj[xid] = 0.0;
	}

	if (print_debug_info)
	{
		OutputFile = fopen("/Users/cangtao/Desktop/tmp/tmp_test_Radio_Heating/dTff_dz_tmp.txt", "w");
	}

	// Start Evolving
	dT_Radio = 0.0;
	for (zid = 1; zid < nz; zid++)
	{
		z = zax[zid];
		dz = zax[zid] - zax[zid - 1];
		H = Hax[zid];
		xe = xe_ax[zid];
		Tk = Tkax[zid];
		Tcmb = 2.728 * (1.0 + z);
		SFRD_II = SFRD_II_ax[zid];
		SFRD_III = SFRD_III_ax[zid];
		for (xid = 0; xid < xax_len; xid++)
		{
			xi = xax[xid] * Tcmb / Tk;
			LBR_ax[xid] = Lambda_BR(z, xe, xax[xid], Tk, OmBh2, YHe);
			dtff = dtauff_dz(z, H, xe, LBR_ax[xid], xi, OmBh2, YHe) * dz;

			// Updating n_inj and DeltaN_inj, not very efficient because this needs to call SoftPhoton_inj which has already been called in DeltaS_tilte, can be optimized if this adds too much time
			S_inj = SoftPhoton_inj(xax[xid], xe, z, SFRD_II, SFRD_III, fR_II, fR_III, aR_II, aR_III, n_inj[xid], H, YHe, OmBh2);
			dninj_dtff = S_inj * pow(xi, 3.0) / (LBR_ax[xid] * (1.0 - exp(-xi)));
			n_inj[xid] = n_inj[xid] + dninj_dtff * dtff;
			// Finished updating n_inj

			exp_dtff_inv = exp(-dtff);
			DST_ax[xid] = DeltaS_tilte(xax[xid], xe, z, SFRD_II, SFRD_III, fR_II, fR_III, aR_II, aR_III, n_inj[xid], H, YHe, OmBh2, Tk, LBR_ax[xid]);
			DeltaN_ax[xid] = DeltaN_ax[xid] * exp_dtff_inv + DST_ax[xid] * (1.0 - exp_dtff_inv);
			DeltaN_inj[xid] = DeltaN_inj[xid] * exp_dtff_inv + dninj_dtff * (1.0 - exp_dtff_inv);
			// DeltaN_ax[xid] = DeltaN_ax[xid]*(exp_dtff_inv) + DST_ax[xid]*dtff;
			/*
			if (isnan(DeltaN_inj[xid]))
			{
				printf("Var3 = %15E\n", S_inj);
			}
			*/
		}
		dTdz[zid] = Build_dTffdz_Kernel(z, xe, Tk, H, YHe, OmBh2, xax, LBR_ax, DeltaN_ax, DeltaN_inj, EMS_inj, xax_len);
		
		// ======== Computing Radio attenuation ========
		if (fabs(z - zax[nz-1]) > delta_z_TRadio)
		{
			EMS21 = 0.0;
		}
		else
		{
			EMS21 = Interp_1D(x21, xax, EMS_inj, xax_len, 1, 0, 0);
		}
		dT_Radio -= dT_Radio_PreFix * EMS21 * dz /(H * (1.0+z));// remember that dz is negative!
		/*
		if (isnan(dT_Radio))
		{
			printf("Var = %15E	%15E\n", EMS21, z);
		}
		*/
		if (print_debug_info)
		{
			fprintf(OutputFile, "%7f   %7E\n", z, dTdz[zid]);
		}
	}
	if (print_debug_info)
	{
		fclose(OutputFile);
	}
	// Still need to figure out how to pass dT_Radio outside
	*dT_Radio_out = dT_Radio;
	return dTdz[nz - 1];
}
