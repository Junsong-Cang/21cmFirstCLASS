// Things needed for Radio excess

// nu0 is degenerate with fR so no reason to leave this as a param
#define astro_nu0 0.15	   // in GHz
#define History_box_DIM 20 // number of quantities to be saved in History_box
#define zax_len_FF 1000 // axis length for soft-photon heating

int Find_Index(double *x_axis, double x, int nx)
{
	/*
	Find the index of closest left element, need this for interpolation
	range handle:
		if x is on the LEFT of x_axis[0] : return -1
		if x is on the RIGHT of x_axis[nx-1] : return nx
	*/
	double x1, x2, x3;
	int id1, id2, id3, Stop, s1, s2, s3, idx, count, reversed;
	id1 = 0;
	id3 = nx - 1;
	Stop = 0;
	x1 = x_axis[id1];
	x3 = x_axis[id3];
	reversed = x1 < x3 ? 0 : 1;
	if (!reversed)
	{
		if (x <= x1)
		{
			Stop = 1;
			idx = -1;
		}
		if (x >= x3)
		{
			Stop = 1;
			idx = nx - 1;
		}
	}
	else
	{
		// printf("x1 = %f, x3 = %f\n", x1, x3);
		if (x >= x1)
		{
			Stop = 1;
			idx = -1;
		}
		if (x <= x3)
		{
			Stop = 1;
			idx = nx - 1;
		}
	}

	count = 0;
	while (Stop == 0)
	{
		count = count + 1;
		id2 = (int)round((((double)(id1 + id3))) / 2.0);
		if (id3 == id1 + 1)
		{
			idx = id1;
			Stop = 1;
		}

		x1 = x_axis[id1];
		x2 = x_axis[id2];
		x3 = x_axis[id3];

		if (!reversed)
		{
			if (x < x2)
			{
				id3 = id2;
			}
			else
			{
				id1 = id2;
			}
		}
		else
		{
			if (x < x2)
			{
				id1 = id2;
			}
			else
			{
				id3 = id2;
			}
		}
		if (count > 100)
		{
			fprintf(stderr, "Error @ Find_Index: solution not found after 100 iterations, x_axis[0] = %E, x = %E, x_axis[-1] = %E.\n", x_axis[0], x, x_axis[nx - 1]);
			exit(1);
		}
	}

	// printf("Stopping, id1 = %d, id3 = %d, x1 = %f, x = %f, x3 = %f, idx = %d\n", id1, id3, x1, x, x3, idx);

	return idx;
}

double Interp_1D(double x, double *x_axis, double *y_axis, int nx, int Use_LogX, int Use_LogY, int Overflow_Handle)
{
	/* Find value of y at x
	Use_LogX : whether to use log axis for x
	Use_LogY : whether to use log axis for y
	Overflow_Handle : what to do if x is not in x_axis
					  0 : raise error and exit
					  1 : give nearest value
					  2 : give -100
	*/
	int id1, id2;
	double x1, x2, y1, y2, x_, r, Small;
	id1 = Find_Index(x_axis, x, nx);
	Small = 1e-280;

	if (id1 == -1)
	{
		if (Overflow_Handle == 1)
		{
			r = y_axis[0];
		}
		else if (Overflow_Handle == 2)
		{
			r = -100.;
		}
		else
		{
			fprintf(stderr, "Error from Interp_1D: x is not in range, axis range: [%E   %E], x = %E\n", x_axis[0], x_axis[nx - 1], x);
			exit(1);
		}
	}
	else if (id1 == nx - 1)
	{
		if (Overflow_Handle == 1)
		{
			r = y_axis[nx - 1];
		}
		else if (Overflow_Handle == 2)
		{
			r = -100.;
		}
		else
		{
			fprintf(stderr, "Error from Interp_1D: x is not in range, axis range: [%E   %E], x = %E\n", x_axis[0], x_axis[nx - 1], x);
			exit(1);
		}
	}
	else
	{
		id2 = id1 + 1;
		if (!Use_LogX)
		{
			x1 = x_axis[id1];
			x2 = x_axis[id2];
			x_ = x;
		}
		else
		{
			// Detect negative element
			x1 = x_axis[id1];
			x2 = x_axis[id2];
			if (((x1 < 0) || (x2 < 0)) || (x < 0))
			{
				fprintf(stderr, "cannot use LogX for axis or x with negative element\n");
				exit(1);
			}

			x1 = log(x1);
			x2 = log(x2);
			x_ = log(x);
		}
		y1 = y_axis[id1];
		y2 = y_axis[id2];

		if (Use_LogY)
		{
			// This is to avoid nan at log
			if ((y1 < 0) || (y2 < 0))
			{
				fprintf(stderr, "Cannot use LogY for axis with negative element. Info: x1 = %E, x =  %E, x2 = %E, y1 = %E, y2 = %E\n", x1, x, x1, y1, y2);
				exit(1);
			}

			y1 = y1 > Small ? y1 : Small;
			y2 = y2 > Small ? y2 : Small;
			y1 = log(y1);
			y2 = log(y2);
		}

		r = (y2 - y1) * (x_ - x1) / (x2 - x1) + y1;

		if (Use_LogY)
		{
			r = exp(r);
		}
		// printf("x_ = %f, x1 = %f, x2 = %f, y1 = %f, y2 = %f\n", x_, x1, x2, y1, y2);
	}

	return r;
}

double History_box_Interp(struct TsBox *previous_spin_temp, double z, int Type, int Overflow_Handle)
{
	/*
	Interpolate to find quantities archived in History_box
	Initial test shows very good (0.3%) consistency
	---- inputs ----
	z: redshift
	Type: what do you want from the box
		1 - Phi_II
		2 - Phi_III
		3 - Tk
		4 - mturn_II
		5 - mturn_III
		6 - SFRD_EoR_MINI
		7 - xH
	*/
	int ArchiveSize, idx, head, zid, fid;
	// Very generous with memory, nobody is gonna run lightcones with 1000 timesteps (right?)
	double z_axis[1000], f_axis[1000], r;

	ArchiveSize = (int)round(previous_spin_temp->History_box[0]);
	if (previous_spin_temp->first_box || ArchiveSize < 2)
	{
		// Too early to get reliable results?
		if ((Type == 1 || Type == 2) || Type == 6)
		{// Set Phi_II, Phi_III, SFRD_EoR_MINI to 0
			return 0.0;
		}
		else if (Type == 3)
		{
			return global_params.TK_at_Z_HEAT_MAX;
		}
		else if (Type == 4 || Type == 5)
		{
			return 1.0E20;
		}
		else
		{
			fprintf(stderr, "Error in History_box_Interp: Exception not set for Type = %d.\n", Type);
		}
	}

	if (ArchiveSize > 800)
	{
		fprintf(stderr, "Error: ArchiveSize exceeds z_axis size.\n");
		Throw(ValueError);
	}

	// Fill axis
	for (idx = 0; idx < ArchiveSize; idx++)
	{
		head = idx * History_box_DIM + 1;
		if ((Type == 1) || (Type == 2))
		{
			// for Phi z_axis should be zpp
			zid = head + 4;
		}
		else
		{
			zid = head;
		}
		if (Type == 1)
		{ // Phi
			fid = head + 1;
		}
		else if (Type == 2)
		{ // Phi3
			fid = head + 3;
		}
		else if (Type == 3)
		{ // Tk
			fid = head + 2;
		}
		else if (Type == 4)
		{ // mturn
			fid = head + 5;
		}
		else if (Type == 5)
		{ // mturn_III
			fid = head + 6;
		}
		else if (Type == 6)
		{ // SFRD_MINI_EoR
			fid = head + 7;
		}
		else if (Type == 7)
		{ // xH
			fid = head + 8;
		}
		else
		{
			LOG_ERROR("Wrong Type setting, must be in [1, 7].\n");
			Throw(ValueError);
		}
		z_axis[idx] = previous_spin_temp->History_box[zid];
		f_axis[idx] = previous_spin_temp->History_box[fid];
	}
	// Use_LogY for all except xH, Interp_1D already sets a floor for very small Y
	if (Type == 7)
	{
		r = Interp_1D(z, z_axis, f_axis, ArchiveSize, 0, 0, Overflow_Handle);
	}
	else
	{
		r = Interp_1D(z, z_axis, f_axis, ArchiveSize, 0, 1, Overflow_Handle);
	}

	return r;
}

double Get_Radio_Temp_HMG(struct TsBox *previous_spin_temp, struct TsBox *this_spin_temp, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, struct FlagOptions *flag_options, double zpp_max, double redshift, double Z_HEAT_MAX)
{

	/* Find Radio Temp from sources in redshifts [zpp_max, Z_Heat_max]
	---- inputs ----
	zpp_max: maximum zpp
	redshift: redshift at which you want to compute radio temp
	*/

	double z1, z2, dz, Phi, Phi_mini, z, fun_ACG, fun_MCG, Radio_Temp, Radio_Prefix_ACG, Radio_Prefix_MCG;
	int nz, zid, RadioSilent;

	nz = 1000;

	if (flag_options->USE_RADIO_ACG)
	{
		Radio_Prefix_ACG = 113.6161 * astro_params->fR * cosmo_params->OMb * (pow(cosmo_params->hlittle, 2)) * (astro_params->F_STAR10) * pow(astro_nu0 / 1.4276, astro_params->aR) * pow(1 + redshift, 3 + astro_params->aR);
	}
	else
	{
		Radio_Prefix_ACG = 0.0;
	}

	if (flag_options->USE_RADIO_MCG)
	{
		Radio_Prefix_MCG = 113.6161 * astro_params->fR_mini * cosmo_params->OMb * (pow(cosmo_params->hlittle, 2)) * (astro_params->F_STAR7_MINI) * pow(astro_nu0 / 1.4276, astro_params->aR_mini) * pow(1 + redshift, 3 + astro_params->aR_mini);
	}
	else
	{
		Radio_Prefix_MCG = 0.0;
	}

	if (flag_options->USE_RADIO_ACG || flag_options->USE_RADIO_MCG)
	{
		RadioSilent = 0;
	}
	else
	{
		RadioSilent = 1;
	}

	if ((RadioSilent || redshift > Z_HEAT_MAX - 0.8) || this_spin_temp->first_box)
	{
		Radio_Temp = 0.0;
	}
	else
	{
		z2 = previous_spin_temp->History_box[5] - 0.01;
		z1 = zpp_max;
		if (z1 > z2)
		{
			Radio_Temp = 0.0;
		}
		else
		{
			dz = (z2 - z1) / (((double)nz) - 1);

			z = z1;
			Radio_Temp = 0.0;

			for (zid = 1; zid <= nz; zid++)
			{
				Phi = History_box_Interp(previous_spin_temp, z, 1, 1);
				Phi_mini = History_box_Interp(previous_spin_temp, z, 2, 1);
				Phi = Phi > 1e-50 ? Phi : 0.;
				Phi_mini = Phi_mini > 1e-50 ? Phi_mini : 0.;
				fun_ACG = Radio_Prefix_ACG * Phi * pow(1 + z, astro_params->X_RAY_SPEC_INDEX - astro_params->aR) * dz;
				fun_MCG = Radio_Prefix_MCG * Phi_mini * pow(1 + z, astro_params->X_RAY_SPEC_INDEX - astro_params->aR_mini) * dz;
				if (z > astro_params->Radio_Zmin)
				{
					Radio_Temp += fun_ACG + fun_MCG;
				}
				z += dz;
			}
		}
	}
	if (isfinite(Radio_Temp) == 0)
	{
		fprintf(stderr, "Error @ Get_Radio_Temp_HMG :  Radio_Temp is NaN! Crash imminent\n");
		Throw(ValueError);
	}
	return Radio_Temp;
}

void Refine_T_Radio(struct TsBox *previous_spin_temp, struct TsBox *this_spin_temp, float prev_redshift, float redshift, struct AstroParams *astro_params, struct FlagOptions *flag_options)
{
	/* An analytic formula to eliminate numerical kinks from Radio_Zmin
	Need to be careful when executed wihin a mpi loop
	This has a number of issues:
	1. Only applicapable to sources with same spectra shape
	2. This is called within a box_ct loop, but I am doing another r_ct looop here. Results are the same but this wastes memory and cpu
	*/
	int box_ct;
	float T_prev, T_now, Conversion_Factor;

	if (redshift < astro_params->Radio_Zmin && (flag_options->USE_RADIO_MCG || flag_options->USE_RADIO_ACG))
	{

		if (flag_options->USE_RADIO_ACG && (!flag_options->USE_RADIO_MCG))
		{ // Only ACG
			Conversion_Factor = pow((1 + redshift) / (1 + prev_redshift), 3 + astro_params->aR);
		}
		else if (flag_options->USE_RADIO_MCG && (!flag_options->USE_RADIO_ACG))
		{ // Only MCG
			Conversion_Factor = pow((1 + redshift) / (1 + prev_redshift), 3 + astro_params->aR_mini);
		}
		else if (flag_options->USE_RADIO_ACG && flag_options->USE_RADIO_MCG)
		{ // co-exist
			if (fabs(astro_params->aR - astro_params->aR_mini) < 0.0001)
			{ // same spectra
				Conversion_Factor = pow((1 + redshift) / (1 + prev_redshift), 3 + astro_params->aR_mini);
			}
			else
			{ // different spectra, raise error
				LOG_ERROR("Using multiple radio sources with different spectra");
				Throw(ValueError);
			}
		}
		for (box_ct = 0; box_ct < HII_TOT_NUM_PIXELS; box_ct++)
		{
			this_spin_temp->Trad_box[box_ct] = Conversion_Factor * previous_spin_temp->Trad_box[box_ct];
		}
	}
}

float Phi_2_SFRD(double Phi, double z, double H, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, int Use_MINI)
{
	/*
	Convert Phi to SFRD in msun/Mpc^3/yr
	*/

	double f710, SFRD;

	if (Use_MINI)
	{
		f710 = astro_params->F_STAR7_MINI;
	}
	else
	{
		f710 = astro_params->F_STAR10;
	}

	SFRD = Phi * cosmo_params->OMb * RHOcrit * f710 * pow(1.0 + z, astro_params->X_RAY_SPEC_INDEX + 1.0) * H * SperYR;

	return SFRD;
}

float find_redshift_step(int idx)
{
	float dx, x1, x, z;
	dx = log(global_params.ZPRIME_STEP_FACTOR);
	x1 = log(1.0 + global_params.Z_HEAT_MAX);
	x = x1 - ((float) idx - 1.0) * dx;
	z = exp(x) - 1.0;
	return z;
}

double Get_EoR_Radio_mini(struct TsBox *this_spin_temp, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, float redshift)
{
	int idx, nz, ArchiveSize, head, terminate;
	double nion, dz, fun, dT, T, Prefix, Phi, z, z_prev, mt, mc, Mlim_Fstar_MINI, z_axis[400], nion_axis[400], zmin, zmax;
	nz = 400;
	terminate = 0; // sometimes in mpi or parralel loops python might proceed even with error, use this to give NaN which will terminate the simulation by various NaN checkpoints
	
	if ((this_spin_temp->first_box) || (redshift > find_redshift_step(2) - 0.5))
	{
		T = 0;
	}
	else
	{
		ArchiveSize = (int)round(this_spin_temp->History_box[0]);
		if (ArchiveSize > 3) // U need to have sufficiently large interp table
		{
			Mlim_Fstar_MINI = Mass_limit_bisection(global_params.M_MIN_INTEGRAL, global_params.M_MAX_INTEGRAL, astro_params->ALPHA_STAR_MINI,
												   astro_params->F_STAR7_MINI * pow(1e3, astro_params->ALPHA_STAR_MINI));

			// First fill z_axis and nion_axis
			if (ArchiveSize > 390)
			{
				fprintf(stderr, "Error @ Get_EoR_Radio_mini: Running with very fine z time steps, z and nion axis is not large enough.\n");
				Throw(ValueError);
			}
			for (idx = 0; idx < ArchiveSize; idx++)
			{
				head = idx * History_box_DIM + 1;
				z = this_spin_temp->History_box[head];
				z_axis[idx] = z;
				mc = atomic_cooling_threshold(z);
				mt = this_spin_temp->History_box[head + 6];
				if (mt < 1.0e2)
				{
					fprintf(stderr, "Error @ Get_EoR_Radio_mini (p21f): mturn is smaller than 100. mturn = %.3E, redshift = %.3f, contaminated z = %.3f\n", mt, redshift, z);
					Throw(ValueError);
					terminate = 1;
				}
				if (mt > 1.0E15)
				{
					nion_axis[idx] = 0.0;
				}
				else
				{
					nion_axis[idx] = Nion_General_MINI(z, global_params.M_MIN_INTEGRAL, mt, mc, astro_params->ALPHA_STAR_MINI, 0., astro_params->F_STAR7_MINI, 1., Mlim_Fstar_MINI, 0.);
				}
			}

			// Now interpolate for finer z
			Prefix = 113.6161 * astro_params->fR_mini * cosmo_params->OMb * (pow(cosmo_params->hlittle, 2)) * astro_params->F_STAR7_MINI * pow(astro_nu0 / 1.4276, astro_params->aR_mini) * pow(1 + redshift, 3.0 + astro_params->aR_mini);
			zmax = z_axis[0];				// z_heat_max
			zmin = z_axis[ArchiveSize - 1]; // this is current redshift
			dz = (zmax - zmin) / (((double)nz) - 1.0);
			T = 0.;
			for (idx = 0; idx < nz; idx++)
			{
				z = zmin + ((double)idx) * dz;
				nion = Interp_1D(z, z_axis, nion_axis, ArchiveSize, 0, 1, 1);
				fun = Prefix * nion / astro_params->t_STAR / pow(1.0 + z, astro_params->aR_mini + 1);
				dT = fun * dz;
				T = T + dT;
			}
		}
		else
		{
			if (terminate == 1)
			{
				T = 0.0/0.0;
			}
			else
			{
				T = 0;
			}
		}
	}
	return T;
}

double Get_SFRD_EoR_MINI(struct TsBox *previous_spin_temp, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, double redshift)
{
	// This is abit buggy cause we are actually using the mturn info from the previous box
	// However at lowz the z timestep is small enough so that prev box gives good enough results
	
	double Phi_EoR, H, SFRD, mturn, mc, Mlim_Fstar_MINI;
	if (redshift > find_redshift_step(2) - 0.5)
	{
		mturn = 1.0E20;
	}
	else
	{
		if (fabs(previous_spin_temp->IonBox_cache[3] - 1.0) > 1E-2)
		{
			printf("mturn is unset, setting now to inf\n");
			mturn = 1.0E20;
		}
		else
		{
			mturn = previous_spin_temp->IonBox_cache[1];
		}
	}
	mc = atomic_cooling_threshold(redshift);
	Mlim_Fstar_MINI = Mass_limit_bisection(global_params.M_MIN_INTEGRAL, global_params.M_MAX_INTEGRAL, astro_params->ALPHA_STAR_MINI, astro_params->F_STAR7_MINI * pow(1e3, astro_params->ALPHA_STAR_MINI));
	if (mturn > 1.0E18)
	{
		Phi_EoR = 0.0;
	}
	else
	{
		Phi_EoR = Nion_General_MINI(redshift, global_params.M_MIN_INTEGRAL, mturn, mc, astro_params->ALPHA_STAR_MINI, 0., astro_params->F_STAR7_MINI, 1., Mlim_Fstar_MINI, 0.);
	}
	Phi_EoR = Phi_EoR / (astro_params->t_STAR * pow(1. + redshift, astro_params->X_RAY_SPEC_INDEX + 1.0));
	H = hubble(redshift);
	SFRD = Phi_2_SFRD(Phi_EoR, redshift, H, astro_params, cosmo_params, 1);
	return SFRD;
}

/*---- Soft Photon Heating Module ----
To test this module, copy entire following codes to another file, define SOFT_PHOTON_TEST_MODE and then include the file
Some functions are already defined in RadioExcess.h, though in 21cmFAST this file has to be wrtten inside RadioExcess.h to avoid repeated includes
*/

#define xax_len 5000
#ifdef SOFT_PHOTON_TEST_MODE
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "RadioExcess.h"
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
    alpha = 1.0/137.03599976;
    
    // Gaunt factor following 10.1093/mnras/stad1539 (or 10.1051/0004-6361/201424643), mathematically identical to Chluba's paper
    vGHz = x * kB * 2.728 * (1.0+z)/(h * 1.0E9);
    T4 = Tk/1.0E4;
    gaunt = 5.9598540023339722 - sqrt(3.0)/M_PI * log(vGHz * pow(T4, -1.5));
    gaunt = log(exp(gaunt) + exp(1.0));

    mec2 = 510998.9461 * 1.602176634E-19;
    theta = Tk*kB/mec2;
    Nb = 0.25186318588882056 * OmBh2/0.02242 * pow(1.0+z, 3.0); // rho_cr*OmB/mp, number density of protons assuming baryons are made entirely of Hydrogen
    Np = xe * (1.0 - YHe) * Nb;
    lambda_e = h*c/mec2;
    r = alpha * pow(lambda_e, 3.0) / (2.0 * M_PI * sqrt(6.0*M_PI)) * pow(theta, -3.5) * Np * gaunt;
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
    
    Tcmb = 2.728 * (1.0+z);
    fHe = YHe/(4.0* (1.0 - YHe));
    nH0 = 0.1901567053460595 * OmBh2/0.02242; // Number density of H neclei today, in m^-3
    ne = nH0 * (1.0+fHe) * xe * pow(1.0+z, 3.0); // electron number density
    f = kB * Tcmb * x / h; // Photon frequency in Hz
    sII = fR_II * pow(f/0.15E9, -aR_II) * SFRD_II;
    sIII = fR_III * pow(f/0.15E9, -aR_III) * SFRD_III;
    src = sII + sIII;
    r1 = 5.507E11 * src / pow(f * (1.0+z), 3.0);
    r2 = 4.0 * n_inj * H;
    r = (r1 + r2) / (sigmaT * ne * c);
    return r;
    // return 0.0;
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
    Tcmb = 2.728 * (1.0+z);
    ex = exp(x);
    Y = x * ex / pow(ex - 1.0, 2.0) * (x * (ex+1.0)/(ex-1.0) - 4.0);
    Smid = kB / mec2 * (Tk - Tcmb) * Y;

    // Sff term, Eq(9)
    // LBR = Lambda_BR(z, xe, x, Tk, OmBh2, YHe);
    xi = x*Tcmb/Tk;
    exi = exp(xi);
    f1 = LBR * (1.0-1.0/exi)/pow(xi, 3.0);
    f2 = 1.0/(exi - 1.0) - 1.0/(ex - 1.0);
    Sff = f1 * f2;

    //DeltaS
    DeltaS = Sinj + Smid + Sff;
    
    r = DeltaS * pow(xi, 3.0) / (LBR * (1.0 - 1.0/exi));
    return r;
}

double dtauff_dz(double z, double H, double xe, double LBR, double xi, double OmBh2, double YHe)
{
    //\frac{d \tau_{ff}}{dz} = - \frac{\sigma_T N_e c \Lambda_{BR} (1 - e^{-x_i})}{(1+z) H x_i^3}
    double sigmaT, ne, c, nH0, fHe, r;
    c = 2.99792458E8;
    sigmaT = 0.665245854E-28;
    nH0 = 0.1901567053460595 * OmBh2/0.02242; // Number density of H neclei today, in m^-3
    fHe = YHe/(4.0* (1.0 - YHe));
    ne = nH0 * xe * (1.0+fHe) * pow(1.0+z, 3.0); // electron number density
    r = - sigmaT * ne * c * LBR * (1 - exp(-xi)) / ((1.0+z) * H * pow(xi, 3.0));
    return r;
}

double Build_dTffdz_Kernel(double z, double xe, double Tk, double H, double *x, double *LBR, double *DeltaN, int nx)
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
    nx: size of pointers x, LBR, DeltaN
    */
    double rho_cmb, Tcmb, kB, hb, c, dtdz, ne_ntot, Prefix, sigmaT, Integrand[xax_len], xi, ex, exi, dTkdz;
    int idx;

    if (fabs((double)xax_len - (double)nx) > 1E-10)
    {
       fprintf(stderr, "Error: unexpected x length nx\n");
	    exit(1);
    }

    kB = 1.38064852E-23;
    hb = 6.626070040818181818E-34/(2.0*M_PI);
    c = 2.99792458E8;
    sigmaT = 0.665245854E-28;

    Tcmb = 2.728*(1.0+z);
    rho_cmb = pow(M_PI, 2.0) * pow(kB*Tcmb, 4.0)/(15.0*pow(hb*c, 3.0));
    dtdz = -1/((1+z)*H);
    ne_ntot = xe/(1.0+xe); // This is ne/(nH(1+xe+fHe+xe*fHe))
    Prefix = 2 * dtdz*ne_ntot/(3*kB);
    Prefix = Prefix * 15 * rho_cmb * sigmaT * c * pow(Tk/Tcmb, 3.0)/pow(M_PI, 4.0); // Part not dependent on x integrand

    // Defining the integrand
    for (idx=0; idx<nx; idx++)
    {
        ex = exp(x[idx]);
        xi = x[idx] * Tcmb/Tk;
        exi = exp(xi);
        Integrand[idx] = LBR[idx] * (1.0 - 1.0/exi) * (DeltaN[idx] + 1.0/(ex - 1.0) - 1.0/(exi - 1));
    }
    dTkdz = Prefix * Integrate(x, Integrand, nx, 1);
    return dTkdz;
}

double Compute_dTffdz(double *zax, double *dTdz, double *Hax, double *xe_ax, double *Tkax, double *SFRD_II_ax, double *SFRD_III_ax, double fR_II, double fR_III, double aR_II, double aR_III, double OmBh2, int nz)
{
    /*
    Compute Soft-Photon heating rate dT/dz for an array of z, many inputs must be passed externally, such that the code can be easily tested outside of 21cmFAST (e.g., with python) and it can 
    merge streamelessly to 21cmFAST. Code computes dT/dz at all redshifts, which can be used by e.g. python for testing. Time used for Python to initialize input arrays seems subdominant compared 
    to the min solver
    -- inputs --
    zax: z axis
    dTdz: axis for dT/dz output
    */
    double xax[xax_len], LBR_ax[xax_len], DeltaN_ax[xax_len], DST_ax[xax_len], n_inj[xax_len];
    double z, xe, Tk, dz, dtff, xi, Tcmb, SFRD_II, SFRD_III, exp_dtff_inv, H, dninj_dtff, S_inj, YHe;
    int xid, zid;

    // Setting params
    YHe = 0.245;

    // Initializing
    logspace(-8.0, 2.0, xax, xax_len);
    z = zax[0];
    if (zax[0] < zax[nz-1])
    {
        fprintf(stderr, "Error: z needs to be descending\n");
	    exit(1);
    }
    dTdz[0] = 0.0;
    for (xid=0; xid<xax_len; xid++)
    {
        DeltaN_ax[xid] = 0.0;
        n_inj[xid] = 0.0;
    }
    
    // Start Evolving
    for (zid=1; zid<nz; zid++)
    {
        z = zax[zid];
        dz = zax[zid] - zax[zid-1];
        H = Hax[zid];
        xe = xe_ax[zid];
        Tk = Tkax[zid];
        Tcmb = 2.728 * (1.0+z);
        SFRD_II = SFRD_II_ax[zid];
        SFRD_III = SFRD_III_ax[zid];
        for (xid=0; xid<xax_len; xid++)
        {
            xi = xax[xid]*Tcmb/Tk;
            LBR_ax[xid] = Lambda_BR(z, xe, xax[xid], Tk, OmBh2, YHe);
            dtff = dtauff_dz(z, H, xe, LBR_ax[xid], xi, OmBh2, YHe)*dz;
            
            // Updating n_inj, not very efficient because this needs to call SoftPhoton_inj which has already been called in DeltaS_tilte, can be optimized if this adds too much time
            S_inj = SoftPhoton_inj(xax[xid], xe, z, SFRD_II, SFRD_III, fR_II, fR_III, aR_II, aR_III, n_inj[xid], H, YHe, OmBh2);
            dninj_dtff = S_inj * pow(xi, 3.0)/(LBR_ax[xid]*(1.0 - exp(-xi)));
            n_inj[xid] = n_inj[xid] + dninj_dtff * dtff;
            // Finished updating n_inj

            exp_dtff_inv = exp(-dtff);
            DST_ax[xid] = DeltaS_tilte(xax[xid], xe, z, SFRD_II, SFRD_III, fR_II, fR_III, aR_II, aR_III, n_inj[xid], H, YHe, OmBh2, Tk, LBR_ax[xid]);
            DeltaN_ax[xid] = DeltaN_ax[xid]*exp_dtff_inv + DST_ax[xid]*(1.0 - exp_dtff_inv);
            // DeltaN_ax[xid] = DeltaN_ax[xid]*(exp_dtff_inv) + DST_ax[xid]*dtff;
        }
        dTdz[zid] = Build_dTffdz_Kernel(z, xe, Tk, H, xax, LBR_ax, DeltaN_ax, xax_len);
    }

    return dTdz[nz-1];
}

double Find_dTff_dz(struct TsBox *previous_spin_temp, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, struct FlagOptions *flag_options)
{
	// TODO:
	// 1 - do Pop II radio
	// 2 - need to do last step?
	double zax[zax_len_FF], dTdz_ax[zax_len_FF], Hax[zax_len_FF], xe_ax[zax_len_FF], Tk_ax[zax_len_FF], SFRD_II_ax[zax_len_FF], SFRD_III_ax[zax_len_FF];
	double OmBh2, YHe, z1, z2, z, r;
	int idx, ArchiveSize, head;
	if (!flag_options->USE_RADIO_HEATING)
	{
		return 0.0;
	}
	else
	{
		// Do some checks, can be done outside in python
		if (!flag_options->USE_RADIO_MCG)
		{
			fprintf(stderr, "Error: You must set USE_RADIO_MCG to True if USE_RADIO_HEATING.\n");
            Throw(ValueError);
		}
		
	}	

	OmBh2 = cosmo_params->OMb * (pow(cosmo_params->hlittle, 2));
	YHe = 0.245;
	
	ArchiveSize = (int)round(previous_spin_temp->History_box[0]);
	if (ArchiveSize < 3)
	{
		return 0.0;
	}
	head = 3 * History_box_DIM + 1;
	z1 = previous_spin_temp->History_box[head];
	z2 = previous_spin_temp->History_box[(ArchiveSize - 1) * History_box_DIM + 1]+1.0E-3;

	logspace(log10(1.0+z1), log10(1.0+z2), zax, zax_len_FF);
	for (idx=0; idx<zax_len_FF; idx++)
	{
		zax[idx] = zax[idx] - 1; // logspace gives 1+z
		z = zax[idx];
		Hax[idx] = hubble(z);
		xe_ax[idx] = 1.0 - History_box_Interp(previous_spin_temp, z, 7, 1);
		Tk_ax[idx] = History_box_Interp(previous_spin_temp, z, 3, 1);
		SFRD_II_ax[idx] = 0.0;
		SFRD_III_ax[idx] = History_box_Interp(previous_spin_temp, z, 6, 1);
	}
	
	r = Compute_dTffdz(zax, dTdz_ax, Hax, xe_ax, Tk_ax, SFRD_II_ax, SFRD_III_ax, astro_params->fR, astro_params->fR_mini, astro_params->aR, astro_params->aR_mini, OmBh2, zax_len_FF);
	
	return r;

}
