// Things needed for Radio excess

// nu0 is degenerate with fR so no reason to leave this as a param
#define astro_nu0 0.15	   // in GHz
#define History_box_DIM 20 // number of quantities to be saved in History_box
#define zax_len_FF 1000	   // axis length for soft-photon heating
#define Clump_Factor_nm 1000
#define Clump_Factor_nx 1000
#define print_debug_info 0
#include "HaloProfile.c"

int Check_Astro_Call_Status(struct TsBox *spin_temp, int Check_Spin)
{
	int idx, status;
	idx = Check_Spin ? 2 : 3;
	if (fabs(spin_temp->IonBox_cache[idx] - 1994.0) < 0.001)
	{
		status = 0; // Astro not called
	}
	else if (fabs(spin_temp->IonBox_cache[idx] - 2026.0) < 0.001)
	{
		status = 1; // Astro called
	}
	else
	{
		status = -1; // Unknown
	}
	if (status == -1)
	{
		fprintf(stderr, "Cannot determine run status from IonBox_cache, v = %.5E, Check_Spin = %d.\n", spin_temp->IonBox_cache[idx], Check_Spin);
		Throw(ValueError);
	}
	return status;
}

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
			Throw(ValueError);
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
			Throw(ValueError);
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
			Throw(ValueError);
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
				Throw(ValueError);
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
				Throw(ValueError);
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
	{// Too early to get reliable results?
		if ((Type == 1 || Type == 2) || Type == 6)
		{ // Set Phi_II, Phi_III, SFRD_EoR_MINI to 0
			return 0.0;
		}
		else if (Type == 3)
		{
			fprintf(stderr, "Error: Too early to get reliable Tk from HistoryBox\n");
			Throw(ValueError);
		}
		else if (Type == 4 || Type == 5)
		{// Give large Mturns so that SFRD is 0 at high z
			return 1.0E20;
		}
		else
		{
			fprintf(stderr, "Error in History_box_Interp: Exception not set for Type = %d.\n", Type);
		}
	}

	if (ArchiveSize > 980)
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
		// Ensure xH is between [0, 1] while leaving for some redundencies for possible numerical precision issues
		if (r > 1.0-1.0E-8)
		{
			if (r < 1.0+1.0E-4)
			{
				r = 1.0-1.0E-8;
			}
			else
			{
				fprintf(stderr, "Interpolated xH is >1\n");
				Throw(ValueError);
			}
		}
		if (r < 1.0E-10)
		{
			if (r > -1.0E-15)
			{
				r = 1.0E-15;
			}
			else
			{
				fprintf(stderr, "Interpolated xH is negative\n");
				Throw(ValueError);
			}
		}
	}
	else
	{
		r = Interp_1D(z, z_axis, f_axis, ArchiveSize, 0, 1, Overflow_Handle);
	}
	return r;
}

double Get_Radio_Temp_HMG(struct TsBox *previous_spin_temp, struct TsBox *this_spin_temp, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, struct FlagOptions *flag_options, double zpp_max, double redshift)
{

	/* 
	Find Radio Temp from sources in redshifts [zpp_max, Z_Heat_max]
	---- inputs ----
	zpp_max: maximum zpp
	redshift: redshift at which you want to compute radio temp
	*/

	double z1, z2, dz, Phi, Phi_mini, z, fun_ACG, fun_MCG, Radio_Temp, Radio_Prefix_ACG, Radio_Prefix_MCG;
	int nz, zid, RadioSilent, AstroCalled, ArchiveSize;
	AstroCalled = Check_Astro_Call_Status(previous_spin_temp, 1);
	ArchiveSize = (int)round(previous_spin_temp->History_box[0]);
	
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

	if (((RadioSilent || (AstroCalled==0)) || this_spin_temp->first_box) || ArchiveSize < 3)
	{
		// Don't proceed if: radio excess is off, astro not called or this is first box, don't have enough samples in HistoryBox
		Radio_Temp = 0.0;
	}
	else
	{
		// Start from the first redshift in History_box, abort if it's below zpp_max
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
	x = x1 - ((float)idx - 1.0) * dx;
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
		// Need to do this more carefully, find_redshift_step does not work when USE_MANY_Z_xxx = T
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
				T = 0.0 / 0.0;
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
		if (fabs(previous_spin_temp->IonBox_cache[3] - 2026.0) > 1E-4)
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

#include "SoftPhoton.h"

void Print_debug_info_HistoryBox(struct TsBox *this_spin_temp)
{
	FILE *OutputFile;
	int ArchiveSize, head, idx;
	// Yell to ensure that the user does not forget this, e.g. when running mcmc
	printf("------------------------------------------------ print_debug_info activated------------------------------------------------\n");
	OutputFile = fopen("/Users/cangtao/Desktop/tmp_History_box.txt", "w");
	ArchiveSize = (int)round(this_spin_temp->History_box[0]);
	fprintf(OutputFile, "   z        Tk          xH        SFRD3		Mturn_III\n");
	for (idx = 1; idx <= ArchiveSize; idx++)
	{
		head = (idx - 1) * History_box_DIM + 1;
		fprintf(OutputFile, "%.3f   ", this_spin_temp->History_box[head]);
		fprintf(OutputFile, "%.3E   ", this_spin_temp->History_box[head + 2]);
		fprintf(OutputFile, "%.3E   ", this_spin_temp->History_box[head + 8]);
		fprintf(OutputFile, "%.3E	", this_spin_temp->History_box[head + 7]);
		fprintf(OutputFile, "%.3E\n", this_spin_temp->History_box[head + 6]);
	}
	fclose(OutputFile);
}

double Find_dTff_dz(struct TsBox *previous_spin_temp, struct AstroParams *astro_params, struct CosmoParams *cosmo_params, struct FlagOptions *flag_options, double *dT_Radio, double redshift)
{
	/* Limitations:
	1 - Not doing Pop II radio
	2 - Need to do last step?
	*/
	double zax[zax_len_FF], dTdz_ax[zax_len_FF], Hax[zax_len_FF], xe_ax[zax_len_FF], Tk_ax[zax_len_FF], SFRD_II_ax[zax_len_FF], SFRD_III_ax[zax_len_FF];
	double OmBh2, YHe, z1, z2, z, r;
	FILE *OutputFile;
	int idx, ArchiveSize, head;
	if (!flag_options->USE_RADIO_HEATING)
	{
		*dT_Radio = 0.0;
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
	
	// Check whether astro has been called in Spin.c previously
	if (Check_Astro_Call_Status(previous_spin_temp,1)==0)
	{
		*dT_Radio = 0.0;
		return 0.0;
	}

	OmBh2 = cosmo_params->OMb * (pow(cosmo_params->hlittle, 2));
	YHe = 0.245;

	ArchiveSize = (int)round(previous_spin_temp->History_box[0]);
	if ((ArchiveSize > 1) && print_debug_info)
	{
		Print_debug_info_HistoryBox(previous_spin_temp);
	}

	if (ArchiveSize < 3)
	{
		*dT_Radio = 0.0;
		return 0.0;
	}
	z1 = previous_spin_temp->History_box[1];
	z2 = previous_spin_temp->History_box[(ArchiveSize - 1) * History_box_DIM + 1] + 1.0E-3;
	if (print_debug_info)
	{
		printf("Running Radio heating, ArchiveSize = %d, z1 = %7f, z2 = %7f\n", ArchiveSize, z1, z2);
		OutputFile = fopen("/Users/cangtao/Desktop/tmp_Radio_Heating.txt", "w");
		fprintf(OutputFile, "   z       Tk        xe       SFRD3\n");
	}

	logspace(log10(1.0 + z1), log10(1.0 + z2), zax, zax_len_FF);

	for (idx = 0; idx < zax_len_FF; idx++)
	{
		zax[idx] = zax[idx] - 1; // logspace gives 1+z
		z = zax[idx];
		Hax[idx] = hubble(z);
		xe_ax[idx] = 1.0 - History_box_Interp(previous_spin_temp, z, 7, 1);
		Tk_ax[idx] = History_box_Interp(previous_spin_temp, z, 3, 1);
		SFRD_II_ax[idx] = 0.0;
		SFRD_III_ax[idx] = History_box_Interp(previous_spin_temp, z, 6, 1);
		if (print_debug_info)
		{
			fprintf(OutputFile, "%7f	%7E    %7E    %7E\n", z, Tk_ax[idx], xe_ax[idx], SFRD_III_ax[idx]);
		}
	}
	if (print_debug_info)
	{
		fclose(OutputFile);
	}

	r = Compute_dTffdz(zax, dTdz_ax, Hax, xe_ax, Tk_ax, SFRD_II_ax, SFRD_III_ax, dT_Radio, astro_params->fR, astro_params->fR_mini, astro_params->aR, astro_params->aR_mini, OmBh2, YHe, redshift, 1.0E10, zax_len_FF);
	if (print_debug_info)
	{
		OutputFile = fopen("/Users/cangtao/Desktop/tmp_dT_Radio.txt", "a");
		fprintf(OutputFile, "%7f	%7E\n", redshift, *dT_Radio);
		fclose(OutputFile);
	}
	return r;
}

double Collisional_Ionization_SigmaV(double T)
{
	/*
	Collisional ionization <sigma v> in m^3/s
	Follows Scholz91_APJ & Voronov97
	Scholz91_APJ fig1 data for log10(T) and Gamma axis are taken from their Fig.1, but their y axis was mislabeled: the unit
	should be 1E-9 cm^3/s rather than 1E-19cm^3/s, otherwise the figure does not match their Table 1 or VORONOV results!
	*/
	double a0, a1, a2, a3, a4, a5, a6, kB, r, y, Q, Gamma, U, A, P, X, K;

	a0 = -9.61443E1;
	a1 = 3.79523E1;
	a2 = -7.96885;
	a3 = 8.83922E-1;
	a4 = -5.34513E-2;
	a5 = 1.66344E-3;
	a6 = -2.08888E-5;
	kB = 1.38064852E-23;
	Q = 1.602176634E-19;

	if (T < 2000.0)
	{
		r = 0.0;
	}
	else if (T <= 1.0E8)
	{
		// Scholz91 fit
		y = log(T);
		Gamma = a0 + a1 * y + a2 * pow(y, 2.0) + a3 * pow(y, 3.0) + a4 * pow(y, 4.0) + a5 * pow(y, 5.0) + a6 * pow(y, 6.0);
		Gamma = exp(Gamma);
		r = Gamma * exp(-13.6 * Q / (kB * T)) * 1E-6;
	}
	else if (kB * T / Q <= 2.0E4)
	{
		// Voronov97 fit, valid for [2eV, 20KeV]
		U = 13.6 * Q / (T * kB);
		P = 0.0;
		A = 0.291E-7;
		X = 0.232;
		K = 0.39;
		r = A * (1 + P * pow(U, 0.5)) * pow(U, K) * exp(-U) / (X + U) * 1E-6;
	}
	else
	{
		// Return NaN&exit at higher energy, at such high T CI would kill 21cm signal
		fprintf(stderr, "T is too high to get reliable SigmaV, T = %.3E\n", T);
		Throw(ValueError);
		r = NAN;
	}
	return r;
}

double Find_dxe_dz_Collisional(double z, double xe, double T, double H, double delta, double ClumpingFactor)
{
	/*
	Just a quick test, many things can be optimized, e.g., repeated calculations of C
	*/
	double r, sigmav, nH0, nb, fHe, dtdz;
	nH0 = 0.1901567053460595;
	fHe = 0.08112582781456953;

	sigmav = Collisional_Ionization_SigmaV(T);
	nb = nH0 * pow(1.0 + z, 3.0) * (1.0 + fHe) * (1.0 + delta);
	dtdz = -1.0 / (H * (1.0 + z));
	r = ClumpingFactor * nb * xe * sigmav * (1.0 - xe) * dtdz / (1.0 + fHe);
	return r;
}

double Find_Clumping_Factor(double z, double T, double HaloTab_Mmin, struct UserParams *user_params, double growthf, struct CosmoParams *cosmo_params)
{
	/*
	TODO:
	check F(cx) limit at x->0
	*/
	double m_ax[Clump_Factor_nm], m_ax_SI[Clump_Factor_nm], x_ax[Clump_Factor_nx], Mmin, m, dndm, RhoSQ, Rho2_Integrand[Clump_Factor_nm], Fcoll_Integrand[Clump_Factor_nm];
	double msun, Mpc, fcoll, RhoSQ_Halo, RhoB_ave, result, RhoCr;
	int idx;
	double Rvir, rho_ax[Clump_Factor_nx]; // dummy variables, we don't really need them here
	Mpc = 3.086E22;
	msun = 1.98847E30;
	RhoCr = 1.879E-26 * pow(cosmo_params->hlittle, 2.0);

	Mmin = 1.3E3 * pow(10.0 * T / (1.0 + z), 1.5); // Above this mass baryons won't varialize and form halos
	if (HaloTab_Mmin>Mmin && !user_params->EVOLVE_MATTER)
	{// if Mmin is too low, it's likely that gas is too cold, then collisional ionization won't matter anyway
		Mmin = HaloTab_Mmin;
	}
	// This is the usual interpolation lower limit for sigma
	if (Mmin < 1.1E4)
	{
		Mmin = 1.1E4;
	}
	
	logspace(log10(Mmin), 19.0, m_ax, Clump_Factor_nm);
	logspace(-4.0, 0.0, x_ax, Clump_Factor_nx);
	
	for (idx = 0; idx < Clump_Factor_nm; idx++)
	{
		m = m_ax[idx];
		m_ax_SI[idx] = m*msun;
		// Get HMF
		if (user_params->HMF == 0)
		{
			dndm = dNdM(growthf, m, z);
		}
		else if (user_params->HMF == 1)
		{
			dndm = dNdM_st(growthf, m, z);
		}
		else if (user_params->HMF == 2)
		{
			dndm = dNdM_WatsonFOF(growthf, m, z);
		}
		else if (user_params->HMF == 3)
		{
			dndm = dNdM_WatsonFOF_z(z, growthf, m);
		}
		else
		{
			LOG_ERROR("Wrong choice of hmf_model!");
			Throw(ValueError);
		}
		Halo_Baryon_Profile(z, m, &Rvir, &RhoSQ, x_ax, rho_ax, Clump_Factor_nx, cosmo_params->hlittle, cosmo_params->OMm, cosmo_params->OMb);
		Rho2_Integrand[idx] = dndm * RhoSQ /(msun* pow(Mpc, 3.0));
		Fcoll_Integrand[idx] = m*dndm /pow(Mpc, 3.0);
	}
	
	fcoll = Integrate(m_ax_SI, Fcoll_Integrand, Clump_Factor_nm, 1);
	RhoSQ_Halo = Integrate(m_ax_SI, Rho2_Integrand, Clump_Factor_nm, 1);
	RhoB_ave = cosmo_params->OMb * RhoCr * pow(1.0 + z, 3.0);
	result = pow(1.0 + fcoll, 2.0) + pow(1.0+z, 3.0) *RhoSQ_Halo / pow(RhoB_ave, 2.0);
	return result;
}

void Print_HMF(double z, int hmf_model)
{
	int nm = 200;
	int idx;
	double m_ax[200], dndm, m, growthf, Mmin, mc;
	FILE *OutputFile;

	printf("-------- Printing HMF ----\n");
	Mmin = 1.001E4;
	logspace(log10(Mmin), 19.0, m_ax, nm);
	growthf = dicke(z);
	OutputFile = fopen("/Users/cangtao/FileVault/Projects/SDM/data/HMF/HMF_Table.txt", "a");
	fprintf(OutputFile, "%.4E  ", z);
	
	for (idx = 0; idx < nm; idx++)
	{
		m = m_ax[idx];
		if (hmf_model == 0)
		{
			dndm = dNdM(growthf, m, z);
		}
		else if (hmf_model == 1)
		{
			dndm = dNdM_st(growthf, m, z);
		}
		else if (hmf_model == 2)
		{
			dndm = dNdM_WatsonFOF(growthf, m, z);
		}
		else if (hmf_model == 3)
		{
			dndm = dNdM_WatsonFOF_z(z, growthf, m);
		}
		else
		{
			LOG_ERROR("Wrong choice of hmf_model!");
			Throw(ValueError);
		}
		fprintf(OutputFile, "%.4E  ", dndm);
	}
	fprintf(OutputFile, "\n");
	fclose(OutputFile);

	mc = atomic_cooling_threshold(z);
	OutputFile = fopen("/Users/cangtao/FileVault/Projects/SDM/data/HMF/mc.txt", "a");
	fprintf(OutputFile, "%.4f  %.4E\n", z, mc);
	fclose(OutputFile);	
}
