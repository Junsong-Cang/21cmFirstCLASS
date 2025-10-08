// Things needed for Radio excess

// nu0 is degenerate with fR so no reason to leave this as a param
#define astro_nu0 0.15	   // in GHz
#define History_box_DIM 20 // number of quantities to be saved in History_box

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
		6 - Phi_III_EoR
	*/
	int ArchiveSize, idx, head, zid, fid;
	// Very generous with memory, nobody is gonna run lightcones with 1000 timesteps (?)
	double z_axis[1000], f_axis[1000], r;

	ArchiveSize = (int)round(previous_spin_temp->History_box[0]);
	if (previous_spin_temp->first_box || ArchiveSize < 2)
	{
		if ((Type == 1 || Type == 2) || Type == 6)
		{
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
		{ // Phi3_EoR
			fid = head + 7;
		}
		else
		{
			LOG_ERROR("Wrong Type setting, must be in [1, 6].\n");
			Throw(ValueError);
		}
		z_axis[idx] = previous_spin_temp->History_box[zid];
		f_axis[idx] = previous_spin_temp->History_box[fid];
	}
	// Use_LogY for all
	r = Interp_1D(z, z_axis, f_axis, ArchiveSize, 0, 1, Overflow_Handle);

	// actually let's not use mturn interp for SFRD

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
		if (fabs(previous_spin_temp->mturns_EoR[3] - 1.0) > 1E-2)
		{
			printf("mturn is unset, setting now to inf\n");
			mturn = 1.0E20;
		}
		else
		{
			mturn = previous_spin_temp->mturns_EoR[1];
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
