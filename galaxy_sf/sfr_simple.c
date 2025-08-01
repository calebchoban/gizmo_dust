#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "../allvars.h"
#include "../proto.h"

/*!
 *  routines for star formation in cosmological simulations
 */
/*
 * This file is largely written by Phil Hopkins (phopkins@caltech.edu) for GIZMO.
 *   It was based on a similar file in GADGET3 by Volker Springel (volker.springel@h-its.org),
 *   but the physical modules for star formation and feedback have been
 *   replaced, and the algorithm is mostly new to GIZMO.
 */


#if defined(GALSF) && (defined(STOCHASTIC_IMF) || defined(SAMPLE_IMF_FROM_GAS))

#ifdef OUTPUT_STARFORMATION_INFO
void save_sf_info_to_file(){
  int i, j;
  int nstars_per_proc = 0;  //for output purpose

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type==807)
      nstars_per_proc++;
  
  //Begin the painful IO stuff...
  // total_nstars = 2+5+0+1 = 8
  int total_nstars;
  MPI_Allreduce(&nstars_per_proc, &total_nstars, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(total_nstars > 0)
    {

      //Each proc has its only array
      MyIDType *id_star = (MyIDType*) mymalloc("id_star",  nstars_per_proc * sizeof(MyIDType));
      double *x_sf = (double*) mymalloc("x_sf",  nstars_per_proc * sizeof(double));
      double *y_sf = (double*) mymalloc("y_sf",  nstars_per_proc * sizeof(double));
      double *z_sf = (double*) mymalloc("z_sf",  nstars_per_proc * sizeof(double));
      double *u_sf = (double*) mymalloc("u_sf",  nstars_per_proc * sizeof(double));
      double *rho_sf = (double*) mymalloc("rho_sf",  nstars_per_proc * sizeof(double));
      double *m_sf = (double*) mymalloc("m_sf",  nstars_per_proc * sizeof(double));
      int ii=0;
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])                                                   
	if(P[i].Type==807)
	  {
	    P[i].Type = 4;
	    id_star[ii] = P[i].ID;
	    x_sf[ii] = P[i].Pos[0];
	    y_sf[ii] = P[i].Pos[1];
	    z_sf[ii] = P[i].Pos[2];
	    u_sf[ii] = SphP[i].InternalEnergy;
	    rho_sf[ii] = SphP[i].Density;
#ifdef STOCHASTIC_IMF
	    m_sf[ii] = P[i].MassMassiveStar;
#else
	    m_sf[ii] = 0.;
#endif
	    ii++;
	  }
      
      /* Only root has the received data */
      int *nstars_per_proc_global_list = NULL;
      if (ThisTask == 0)
	nstars_per_proc_global_list = (int*) mymalloc("nstars_per_proc_global_list", NTask * sizeof(int)) ;
      
      // nstars_per_proc_global_list holds info for how many SN occur in each proc (e.g. [2, 5, 0, 1] means 2 SNe in proc=0, 5 SNe in proc=1, 0 SNe in proc=2, 1 SNe in proc=3)
      MPI_Gather(&nstars_per_proc, 1, MPI_INT, nstars_per_proc_global_list, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      
      int* displs = NULL;
      MyIDType* global_id_star = NULL;
      double* global_x_sf = NULL;
      double* global_y_sf = NULL;
      double* global_z_sf = NULL;
      double* global_u_sf = NULL;
      double* global_rho_sf = NULL;
      double* global_m_sf = NULL;
      //nstars_per_proc_global_list & displs are only allocated in the root node 
      if (ThisTask == 0)
	{
	  //displs is the displacement array for MPI_Gatherv
	  displs = (int*) mymalloc("displs", NTask * sizeof(int) );
	  for(i=1, displs[0]=0; i<NTask; i++)
	    displs[i] = displs[i-1] + nstars_per_proc_global_list[i-1];
	  
	  global_id_star = (MyIDType*) mymalloc("global_id_star", total_nstars * sizeof(MyIDType) );
	  global_x_sf = (double*) mymalloc("global_x_sf", total_nstars * sizeof(double) );
	  global_y_sf = (double*) mymalloc("global_y_sf", total_nstars * sizeof(double) );
	  global_z_sf = (double*) mymalloc("global_z_sf", total_nstars * sizeof(double) );
	  global_u_sf = (double*) mymalloc("global_u_sf", total_nstars * sizeof(double) );
	  global_rho_sf = (double*) mymalloc("global_rho_sf", total_nstars * sizeof(double) );
	  global_m_sf = (double*) mymalloc("global_m_sf", total_nstars * sizeof(double) );
	}
      
      /* Now we have the receive buffer, counts, and displacements, we are ready to gather the data to the root node */
      /* Don't put MPI_Gatherv inside if(ThisTask==0) !!! That would lead to deadlock... */
#ifdef LONGIDS
      MPI_Gatherv(id_star, nstars_per_proc, MPI_UNSIGNED_LONG_LONG, global_id_star, nstars_per_proc_global_list, displs, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#else
      MPI_Gatherv(id_star, nstars_per_proc, MPI_UNSIGNED, global_id_star, nstars_per_proc_global_list, displs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif
      MPI_Gatherv(x_sf, nstars_per_proc, MPI_DOUBLE, global_x_sf, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(y_sf, nstars_per_proc, MPI_DOUBLE, global_y_sf, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(z_sf, nstars_per_proc, MPI_DOUBLE, global_z_sf, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(u_sf, nstars_per_proc, MPI_DOUBLE, global_u_sf, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(rho_sf, nstars_per_proc, MPI_DOUBLE, global_rho_sf, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(m_sf, nstars_per_proc, MPI_DOUBLE, global_m_sf, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      //Finally we can do the IO here!!!
      if (ThisTask == 0)
	for(j=0; j<total_nstars; j++)
	  {
#ifdef LONGIDS
            printf("All.Time=%g, x_sf=%g, y_sf=%g, z_sf=%g, u_sf=%g, rho_sf=%g, id_star=%llu, m_sf=%g\n",
	      All.Time, global_x_sf[j], global_y_sf[j], global_z_sf[j], global_u_sf[j], global_rho_sf[j], global_id_star[j], global_m_sf[j]);
            fprintf(FdSFinfo, "%g %g %g %g %g %g %llu %g\n",
	      All.Time, global_x_sf[j], global_y_sf[j], global_z_sf[j], global_u_sf[j], global_rho_sf[j], global_id_star[j], global_m_sf[j]);
#else
            printf("All.Time=%g, x_sf=%g, y_sf=%g, z_sf=%g, u_sf=%g, rho_sf=%g, id_star=%u, m_sf=%g\n",
	      All.Time, global_x_sf[j], global_y_sf[j], global_z_sf[j], global_u_sf[j], global_rho_sf[j], global_id_star[j], global_m_sf[j]);
            fprintf(FdSFinfo, "%g %g %g %g %g %g %u %g\n",
	      All.Time, global_x_sf[j], global_y_sf[j], global_z_sf[j], global_u_sf[j], global_rho_sf[j], global_id_star[j], global_m_sf[j]);
#endif
	    fflush(FdSFinfo); // can flush it, because only occuring on master steps anyways
	  }
      
      
      //free the arrays
      if (ThisTask == 0)
	{
	  myfree(global_m_sf);
	  myfree(global_rho_sf);
	  myfree(global_u_sf);
	  myfree(global_z_sf);
	  myfree(global_y_sf);
	  myfree(global_x_sf);
	  myfree(global_id_star);
	  myfree(displs);
	  myfree(nstars_per_proc_global_list);
	}
      myfree(m_sf);
      myfree(rho_sf);
      myfree(u_sf);
      myfree(z_sf);
      myfree(y_sf);
      myfree(x_sf);
      myfree(id_star);
    }//total_nstars > 0
  //Done with the painful IO stuff...
}
#endif

#define WindInitialVelocityBoost 1.0 // (optional) boost velocity coupled (fixed momentum)

/* return the stellar age in Gyr for a given labeled age, needed throughout for stellar feedback */
double evaluate_stellar_age_Gyr_simple(double stellar_tform)
{
    double age,a0,a1,a2,x0,x1,x2;
    if(All.ComovingIntegrationOn)
    {
        a0 = stellar_tform;
        a2 = All.Time;
        if(fabs(1-(All.OmegaMatter+All.OmegaLambda))<=0.01)
        {
            /* use exact solution for flat universe */
            x0 = (All.OmegaMatter/(1-All.OmegaMatter))/(a0*a0*a0);
            x2 = (All.OmegaMatter/(1-All.OmegaMatter))/(a2*a2*a2);
            age = (2./(3.*sqrt(1-All.OmegaMatter)))*log(sqrt(x0*x2)/((sqrt(1+x2)-1)*(sqrt(1+x0)+1)));
            age *= 1./All.Hubble_H0_CodeUnits;
        } else {
            /* use simple trap rule integration */
            a1 = 0.5*(a0+a2);
            x0 = 1./(a0*hubble_function(a0));
            x1 = 1./(a1*hubble_function(a1));
            x2 = 1./(a2*hubble_function(a2));
            age = (a2-a0)*(x0+4.*x1+x2)/6.;
        }
    } else {
        /* time variable is simple time, when not in comoving coordinates */
        age=All.Time-stellar_tform;
    }
    age *= 0.001*UNIT_TIME_IN_MYR/All.HubbleParam; // convert to absolute Gyr
    //if((age<=1.e-5)||(isnan(age))) {age=1.e-5;}
    return age;
}


/* simple routine to determine density thresholds and other common units for SF routines */
void set_units_sfr(void)
{
    All.OverDensThresh = All.CritOverDensity * All.OmegaBaryon * 3 * All.Hubble_H0_CodeUnits * All.Hubble_H0_CodeUnits / (8 * M_PI * All.G);
    All.PhysDensThresh = All.CritPhysDensity * PROTONMASS / (HYDROGEN_MASSFRAC * UNIT_DENSITY_IN_CGS * All.HubbleParam*All.HubbleParam);
}


/* Routine to actually determine the SFR assigned to an individual gas particle at each time */
double get_starformation_rate_simple(int i)
{
  double rateOfSF,tsfr;

  tsfr = sqrt( 3.*M_PI / ( 32. * All.G * SphP[i].Density * All.cf_a3inv) );
  if(tsfr<=0) return 0;
    
  rateOfSF = All.SfEffPerFreeFall * P[i].Mass / tsfr;
    
  return rateOfSF;
}



/* master routine for star formation. for 'effective equation of state' models for star-forming gas, this also updates their effective EOS parameters */
void star_formation_parent_routine(void)
{
  int i, bin, flag, stars_spawned, tot_spawned, stars_converted, tot_converted, number_of_stars_generated;
  unsigned int bits;
  double dtime, mass_of_star, p, prob, rate_in_msunperyear, sfrrate, totsfrrate;
  double sum_sm, total_sm, sm=0, rate, sum_mass_stars, total_sum_mass_stars;
    
  for(bin = 0; bin < TIMEBINS; bin++) {if(TimeBinActive[bin]) {TimeBinSfr[bin] = 0;}}
  
  stars_spawned = stars_converted = 0; sum_sm = sum_mass_stars = 0;
  
  for(bits = 0; GALSF_GENERATIONS > (1 << bits); bits++);
  
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if((P[i].Type == 0)&&(P[i].Mass>0))
	{
	  SphP[i].Sfr = 0; flag = 1; /* will be reset below if flag==0, but default to flag = 1 (non-eligible) */
	  dtime = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval / All.cf_hubble_a; /*  the actual time-step */
	  
          double gamma = 5./3.;
          double gamma_minus1 = 5./3. -1.;
#ifdef JEANS_LENGTH_THRESHOLD
	  double soundspeed = sqrt(gamma * gamma_minus1 * SphP[i].InternalEnergy);
	  double L_J = 1.77245 * soundspeed / sqrt(All.G * SphP[i].Density * All.cf_a3inv);  //sqrt(pi) = 1.77245
	  if(L_J < All.SfThreshJeansLength)
	    flag = 0;
#elif defined JEANS_MASS_THRESHOLD
#ifdef JEANS_MASS_THRESHOLD_VTURB
	  double vgrad2 = ( pow(SphP[i].Gradients.Velocity[0][0],2) + pow(SphP[i].Gradients.Velocity[0][1],2) + pow(SphP[i].Gradients.Velocity[0][2],2) + 
			    pow(SphP[i].Gradients.Velocity[1][0],2) + pow(SphP[i].Gradients.Velocity[1][1],2) + pow(SphP[i].Gradients.Velocity[1][2],2) + 
			    pow(SphP[i].Gradients.Velocity[2][0],2) + pow(SphP[i].Gradients.Velocity[2][1],2) + pow(SphP[i].Gradients.Velocity[2][2],2) );
	  double vturb2 = 0.2 * pow(PPP[i].Hsml,2) * vgrad2;
	  double soundspeed = sqrt(gamma * gamma_minus1 * SphP[i].InternalEnergy + vturb2);
#else
	  double soundspeed = sqrt(gamma * gamma_minus1 * SphP[i].InternalEnergy);
#endif	  
	  double M_J = 2.9155697 * pow(soundspeed, 3) / pow(All.G, 1.5) / sqrt(SphP[i].Density * All.cf_a3inv);  // pi^2.5 / 6 = 2.9155697
	  double M_th = All.FacSfThreshMJ * All.DesNumNgb * P[i].Mass;  //code unit
	  if(M_J < M_th)
	    flag = 0;
#else
	  if(SphP[i].Density*All.cf_a3inv >= All.PhysDensThresh)
	    flag = 0;
#endif

	  
	  if(All.ComovingIntegrationOn) {if(SphP[i].Density < All.OverDensThresh) flag = 1;} // (additional density check for cosmological runs) //

#ifdef SF_ONLY_IN_CONV_FLOWS
	  if((SphP[i].Gradients.Velocity[0][0] >= 0.) || (SphP[i].Gradients.Velocity[1][1] >= 0.) || (SphP[i].Gradients.Velocity[2][2] >= 0.))
	    flag = 1; //don't form stars
#endif

#ifdef SF_INSTANT_CUTOFF
	  double nH = SphP[i].Density * UNIT_DENSITY_IN_CGS * HYDROGEN_MASSFRAC / PROTONMASS;
	  if(nH > All.nHcutoffSF)
	    flag = 0; //form stars
#endif

#ifdef PHOTO_IONIZATION
	  //only needed for density threshold; for Jeans threshold, the photo-heated gas is naturally not star forming
	  if(SphP[i].Ionized == 1)
	    flag = 1; //don't form stars
#endif

#ifdef DETERMINISTIC_SF
	  if(flag == 1){
	    SphP[i].TimeBeginSF = -1.; //reset to -1
	    SphP[i].TimeSF = 0.;
	  }
#endif

	  if((flag == 0)&&(dtime>0)&&(P[i].TimeBin))		/* active star formation (upon start-up, we need to protect against dt==0) */
	    {
	      /* the upper bits of the gas particle ID store how many stars this gas particle gas already generated */
	      if(bits == 0)
		number_of_stars_generated = 0;
	      else
		number_of_stars_generated = (P[i].ID >> (sizeof(MyIDType) * 8 - bits));
	      
#ifdef DETERMINISTIC_SF
	      //we don't need mass_of_star as we always assume GALSF_GENERATIONS==1
	      double rateOfSF,tff, tSF, age;
	      if(SphP[i].TimeBeginSF <= 0.) //newly crossing the SF threshold
		{
		  SphP[i].TimeBeginSF = All.Time;
		  tff = sqrt( 3.*M_PI / ( 32. * All.G * SphP[i].Density * All.cf_a3inv) ); //free-fall time
		  SphP[i].TimeSF = tff / All.SfEffPerFreeFall;
		}
	      age = All.Time - SphP[i].TimeBeginSF;
	      if( (age > SphP[i].TimeSF) && (SphP[i].TimeSF > 0.) )
#else
	      sm = get_starformation_rate_simple(i) * dtime; // expected stellar mass formed this timestep
	      // (this also updates entropies for the effective equation-of-state model) //
	      p = sm / P[i].Mass;
	      sum_sm += P[i].Mass * (1 - exp(-p));
	      
	      /* Alright, now we consider the actual gas-to-star particle conversion and associated steps */
	      mass_of_star = P[i].Mass / (GALSF_GENERATIONS - number_of_stars_generated);
	      if(number_of_stars_generated >= GALSF_GENERATIONS-1) mass_of_star=P[i].Mass;
	      
	      SphP[i].Sfr = sm / dtime *
		(All.UnitMass_in_g / SOLAR_MASS_CGS) / (UNIT_TIME_IN_CGS / SECONDS_PER_YEAR);
	      if(dtime>0) TimeBinSfr[P[i].TimeBin] += SphP[i].Sfr;
	      
	      prob = P[i].Mass / mass_of_star * (1 - exp(-p));
#ifdef SF_INSTANT_CUTOFF
	      if(nH > All.nHcutoffSF)
		prob = 1.0; //form stars instantaneously
#endif
	      //if(get_random_number(P[i].ID + 1) < prob)	/* ok, make a star */
	      if(gsl_rng_uniform(random_generator) < prob)	/* ok, make a star */
#endif
		{
		  
		  
		  /* ok, we're going to make a star! */
		  
		  if(number_of_stars_generated == (GALSF_GENERATIONS - 1))
		    {
		      /* here we turn the gas particle itself into a star */
		      Stars_converted++;
		      stars_converted++;
		      sum_mass_stars += P[i].Mass;
		      
		      P[i].Type = 4;
#ifdef OUTPUT_STARFORMATION_INFO
		      P[i].Type = 807;
#endif
		      TimeBinCountGas[P[i].TimeBin]--;
		      TimeBinSfr[P[i].TimeBin] -= SphP[i].Sfr;
		      
		      P[i].StellarAge = All.Time;
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
		      P[i].DensAroundStar = SphP[i].Density;
#endif
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
		      P[i].Mass = SphP[i].MassTrue + SphP[i].dMass;
#endif
		      
#ifdef STELLAR_FEEDBACK
		      double mass, prob, random, m_imf;
		      mass = P[i].Mass * All.UnitMass_in_g / SOLAR_MASS_CGS;
		      prob = mass / All.MassPerStarIMF;
		      //random = get_random_number(P[i].ID+2);
		      random = gsl_rng_uniform(random_generator);
		      //printf("prob=%g, random=%g, P[i].ID=%d\n", prob, random, P[i].ID);
		      if(random < prob)
			{
			  P[i].flag_SNII = 0;
#ifdef STOCHASTIC_IMF
			  double xmin = All.MinMassIMF;
			  double xmax = All.MaxMassIMF;
			  //double y = get_random_number(P[i].ID+3);
			  double y = gsl_rng_uniform(random_generator);
			  double nor = pow(xmax, -1.3) - pow(xmin, -1.3);
			  P[i].MassMassiveStar = pow( (pow(xmin, -1.3) + nor * y ), -1/1.3);
			  printf("sampling a star from IMF: random=%g,  y=%g,  MassMassiveStar=%g\n", random, y, P[i].MassMassiveStar);
			  if(P[i].MassMassiveStar >= 8.0)
			    P[i].flag_SNII = 0;
			  else
			    P[i].flag_SNII = -1;
#ifdef G0_VARIABLE			  
			  P[i].UV_luminosity = pow(10., get_logL_pe( P[i].MassMassiveStar ) );
#endif
#ifdef DUST_IN_AGB
			  if(P[i].MassMassiveStar < 8.0 && P[i].MassMassiveStar > 1.0)
			    {
			      P[i].MdustCarbonAGB   = 1e-10*pow(10., get_logMdustC_agb( P[i].MassMassiveStar ) ); //code unit
			      //P[i].MdustSilicateAGB = 1e-10*pow(10., get_logMdustSi_agb( P[i].MassMassiveStar ) ); //code unit
			      printf("DUST_IN_AGB... i=%d, MassMassiveStar=%g, MdustCarbonAGB=%g\n", i, P[i].MassMassiveStar, P[i].MdustCarbonAGB);
			    }
#endif
#endif //STOCHASTIC_IMF
			}
		      else
			{
			  P[i].flag_SNII = -1; //mark as exploded
#ifdef G0_VARIABLE
			  P[i].UV_luminosity = 0.0;
#endif
			}
#endif
		      
		    } /* closes final generation from original gas particle */
		  else
		    {
		      /* here we spawn a new star particle */
		      
		      if(NumPart + stars_spawned >= All.MaxPart)
			{
			  printf
			    ("On Task=%d with NumPart=%d we try to spawn %d particles. Sorry, no space left...(All.MaxPart=%d)\n",
			     ThisTask, NumPart, stars_spawned, All.MaxPart);
			  fflush(stdout);
			  endrun(8888);
			}
		      
		      P[NumPart + stars_spawned] = P[i];
		      P[NumPart + stars_spawned].Type = 4;
#ifdef DO_DENSITY_AROUND_STAR_PARTICLES
		      P[NumPart + stars_spawned].DensAroundStar = SphP[i].Density;
#endif
		      NextActiveParticle[NumPart + stars_spawned] = FirstActiveParticle;
		      FirstActiveParticle = NumPart + stars_spawned;
		      NumForceUpdate++;
		      
		      TimeBinCount[P[NumPart + stars_spawned].TimeBin]++;
		      PrevInTimeBin[NumPart + stars_spawned] = i;
		      NextInTimeBin[NumPart + stars_spawned] = NextInTimeBin[i];
		      if(NextInTimeBin[i] >= 0)
			PrevInTimeBin[NextInTimeBin[i]] = NumPart + stars_spawned;
		      NextInTimeBin[i] = NumPart + stars_spawned;
		      if(LastInTimeBin[P[i].TimeBin] == i)
			LastInTimeBin[P[i].TimeBin] = NumPart + stars_spawned;
		      
		      P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - bits));
		      
		      P[NumPart + stars_spawned].Mass = mass_of_star;
		      P[i].Mass -= P[NumPart + stars_spawned].Mass;
		      if(P[i].Mass<0) P[i].Mass=0;
#ifdef HYDRO_MESHLESS_FINITE_VOLUME
		      SphP[i].MassTrue -= P[NumPart + stars_spawned].Mass;
		      if(SphP[i].MassTrue<0) SphP[i].MassTrue=0;
#endif
		      sum_mass_stars += P[NumPart + stars_spawned].Mass;
		      P[NumPart + stars_spawned].StellarAge = All.Time;
		      force_add_star_to_tree(i, NumPart + stars_spawned);
		      
		      stars_spawned++;
		    }
		}
	    } // closes check of flag==0 for star-formation operation
	} /* End of If Type = 0 */
    } /* end of main loop over active particles, huzzah! */

#ifdef OUTPUT_STARFORMATION_INFO
  save_sf_info_to_file();
#endif

#ifdef G0_SCALE_WITH_TOTAL_SFR
  double sfr_sc = 0.; 
  double tot_sfr_sc;
  double age;
  double time_define_SFR = 3e7; //in years

  for(i=0; i<NumPart; i++)
    {
      if(P[i].Type == 4)
	{
          if(All.ComovingIntegrationOn)
            age = 0.;//TODO
          else
            age = All.Time - P[i].StellarAge; //in code units

	  age *= (UNIT_TIME_IN_CGS / SECONDS_PER_YEAR / All.HubbleParam); //in years
	  if(age < time_define_SFR)
	    sfr_sc += P[i].Mass;
	}
    }
  MPI_Allreduce(&sfr_sc, &tot_sfr_sc, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  tot_sfr_sc = tot_sfr_sc / time_define_SFR; //mass is code unit, but it's okay since All.Mgas0 is also in code unit
  double depletion_time_0 = 4e9; //in years, for solar neighborhood conditions
  double SFR_0 = All.Mgas0 / depletion_time_0; 
  All.FactorG0 = tot_sfr_sc / SFR_0;
  if(ThisTask == 0)
    printf("All.Mgas0 = %g   SFR_0 = %g   All.FactorG0 = %g\n", All.Mgas0, SFR_0, All.FactorG0);
#endif


  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  MPI_Allreduce(&stars_converted, &tot_converted, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0 || tot_converted > 0)
    {
      if(ThisTask == 0)
	{
	  //#ifndef IO_REDUCED_MODE
	  printf("SFR: spawned %d stars, converted %d gas particles into stars\n",
		 tot_spawned, tot_converted);
	  //#endif
    }
      All.TotNumPart += tot_spawned;
      All.TotN_gas -= tot_converted;
      NumPart += stars_spawned;
      /* Note: N_gas is only reduced once rearrange_particle_sequence is called */
      /* Note: New tree construction can be avoided because of  `force_add_star_to_tree()' */
    } //(tot_spawned > 0 || tot_converted > 0)


  for(bin = 0, sfrrate = 0; bin < TIMEBINS; bin++)
    if(TimeBinCount[bin])
      sfrrate += TimeBinSfr[bin];

#ifdef IO_REDUCED_MODE
    if(All.HighestActiveTimeBin == All.HighestOccupiedTimeBin)
#endif
    {
        MPI_Allreduce(&sfrrate, &totsfrrate, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        MPI_Reduce(&sum_sm, &total_sm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        MPI_Reduce(&sum_mass_stars, &total_sum_mass_stars, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
        if(ThisTask == 0)
        {
            if(All.TimeStep > 0)
                rate = total_sm / (All.TimeStep / (All.cf_atime*All.cf_hubble_a));
            else
                rate = 0;
            /* convert to solar masses per yr */
            rate_in_msunperyear = rate * (All.UnitMass_in_g / SOLAR_MASS_CGS) / (UNIT_TIME_IN_CGS / SECONDS_PER_YEAR);
            fprintf(FdSfr, "%g %g %g %g %g\n", All.Time, total_sm, totsfrrate, rate_in_msunperyear, total_sum_mass_stars);
            fflush(FdSfr); // can flush it, because only occuring on master steps anyways
        } // thistask==0
    }

    if(tot_converted+tot_spawned > 0) {rearrange_particle_sequence();}

    CPU_Step[CPU_COOLINGSFR] += measure_time();
} /* end of main sfr_cooling routine!!! */


#endif // GALSF


