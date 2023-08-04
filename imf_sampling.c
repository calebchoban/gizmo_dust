#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include "allvars.h"
#include "proto.h"
//#include "forcetree.h"


#ifdef SAMPLE_IMF

static inline double drawMassFromIMF(double x)
{
  double A = 0.126512;
  double y = 0.0;
  
  if( x < 0.5)
    y = 2. * A * pow(x, -1.3);
  else if(x >= 0.5)
    y = A * pow(x, -2.3);
  return y;
}


static inline double envelope_function(double x)
{
  double A = 0.126512;
  return 2. * A * pow(x, -1.3);
}


const double M_max = 50.;
const double M_min = 0.08;

#ifdef METALSSSS
static void IMFsampling(double* Mass, double* Pos_x, double* Pos_y, double* Pos_z, double* MstarSampleIMF, double* ZMass, int* flag, int NumNewStars);
#else
static void IMFsampling(double* Mass, double* Pos_x, double* Pos_y, double* Pos_z, double* MstarSampleIMF, int* flag, int NumNewStars);
#endif

void assign_stellar_masses() {

  int i;
  int nstars_per_proc = 0;

  if(ThisTask == 0)
    printf("assigning stellar masses to star particles...\n");

  for(i = 0; i < NumPart; i++)
    if(P[i].Type == 4 && P[i].sampled == 0)
      nstars_per_proc++;

  double *mass = (double*) mymalloc("mass",  nstars_per_proc * sizeof(double));
  int *flag = (int*) mymalloc("flag",  nstars_per_proc * sizeof(int));
  double *xstars = (double*) mymalloc("xstars",  nstars_per_proc * sizeof(double));
  double *ystars = (double*) mymalloc("ystars",  nstars_per_proc * sizeof(double));
  double *zstars = (double*) mymalloc("zstars",  nstars_per_proc * sizeof(double));
  double* stellar_masses = (double*) mymalloc("stellar_masses",  nstars_per_proc*N_STELLAR_MASS * sizeof(double));
#ifdef METALSSSS
  double *Zmass = (double*) mymalloc("Zmass",  12 * nstars_per_proc * sizeof(double));
#endif

  int j=0;
  int ik;
  for(i = 0; i < NumPart; i++)
    {
      //if(P[i].Type == 4)
      if(P[i].Type == 4 && P[i].sampled == 0)
	{
	  mass[j] = P[i].Mass * All.UnitMass_in_g / SOLAR_MASS;
#ifdef METALSSSS
	  for(ik = 0; ik < 12; ik++)
	    Zmass[j*12+ik] = P[i].Zm[ik];
#endif
	  xstars[j] = P[i].Pos[0];
	  ystars[j] = P[i].Pos[1];
	  zstars[j] = P[i].Pos[2];
	  //printf("ThisTask=%d  j=%d  mass=%g  xstars=%g  ystars=%g  zstars=%g\n", ThisTask, j, mass[j], xstars[j], ystars[j], zstars[j]);
	  flag[j] = 0;
	  j++;
	}
    }

  //memset(stellar_masses, 0, nstars_per_proc*N_STELLAR_MASS);
  for(i=0; i<nstars_per_proc*N_STELLAR_MASS; i++)
    stellar_masses[i] = 0.;

  /*
   * Now, we Gather the mass and positions to the root process, 
   * so we can create the buffer into which we'll receive the strings
   */

  int *nstars_per_proc_global_list = NULL;
  /* Only root has the received data */
  if (ThisTask == 0)
    nstars_per_proc_global_list = (int*) mymalloc("nstars_per_proc_global_list", NTask * sizeof(int)) ;

  MPI_Gather(&nstars_per_proc, 1, MPI_INT, nstars_per_proc_global_list, 1, MPI_INT, 0, MPI_COMM_WORLD);




  int total_nstars;
  int* displs = NULL;
  double* global_mass = NULL;
#ifdef METALSSSS
  double* global_Zmass = NULL;
#endif
  int* global_flag = NULL;
  double* global_xstars = NULL;
  double* global_ystars = NULL;
  double* global_zstars = NULL;
  double* global_stellar_masses = NULL;

  
  MPI_Allreduce(&nstars_per_proc, &total_nstars, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  //printf("ThisTask=%d  total_nstars=%d\n", ThisTask, total_nstars);

  if (ThisTask == 0)
    {
      displs = (int*) mymalloc("displs", NTask * sizeof(int) );
      
      for(i=1, displs[0]=0; i<NTask; i++)
	{
	  //printf("nstars_per_proc_global_list[%d]=%d\n", i, nstars_per_proc_global_list[i]);
	  displs[i] = displs[i-1] + nstars_per_proc_global_list[i-1];
	  //printf("displs[%d]=%d\n", i, displs[i]);
	}
      
      global_mass = (double*) mymalloc("global_mass", total_nstars * sizeof(double) );
      global_flag = (int*) mymalloc("global_flag", total_nstars * sizeof(int) );
      global_xstars = (double*) mymalloc("global_xstars", total_nstars * sizeof(double) );
      global_ystars = (double*) mymalloc("global_ystars", total_nstars * sizeof(double) );
      global_zstars = (double*) mymalloc("global_zstars", total_nstars * sizeof(double) );
      global_stellar_masses = (double*) mymalloc("global_zstars", total_nstars*N_STELLAR_MASS * sizeof(double) );
#ifdef METALSSSS
      global_Zmass = (double*) mymalloc("global_Zmass", 12 * total_nstars * sizeof(double) );
#endif
      //memset(global_stellar_masses, 0, total_nstars*N_STELLAR_MASS);
      for(j=0; j<total_nstars*N_STELLAR_MASS; j++)
	global_stellar_masses[j] = 0.;
    }
  
  
  /* Now we have the receive buffer, counts, and displacements, we are ready to gather the data to the root node */

  if(total_nstars > 0)
    {
      /* Don't put MPI_Gatherv inside if(ThisTask==0) !!! That would lead to deadlock... */
      MPI_Gatherv(mass, nstars_per_proc, MPI_DOUBLE, global_mass, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(flag, nstars_per_proc, MPI_INT, global_flag, nstars_per_proc_global_list, displs, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Gatherv(xstars, nstars_per_proc, MPI_DOUBLE, global_xstars, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(ystars, nstars_per_proc, MPI_DOUBLE, global_ystars, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(zstars, nstars_per_proc, MPI_DOUBLE, global_zstars, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

#ifdef METALSSSS
      nstars_per_proc *= 12;
      if(ThisTask == 0) //nstars_per_proc_global_list & displs are only allocated in the root node
	{
	  for(i=0; i<NTask; i++)
	    {
	      nstars_per_proc_global_list[i] *= 12;
	      displs[i] *= 12;
	    }
	}
      MPI_Gatherv(Zmass, nstars_per_proc, MPI_DOUBLE, global_Zmass, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      nstars_per_proc /= 12;
      if(ThisTask == 0) //nstars_per_proc_global_list & displs are only allocated in the root node
	{
	  for(i=0; i<NTask; i++)
	    {
	      nstars_per_proc_global_list[i] /= 12;
	      displs[i] /= 12;
	    }
	}
#endif      
      
      if(ThisTask == 0)
	{
	  /* do the sampling only on the root node */
#ifdef METALSSSS
	  IMFsampling(global_mass, global_xstars, global_ystars, global_zstars, global_stellar_masses, global_Zmass, global_flag, total_nstars);
#else
	  IMFsampling(global_mass, global_xstars, global_ystars, global_zstars, global_stellar_masses, global_flag, total_nstars);
#endif
	  printf("sampling done.\n");
	  
	  
	  //for(i=0; i<total_nstars; i++){
	    //printf("global_flag[%d]=%d\n", i, global_flag[i]);
	    //for(j=0; j<12; j++)
	      //printf("global_Zmass(%d,%d) = %g\n", i, j, global_Zmass[i*12+j] );
	  //}
	  
	}

      MPI_Scatterv(global_mass, nstars_per_proc_global_list, displs, MPI_DOUBLE, mass, nstars_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Scatterv(global_flag, nstars_per_proc_global_list, displs, MPI_INT, flag, nstars_per_proc, MPI_INT, 0, MPI_COMM_WORLD);


#ifdef METALSSSS
      nstars_per_proc *= 12;
      if(ThisTask == 0) //nstars_per_proc_global_list & displs are only allocated in the root node
	{
	  for(i=0; i<NTask; i++)
	    {
	      nstars_per_proc_global_list[i] *= 12;
	      displs[i] *= 12;
	    }
	}
      MPI_Scatterv(global_Zmass, nstars_per_proc_global_list, displs, MPI_DOUBLE, Zmass, nstars_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      nstars_per_proc /= 12;
      if(ThisTask == 0) //nstars_per_proc_global_list & displs are only allocated in the root node
	{
	  for(i=0; i<NTask; i++)
	    {
	      nstars_per_proc_global_list[i] /= 12;
	      displs[i] /= 12;
	    }
	}
      //for(i=0; i<nstars_per_proc; i++){
      //for(j=0; j<12; j++)
      //printf("Zmass(%d,%d) = %g\n", i, j, Zmass[i*12+j] );
      //}
#endif      

      
      /* Hack: reuse the arrays of counts and offset for stellar_masses without allocating new arrays */
      nstars_per_proc *= N_STELLAR_MASS;  
      
      if(ThisTask == 0) //nstars_per_proc_global_list & displs are only allocated in the root node
	{
	  for(i=0; i<NTask; i++)
	    {
	      nstars_per_proc_global_list[i] *= N_STELLAR_MASS;
	      displs[i] *= N_STELLAR_MASS;
	    }
	}
      
      //MPI_Barrier(MPI_COMM_WORLD);
      MPI_Scatterv(global_stellar_masses, nstars_per_proc_global_list, displs, MPI_DOUBLE, stellar_masses, nstars_per_proc, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
    }

  // Undo the hack. We don't need to do this for nstars_per_proc_global_list & displs since they will no longer be used 
  nstars_per_proc /= N_STELLAR_MASS;  
  
  /*
  for(i=0; i<nstars_per_proc; i++){
    for(j=0; j<N_STELLAR_MASS; j++)
      if(stellar_masses[i*N_STELLAR_MASS+j] > 0.)
	printf("ThisTask=%d  stellar_masses(%d,%d) = %g\n", ThisTask, i, j, stellar_masses[i*N_STELLAR_MASS+j] );
    printf("ThisTask=%d  i=%d  flag[i]=%d\n", ThisTask, i, flag[i]);
  }
  */

  //assign the results to the actual variables
  
  int iz;
  int k=0;
  int tot_spawned;
  int stars_spawned=0;
  for(i = 0; i < NumPart; i++)
    {
      //if(P[i].Type == 4)
      if(P[i].Type == 4 && P[i].sampled == 0)
        {
	  if(flag[k] == 1) //if flag is not 1 we don't assign (wait for the round)
	    {
	      P[i].sampled = 1;

	      
	      //change the dynamical mass of the star particle
	      double ratio = (mass[k] * SOLAR_MASS / All.UnitMass_in_g) / P[i].Mass;
	      //for(iz=0; iz<12; iz++)
	      //P[i].Zm[iz] *= ratio;
	      
#ifdef METALSSSS
	      //merged metals
	      for(iz=0; iz<12; iz++){
		//printf("!!!Zmass(%d,%d) = %g\n", k, iz, Zmass[k*12+iz] );
		//printf("before: ID=%llu  Zm[%d] = %g   ", P[i].ID, iz, P[i].Zm[iz] );
		P[i].Zm[iz] = Zmass[k*12 + iz];
		//printf("after:  ID=%llu  Zm[%d] = %g\n", P[i].ID, iz, P[i].Zm[iz] );
	      }
#endif

	      P[i].Mass = mass[k] * SOLAR_MASS / All.UnitMass_in_g;
	      //printf("P[i].Mass = %g, P[i].ID = %u\n", P[i].Mass, P[i].ID);                                        
	      for(iz=0; iz<3; iz++)
		P[i].dp[iz] *= ratio;

#ifdef G0_VARIABLE
	      P[i].UV_luminosity = 0.;
#ifdef CS_PHOTO_IONIZE
	      P[i].Lyman_photons_per_sec = 0.;
#endif
#endif	  


#ifdef INDIVIDUAL_STARS_SPLIT
	      int counter=0;
	      for(j=0; j<N_STELLAR_MASS; j++)
		{
		  if(stellar_masses[k*N_STELLAR_MASS+j] > 0.0)
		    {
		      counter++;
		      if(counter > 1) //more than one star, so split!
			{
			  
			  if(counter==2){
			    printf("Initial P[i].Mass = %g\n", P[i].Mass*1e10);
			    for(int jj=0; jj<N_STELLAR_MASS; jj++)
			      printf("stellar_masses[%d] = %g\n", jj, stellar_masses[k*N_STELLAR_MASS+jj]);
			  }

			  //spawn a new star particle
			  P[NumPart + stars_spawned] = P[i];
			  //P[NumPart + stars_spawned].Type = 4;
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
			  //P[i].ID += ((MyIDType) 1 << (sizeof(MyIDType) * 8 - bits));
			  //P[i].ID = TODO
#ifdef METALSSSS
			  //TODO
#endif
			  P[NumPart + stars_spawned].Mass = stellar_masses[k*N_STELLAR_MASS+j] * SOLAR_MASS / All.UnitMass_in_g;
			  P[i].Mass -= P[NumPart + stars_spawned].Mass;
			  //if(P[i].Mass<0) P[i].Mass=0;
			  //if(P[i].Mass<0) 
			  printf("ID=%d, P[i].Mass = %g, P[NumPart + stars_spawned].Mass = %g\n", P[i].ID, P[i].Mass*1e10, P[NumPart + stars_spawned].Mass*1e10);
			  
			  //sum_mass_stars += P[NumPart + stars_spawned].Mass;
			  //P[NumPart + stars_spawned].StellarAge = All.Time;
			  force_add_star_to_tree(i, NumPart + stars_spawned);
			  
			  stars_spawned++;
			}
		    }
		}
#else //INDIVIDUAL_STARS_SPLIT
	      for(j=0; j<N_STELLAR_MASS; j++)
		{
		  P[i].MstarSampleIMF[j] = stellar_masses[k*N_STELLAR_MASS+j];
		  //printf("ThisTask=%d  k=%d  j=%d  i=%d  P[i].Mass=%g  P[i].MstarSampleIMF[j] = %g\n", ThisTask, k, j, i, P[i].Mass, P[i].MstarSampleIMF[j]);
#ifdef G0_VARIABLE
		  double mass_j = P[i].MstarSampleIMF[j];
		  if( (mass_j > 0.) ) //those already exploded have negative masses and so would be filtered out here
		    {
#ifdef DEBUG_NOHYDRO_CHEM
		      double time = All.TimeFreeze;
#else
		      double time = All.Time;
#endif
		      if( (All.UnitTime_in_s/SEC_PER_YEAR/All.HubbleParam*(time - P[i].StellarAge)) < get_lifetime(mass_j) )
			{
			  P[i].UV_luminosity += pow(10., get_logL_pe( mass_j ) ); //test   
#ifdef CS_PHOTO_IONIZE
			  P[i].Lyman_photons_per_sec += pow(10., get_logS_ly( mass_j ) );
#endif
			}
		    }
#endif //G0_VARIABLE
		}
#endif //INDIVIDUAL_STARS_SPLIT
	    }//flag[k]==1
	  //printf("ThisTask=%d  k=%d  i=%d  P[i].Mass=%g\n", ThisTask, k, i, P[i].Mass);
          k++;
        }
    }
  MPI_Allreduce(&stars_spawned, &tot_spawned, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(tot_spawned > 0){
    All.TotNumPart += tot_spawned;
    NumPart += stars_spawned;
    rearrange_particle_sequence();
  }

  if(ThisTask == 0)
    {
#ifdef METALSSSS
      myfree(global_Zmass);
#endif      
      myfree(global_stellar_masses);
      myfree(global_zstars);
      myfree(global_ystars);
      myfree(global_xstars);
      myfree(global_flag);
      myfree(global_mass);
      myfree(displs);
      myfree(nstars_per_proc_global_list);
    }
#ifdef METALSSSS
  myfree(Zmass);
#endif
  myfree(stellar_masses);
  myfree(zstars);
  myfree(ystars);
  myfree(xstars);
  myfree(flag);
  myfree(mass);
  
  /*
  if (ThisTask == 0) {
    for(i=0; i<total_nstars; i++){
      printf("i=%d  global_mass=%g  global_xstars=%g  global_ystars=%g  global_zstars=%g\n", i, global_mass[i], global_xstars[i], global_ystars[i], global_zstars[i]);
    }
  }
  */
  
}

/*
 * This function is sequential and should only be called by the root node, after the MPI_Gatherv is done.
 * Input: masses and positions of the newly formed stars
 * Output: the sampled stellar masses MstarSampleIMF
 */
#ifdef METALSSSS
static void IMFsampling(double* Mass, double* Pos_x, double* Pos_y, double* Pos_z, double* MstarSampleIMF, double* ZMass, int* flag, int NumNewStars)
#else
static void IMFsampling(double* Mass, double* Pos_x, double* Pos_y, double* Pos_z, double* MstarSampleIMF, int* flag, int NumNewStars)
#endif 
      {
  int i, j, k;

  double m_proposed = 0.;
  
  double randU = 0.;

  //#define MstarSampleIMF(a,b) MstarSampleIMF[a*N_STELLAR_MASS+b]

  double r_search_2 = All.SearchingRadius * All.SearchingRadius;

  //srand(10);

  double M_clus;
  
  //int flag[NumNewStars];   // 0 = not yet flagged; 1 = flagged; -1 = flag for mass tranfer

  double totalSampledMass=0.;
  int iz;
  printf("IMF sampling started...\n");
  for(i=0; i<NumNewStars; i++)
    {  
      //printf("i=%d\n", i);
      M_clus = 0.;  
      j = 0;
      if( fabs(Mass[i]) < 1e-30) //if Mass[i] == 0.; this can happen if previous particles have borrowed mass from the current particle
	continue;

      while(1)
	{
	  //put the sampling operation before the break condition to make sure that we do the sampling at least once no matter how small Mass[i] is
	  //otherwise a small Mass[i] (smaller than "All.MassTolerance") will pass the break condition even if M_clus=0 (haven't even done one sampling yet)
	  do
	    {
	      m_proposed = pow( ( (pow(M_max, -0.3) - pow(M_min, -0.3)) * gsl_rng_uniform(random_generator) + pow(M_min, -0.3) ), (-1.0/0.3) );  //proposed m, power-law distribution                                              
	      randU = gsl_rng_uniform(random_generator);
	      //m_proposed = pow( ( (pow(M_max, -0.3) - pow(M_min, -0.3)) * ((double)rand() / RAND_MAX) + pow(M_min, -0.3) ), (-1.0/0.3) );  //proposed m, power-law distribution     
	      //randU = (double)rand() / RAND_MAX;
	    } while(randU > ( drawMassFromIMF(m_proposed) / envelope_function(m_proposed) ) ); //reject                                                                                
	  
	  //if we pass the do-loop, it means that m_proposed is accepted
	  M_clus += m_proposed;

	  totalSampledMass += m_proposed;
	  //printf("m_proposed=%g  M_clus=%g  i=%d  Mass[i]=%g\n", m_proposed, M_clus, i, Mass[i]);

	  if(m_proposed > All.IMFSampleStellarMassCut){
	    MstarSampleIMF[i*N_STELLAR_MASS+j] = m_proposed;
	    //printf("MstarSampleIMF(%d,%d) = %g\n", i, j, MstarSampleIMF[i*N_STELLAR_MASS+j]);
	    if(j>=N_STELLAR_MASS){
	      printf("Not enough space for IMF sampling, M_clus=%g\n", M_clus);
	      fflush(stdout);
	      endrun(807);
	    }
	    j++;
	  }


	  if(  fabs(M_clus - Mass[i]) <= All.MassTolerance )
	    {
	      flag[i] = 1; //perfect match! move on to the next particle
	      //printf("perfect match! move on to the next particle...\n");
	      //MstarSampleIMF[i][j-1] -= (M_clus - Mass[i]); //adjust the last sample to match the dynamical mass (do we really need this?)
	      break;
	    }
	  if(M_clus - Mass[i] > All.MassTolerance)
	    {
	      double mass_sum = Mass[i];
	      flag[i] = -1; //overflow... need to see if there are ngb stars that can kindly borrow me some mass
#ifdef METALSSSS
	      int iz;
	      double Zmass_sum[12];
	      for(iz=0; iz<12; iz++)
		Zmass_sum[iz] = ZMass[i*12+iz];
	      double frac_mass_overshoot;
#endif	      
	      int idx=0; //index of the last ngb star particle
	      for(k=i+1; k<NumNewStars; k++)
		{
		  double distance_2 = pow(Pos_x[i]-Pos_x[k], 2) + pow(Pos_y[i]-Pos_y[k], 2) + pow(Pos_z[i]-Pos_z[k], 2);
		  if(distance_2 < r_search_2) //found a ngb star within the searching radius, so borrow me some mass!
		    {
		      mass_sum += Mass[k]; 
		      flag[k] = -1;
#ifdef METALSSSS
		      for(iz=0; iz<12; iz++)
			Zmass_sum[iz] += ZMass[k*12+iz];
#endif
		    }
		  
		  if(mass_sum > M_clus) //found enough ngbs in this timebin, so the assigned mass is allowed!
		    {
#ifdef METALSSSS
		      frac_mass_overshoot = (mass_sum - M_clus) / Mass[k];
		      for(iz=0; iz<12; iz++)
			Zmass_sum[iz] -= (ZMass[k*12+iz] * frac_mass_overshoot);
#endif
		      idx = k;
		      break;
		    }
		}
	      if(idx > 0) //found enough ngbs
		{
		  //printf("found enough ngbs!!! borrowing mass...\n");
		  Mass[i] = M_clus; //change the dynamical mass here as we have enough ngbs 
#ifdef METALSSSS
		  for(iz=0; iz<12; iz++)
		    ZMass[i*12+iz] = Zmass_sum[iz];
#endif
		  flag[i] = 1;
		  for(k=i+1; k<idx; k++)
		    if (flag[k] == -1)
		      {
			//printf("setting Mass[%d] from %g to zero...\n",  k, Mass[k]);
			Mass[k] = 0.;
			flag[k] = 1;
		      }
		  //printf("setting Mass[%d] from %g to %g...\n",  idx, Mass[idx], (mass_sum - M_clus));
		  Mass[idx] = (mass_sum - M_clus); //there is still this amount of mass left, so let flag[idx] remains zero for later sampling (borrowing)
#ifdef METALSSSS
		  for(iz=0; iz<12; iz++)
		    ZMass[idx*12+iz] *= frac_mass_overshoot;
#endif		  
		}
	      else //didn't find enough ngbs
		{
		  printf("P[i].ID = %llu, didn't find enough ngbs, I'll wait for the next round...\n", P[i].ID);
		  //printf("MstarSampleIMF(i,j-1) = %g,  m_proposed = %g\n", MstarSampleIMF[i*N_STELLAR_MASS+j-1] , m_proposed);
		  //totalSampledMass -= MstarSampleIMF[i][j-1];	    //the last m_proposed is not necessarily massive enough to be recorded in MstarSampleIMF[i][j-1]
		  totalSampledMass -= m_proposed;		  
		  if(MstarSampleIMF[i*N_STELLAR_MASS+j-1] > 0) 
		    MstarSampleIMF[i*N_STELLAR_MASS+j-1] = 0.; //undo this sampling
		  for(k=i+1; k<NumNewStars; k++)
		    if (flag[k] == -1)
		      flag[k] = 0;
		}
	      break;
	    } //if(M_clus - Mass[i] > All.MassTolerance)
	}//while loop
    }
  printf("totalSampledMass=%g\n", totalSampledMass);
  printf("done.\n");


  /*
  for(i=0; i<NumNewStars; i++){
    printf("i=%d\n", i);
    for(j=0; j<N_STELLAR_MASS; j++)
      printf("%g\n", MstarSampleIMF[i][j]);
  }
  */

}

#endif
