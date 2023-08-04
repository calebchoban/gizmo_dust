#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"

#ifdef PHOTO_IONIZATION

//MyFloat *N_tbi;
//MyFloat *N_ionized;
//MyFloat *Rs_old;    
//MyFloat *N_unmark;
//int *OverIonized;
//int *NumNgb;

static struct photoionize_data_in
{
  MyDouble Pos[3];
  MyFloat StroemgrenRadius;
  MyFloat OverIonized;
  //MyFloat Rs_old;
  MyFloat Right;
  int Type;
  int NodeList[NODELISTLENGTH];
  MyIDType ID;
}
  *PhotoIonizeDataIn, *PhotoIonizeDataGet;

static struct photoionize_data_out
{
  double N_ionized;
  double unmark;
  int Ngb;
}
  *PhotoIonizeDataResult, *PhotoIonizeDataOut;

void particle2in_photoionize(struct photoionize_data_in *in, int i);
void out2particle_photoionize(struct photoionize_data_out *out, int i, int mode);

MyFloat *Left, *Right;

void particle2in_photoionize(struct photoionize_data_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    in->Pos[k] = P[i].Pos[k];
  
  in->StroemgrenRadius = P[i].StroemgrenRadius;  
  in->OverIonized = P[i].OverIonized;
  //in->Rs_old = P[i].Rs_old;
  in->Right = Right[i];
  //in->Type = P[i].Type;
  in->ID = P[i].ID;
}


void out2particle_photoionize(struct photoionize_data_out *out, int i, int mode)
{
  int k, j;
  
  //printf("!!!! at out2particle, it's a type %d !!!\n", P[i].Type);

  //ASSIGN_ADD(PPP[i].n.NumNgb, out->Ngb, mode);
  if(P[i].OverIonized == 0)
    ASSIGN_ADD(P[i].N_ionized, out->N_ionized, mode);
  
  ASSIGN_ADD(P[i].NumNgbAll, out->Ngb, mode);
  ASSIGN_ADD(P[i].N_unmark, out->unmark, mode);
}



/*
 * Based on density.c
 */

void photoionize(void)
{
  //MyFloat *Left, *Right;
  int i, j, k, ii, ndone, ndone_flag, npleft, iter = 0;
  int ngrp, recvTask, place;
  long long ntot;
  double fac;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0;
  double timecomp, timecomm, timewait;
  double tstart, tend, t0, t1;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported = 0;
  int redo_particle;
  int NTaskTimesNumPart;
  int numPI=0, ntotPI=0;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type == 4)
      {
	if(P[i].flag_SNII == 0)
	  {
	    P[i].flag_PI = 1; //do the heating as long as SNII hasn't happened
	    numPI += 1;
	  }
	else
	  P[i].flag_PI = 0;
      }
  
  MPI_Allreduce(&numPI, &ntotPI, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(ntotPI == 0){
    if(ThisTask == 0)
      printf("No active ionizing stars in this timestep, so we skip photoionize...\n" );
    return;
  }
  
  if(ThisTask == 0)
    printf("....... start photoionize......................\n");

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  //N_tbi       = (MyFloat *) mymalloc("N_tbi", NumPart * sizeof(MyFloat) );  //number of ngbs to be ionized
  //N_ionized   = (MyFloat *) mymalloc("N_ionized", NumPart * sizeof(MyFloat) ); //number of ngbs get ionized
  //Rs_old      = (MyFloat *) mymalloc("Rs_old", NumPart * sizeof(MyFloat) );  
  //N_unmark    = (MyFloat *) mymalloc("N_unmark", NumPart * sizeof(MyFloat) );  
  //OverIonized = (int *) mymalloc("OverIonized", NumPart * sizeof(int) ); 
  //NumNgb      = (int *) mymalloc("NumNgb", NumPart * sizeof(int) ); //number of ngbs



  //memset(N_tbi, 0, sizeof(MyFloat));
  //memset(N_ionized, 0, sizeof(MyFloat));

  

  
  //for(i = 0; i < N_gas; i++){
  //SphP[i].Ionized = 0;
  //}
  
  for(i = 0; i < NumPart; i++){
    //P[i].Rs_old = 0.;
    P[i].OverIonized = 0;
    P[i].oldNumNgbAll = 0;
  }


  double sum_tbi, total_tbi;
  sum_tbi = total_tbi = 0.;
  double photons_per_sec = 1e48;

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if( (P[i].flag_PI == 1) && P[i].TimeBin >= 0)
	{
	  Left[i] = Right[i] = 0;	  
	  double a3inv = 1.;
	  double ascale = 1.;
	  double temp = 1.e4;
	  double beta = 2.56 * 1e-13 * pow(temp*1e-4, -0.83); //from Draine's book, p.163

#ifdef STOCHASTIC_IMF
          photons_per_sec = pow(10., get_logS_ly(P[i].MassMassiveStar) );
	  //printf("P[i].MassMassiveStar=%g,  photons_per_sec=%g\n", P[i].MassMassiveStar, photons_per_sec);
#endif
#ifdef SAMPLE_IMF 
	  photons_per_sec = P[i].Lyman_photons_per_sec;
#endif

	  /*
	  double age, hubble_a, time_hubble_a;
	  if(All.ComovingIntegrationOn){
	    hubble_a = hubble_function(All.Time);
	    time_hubble_a = All.Time * hubble_a;
	    age = integrated_time(i, time_hubble_a);
	  }
	  else{
	    age = All.Time - P[i].StellarAge;
	  }
	  int ik;
	  double metal = 0.;
	  for(ik = 1; ik < 12; ik++)	// all chemical elements but H & He 
	    if(ik != 6)
	      metal += P[i].Zm[ik];	  
	  metal /= P[i].Mass;	//metallicity in absolute units 	  

	  double dt_table = 1e5 / (All.UnitTime_in_s / SEC_PER_YEAR / All.HubbleParam); //internal unit
	  int idx = (int)(age / dt_table);

	  photons_per_sec = get_L_photo_ionize(idx, metal);
	  */


	  //photons_per_sec = 1e52;
	  //if(ThisTask == 0)
	  //printf("photons_per_sec = %g, age = %g\n", photons_per_sec, age);
	  //photons_per_sec *=  (P[i].Mass * All.UnitMass_in_g / All.HubbleParam / SOLAR_MASS / 1e6); //normalize to P[i].Mass


	  double nd = P[i].DensAroundStar_new * UNIT_DENSITY_IN_CGS * HYDROGEN_MASSFRAC / PROTONMASS; // hydrogen number density [cm^-3]

	  P[i].StroemgrenRadius = pow( 3./(4. * M_PI) * ( photons_per_sec /(nd*nd * beta ) ) , 1./3.); //cgs units
	  P[i].StroemgrenRadius /= ascale / All.HubbleParam * All.UnitLength_in_cm;         //code units
	  //!!!P[i].N_tbi = photons_per_sec /(nd*nd * beta ) * P[i].DensAroundStar_new* UNIT_DENSITY_IN_CGS / (P[i].Mass* All.UnitMass_in_g) ;

	  P[i].StroemgrenRadius  = 0.33 * P[i].Hsml; //start from a searching radius of 0.33 hsml

	  P[i].N_tbi = photons_per_sec; //photon budget
	  
	  //printf("Stroemgren radius for %u: %g, density = %g, N_tbi = %g\n", P[i].ID, P[i].StroemgrenRadius, P[i].DensAroundStar_new, P[i].N_tbi);
	  sum_tbi += P[i].N_tbi;
	  //P[i].StroemgrenRadius /= 3.;
	}
    }
  MPI_Reduce(&sum_tbi, &total_tbi, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    printf("---------------------- Aim to ionize %g 10^48 hydrogen particles -------------------------\n", total_tbi*1e-48);




  /* allocate buffers to arrange communication */

  size_t MyBufferSize = All.BufferSize;

  All.BunchSize =
    (int) ((MyBufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					   sizeof(struct photoionize_data_in) + sizeof(struct photoionize_data_out) +
					   sizemax(sizeof(struct photoionize_data_in),
						   sizeof(struct photoionize_data_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  t0 = my_second();




  if(ThisTask == 0)
    printf("Initialization of photoionize finished.......\n");


  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
    {

      NextParticle = FirstActiveParticle;	/* beginn with this index */

      do
	{
	  BufferFullFlag = 0;
	  Nexport = 0;
	  save_NextParticle = NextParticle;

	  tstart = my_second();

	  {
	    int mainthreadid = 0;
	    //if(ThisTask == 0)
	    //printf("start photoionize_evaluate_primary()...\n");
	    photoionize_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
	    //if(ThisTask == 0)
	    //printf("end photoionize_evaluate_primary()...\n");
	  }

	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);

	  if(BufferFullFlag)
	    {
	      int last_nextparticle = NextParticle;

	      NextParticle = save_NextParticle;

	      while(NextParticle >= 0)
		{
		  if(NextParticle == last_nextparticle)
		    break;

		  if(ProcessedFlag[NextParticle] != 1)
		    break;

		  ProcessedFlag[NextParticle] = 2;

		  NextParticle = NextActiveParticle[NextParticle];
		}

	      if(NextParticle == save_NextParticle)
		{
		  /* in this case, the buffer is too small to process even a single particle */
		  printf("Task %d: Type=%d pos=(%g,%g,%g) mass=%g\n",ThisTask,P[NextParticle].Type,
			 P[NextParticle].Pos[0],P[NextParticle].Pos[1],P[NextParticle].Pos[2],P[NextParticle].Mass);
		  if(P[NextParticle].Type == 0)
		    printf("   rho=%g hsml=%g\n",SphP[NextParticle].Density,PPP[NextParticle].Hsml);

		  endrun(12998);
		}


	      int new_export = 0;

	      for(j = 0, k = 0; j < Nexport; j++)
		if(ProcessedFlag[DataIndexTable[j].Index] != 2)
		  {
		    if(k < j + 1)
		      k = j + 1;

		    for(; k < Nexport; k++)
		      if(ProcessedFlag[DataIndexTable[k].Index] == 2)
			{
			  int old_index = DataIndexTable[j].Index;

			  DataIndexTable[j] = DataIndexTable[k];
			  DataNodeList[j] = DataNodeList[k];
			  DataIndexTable[j].IndexGet = j;
			  new_export++;

			  DataIndexTable[k].Index = old_index;
			  k++;
			  break;
			}
		  }
		else
		  new_export++;

	      Nexport = new_export;

	    }


	  n_exported += Nexport;

	  for(j = 0; j < NTask; j++)
	    Send_count[j] = 0;
	  for(j = 0; j < Nexport; j++)
	    Send_count[DataIndexTable[j].Task]++;

	  MYSORT_DATAINDEX(DataIndexTable, Nexport, sizeof(struct data_index), data_index_compare);

	  tstart = my_second();

	  MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);

	  tend = my_second();
	  timewait1 += timediff(tstart, tend);

	  for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
	    {
	      Nimport += Recv_count[j];

	      if(j > 0)
		{
		  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
		  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
		}
	    }

	  PhotoIonizeDataGet = (struct photoionize_data_in *) mymalloc("PhotoIonizeDataGet", Nimport * sizeof(struct photoionize_data_in));
	  PhotoIonizeDataIn = (struct photoionize_data_in *) mymalloc("PhotoIonizeDataIn", Nexport * sizeof(struct photoionize_data_in));

	  /* prepare particle data for export */
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;

	      particle2in_photoionize(&PhotoIonizeDataIn[j], place);

	      memcpy(PhotoIonizeDataIn[j].NodeList,
		     DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
	    }
	  /* exchange particle data */
	  tstart = my_second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;

	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* get the particles */
		      MPI_Sendrecv(&PhotoIonizeDataIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct photoionize_data_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &PhotoIonizeDataGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct photoionize_data_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}
	    }
	  tend = my_second();
	  timecommsumm1 += timediff(tstart, tend);

	  myfree(PhotoIonizeDataIn);
	  PhotoIonizeDataResult =
	    (struct photoionize_data_out *) mymalloc("PhotoIonizeDataResult", Nimport * sizeof(struct photoionize_data_out));
	  PhotoIonizeDataOut =
	    (struct photoionize_data_out *) mymalloc("PhotoIonizeDataOut", Nexport * sizeof(struct photoionize_data_out));

	  //report_memory_usage(&HighMark_sphdensity, "SPH_DENSITY");

	  /* now do the particles that were sent to us */

	  tstart = my_second();

	  NextJ = 0;

	  {
	    int mainthreadid = 0;
	    //if(ThisTask == 0)
	    //printf("start photoionize_evaluate_secondary()...\n");
	    photoionize_evaluate_secondary(&mainthreadid);
	    //if(ThisTask == 0)
	    //printf("end photoionize_evaluate_secondary()...\n");
	  }

	  tend = my_second();
	  timecomp2 += timediff(tstart, tend);

	  if(NextParticle < 0)
	    ndone_flag = 1;
	  else
	    ndone_flag = 0;

	  tstart = my_second();
	  MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	  tend = my_second();
	  timewait2 += timediff(tstart, tend);


	  /* get the result */
	  tstart = my_second();
	  for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
	    {
	      recvTask = ThisTask ^ ngrp;
	      if(recvTask < NTask)
		{
		  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
		    {
		      /* send the results */
		      MPI_Sendrecv(&PhotoIonizeDataResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct photoionize_data_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &PhotoIonizeDataOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct photoionize_data_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		    }
		}

	    }
	  tend = my_second();
	  timecommsumm2 += timediff(tstart, tend);


	  /* add the result to the local particles */
	  tstart = my_second();
	  for(j = 0; j < Nexport; j++)
	    {
	      place = DataIndexTable[j].Index;
	      out2particle_photoionize(&PhotoIonizeDataOut[j], place, 1);
	    }
	  tend = my_second();
	  timecomp1 += timediff(tstart, tend);


	  myfree(PhotoIonizeDataOut);
	  myfree(PhotoIonizeDataResult);
	  myfree(PhotoIonizeDataGet);
	}
      while(ndone < NTask);


      /* do final operations on results */
      tstart = my_second();

      npleft = 0;



      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	{
	  if( (P[i].flag_PI == 1) && P[i].TimeBin >= 0)
	    {
	     
	      /* now check whether we had enough neighbours */

	      if(P[i].OverIonized == 1 && P[i].N_unmark > 0)
		P[i].N_ionized = P[i].N_ionized - P[i].N_unmark;

#ifdef STOCHASTIC_IMF
	      photons_per_sec = pow(10., get_logS_ly(P[i].MassMassiveStar) );
#endif
#ifdef SAMPLE_IMF 
	      photons_per_sec = P[i].Lyman_photons_per_sec;
#endif

	      double temp = 1.e4;
	      double beta = 2.56 * 1e-13 * pow(temp*1e-4, -0.83); //from Draine's book, p.163
	      double nd = (P[i].DensAroundStar_new * UNIT_DENSITY_IN_CGS * HYDROGEN_MASSFRAC / PROTONMASS); //use DensAroundStar to estimate typical density of ngbs
	      //if(nd < 1.0)
	      //nd = 1.0; // ngbs are unlikely to be all < nd 	      
	      double tolerance = beta * 1000. * (All.InitialGasMass * All.UnitMass_in_g * HYDROGEN_MASSFRAC / PROTONMASS); // for SPH particle of 4 solar mass and nH = 10 cm^-3
	      //double tolerance = photons_per_sec * 0.01;
	      //if(tolerance < 1e46)
	      //tolerance = 1e46;
	      //if(ThisTask == 0)
	      //printf("tolerance=%g\n", tolerance);

	      //double tolerance = 1.;


	      
	      if(iter >= MAXITER-10)
	      //if(iter >= 0)
		{
		  printf
		    ("i=%d task=%d, ID=%llu, Left=%g, Right=%g, Rs=%g, Ngbs=%d, old_Ngbs=%d, N_tbi=%g, N_ionized=%g, tol=%g, nd=%g\n",
		     i, ThisTask, (unsigned long long) P[i].ID, Left[i], Right[i], P[i].StroemgrenRadius, P[i].NumNgbAll, P[i].oldNumNgbAll, 
		     P[i].N_tbi, P[i].N_ionized, tolerance, nd );
		  fflush(stdout);
		}


	      redo_particle = 1;
	      if( fabs(P[i].N_ionized - P[i].N_tbi) > tolerance ) //need to satisfy this condition
		{
		  //increase by only one particle and we go from underIonized to overIonized, meaning that we can never satisfy the tolerance
		  if( ((P[i].NumNgbAll - P[i].oldNumNgbAll) == 1) && (P[i].OverIonized == 0) && (P[i].N_ionized > P[i].N_tbi) )
		    redo_particle = 0;
		  
		  //decrease by only one particle and we go from overIonized to underIonized, meaning that we can never satisfy the tolerance
		  if( ((P[i].NumNgbAll - P[i].oldNumNgbAll) == -1) && (P[i].OverIonized == 1) && (P[i].N_ionized < P[i].N_tbi) )
		    redo_particle = 0;
		  
		}
	      else
		redo_particle = 0;
	      
	      //no larger than 50 pc, as this method is biased toward dense clouds (angular resolution = 4 pi)
	      if(P[i].StroemgrenRadius > 0.05){
		redo_particle = 0;
		//printf("ID=%llu, HII region terminating at 50 pc....\n", (unsigned long long) P[i].ID);
	      }

	      if(iter >= MAXITER)
		{
		  printf
		    ("!!!!!WARNING!!!!! i=%d task=%d, ID=%llu, Left=%g, Right=%g, Rs=%g, Ngbs=%d, N_tbi=%g, N_ionized=%g, tolerance=%g\n",
		     i, ThisTask, (unsigned long long) P[i].ID, Left[i], Right[i], P[i].StroemgrenRadius, P[i].NumNgbAll, P[i].N_tbi, P[i].N_ionized, tolerance );
		  fflush(stdout);		  
		  redo_particle = 0;
		}

	      if(redo_particle)
		{
		  npleft++;
		  //printf("StroemgrenRadius = %g\n", P[i].StroemgrenRadius);

		  if(P[i].N_ionized < P[i].N_tbi) //need more ngbs
		    {
		      P[i].OverIonized = 0;
		      //Left[i] = DMAX(P[i].StroemgrenRadius, Left[i]);
		      Left[i] = P[i].StroemgrenRadius;

		      P[i].N_tbi = P[i].N_tbi - P[i].N_ionized; //sucessfully used up this amount of photons!!!
		    }
		  else
		    {
		      P[i].OverIonized = 1;
		      Right[i] = P[i].StroemgrenRadius;
		    }

		  //-------- modify StroemgrenRadius -------
		  if(Left[i] > 0. && Right[i] > 0.)
		    {
		      if(Left[i] > Right[i])
			{
			  printf("Left > Right, this shouldn't happen......\n");
			  fflush(stdout);
			  endrun(517);
			}
		      
		      if( ( fabs(Right[i] - Left[i]) < 1.e-6 * Left[i])  ) //to prevent from trapping
			{
			  if(P[i].OverIonized == 0)
			    {
			      fac = 1. + fabs( P[i].N_tbi - P[i].N_ionized ) / (NUMDIMS * photons_per_sec);
			      if(fac > 1.1)
				fac = 1.1;

			      Left[i] = P[i].StroemgrenRadius;
			      P[i].StroemgrenRadius *= fac;
			      Right[i] = P[i].StroemgrenRadius * fac;
			      //printf("need %g particles, increasing Stroemgren Radius......\n", (P[i].N_tbi - P[i].N_ionized));
			    }
			  else if(P[i].OverIonized == 1)
			    {
			      fac = 1. + fabs( P[i].N_tbi - P[i].N_ionized ) / (NUMDIMS * photons_per_sec);
			      if(fac > 1.1)
				fac = 1.1;
			      
			      Right[i] = P[i].StroemgrenRadius;
			      P[i].StroemgrenRadius /= fac;
			      Left[i] = P[i].StroemgrenRadius / fac;	
			      //printf("need %g particles, decreasing Stroemgren Radius......\n", (P[i].N_tbi - P[i].N_ionized));
			    }
			}
		      else //squeeze
			P[i].StroemgrenRadius = pow( 0.5 * ( pow(Left[i], 3) + pow(Right[i], 3) ), 1./3 );
		    }
		  else if(Left[i] > 0. && Right[i] == 0. )//need more ngbs
		    {
		      fac = 1. + fabs( P[i].N_tbi - P[i].N_ionized ) / (NUMDIMS * photons_per_sec);
		      if(fac < 1.1)
			P[i].StroemgrenRadius *= fac;
		      else
			P[i].StroemgrenRadius *= 1.1;
		      //printf("right = 0, need %g particles, increasing Stroemgren Radius......\n", (P[i].N_tbi - P[i].N_ionized) );
		    }
		  else if(Left[i] == 0. && Right[i] > 0. )//too many ngbs
		    {
		      fac = 1. + fabs( P[i].N_tbi - P[i].N_ionized ) / (NUMDIMS * photons_per_sec);
		      if(fac < 1.1)
			P[i].StroemgrenRadius /= fac;
		      else
			P[i].StroemgrenRadius /= 1.1;
		      //printf("left = 0, need %g particles, decreasing Stroemgren Radius......\n", (P[i].N_tbi - P[i].N_ionized));
		    }
		  P[i].oldNumNgbAll = P[i].NumNgbAll;
		}
	      else
		{
		  P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
		  //printf("ID=%d converged, marked as inactive!!!\n", P[i].ID );
		}
	    }//P[i].flag_PI == 1 && P[i].TimeBin >= 0
	}//loop over active particles
      

      tend = my_second();
      timecomp1 += timediff(tstart, tend);

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
	{
	  iter++;

	  if(iter > 0 && iter%10 == 0 && ThisTask == 0)
	    {
	      printf("photoionize ngb iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

	  if(iter > MAXITER*2)
	    {
	      printf("failed to converge in neighbour iteration in photoionize()\n");
	      fflush(stdout);
	      endrun(1155);
	    }
	}


    }
  while(ntot > 0);

  myfree(DataNodeList);
  myfree(DataIndexTable);
  myfree(Right);
  myfree(Left);
  
  //myfree(NumNgb);
  //myfree(OverIonized);
  //myfree(N_unmark);
  //myfree(Rs_old);
  //myfree(N_ionized);
  //myfree(N_tbi);

  myfree(Ngblist);
  
  /* mark as active again */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].TimeBin < 0){
	P[i].TimeBin = -P[i].TimeBin - 1;
	//printf("type= %d\n", P[i].Type);
      }
      //if(P[i].Type == 9)
	//printf("After photoionize... Stroemgren radius for %u: %g, density = %g, N_tbi = %g\n", P[i].ID, P[i].StroemgrenRadius, P[i].DensAroundStar_new, P[i].N_tbi);

    }
  
  

  int sum_ionized, total_ionized;
  sum_ionized = total_ionized = 0;

  for(i = 0; i < N_gas; i++)
    {

      //unmark those that are now beyond HII regions
      if(SphP[i].Ionized == 1)
	SphP[i].Ionized = 0;

      //mark the newly ionized ngbs
      if(SphP[i].Ionized == 2){
	SphP[i].Ionized = 1;
	SphP[j].wakeup = 1;
	SphP[j].flagFBinj = 1;
      }
      //counting
      if(SphP[i].Ionized == 1)
	{   
	  //We don't do heating here. Instead, we'll do it in the cooling routine	  
	  sum_ionized++;
	}
    }

  
  MPI_Reduce(&sum_ionized, &total_ionized, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  if(ThisTask == 0)
    printf("---------------- Ionized %d particles ----------------\n", total_ionized);

  
  /* collect some timing information */

  t1 = WallclockTime = my_second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  //CPU_Step[CPU_DENSCOMPUTE] += timecomp;
  //CPU_Step[CPU_DENSWAIT] += timewait;
  //CPU_Step[CPU_DENSCOMM] += timecomm;
  //CPU_Step[CPU_DENSMISC] += timeall - (timecomp + timewait + timecomm);
}


/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int photoionize_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		     int *ngblist)
{
  int k, j, n;
  int startnode, numngb_inbox, listindex = 0;
  double r2, h2, u, mass_j, r_search;
  double dx, dy, dz;

  //struct kernel_density kernel;
  struct photoionize_data_in local;
  struct photoionize_data_out out;
  memset(&out, 0, sizeof(struct photoionize_data_out));



  if(mode == 0)
    particle2in_photoionize(&local, target);
  else
    local = PhotoIonizeDataGet[target];

  //out.MinTimeStep = local.TimeStep; 
  //out.MinTimeStep = 1<<28; //start from a large number

  //h2 = local.StroemgrenRadius * local.StroemgrenRadius;

  if(local.OverIonized == 0)
    r_search = local.StroemgrenRadius;
  if(local.OverIonized == 1)
    //r_search = pow( (2.*pow(local.StroemgrenRadius, 3) - pow(local.Rs_old, 3)), 1./3);
    r_search = local.Right;

  //kernel_hinv(local.Hsml, &kernel.hinv, &kernel.hinv3, &kernel.hinv4);

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = PhotoIonizeDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{
	  
	  numngb_inbox =
	    ngb_treefind_variable_threads(local.Pos, r_search, target, &startnode, mode, exportflag,
					  exportnodecount, exportindex, ngblist);
	  
	  if(numngb_inbox < 0)
	    return -1;
	  
	  for(n = 0; n < numngb_inbox; n++)
	    {
	      j = ngblist[n];
	      
	      dx = local.Pos[0] - P[j].Pos[0];
	      dy = local.Pos[1] - P[j].Pos[1];
	      dz = local.Pos[2] - P[j].Pos[2];
#if defined(BOX_PERIODIC) || defined (TALLBOX)
              NEAREST_XYZ(dx, dy, dz, 1);
#endif
	      
	      r2 = dx * dx + dy * dy + dz * dz;
	      
	      if(r2 < pow(r_search, 2) )
		{
		  double nH_j = SphP[j].Density * UNIT_DENSITY_IN_CGS * HYDROGEN_MASSFRAC / PROTONMASS; // hydrogen number density [cm^-3]
		  double NH_j = P[j].Mass * All.UnitMass_in_g * HYDROGEN_MASSFRAC / PROTONMASS; //number of hydrogen particle in m_j [#]
		  
		  double temp = 1.e4;
		  double beta = 2.56 * 1e-13 * pow(temp*1e-4, -0.83); //from Draine's book, p.163
		  
		  if(local.OverIonized == 0) 
		    {
		      out.Ngb += 1; //all ngbs regardless whether they are ionized or not
		      
		      if(SphP[j].Ionized != 2) //Ionized == 0 or 1 (2 means this ngb is taken already by another star)
			{
			  SphP[j].Ionized = 2; //mark the nearest ngbs as Ionized=2
			  SphP[j].IonizedBy = local.ID;
			  //!!!out.N_ionized += 1.;
			  out.N_ionized += NH_j * nH_j * beta;
			  //printf("ID=%d got ionized by ID=%d\n", P[j].ID, local.ID);
			}
		    }
		  else if(local.OverIonized == 1) //ionized too many ngbs during last iteration, so unmark the farthest ngbs
		    {
		      if( r2 < pow(local.StroemgrenRadius, 2) )
			out.Ngb += 1; //all ngbs regardless whether they are ionized or not
		      else if( (SphP[j].Ionized == 2) && (SphP[j].IonizedBy == local.ID) ) //r2 > pow(local.StroemgrenRadius, 2)
			{
			  //gas particles can only be unmarked by the star particle that ionized them (otherwise we get infinite loops) 
			  SphP[j].Ionized = 0; //unmark
			  out.unmark += NH_j * nH_j * beta;
			  //printf("ID=%d got unmarked by ID=%d\n", P[j].ID, local.ID);
			}
		    }
		}
	      
	    }
	}
      
      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = PhotoIonizeDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }

  if(mode == 0)
    out2particle_photoionize(&out, target, 0);
  else
    PhotoIonizeDataResult[target] = out;

  return 0;
}

void *photoionize_evaluate_primary(void *p)
{
  int thread_id = *(int *) p;
  int i, j;
  int *exportflag, *exportnodecount, *exportindex, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;
  exportflag = Exportflag + thread_id * NTask;
  exportnodecount = Exportnodecount + thread_id * NTask;
  exportindex = Exportindex + thread_id * NTask;

  /* Note: exportflag is local to each thread */
  for(j = 0; j < NTask; j++)
    exportflag[j] = -1;

  while(1)
    {
      int exitFlag = 0;
      {
	if(BufferFullFlag != 0 || NextParticle < 0)
	  {
	    exitFlag = 1;
	  }
	else
	  {
	    i = NextParticle;
	    ProcessedFlag[i] = 0;
	    NextParticle = NextActiveParticle[NextParticle];
	  }
      }
      if(exitFlag)
	break;

      //if( (P[i].Type == 9 || (P[i].Type & 7) == 7 ) && P[i].TimeBin >= 0)
      if( (P[i].flag_PI == 1) && P[i].TimeBin >= 0)
	{
	  if(photoionize_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */
    }

  return NULL;
}

void *photoionize_evaluate_secondary(void *p)
{
  int thread_id = *(int *) p;

  int j, dummy, *ngblist;

  ngblist = Ngblist + thread_id * NumPart;


  while(1)
    {
      {
	j = NextJ;
	NextJ++;
      }

      if(j >= Nimport)
	break;

      photoionize_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;
}
#endif
