#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>

#include "allvars.h"
#include "proto.h"
#include "kernel.h"


/* update_weights ---  Based on density.c */
//#if defined SINGLE_SN_INJECT || defined MULTI_SN_INJECT

struct kernel_hydra
{
  double dp[3], dv[3];
  double r;
  double wk_i, wk_j, dwk_i, dwk_j;
  double h_i, h_j;
};


static struct updateweight_in
{
  MyDouble Pos[3];
  MyFloat Hsml;
  int NodeList[NODELISTLENGTH];
}
*UpdateweightIn, *UpdateweightGet;

static struct updateweight_out
{
  MyLongDouble Ngb;
#if defined STELLAR_FEEDBACK || defined CS_PHOTO_IONIZE
  MyLongDouble Rho;
  MyLongDouble U;
#endif
#ifdef EQUAL_WEIGHT_SN_INJECT
  int NumNgb;
#endif
}
*UpdateweightResult, *UpdateweightOut;


void particle2in_update(struct updateweight_in *in, int i);
void out2particle_update(struct updateweight_out *out, int i, int mode);


/* ------------------------------------------------------------------------- */
void particle2in_update(struct updateweight_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++){
#if defined MULTI_SN_INJECT
    in->Pos[k] = All.BoxSize * 0.5; //explode at the center of the box
#else
    in->Pos[k] = P[i].Pos[k];
#endif
  }

  in->Hsml = P[i].HsmlSN; //initial guess
}


/* ------------------------------------------------------------------------- */
void out2particle_update(struct updateweight_out *out, int i, int mode)
{
#ifdef EQUAL_WEIGHT_SN_INJECT
  ASSIGN_ADD(P[i].NumNgbSN, out->NumNgb, mode);
#else
  ASSIGN_ADD(P[i].NumNgbSN, out->Ngb, mode);
#endif


#if defined STELLAR_FEEDBACK || defined CS_PHOTO_IONIZE
  //if (P[i].Type >= 6)
  ASSIGN_ADD(P[i].DensAroundStar_new, out->Rho, mode);
  ASSIGN_ADD(P[i].InternalEnergyAroundStar, out->U, mode);
#endif

}


/* ------------------------------------------------------------------------- */
/* This function updates the weights for SN before exploding and calculates the
 * gas density around them (the latter for kinetic feedback purposes). Necessary
 * because gas particles neighbours of a given star could have turned into
 * stars and they need to be taken off the neighbour list for the exploding star. */
void update_weights(void)
{
  MyFloat *Left, *Right;
  int i, j, k, ndone, ndone_flag, npleft, iter=0;
  int ngrp, recvTask, place;
  long long ntot;
  //double dmax1, dmax2;
  double desnumngb, desnumngbdev;
  int save_NextParticle;
  long long n_exported=0;
  int redo_particle;


#ifndef STELLAR_FEEDBACK
#ifdef MULTI_SN_INJECT
    if(All.Time >= All.TimeToExplode && All.Time < (All.NumOfSN) * All.TimeIntervalSN){
      if(ThisTask == 0)
	printf("All.Time = %g... time to explode an SN!\n", All.Time);

      int i;
      for(i = 0; i < N_gas; i++){
	if( P[i].ID == 1 ){
          printf("marked this particle\n");
	  P[i].flag_SNII = 1;
	}
      }

      All.TimeToExplode += All.TimeIntervalSN;
      if(ThisTask == 0)
	printf("Exploded... All.TimeToExplode = %g\n", All.TimeToExplode);
    }
    //this has to be outside the for-loop above!
    //so that the inactive particle has the chance to explode when it's active later on
    //update_weights();
    //random_SN_inject();

#else
    //int i;
    int nstars_per_proc = 0;
    double tSN = All.TotN_gas * (1e6 / All.RateSN / 977813106.);
    for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
      if(P[i].Type==0)
	{
	  //SphP[i].flagSN = 0;
	  SphP[i].flagFBinj = 0;
	  double dt = (P[i].TimeBin ? (1 << P[i].TimeBin) : 0) * All.Timebase_interval;
	  double prob = dt / (tSN) * pow(SphP[i].Density / All.InitDensity, All.DensScaleIndex - 1.0);  //random=0, linear=1, Schmidt=1.5
	  //if(P[i].ID == 100 && All.Time == 0.){
	  //if(get_random_number(P[i].ID + 1) < prob){
	  double rand_number = gsl_rng_uniform(random_generator);
	  if(rand_number < prob)
	    {
	      printf("ID=%d   random number=%g   prob=%g   dt=%g   tSN=%g   SphP[i].Density=%g\n", P[i].ID, rand_number, prob, dt, tSN, SphP[i].Density);
	      P[i].flag_SNII = 1; //for RANDOM_SN_INJECT we mark it here
	      printf("Mark an explosion site at %g|%g|%g\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
	      nstars_per_proc++;
	      printf("P[i].NumNgbSN = %g, P[i].HsmlSN = %g\n", P[i].NumNgbSN, P[i].HsmlSN);
	    }
	}

    //update_weights();
    //random_SN_inject();
#endif
#endif


#ifdef STELLAR_FEEDBACK
  double age;
  double prob;
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].Type==0)
	SphP[i].flagFBinj = 0; //for timestep limiter

      if(P[i].Type==4) //check if any unexploded stars should explode (if so, set flag_SNII from 0 to 1)
	{
	  if(All.ComovingIntegrationOn)
	    age = 0.;//TODO
	  else
	    age = All.Time - P[i].StellarAge; //in code units
	  
	  age *= (UNIT_TIME_IN_CGS / SECONDS_PER_YEAR / All.HubbleParam); //in years
	  
	  if(P[i].flag_SNII == 0)
	    {
	      double lifetime, mass;
	      lifetime = All.LifeTimeSNII;
#ifdef STOCHASTIC_IMF
	      mass = P[i].MassMassiveStar;
	      lifetime = get_lifetime(mass);
	      //lifetime = All.LifeTimeSNII;
#endif
#ifdef INDIVIDUAL_STARS_SPLIT
	      mass = P[i].Mass * All.UnitMass_in_g / SOLAR_MASS;
	      lifetime = get_lifetime(mass);
#endif
#ifdef INSTANT_SN_FEEDBACK_FOR_SOME_TIME
	      if(All.Time < All.TimeInstantSN)
		lifetime = 0.;
#endif

	      if(age >= lifetime && mass >= 8.0 )
		{
		  printf("marked as SNII!  mass=%g,  lifetime=%g,  age=%g,  ID=%g\n", mass, lifetime, age, P[i].ID);
		  P[i].flag_SNII = 1; //explode in this timestep, need to find its hsml (and will be set to -1 later in stellarfeedback.c)
		}
	    }//P[i].flag_SNII == 0
#if defined STOCHASTIC_IMF && defined G0_VARIABLE
	  if(P[i].UV_luminosity > 0.0)
	    if( age > get_lifetime(P[i].MassMassiveStar) ){
	      P[i].UV_luminosity = 0.0;
	      printf("unmarked a FUV source!  mass=%g,  lifetime=%g,  age=%g,  ID=%g\n", 
		     P[i].MassMassiveStar, get_lifetime(P[i].MassMassiveStar), age, P[i].ID);
	    }
#endif
#ifdef DUST_IN_AGB
	  if(P[i].flag_AGB == 0)
	    {
	      double lifetime, mass;
#ifdef STOCHASTIC_IMF
              mass = P[i].MassMassiveStar;
              lifetime = get_lifetime(mass);
#endif
#ifdef INDIVIDUAL_STARS_SPLIT
              mass = P[i].Mass * All.UnitMass_in_g / SOLAR_MASS;
              lifetime = get_lifetime(mass);
#endif
              if(age >= lifetime && mass < 8.0 && mass > 1.0)
                {
                  printf("marked as AGB!  mass=%g,  lifetime=%g,  age=%g,  ID=%g\n", mass, lifetime, age, P[i].ID);
                  P[i].flag_AGB = 1; //enrich in this timestep, need to find its hsml (and will be set to -1 later in agb_enrich.c)                          
                }
	    }//P[i].flag_AGB == 0
#endif
	}//P[i].Type==4
    }
#endif


  if(ThisTask == 0)
    {
      printf("... start update weights ...\n");
      fflush(stdout);
    }

  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  /* don't need this anymore R2ngblist = (double *) mymalloc("R2ngblist", NTaskTimesNumPart * sizeof(double));*/

  Left = (MyFloat *) mymalloc("Left", NumPart * sizeof(MyFloat));
  Right = (MyFloat *) mymalloc("Right", NumPart * sizeof(MyFloat));

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      //#ifdef MULTI_SN_INJECT
#ifdef DUST_IN_AGB
      if(P[i].flag_SNII == 1 || P[i].flag_AGB == 1)
#else
      if(P[i].flag_SNII == 1)
#endif
	//if(i == FirstActiveParticle && ThisTask == 0)
	//#endif
	{
          Left[i] = Right[i] = 0;
	  P[i].HsmlSN = pow(All.FacNumNgbSN, 1./3) * PPP[i].Hsml; //initial guess
	}
    }

  /* allocate buffers to arrange communication */
  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct updateweight_in) + sizeof(struct updateweight_out) +
					     sizemax(sizeof(struct updateweight_in),
                         sizeof(struct updateweight_out))));

  DataIndexTable =
      (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
      (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));

  desnumngb = All.FacNumNgbSN * All.DesNumNgb;
  desnumngbdev = 1.0;

  /* we will repeat the whole thing for those particles where we didn't find enough neighbours */
  do
  {

      NextParticle = FirstActiveParticle;	/* begin with this index */

      do
      {
          BufferFullFlag = 0;
          Nexport = 0;
          save_NextParticle = NextParticle;

          {
              int mainthreadid = 0;
	      update_weight_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
          }

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


          MPI_Alltoall(Send_count, 1, MPI_INT, Recv_count, 1, MPI_INT, MPI_COMM_WORLD);


          for(j = 0, Nimport = 0, Recv_offset[0] = 0, Send_offset[0] = 0; j < NTask; j++)
          {
              Nimport += Recv_count[j];

              if(j > 0)
              {
                  Send_offset[j] = Send_offset[j - 1] + Send_count[j - 1];
                  Recv_offset[j] = Recv_offset[j - 1] + Recv_count[j - 1];
              }
          }

	  UpdateweightGet = (struct updateweight_in *) mymalloc("UpdateweightGet", Nimport * sizeof(struct updateweight_in));
	  UpdateweightIn = (struct updateweight_in *) mymalloc("UpdateweightIn", Nexport * sizeof(struct updateweight_in));


          /* prepare particle data for export */
          for(j = 0; j < Nexport; j++)
          {
              place = DataIndexTable[j].Index;

              particle2in_update(&UpdateweightIn[j], place);

              memcpy(UpdateweightIn[j].NodeList,
                  DataNodeList[DataIndexTable[j].IndexGet].NodeList, NODELISTLENGTH * sizeof(int));
          }
          /* exchange particle data */
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
          {
              recvTask = ThisTask ^ ngrp;

              if(recvTask < NTask)
              {
                  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                  {
                      /* get the particles */
		      MPI_Sendrecv(&UpdateweightIn[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct updateweight_in), MPI_BYTE,
				   recvTask, TAG_DENS_A,
				   &UpdateweightGet[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct updateweight_in), MPI_BYTE,
				   recvTask, TAG_DENS_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
              }
          }

          myfree(UpdateweightIn);

          UpdateweightResult =
	    (struct updateweight_out *) mymalloc("UpdateweightResult",
						 Nimport * sizeof(struct updateweight_out));
          UpdateweightOut =
	    (struct updateweight_out *) mymalloc("UpdateweightOut",
						 Nexport * sizeof(struct updateweight_out));


          /* now do the particles that were sent to us */

          NextJ = 0;

          {
              int mainthreadid = 0;
              update_weight_evaluate_secondary(&mainthreadid);
          }

          /* Check if this is last iteration */
          if(NextParticle < 0)
              ndone_flag = 1;
          else
              ndone_flag = 0;

          MPI_Allreduce(&ndone_flag, &ndone, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

          /* get the result */
          for(ngrp = 1; ngrp < (1 << PTask); ngrp++)
          {
              recvTask = ThisTask ^ ngrp;
              if(recvTask < NTask)
              {
                  if(Send_count[recvTask] > 0 || Recv_count[recvTask] > 0)
                  {
                      /* send the results */
		      MPI_Sendrecv(&UpdateweightResult[Recv_offset[recvTask]],
				   Recv_count[recvTask] * sizeof(struct updateweight_out),
				   MPI_BYTE, recvTask, TAG_DENS_B,
				   &UpdateweightOut[Send_offset[recvTask]],
				   Send_count[recvTask] * sizeof(struct updateweight_out),
				   MPI_BYTE, recvTask, TAG_DENS_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                  }
              }

          }

          /* add the result to the local particles */
          for(j = 0; j < Nexport; j++)
          {
              place = DataIndexTable[j].Index;
              out2particle_update(&UpdateweightOut[j], place, 1);
          }

          myfree(UpdateweightOut);
          myfree(UpdateweightResult);
          myfree(UpdateweightGet);
      }
      while(ndone < NTask);


      /* do final operations on results */
      for(i = FirstActiveParticle, npleft = 0; i >= 0; i = NextActiveParticle[i])
      {
	//#ifdef MULTI_SN_INJECT
#ifdef DUST_IN_AGB
      if(P[i].flag_SNII == 1 || P[i].flag_AGB == 1)
#else
      if(P[i].flag_SNII == 1)
#endif
	  //if(i == FirstActiveParticle && ThisTask == 0)
	  //#endif
          {

              /* now check whether we had enough neighbours */

	      desnumngb = All.FacNumNgbSN * All.DesNumNgb;
	      desnumngbdev = 1.0;

              redo_particle = 0;

	      if(P[i].NumNgbSN < (desnumngb - desnumngbdev) || (P[i].NumNgbSN > (desnumngb + desnumngbdev)))
		redo_particle = 1; /* CECILIA COMMENTED OUT THE FOLLOWING CONDITION && PPP[i].Hsml > (1.01 * All.MinGasHsml)))*/


              if(redo_particle)
              {
                  /* need to redo this particle */
                  npleft++;

                  if(Left[i] > 0 && Right[i] > 0)
                      if((Right[i] - Left[i]) < 1.0e-3 * Left[i])
                      {
                          /* this one should be ok */
                          npleft--;
			P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */
                          continue;
                      }
                  if(P[i].NumNgbSN < (desnumngb - desnumngbdev))
                      Left[i] = DMAX(P[i].HsmlSN, Left[i]);
                  else
                  {
                      if(Right[i] != 0)
                      {
                          if(P[i].HsmlSN < Right[i])
                              Right[i] = P[i].HsmlSN;
                      }
                      else
                          Right[i] = P[i].HsmlSN;
                  }

                  if(iter >= MAXITER - 10)
                  {
		      printf
			("i=%d task=%d ID=%llu Hsml=%g Left=%g Right=%g Ngbs=%g Right-Left=%g\n   pos=(%g|%g|%g)\n",
			 i, ThisTask, (unsigned long long) P[i].ID, P[i].HsmlSN, Left[i], Right[i],
			 (float) P[i].NumNgbSN, Right[i] - Left[i], P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
                      fflush(stdout);
                  }


                  if(Right[i] > 0 && Left[i] > 0)
                      P[i].HsmlSN = pow(0.5 * (pow(Left[i], 3) + pow(Right[i], 3)), 1.0 / 3);
                  else
                  {
                      if(Right[i] == 0 && Left[i] == 0)
                      {
                          char buf[1000];
                          sprintf(buf, "Right[i] == 0 && Left[i] == 0, P[i].HsmlSN=%g, P[i].Type=%d, P[i].ID=%d\n", P[i].HsmlSN, P[i].Type, P[i].ID);
                          terminate(buf);
                      }

                      if(Right[i] == 0 && Left[i] > 0)
                          P[i].HsmlSN *= 1.26;

                      if(Right[i] > 0 && Left[i] == 0)
                          P[i].HsmlSN /= 1.26;

                  }

                  /*CECILIA CHECK: why do we restrict Hsml to MinGasHsml here, since this is for stars? I COMMENTED THIS OUT */
                  /* if(PPP[i].Hsml < All.MinGasHsml)
                         PPP[i].Hsml = All.MinGasHsml; */

              }
              else
		P[i].TimeBin = -P[i].TimeBin - 1;	/* Mark as inactive */

              //if(iter == MAXITER)
              //{
              //    int old_type;
              //    old_type = P[i].Type;
	      //P[i].Type = 4;	/* no SN mark any more */
	      //SphP[i].NumNgb = 0;
	      //printf("part=%d of type=%d was assigned NumNgb=%g and type=%d\n", i, old_type,
	      //	 PPP[i].n.NumNgb, P[i].Type);
              //}
          }
      }

      sumup_large_ints(1, &npleft, &ntot);

      if(ntot > 0)
      {
          iter++;

	  if(iter > 0 && ThisTask == 0)
	    {
	      printf("update weights iteration %d: need to repeat for %d%09d particles.\n", iter,
		     (int) (ntot / 1000000000), (int) (ntot % 1000000000));
	      fflush(stdout);
	    }

          if(iter > MAXITER)
          {
              printf("failed to converge in neighbour iteration in update_weights\n");
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
  /*  myfree(R2ngblist); */
  myfree(Ngblist);


  /* mark as active again */
  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    {
      if(P[i].TimeBin < 0)
          P[i].TimeBin = -P[i].TimeBin - 1;
    }

  if(ThisTask == 0)
  {
      printf("... update weights finished \n");
      fflush(stdout);
  }

}


/* ------------------------------------------------------------------------- */
/*! This function represents the core of the SPH density computation. The
 *  target particle may either be local, or reside in the communication
 *  buffer.
 */
int update_weight_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		     int *ngblist)
{

  int j, n;
  int startnode, numngb_inbox, listindex=0;
  double r2, u, hinv, hinv3, hinv4;
  //double dx, dy, dz;

  struct kernel_hydra kernel;
  struct updateweight_in local;
  struct updateweight_out out;
  memset(&out, 0, sizeof(struct updateweight_out));

  /*  CECILIA CHECK MyLongDouble weighted_numngb;
      weighted_numngb = 0;*/

  out.Ngb = 0;

  if(mode == 0)
      particle2in_update(&local, target);
  else
      local = UpdateweightGet[target];

  kernel.h_i = local.Hsml;

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
  {
      startnode = UpdateweightGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
  }

  numngb_inbox = 0;

  while(startnode >= 0)
  {
      while(startnode >= 0)
      {
	  numngb_inbox =
	    ngb_treefind_variable_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
					  exportnodecount, exportindex, ngblist);


          if(numngb_inbox < 0)
              return -1;

          for(n = 0; n < numngb_inbox; n++)
          {
              j = ngblist[n];

              kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
              kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
              kernel.dp[2] = local.Pos[2] - P[j].Pos[2];

#if defined(BOX_PERIODIC) || defined (TALLBOX)
	      NEAREST_XYZ(kernel.dp[0], kernel.dp[1], kernel.dp[2], 1);
#endif

              r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];

	      if (r2 < kernel.h_i * kernel.h_i)
		{
		  kernel.r = sqrt(r2);
                  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
                  u = kernel.r * hinv;
                  kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, -1);

                  //double weight = kernel.wk_i * NORM_COEFF / (hinv3 * local.Ngb);
		}
	      else
		continue;

              /* Before: weighted_numngb += FLT(NORM_COEFF * wk / hinv3); */
	      out.Ngb +=  (NORM_COEFF * kernel.wk_i / hinv3);	/* 4.0/3 * PI = 4.188790204786 */
#ifdef EQUAL_WEIGHT_SN_INJECT
	      out.NumNgb +=  1;
#endif

#if defined STELLAR_FEEDBACK || defined CS_PHOTO_IONIZE
              out.Rho += (P[j].Mass * kernel.wk_i);
              out.U += (P[j].Mass / SphP[j].Density * kernel.wk_i * SphP[j].InternalEnergy);
#endif

          }
      }

      if(mode == 1)
      {
          listindex++;
          if(listindex < NODELISTLENGTH)
          {
              startnode = UpdateweightGet[target].NodeList[listindex];
              if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
          }
      }
  }

   if(mode == 0)
      out2particle_update(&out, target, 0);
  else
      UpdateweightResult[target] = out;

  return 0;
}


/* ------------------------------------------------------------------------- */
void *update_weight_evaluate_primary(void *p)
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

      //#ifdef MULTI_SN_INJECT
#ifdef DUST_IN_AGB
      if((P[i].flag_SNII == 1 || P[i].flag_AGB == 1) && P[i].TimeBin >= 0)
#else
      if(P[i].flag_SNII == 1 && P[i].TimeBin >= 0)
#endif
      //if(i == FirstActiveParticle && ThisTask == 0 && P[i].TimeBin >= 0)
	//#endif
	{
          if(update_weight_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

  }

  return NULL;

}


/* ------------------------------------------------------------------------- */
void *update_weight_evaluate_secondary(void *p)
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

      update_weight_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}

//#endif
