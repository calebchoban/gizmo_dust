#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <gsl/gsl_math.h>


#include "allvars.h"
#include "proto.h"
#include "kernel.h"


#ifdef STELLAR_FEEDBACK


struct kernel_hydra
{
  double dp[3], dv[3];
  //double dx, dy, dz;
  double r;
  double wk_i, wk_j, dwk_i, dwk_j;
  double h_i, h_j;
};

struct stellarFBdata_in
{
  MyDouble Pos[3];
  MyFloat Vel[3];
  MyFloat Hsml;
  MyFloat Mass;
  MyIDType ID;
  int Timestep;
  MyFloat Ngb;

  int NodeList[NODELISTLENGTH];
}
  *StellarFBDataIn, *StellarFBDataGet;


struct stellarFBdata_out
{

}
  *StellarFBDataResult, *StellarFBDataOut;


//double PosSN[3];


static inline void particle2in_stellarFB(struct stellarFBdata_in *in, int i);
static inline void out2particle_stellarFB(struct stellarFBdata_out *out, int i, int mode);

static inline void particle2in_stellarFB(struct stellarFBdata_in *in, int i)
{
  int k;

  for(k = 0; k < 3; k++)
    {
      in->Pos[k] = P[i].Pos[k];
    }

  in->Ngb = P[i].NumNgbSN;
  in->Hsml = P[i].HsmlSN;

  in->Mass = P[i].Mass;
  in->ID = P[i].ID;


}

static inline void out2particle_stellarFB(struct stellarFBdata_out *out, int i, int mode)
{

}


void save_sn_info_to_file(){
  int i, j;
  int nstars_per_proc = 0;  //for output purpose

  for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
    if(P[i].Type==4 && P[i].flag_SNII == 1)
      nstars_per_proc++;

  //Begin the painful IO stuff...
  // total_nstars = 2+5+0+1 = 8
  int total_nstars;
  MPI_Allreduce(&nstars_per_proc, &total_nstars, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
  if(total_nstars > 0)
    {
      //Each proc has its only array
      MyIDType *id_star = (MyIDType*) mymalloc("id_star",  nstars_per_proc * sizeof(MyIDType));
      double *x_sn = (double*) mymalloc("x_sn",  nstars_per_proc * sizeof(double));
      double *y_sn = (double*) mymalloc("y_sn",  nstars_per_proc * sizeof(double));
      double *z_sn = (double*) mymalloc("z_sn",  nstars_per_proc * sizeof(double));
      double *u_sn = (double*) mymalloc("u_sn",  nstars_per_proc * sizeof(double));
      double *rho_sn = (double*) mymalloc("rho_sn",  nstars_per_proc * sizeof(double));
      double *m_sn = (double*) mymalloc("m_sn",  nstars_per_proc * sizeof(double));
      int ii=0;
      for(i = FirstActiveParticle; i >= 0; i = NextActiveParticle[i])
	if(P[i].Type==4 && P[i].flag_SNII == 1)
	  {
	    id_star[ii] = P[i].ID;
	    x_sn[ii] = P[i].Pos[0];
	    y_sn[ii] = P[i].Pos[1];
	    z_sn[ii] = P[i].Pos[2];
	    u_sn[ii] = P[i].InternalEnergyAroundStar;
	    rho_sn[ii] = P[i].DensAroundStar_new;
#ifdef STOCHASTIC_IMF
	    m_sn[ii] = P[i].MassMassiveStar;
#else
	    m_sn[ii] = 0.0;
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
      double* global_x_sn = NULL;
      double* global_y_sn = NULL;
      double* global_z_sn = NULL;
      double* global_u_sn = NULL;
      double* global_rho_sn = NULL;
      double* global_m_sn = NULL;
      //nstars_per_proc_global_list & displs are only allocated in the root node
      if (ThisTask == 0)
	{
	  //displs is the displacement array for MPI_Gatherv
	  displs = (int*) mymalloc("displs", NTask * sizeof(int) );
	  for(i=1, displs[0]=0; i<NTask; i++)
	    displs[i] = displs[i-1] + nstars_per_proc_global_list[i-1];

	  global_id_star = (MyIDType*) mymalloc("global_id_star", total_nstars * sizeof(MyIDType) );
	  global_x_sn = (double*) mymalloc("global_x_sn", total_nstars * sizeof(double) );
	  global_y_sn = (double*) mymalloc("global_y_sn", total_nstars * sizeof(double) );
	  global_z_sn = (double*) mymalloc("global_z_sn", total_nstars * sizeof(double) );
	  global_u_sn = (double*) mymalloc("global_u_sn", total_nstars * sizeof(double) );
	  global_rho_sn = (double*) mymalloc("global_rho_sn", total_nstars * sizeof(double) );
	  global_m_sn = (double*) mymalloc("global_m_sn", total_nstars * sizeof(double) );
	}

      /* Now we have the receive buffer, counts, and displacements, we are ready to gather the data to the root node */
      /* Don't put MPI_Gatherv inside if(ThisTask==0) !!! That would lead to deadlock... */
#ifdef LONGIDS
      MPI_Gatherv(id_star, nstars_per_proc, MPI_UNSIGNED_LONG_LONG, global_id_star, nstars_per_proc_global_list, displs, MPI_UNSIGNED_LONG_LONG, 0, MPI_COMM_WORLD);
#else
      MPI_Gatherv(id_star, nstars_per_proc, MPI_UNSIGNED, global_id_star, nstars_per_proc_global_list, displs, MPI_UNSIGNED, 0, MPI_COMM_WORLD);
#endif
      //MPI_Gatherv(id_star, nstars_per_proc, MPI_INT, global_id_star, nstars_per_proc_global_list, displs, MPI_INT, 0, MPI_COMM_WORLD);
      MPI_Gatherv(x_sn, nstars_per_proc, MPI_DOUBLE, global_x_sn, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(y_sn, nstars_per_proc, MPI_DOUBLE, global_y_sn, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(z_sn, nstars_per_proc, MPI_DOUBLE, global_z_sn, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(u_sn, nstars_per_proc, MPI_DOUBLE, global_u_sn, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(rho_sn, nstars_per_proc, MPI_DOUBLE, global_rho_sn, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Gatherv(m_sn, nstars_per_proc, MPI_DOUBLE, global_m_sn, nstars_per_proc_global_list, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);

      //Finally we can do the IO here!!!
      if (ThisTask == 0)
	for(j=0; j<total_nstars; j++)
	  {
#ifdef LONGIDS
            printf("All.Time=%g, x_sn=%g, y_sn=%g, z_sn=%g, u_sn=%g, rho_sn=%g, id_star=%llu, m_sn=%g\n",
                   All.Time, global_x_sn[j], global_y_sn[j], global_z_sn[j], global_u_sn[j], global_rho_sn[j], global_id_star[j], global_m_sn[j]);
            fprintf(FdSNinfo, "%g %g %g %g %g %g %llu %g\n",
                    All.Time, global_x_sn[j], global_y_sn[j], global_z_sn[j], global_u_sn[j], global_rho_sn[j], global_id_star[j], global_m_sn[j]);
#else
            printf("All.Time=%g, x_sn=%g, y_sn=%g, z_sn=%g, u_sn=%g, rho_sn=%g, id_star=%u, m_sn=%g\n",
                   All.Time, global_x_sn[j], global_y_sn[j], global_z_sn[j], global_u_sn[j], global_rho_sn[j], global_id_star[j], global_m_sn[j]);
            fprintf(FdSNinfo, "%g %g %g %g %g %g %u %g\n",
                    All.Time, global_x_sn[j], global_y_sn[j], global_z_sn[j], global_u_sn[j], global_rho_sn[j], global_id_star[j], global_m_sn[j]);
#endif
	    fflush(FdSNinfo); // can flush it, because only occuring on master steps anyways
	  }


      //free the arrays
      if (ThisTask == 0)
	{
	  myfree(global_m_sn);
	  myfree(global_rho_sn);
	  myfree(global_u_sn);
	  myfree(global_z_sn);
	  myfree(global_y_sn);
	  myfree(global_x_sn);
	  myfree(global_id_star);
	  myfree(displs);
	  myfree(nstars_per_proc_global_list);
	}
      myfree(m_sn);
      myfree(rho_sn);
      myfree(u_sn);
      myfree(z_sn);
      myfree(y_sn);
      myfree(x_sn);
      myfree(id_star);
    }//total_nstars > 0
  //Done with the painful IO stuff...

}

/*! This function is the driver routine for the calculation of hydrodynamical
 *  force and rate of change of entropy due to shock heating for all active
 *  particles .
 */
void stellarfeedback(void)
{
  int i, j, k, ngrp, ndone, ndone_flag;
  int recvTask, place;
  double timeall = 0, timecomp1 = 0, timecomp2 = 0, timecommsumm1 = 0, timecommsumm2 = 0, timewait1 =
    0, timewait2 = 0, timenetwork = 0;
  double timecomp, timecomm, timewait, tstart, tend, t0, t1;

  int save_NextParticle;

  long long n_exported = 0;


  save_sn_info_to_file();


  /* allocate buffers to arrange communication */

  int NTaskTimesNumPart;

  NTaskTimesNumPart = maxThreads * NumPart;

  Ngblist = (int *) mymalloc("Ngblist", NTaskTimesNumPart * sizeof(int));

  All.BunchSize =
    (int) ((All.BufferSize * 1024 * 1024) / (sizeof(struct data_index) + sizeof(struct data_nodelist) +
					     sizeof(struct stellarFBdata_in) +
					     sizeof(struct stellarFBdata_out) +
					     sizemax(sizeof(struct stellarFBdata_in),
						     sizeof(struct stellarFBdata_out))));
  DataIndexTable =
    (struct data_index *) mymalloc("DataIndexTable", All.BunchSize * sizeof(struct data_index));
  DataNodeList =
    (struct data_nodelist *) mymalloc("DataNodeList", All.BunchSize * sizeof(struct data_nodelist));


  CPU_Step[CPU_HYDMISC] += measure_time();
  t0 = my_second();

  NextParticle = FirstActiveParticle;	/* begin with this index */

  do
    {

      BufferFullFlag = 0;
      Nexport = 0;
      save_NextParticle = NextParticle;

      for(j = 0; j < NTask; j++)
	{
	  Send_count[j] = 0;
	  Exportflag[j] = -1;
	}

      /* do local particles and prepare export list */
      tstart = my_second();

      {
	int mainthreadid = 0;
	stellarFB_evaluate_primary(&mainthreadid);	/* do local particles and prepare export list */
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

      StellarFBDataGet = (struct stellarFBdata_in *) mymalloc("StellarFBDataGet", Nimport * sizeof(struct stellarFBdata_in));
      StellarFBDataIn = (struct stellarFBdata_in *) mymalloc("StellarFBDataIn", Nexport * sizeof(struct stellarFBdata_in));

      /* prepare particle data for export */

      for(j = 0; j < Nexport; j++)
	{
	  place = DataIndexTable[j].Index;
	  particle2in_stellarFB(&StellarFBDataIn[j], place);
	  memcpy(StellarFBDataIn[j].NodeList,
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
		  MPI_Sendrecv(&StellarFBDataIn[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct stellarFBdata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A,
			       &StellarFBDataGet[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct stellarFBdata_in), MPI_BYTE,
			       recvTask, TAG_HYDRO_A, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = my_second();
      timecommsumm1 += timediff(tstart, tend);


      myfree(StellarFBDataIn);
      StellarFBDataResult =
	(struct stellarFBdata_out *) mymalloc("StellarFBDataResult", Nimport * sizeof(struct stellarFBdata_out));
      StellarFBDataOut =
	(struct stellarFBdata_out *) mymalloc("StellarFBDataOut", Nexport * sizeof(struct stellarFBdata_out));


      //report_memory_usage(&HighMark_sphhydro, "SPH_HYDRO");

      /* now do the particles that were sent to us */

      tstart = my_second();

      NextJ = 0;

      {
	int mainthreadid = 0;
	stellarFB_evaluate_secondary(&mainthreadid);
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
		  MPI_Sendrecv(&StellarFBDataResult[Recv_offset[recvTask]],
			       Recv_count[recvTask] * sizeof(struct stellarFBdata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B,
			       &StellarFBDataOut[Send_offset[recvTask]],
			       Send_count[recvTask] * sizeof(struct stellarFBdata_out),
			       MPI_BYTE, recvTask, TAG_HYDRO_B, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		}
	    }
	}
      tend = my_second();
      timecommsumm2 += timediff(tstart, tend);



      myfree(StellarFBDataOut);
      myfree(StellarFBDataResult);
      myfree(StellarFBDataGet);
    }
  while(ndone < NTask);


  myfree(DataNodeList);
  myfree(DataIndexTable);

  myfree(Ngblist);


  /* do final operations on results */




  /* collect some timing information */

  t1 = WallclockTime = my_second();
  timeall += timediff(t0, t1);

  timecomp = timecomp1 + timecomp2;
  timewait = timewait1 + timewait2;
  timecomm = timecommsumm1 + timecommsumm2;

  CPU_Step[CPU_HYDCOMPUTE] += timecomp;
  CPU_Step[CPU_HYDWAIT] += timewait;
  CPU_Step[CPU_HYDCOMM] += timecomm;
  CPU_Step[CPU_HYDNETWORK] += timenetwork;
  CPU_Step[CPU_HYDMISC] += timeall - (timecomp + timewait + timecomm + timenetwork);
}




/*! This function is the 'core' of the SPH force computation. A target
 *  particle is specified which may either be local, or reside in the
 *  communication buffer.
 */
int stellarFB_evaluate(int target, int mode, int *exportflag, int *exportnodecount, int *exportindex,
		   int *ngblist)
{
  int startnode, numngb, listindex = 0;
  int j, n;

  double hinv, hinv3, hinv4, r2, u;

  struct kernel_hydra kernel;
  struct stellarFBdata_in local;
  struct stellarFBdata_out out;
  memset(&out, 0, sizeof(struct stellarFBdata_out));



  if(mode == 0)
    particle2in_stellarFB(&local, target);
  else
    local = StellarFBDataGet[target];

  kernel.h_i = local.Hsml;


  /* Now start the actual SPH computation for this particle */

  if(mode == 0)
    {
      startnode = All.MaxPart;	/* root node */
    }
  else
    {
      startnode = StellarFBDataGet[target].NodeList[0];
      startnode = Nodes[startnode].u.d.nextnode;	/* open it */
    }

  while(startnode >= 0)
    {
      while(startnode >= 0)
	{

	  numngb =
	    ngb_treefind_variable_threads(local.Pos, kernel.h_i, target, &startnode, mode, exportflag,
				       exportnodecount, exportindex, ngblist);

	  if(numngb < 0)
	    return -1;

	  for(n = 0; n < numngb; n++)
	    {
	      j = ngblist[n];

	      kernel.dp[0] = local.Pos[0] - P[j].Pos[0];
	      kernel.dp[1] = local.Pos[1] - P[j].Pos[1];
	      kernel.dp[2] = local.Pos[2] - P[j].Pos[2];

#if defined(BOX_PERIODIC) || defined (TALLBOX)
	      NEAREST_XYZ(kernel.dp[0], kernel.dp[1], kernel.dp[2], 1);
#endif

	      r2 = kernel.dp[0] * kernel.dp[0] + kernel.dp[1] * kernel.dp[1] + kernel.dp[2] * kernel.dp[2];
	      //kernel.h_j = P[j].Hsml;

	      if(r2 < kernel.h_i * kernel.h_i)
		{
		  kernel.r = sqrt(r2);

		  kernel_hinv(kernel.h_i, &hinv, &hinv3, &hinv4);
		  u = kernel.r * hinv;
		  kernel_main(u, hinv3, hinv4, &kernel.wk_i, &kernel.dwk_i, -1);

		  //double weight = kernel.wk_i / (hinv3 * local.Ngb);
#ifdef EQUAL_WEIGHT_SN_INJECT
		  double weight = 1.0 / local.Ngb;
#else
		  double weight = kernel.wk_i * NORM_COEFF / (hinv3 * local.Ngb);
#endif
		  double totalEnergy = 1 * 1e51 / 1.989e53;  //1 E_51
		  double energy = totalEnergy * weight;

#ifdef SPH
		  SphP[j].Entropy += energy / P[j].Mass
		    * GAMMA_MINUS1 / pow(SphP[j].d.Density, GAMMA_MINUS1);

		  SphP[j].EntropyPred = SphP[j].Entropy;
#endif

#ifdef SPH_PRESSURE_ENTROPY
		  SphP[j].EntVarPred = pow(SphP[j].EntropyPred, 1/GAMMA);
#endif

#ifdef HYDRO_MESHLESS_FINITE_MASS
		  //printf("ID=%d, Ngb=%g, weight=%g, SphP[j].InternalEnergy=%g, energy=%g, energy/P[j].Mass=%g\n", local.ID, local.Ngb, weight, SphP[j].InternalEnergy, energy, energy / P[j].Mass);
		  SphP[j].InternalEnergy += (energy / P[j].Mass);
		  SphP[j].InternalEnergyPred = SphP[j].InternalEnergy;
#ifdef PRESSURE_FLOOR
		  SphP[j].InternalEnergyTrue += (energy / P[j].Mass);
#endif
#endif
		  SphP[j].Pressure = get_pressure(j);

		  //#ifdef METALS
		  //double totalMetallicity = 1.0; //tracer field with arbitrary unit
		  //P[j].Metallicity[0] += totalMetallicity * weight; //Metallicity is a 1D array with size=1
		  //#endif

#ifdef LES_FILTER
		  double totalScalar = 1.0; //tracer field with arbitrary unit
		  SphP[j].Scalar += totalScalar * weight;
#endif


#ifdef TIMESTEP_LIMITER
		  SphP[j].Skip = 0;  //have to wake up, so can't skip it.
#endif

		  //#ifdef WAKEUP
		  SphP[j].wakeup = 1;
		  SphP[j].flagFBinj = 1;

		  //#endif
		  //printf("at random_SN_inject()... ID=%llu, energy=%g, weight=%g\n", P[j].ID, energy, weight);
		  /*
		      kernel.dvx = local.Vel[0] - SphP[j].VelPred[0];
		      kernel.dvy = local.Vel[1] - SphP[j].VelPred[1];
		      kernel.dvz = local.Vel[2] - SphP[j].VelPred[2];
		      kernel.vdotr2 = kernel.dx * kernel.dvx + kernel.dy * kernel.dvy + kernel.dz * kernel.dvz;
		  */

		}
	    }
	}

      if(mode == 1)
	{
	  listindex++;
	  if(listindex < NODELISTLENGTH)
	    {
	      startnode = StellarFBDataGet[target].NodeList[listindex];
	      if(startnode >= 0)
		startnode = Nodes[startnode].u.d.nextnode;	/* open it */
	    }
	}
    }


  /* Now collect the result at the right place */

  if(mode == 0)
    out2particle_stellarFB(&out, target, 0);
  else
    StellarFBDataResult[target] = out;


  return 0;
}

void *stellarFB_evaluate_primary(void *p)
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

      if(P[i].Type == 4 && P[i].flag_SNII == 1)
	{
	  P[i].flag_SNII = -1; //unmark it!
	  printf("SNII exploded... P[i].ID=%d, mark flag_SNII as -1...\n", P[i].ID);
	  //printf("Got here!!! Pos = %g|%g|%g\n", P[i].Pos[0], P[i].Pos[1], P[i].Pos[2]);
	  if(stellarFB_evaluate(i, 0, exportflag, exportnodecount, exportindex, ngblist) < 0)
	    break;		/* export buffer has filled up */
	}

      ProcessedFlag[i] = 1;	/* particle successfully finished */

    }

  return NULL;

}



void *stellarFB_evaluate_secondary(void *p)
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

      stellarFB_evaluate(j, 1, &dummy, &dummy, &dummy, ngblist);
    }

  return NULL;

}


#endif
