#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>
#include "../allvars.h"
#include "../proto.h"
#include "f2c.h"

#ifdef CHEMCOOL

/*initialization of chemcool*/
void chemcool_init(void)
{
  if(ThisTask==0)
    {
      printf("initialize SG cooling and chemistry...\n");
      fflush(stdout);
    }
  
  COOLINMO();
  CHEMINMO();
  INIT_TOLERANCES();
  LOAD_H2_TABLE();
  INIT_TEMPERATURE_LOOKUP();

  
  if(ThisTask==0)
    {
      printf("initialization of SG cooling and chemistry finished.\n");
      fflush(stdout);
    }
}

/* Compute new entropy and abundances at end of timestep dt. 
 * Mode = 0 ==> update TracAbund, Entropy, Gamma, DustTemp
 * Mode = 1 ==> no update, return cooling rate
 * Mode = 2 ==> no update, return temperature
 * Mode = 3 ==> no update, return final energy
 */
double do_chemcool_step(int target, double dt, double dl, int mode)
{
  double rho, timestep, divv, energy, ekn;
  double temp;
  double yn, abh2, abhd, abco, abe;
  double abundances[TRAC_NUM], column_est;
  //double hubble_a, a3inv;
  double rpar[NRPAR];

  
#if defined TREE_RAD || defined TREE_RAD_H2
  double columni;
#endif
#ifdef TREE_RAD
  double NH;
#endif
#ifdef TREE_RAD_H2
  double NH2;
  double NCO;
#endif
  int i;


  /* XXX: CHECKME - correct units for comoving sims? */
  if (All.ComovingIntegrationOn) { /* comoving variables */
    /*
    a3inv    =  1 / (All.Time*All.Time*All.Time); 
    hubble_a = All.Omega0 / (All.Time * All.Time * All.Time)
	+ (1 - All.Omega0 - All.OmegaLambda) / (All.Time * All.Time) 
        + All.OmegaLambda;
    hubble_a = All.Hubble * All.HubbleParam * sqrt(hubble_a);
    rho      = SphP[target].Density * a3inv;
    dl       *= All.Time;
    timestep = dt / hubble_a;
    divv     = 3.0*hubble_a + P[target].DivVel / sqrt(All.Time); // XXX: check this
    COOLR.redshift = (1.0 / All.Time) - 1.0;
    */
  }
  else {
    rho      = SphP[target].Density * All.cf_a3inv;
    /* We assume that All.Time = 0 at All.InitRedshift, and that
     * all redshifts of interest are >> 1, so we can use a simple
     * approximation for the lookback times; if this isn't the case,
     * then we should probably be using comoving coordinates anyway.
     * Note that we assume that the redshift is fixed for the duration 
     * of the timestep (or in other words that dt << t_Hubble)
     */
    /*
    COOLR.redshift = pow((3. * All.Hubble * All.HubbleParam *
			    sqrt(All.Omega0) * All.Time / 2.)  + 
			    pow(1 + All.InitRedshift, -1.5), -2./3.) - 1;
    */
    COOLR.redshift = 0.;
    //printf("######################### redshift = %g\n", COOLR.redshift);
    timestep = dt;
    divv     = SphP[target].Gradients.Velocity[0][0] + SphP[target].Gradients.Velocity[1][1] + SphP[target].Gradients.Velocity[2][2];
  }  

  //COOLR.abundc  = All.InitialMetallicity * 2.46e-4;
  //COOLR.abundo  = All.InitialMetallicity * 4.90e-4;
  //COOLR.abundsi = All.InitialMetallicity * 3.47e-5;

  //WNM values (Sembach+ 2000)
  COOLR.abundc  = All.InitialMetallicity * 1.4e-4;
  COOLR.abundo  = All.InitialMetallicity * 3.2e-4;
  COOLR.abundsi = All.InitialMetallicity * 1.5e-5;
  //COOLR.abundsi = All.InitialMetallicity * 1.7e-6; //depleted
  COOLR.G0      = All.G0;
  COOLR.cosmic_ray_ion_rate  = All.CosmicRayIonRate;

  COOLI.id_current = P[target].ID;
  
#ifdef DUST_EVOLUTION
  double M_C = 0.;
  double M_Si = 0.;
  int k;
  for(k=0; k<N_BIN_SIZE; k++){
    M_C += SphP[target].CarbonDustMass[k];
    M_Si += SphP[target].SilicateDustMass[k];
  }
  double DGR_solar = 0.01;
  COOLR.dust_to_gas_ratio = ( (M_C + M_Si) / P[target].Mass ) / DGR_solar;
  //if(P[target].ID==6984461)
  //printf("COOLR.dust_to_gas_ratio = %g, M_C = %g, M_Si = %g, P[target].Mass = %g\n", COOLR.dust_to_gas_ratio, M_C, M_Si, P[target].Mass);  
#else
  COOLR.dust_to_gas_ratio = All.DGRnormalized;
#endif

  
#ifdef G0_VARIABLE  
  double u_Habing = 5.29e-14; //Habing field, in erg cm^-3
  double fac_flux2habing = 1.0 / (4.*M_PI* C_LIGHT_CGS * pow(All.UnitLength_in_cm, 2) ) / u_Habing;

  double UV_flux_tot = 0.0;
  double UV_flux_min_pix = 0.324e-2/NPIX / fac_flux2habing;
  for (i = 0; i < NPIX; i++){
    SphP[target].UV_flux[i] = DMAX(SphP[target].UV_flux[i], UV_flux_min_pix);    
    UV_flux_tot += SphP[target].UV_flux[i];
    //printf("SphP[target].UV_flux[i] = %g\n", SphP[target].UV_flux[i]);
  }
  
  double G0_tot = UV_flux_tot * fac_flux2habing * All.G0;
  //if(G0_tot < 0.324e-2)
  //G0_tot = 0.324e-2;
  COOLR.G0 = G0_tot;
  //printf("time = %g   G0 = %g\n", All.Time, COOLR.G0);

  //#ifdef G0_CONSTANT
  //COOLR.G0 = All.UVField; //overwrite with const. G0
  //#endif
#ifdef CR_SCALE_WITH_G0
  COOLR.cosmic_ray_ion_rate = (COOLR.G0/1.7) * All.CosmicRayIonRate;
  //printf("time = %g   CRrate = %g\n", All.Time, COOLR.cosmic_ray_ion_rate);
#endif
#endif


#ifdef G0_SCALE_WITH_TOTAL_SFR
  double g0 = All.FactorG0 * All.G0;
  if(g0 < 0.324e-2)
    g0 = 0.324e-2;
  COOLR.G0 = g0;
  if(P[target].ID==1)
    printf("time = %g   G0 = %g\n", All.Time, COOLR.G0);
#ifdef CR_SCALE_WITH_G0
  COOLR.cosmic_ray_ion_rate = All.FactorG0 * All.CosmicRayIonRate;
  if(P[target].ID==1)
    printf("time = %g   CRrate = %g\n", All.Time, COOLR.cosmic_ray_ion_rate);
#endif
#endif


  /* Set correct dust temperature in coolr common block */
  COOLR.tdust = SphP[target].DustTemp;

  //#ifdef CS_CHEMCOOL
  //COOLR.abundo  = P[target].Zm[3] / 16. / P[target].Zm[6];
  //COOLR.abundc  = P[target].Zm[1] / 12. / P[target].Zm[6];
  //COOLR.abundsi = P[target].Zm[5] / 28. / P[target].Zm[6];
  //#endif

#ifdef WSS_CIE_COOL
    //for(i=0; i<12; i++)
    //COOLR.Zmass[i] = P[target].Zm[i] / P[target].Mass;
#endif




  /* 'energy' is internal energy density, NOT specific internal energy [in code units] */ 
  energy = rho * SphP[target].InternalEnergy;
  if (energy < rho * All.MinEgySpec) {
      energy = rho * All.MinEgySpec;
  }

  //if(P[target].ID==10000)
  //printf("before EVOLVE_ABUNDANCES, rho = %g, u = %g\n", rho, energy/rho);

  /* Convert to cgs units */
  rho      *= UNIT_DENSITY_IN_CGS;
  timestep *= UNIT_TIME_IN_CGS;
  energy   *= UNIT_ENERGY_IN_CGS / pow(All.UnitLength_in_cm, 3);
  dl       *= All.UnitLength_in_cm;
  divv     *= All.UnitVelocity_in_cm_per_s / All.UnitLength_in_cm;
  for (i=0; i<TRAC_NUM; i++) {
    abundances[i] = SphP[target].TracAbund[i];
  }
  yn        = rho / ((1.0 + 4.0 * ABHE) * PROTONMASS);   //number density of hydrogen only

  
  //ID = COOLI.id_current;
    rpar[0] = yn;
    rpar[1] = dl;
    rpar[2] = divv;
    abh2 = abundances[IH2];
    abhd = 0.0;
#if CHEMISTRYNETWORK != 4
    abco = abundances[ICO];
#else
    abco = 0.0;
#endif
    abe  = abundances[IHP];


    ekn  = energy / (BOLTZMANN * (1.0 + ABHE - abh2 + abe) * yn);
    CALC_TEMP(&abh2, &ekn, &temp);
    //if(P[target].ID==1345975)
    //printf("energy=%g, temp=%g, yn=%g, mode=%d\n",
    //energy, temp, yn, mode);    
    if(mode==2)
      return temp;


#ifdef TREE_RAD
  for (i = 0; i < NPIX; i++) {
    columni = SphP[target].Projection[i] * UNIT_DENSITY_IN_CGS * All.UnitLength_in_cm;
    NH     = columni / ((1.0 + 4.0 * ABHE) * PROTONMASS);
    PROJECT.column_density_projection[i] = NH;
#ifdef G0_VARIABLE
    PROJECT.fac_uv[i] = SphP[target].UV_flux[i] / UV_flux_tot; //sum up to 1
    //printf("PROJECT.fac_uv[i]=%g\n", PROJECT.fac_uv[i]);
#endif
  }
#endif

#ifdef TREE_RAD_H2
  for (i = 0; i < NPIX; i++) {
    columni = SphP[target].ProjectionH2[i] * UNIT_DENSITY_IN_CGS * All.UnitLength_in_cm;
    NH2     = columni / (2.0 * PROTONMASS);
    PROJECT.column_density_projection_h2[i] = NH2;
    //printf("projection_h2=%g\n", PROJECT.column_density_projection_h2[i]);
  }

#if CHEMISTRYNETWORK != 1 && CHEMISTRYNETWORK != 4
  for (i = 0; i < NPIX; i++) {
    columni = SphP[target].ProjectionCO[i] * UNIT_DENSITY_IN_CGS * All.UnitLength_in_cm;
    NCO     = columni / (28.0 * PROTONMASS);
    PROJECT.column_density_projection_co[i] = NCO;
    //printf("projection_co=%g\n", PROJECT.column_density_projection_co[i]);
  }
#else
  /* No CO for these networks */
  for (i = 0; i < NPIX; i++) {
    PROJECT.column_density_projection_co[i] = 0.0;
  }
#endif
#endif /* TREE_RAD_H2 */


  
  /* Switch off chemistry for high-density particles */
    COOLI.no_chem = 0;

  /* When NOT using the raytrace, or the local approx for the shielding, use this bit to 
  calculate the column density for the particle and pass it to 'column_est'. 
  Note, this is used for options All.PhotochemApprox = 4, 5, and 6.
  This is then fed through the fortran until it ends up in calc_photo.F */

  column_est = 0.0;
#ifdef H2_FORM_TEST
  column_est = rpar[4];
#endif


  
  COOLR.pdv_term = 0.;
  //ID = COOLI.id_current;


    CALC_PHOTO_WRAPPER(&temp, rpar, &abh2, &abhd, &abco);


    //double temp_HII = 1e4;
    int skip_evolve_abundances = 0;
    
#ifdef PHOTO_IONIZATION
    double temp_HII = 1e4; 
    double energy_HII = temp_HII * 1.5 * BOLTZMANN * yn * (1.0 + ABHE - 0. + 1.); //assume abh2=0 and abe=1
    if(SphP[target].Ionized == 1){

      skip_evolve_abundances = 1;          

      abundances[IH2] = 0.0;
      abundances[IHP] = 0.9998;
      abundances[ICO] = 0.0;
      if(energy < energy_HII) //set a minimum T to 1e4 K (higher than 1e4 K is fine)
	energy = energy_HII;
    }
#endif
    

#ifdef CS_SFR_BUT_NOCOOL
    skip_evolve_abundances = 1;
#endif


    if(mode == 1 || mode == 2) //get cooling rate or temperature, don't evolve
      timestep = 0.0; 

    //if(P[target].ID==1345975)
    //printf("energy=%g, temp=%g, COOLR.dust_to_gas_ratio=%g, COOLR.G0=%g, mode=%d\n",
    //energy, temp, COOLR.dust_to_gas_ratio, COOLR.G0, mode);
    
    /* Evolve abundances */
    if(skip_evolve_abundances==0)
      EVOLVE_ABUNDANCES(&timestep, &dl, &yn, &divv, &energy, abundances, &column_est);

    //This needs to be after EVOLVE_ABUNDANCES s.t. RATE_EQ has already been called and lambda & lambda_chem have been updated for this particle
    double cooling_rate = 0.;  //[erg sec^-1 cm^-3] 
    for(i=0; i<28; i++){
#ifdef OUTPUT_INDIVIDUAL_COOLRATES
      SphP[target].Lambda[i] = COOLR.lambda[i] / yn / yn;
#endif
      cooling_rate += COOLR.lambda[i];
    }
    for(i=0; i<6; i++){
#ifdef OUTPUT_INDIVIDUAL_COOLRATES
      SphP[target].LambdaChem[i] = COOLR.lambda_chem[i] / yn / yn;
#endif
      cooling_rate += COOLR.lambda[i];
    }
    cooling_rate *= (All.UnitLength_in_cm * pow(UNIT_TIME_IN_CGS, 3) / All.UnitMass_in_g); //to code units 


  abh2 = abundances[IH2];
  abe  = abundances[IHP];


  /* Compute final temperature */
  ekn  = energy / (BOLTZMANN * (1.0 + ABHE - abh2 + abe) * yn);
  CALC_TEMP(&abh2, &ekn, &temp);

  /* Convert back to code units from cgs */
  energy *= pow(All.UnitLength_in_cm, 3) / UNIT_ENERGY_IN_CGS;
  rho    /= UNIT_DENSITY_IN_CGS;

  //if(P[target].ID==10000)
  //printf("after EVOLVE_ABUNDANCES, rho = %g, u = %g\n", rho, energy/rho);


  if (mode == 0) {
    for (i=0; i<TRAC_NUM; i++) {
      SphP[target].TracAbund[i] = abundances[i];
    }
    //SphP[target].Gamma   = gamma;
    SphP[target].InternalEnergy = energy / rho;
    SphP[target].InternalEnergyPred = SphP[target].InternalEnergy;
    SphP[target].Temp       = temp;
    SphP[target].DustTemp = COOLR.tdust;
#ifdef OUTPUT_SHIELD_FAC
    SphP[target].Fac_shield_h2 = COOLR.fac_shield_h2;
    SphP[target].Fac_shield_dust = COOLR.fac_shield_dust;
#endif
    return SphP[target].InternalEnergy;
  }
  else if(mode == 1){
    return cooling_rate;
  }
  else if(mode == 2){
    return temp;
  }
  else if(mode == 3){
    //don't update
    return SphP[target].InternalEnergy;
  }
  else {
    printf("Unknown mode: %d!\n", mode);
    endrun(101);
  }
  return SphP[target].InternalEnergy;
}


#endif /* CHEMCOOL */
