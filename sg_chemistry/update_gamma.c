#include <math.h>
#include <stdio.h>
#include "allvars.h"
#include "proto.h"

#ifdef CHEMCOOL

double calc_gamma_from_entropy(double entropy, double density, double abh2, double old_gamma)
{
  double density_cgs, new_gamma, final_gamma, ekn, energy_estimate;
  double dgamma, entropy_error, entropy_error_tolerance, bisection_tolerance;
  double F1, etest1, Ftest1, etest2, Ftest2, gtest1, gtest2;
  double low, high, mid, gmid, Fmid;
  double entropy_new;
  int count;

  bisection_tolerance = entropy_error_tolerance = 1e-3;
  
  density_cgs = density * All.UnitDensity_in_cgs;
  energy_estimate = entropy * pow(density, old_gamma) / (old_gamma - 1.0);
  ekn = calc_ekn(energy_estimate, density_cgs, abh2);
  CALC_GAMMA(&abh2, &ekn, &new_gamma);  
  dgamma = old_gamma - new_gamma;

  if (energy_estimate < 0.0) {
    printf("Error - negative energy!\n");
    printf("Energy = %g, gamma = %g, entropy = %g\n", energy_estimate, old_gamma, entropy);
  }

  /* Actually, we should use the old density here, but since we're just estimating the size
   * of the error and since the density is unlikely to change by a large amount between calls,
   * this should be OK
   */  
  entropy_error = ((1.0 - dgamma / (old_gamma - 1.0)) * pow(density, dgamma)) - 1.0;
  if (fabs(entropy_error) < entropy_error_tolerance) {
    /* Small error, so just return the updated estimate of gamma; this should be the most
     * common case
     */
    return new_gamma;
  }

  /* If the error in the entropy is too large, then we need to determine the new gamma and
   * internal energy in a self-consistent fashion. For now, we use a simple bisection algorithm,
   * but this could be replaced with a better algorithm in the future if this proves to be a
   * hotspot
   */

   F1     = (new_gamma - 1.0) * energy_estimate - entropy * pow(density, new_gamma);
   etest1 = etest2 = energy_estimate;
   Ftest1 = Ftest2 = F1;
 
   count = 0;
   while ((F1 < 0.0 && Ftest1 < 0.0 && Ftest2 < 0.0) || (F1 > 0.0 && Ftest1 > 0.0 && Ftest2 > 0.0)) {
     etest1 *= 1.1;
     ekn = calc_ekn(etest1, density_cgs, abh2);
     CALC_GAMMA(&abh2, &ekn, &gtest1);   
     Ftest1 = (gtest1 - 1.0) * etest1 - entropy * pow(density, gtest1);

     etest2 *= 0.9;
     ekn = calc_ekn(etest2, density_cgs, abh2);
     CALC_GAMMA(&abh2, &ekn, &gtest2);   
     Ftest2 = (gtest2 - 1.0) * etest2 - entropy * pow(density, gtest2);

     count++;
     if (count > 1000) {
       printf("Error: can't find starting points for bisection\n");
       printf("F1 = %g, Ftest1 = %g, Ftest2 = %g\n", F1, Ftest1, Ftest2);
       printf("energy = %g, etest1 = %g, etest2 = %g\n", energy_estimate, etest1, etest2);
       fflush(stdout);
       endrun(1010);
     }
   }

   if (count == 0) {
     printf("Error: we skipped the bounds-setting code!\n");
     printf("F1 = %g, Ftest1 = %g, Ftest2 = %g\n", F1, Ftest1, Ftest2);
     printf("energy = %g, etest1 = %g, etest2 = %g\n", energy_estimate, etest1, etest2);
     printf("energy estimate = %g, entropy error = %g, entropy %g, density %g, new_gamma %g\n", 
            energy_estimate, entropy_error, entropy, density, new_gamma);
     printf("old gamma = %g, H2 fraction = %g\n", old_gamma, abh2);
     endrun(1010);
   }

   /* Find the bounds for our bisection algorithm */

   if (F1 < 0.0) {
     low = energy_estimate;
     if (Ftest1 > 0.0) {
       high = etest1;
     }
     else {
       high = etest2;
     }
   }
   else {
     high = energy_estimate;
     if (Ftest1 < 0.0) {
       low = etest1;
     }
     else {
       low = etest2;
     }
   }

   /* Now do the bisection */
   count = 0;
   mid = low + (high - low) / 2.0;
#if 0
   printf("starting values: low %g, mid %g, high %g\n", low, mid, high);
#endif
   while (fabs(mid - low) > bisection_tolerance * mid && fabs(mid - high) > bisection_tolerance * mid) {
     ekn = calc_ekn(mid, density_cgs, abh2);
     CALC_GAMMA(&abh2, &ekn, &gmid);
     if (gmid < 1.0) {
       printf("Small gmid: %g, ekn = %g, abh2 = %g\n", gmid, ekn, abh2);
     }
     Fmid = (gmid - 1.0) * mid - entropy * pow(density, gmid);
     if (Fmid <= 0.0) {
       low = mid;
     }
     else {
       high = mid;
     }
     mid = low + (high - low) / 2.0;

     count++;
     if (count > 1000) {
       printf("Error: too many iterations in bisection");
       printf("low = %g, high = %g, mid = %g, Fmid = %g\n", low, high, mid, Fmid);
       fflush(stdout);
       endrun(1011);
     }
   }

   entropy_new = (gmid - 1.0) * mid / pow(density, gmid);
#if 0
   printf("gmid = %g, Fmid = %g, density = %g, entropy = %g, abh2 = %g, entropy_new = %g\n", gmid, Fmid, density, entropy, abh2, entropy_new);
#endif
   final_gamma = gmid;
   if (final_gamma < 1.0) {
     printf("Final gamma %g, low %g, mid %g, high %g\n", final_gamma, low, mid, high);
   }
   return final_gamma;
}

double calc_ekn(double energy, double density_cgs, double abh2)
{
  double ekn, energy_cgs, yn;

  energy_cgs = energy * All.UnitEnergy_in_cgs / pow(All.UnitLength_in_cm, 3);
  yn = density_cgs / ((1.0 + 4.0 * ABHE) * PROTONMASS);
  ekn = energy_cgs / (BOLTZMANN * (1.0 + ABHE - abh2) * yn);

  return ekn;
}
#endif
