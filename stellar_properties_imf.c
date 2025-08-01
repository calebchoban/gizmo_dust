#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "allvars.h"
#include "proto.h"

#ifdef SAMPLE_IMF_FROM_GAS

double get_lifetime(double mass){
  double A, B;
  if (mass<3.0){
    A = -2.926;
    B =  9.892;
  }
  else if ( (mass>=3.0) && (mass<7.0) ){
    A = -2.405;
    B =  9.641;
  }
  else if ( (mass>=7.0) && (mass<15.0) ){
    A = -1.765;
    B =  9.105;
  }
  else if (mass>=15.0){
    A = -0.808;
    B =  7.954;
  }

  double logAge = A * log10(mass) + B;
  double age =  pow(10., logAge); // [yr]

  return age;
}



double tbl_Mass[21] = {  0.8     ,   0.9     ,   1.      ,   
			 1.1     ,   1.25    ,   1.35    ,   
			 1.5     ,   1.7     ,   2.      ,   
			 2.5     ,   3.      ,   4.      ,   
			 5.      ,   6.999998,   8.999978,
			 11.999839,  14.999568,  19.999581,  
			 24.999077,  31.998098,  39.996651};

#ifdef DUST_IN_AGB
double tbl_logMdustC_agb[21] = { 
  -20, -20,
  -9.85511, -5.40362, -3.22366, -3.05173, -2.70223, -2.50045, -2.34332, -2.23844, -2.15583, -1.89245, -3.23912, -3.20300,
  -20, -20, -20, -20, -20, -20, -20
};//0.1 solar metallicity
#endif

/*
double tbl_logL_pe[21] = { 
  27.96209322, 28.85707938, 29.59881618, 30.21651755, 
  30.82357317, 31.32654796, 32.12136499, 33.04425441, 
  33.99946029, 34.819052, 35.22220551, 35.84163201, 
  36.319291, 36.98056909, 37.34013383, 37.68326336, 
  37.93895345, 38.17665699, 38.40573288, 38.67544519, 
  38.84857662
}; // solar metallicity

double tbl_logS_ly[21] = { 
  24.41647, 25.683735, 26.734005, 27.608646, 
  28.468212, 29.180407, 30.305836, 32.479298, 
  35.61226, 38.042862, 39.484432, 41.623093, 
  43.022007, 44.548347, 45.449894, 46.48047, 
  47.296062, 48.106155, 48.593216, 49.006676, 
  49.291897
}; // solar metallicity


double tbl_logL_pd[21] = { 
19.16405285, 20.42093579, 21.46260143, 22.33007665, 23.1826014 ,
23.88895918, 25.00516858, 26.75938705, 29.0678275 , 31.59921022,
33.47594486, 34.73365074, 35.43984705, 36.28939792, 36.73162131,
37.13675621, 37.42216435, 37.6834924 , 37.92283649, 38.20300079,
38.38144516
}; // solar metallicity
*/


double tbl_logL_pe[21] = { 30.20996888,  30.93678454,  31.54526414,  32.16387272,
			   32.84424722,  33.36105133,  34.06013471,  34.63134049,
			   35.10873089,  35.57254876,  35.96533313,  36.4787194 ,
			   36.89274178,  37.36706648,  37.60982623,  37.96216217,
			   38.13282568,  38.41724208,  38.6726397 ,  38.88445242,  
			   39.10040492}; // 0.1 solar metallicity


double tbl_logL_pd[21] = { 
  22.25713315, 23.5007376 , 24.54186533, 25.60032405, 27.02670099,
  28.14253881, 29.4652504 , 32.1309753 , 33.76070995, 34.54036975,
  35.06710103, 35.72145087, 36.21294897, 36.774404  , 37.06215449,
  37.44806061, 37.63512952, 37.9344108 , 38.20002426, 38.42028279,
  38.64203753
}; // 0.1 solar metallicity


//Ionizing photon rate [sec^-1]
double tbl_logS_ly[21] = { 27.53371648,  28.68920307,  29.6565598 ,  30.64001949,
			   32.20576521,  33.54968817,  36.6087002 ,  38.18852184,
			   39.57356117,  40.9447155 ,  41.99112754,  43.29239514,
			   44.16185591,  45.28618199,  46.06581943,  47.16068103,
			   47.87237443,  48.57095918,  48.98145146,  49.33857457,  
			   49.66294627}; // 0.1 solar metallicity

#ifdef DUST_IN_AGB
double get_logMdustC_agb(double mass){
  double a, y;
  int idx = get_index(mass);
  if ( mass < tbl_Mass[0] )
    y = tbl_logMdustC_agb[0];
  else if ( mass > tbl_Mass[21-1])
    y = tbl_logMdustC_agb[21-1];
  else{
    a = (mass - tbl_Mass[idx]) / (tbl_Mass[idx+1] - tbl_Mass[idx]);
    y = (1.0 - a) * tbl_logMdustC_agb[idx] + a * tbl_logMdustC_agb[idx+1];
  }
  return y;
}
#endif

double get_logL_pe(double mass){
  double a, y;
  int idx = get_index(mass);
  if ( mass < tbl_Mass[0] )
    y = tbl_logL_pe[0];
  else if ( mass > tbl_Mass[21-1])
    y = tbl_logL_pe[21-1];
  else{
    a = (mass - tbl_Mass[idx]) / (tbl_Mass[idx+1] - tbl_Mass[idx]);
    y = (1.0 - a) * tbl_logL_pe[idx] + a * tbl_logL_pe[idx+1];
  }
  //printf("tbl_logL_pe[idx], tbl_logL_pe[idx+1] = %g  %g  \n", tbl_logL_pe[idx], tbl_logL_pe[idx+1]);
  //printf("mass = %g,  y = %g\n", mass, y);
  return y;
}

double get_logS_ly(double mass){
  double a, y;
  int idx = get_index(mass);
  if ( mass < tbl_Mass[0] )
    y = tbl_logS_ly[0];
  else if ( mass > tbl_Mass[21-1])
    y = tbl_logS_ly[21-1];
  else{
    a = (mass - tbl_Mass[idx]) / (tbl_Mass[idx+1] - tbl_Mass[idx]);
    y = (1.0 - a) * tbl_logS_ly[idx] + a * tbl_logS_ly[idx+1];
  }
  //printf("tbl_logS_ly[idx], tbl_logS_ly[idx+1] = %g  %g  \n", tbl_logS_ly[idx], tbl_logS_ly[idx+1]);
  //printf("mass = %g,  y = %g\n", mass, y);
  return y;
}



/* Get the index of table using binary search (note: mass array has to be sorted) */
int get_index(double search){
  int n = 21;
  int first, last, middle;
  int idx;

  first = 0;
  last = n;
  middle = (first+last)/2;

  while (first < last) {
    middle = first + (last-first)/2;
    if (tbl_Mass[middle] == search) {
      //printf("%g found at location %d.\n", search, middle);
      idx = middle;
      return idx;
      //break;
    }
    else if (tbl_Mass[middle] < search)
      first = middle + 1;
    else
      last = middle;

    //printf("first = %d, last = %d, middle = %d\n", first, last, middle);
  }
  idx = last - 1;

  return idx;
}



#endif //G0_VARIABLE
