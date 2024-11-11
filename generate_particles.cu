#include "cabecera.hu"

/**************************************************************************************************
 This routine generates a distribution of particles radii. This routine runs only in HOST
**************************************************************************************************/

void generate_particles (particle *pp, system_par sys)
{
  int    nn, ii;
  double rad, rmin, rmax;
  
  /* Fetches system parameters */
  nn   = sys.npart;
  rmin = sys.rmin;
  rmax = sys.rmax;
  
  /* Sets density of particles */
  for (ii = 0; ii < nn; ii++)
    pp[ii].density = 2.6;
  
  /* Sets particle radii. */
  for (ii = 0; ii < nn; ii++)
  {
    /* Bidisperse distribution: small particles are 4/5 of total and big ones are 1/5 of total */
    rad = (ii%20) ? rmin : rmax;
    
    /* Saves */
    pp[ii].radius = rad;
  }

  return;
}
