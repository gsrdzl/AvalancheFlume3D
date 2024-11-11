#include "cabecera.hu"

/**************************************************************************************************
   Routine for calculating the masses, centers of mass and moments of inertia of dimers, given the
   density, separation between centers, and component radiuses.
**************************************************************************************************/

/* ============================================ HOST =========================================== */

void calculate_particle_parameters_hst (system_par sys_par, particle *pp_par)
{
  int     nn, ii;
  double  rad, rad3, rad5, rho, mass, aux1, aux2;
  double3 I_mom;

  /* Fetches system parameters */
  nn = sys_par.npart;
  
  aux1 = 4.0*PI/3.0;
  aux2 = 112.0*PI/15.0;
  for (ii = 0; ii < nn; ii++)
  {
    /* Fetches particle parameters. */
    rad = pp_par[ii].radius;
    rho = pp_par[ii].density;

    rad3 = rad*rad*rad;
    rad5 = rad3*rad*rad;
    /* Calculates particle's mass. */
    mass = 4.0*aux1*rad3*rho;

    /* Calculates particle's moment of inertia. */
    I_mom.x = I_mom.y = I_mom.z = aux2*rho*rad5;

    /* Saves */
    pp_par[ii].mass = mass;
    pp_par[ii].moment_of_inertia = I_mom;
  }
  
  return;
}

/* =========================================== DEVICE ========================================== */

__global__ void calculate_particle_parameters_dev (system_par sys_par, particle *pp_par)
{
  int     nn, ii;
  double   rad, rad3, rad5, rho, mass, aux1, aux2;
  double3  I_mom;

  /* Fetches system parameters */
  nn = sys_par.npart;
  
  aux1 = 4.0*PI/3.0;
  aux2 = 112.0*PI/15.0;
  
  ii = threadIdx.x + blockIdx.x*blockDim.x;
  
  for (ii = 0; ii < nn; ii++)
  {
    /* Fetches particle parameters. */
    rad = pp_par[ii].radius;
    rho = pp_par[ii].density;

    rad3 = rad*rad*rad;
    rad5 = rad3*rad*rad;
    /* Calculates particle's mass. */
    mass = aux1*rad3*rho;

    /* Calculates particle's moment of inertia. */
    I_mom.x = I_mom.y = I_mom.z = aux2*rho*rad5;

    /* Saves */
    pp_par[ii].mass = mass;
    pp_par[ii].moment_of_inertia = I_mom;
  }
  
  return;
}
