/**************************************************************************************************
 This routine sets all forces and torques to zero.
**************************************************************************************************/
#include "cabecera.hu"

/* =========================================== HOST ============================================ */

void clear_forces_hst (system_par sys, double3 *force, double3 *torque)
{
  long    ii, nn;
  double3 zero;
  
  /* Fetches system parameters. */
  nn = sys.npart;

  /* Sets forces and torques to zero. */
  zero.x = zero.y = zero.z = 0.0;
  
  for (ii = 0; ii < nn; ii++)
  {
    /* Saves forces and torques. */
    force[ii] = zero;
    torque[ii] = zero;
  }

  return;
}

/* ========================================== DEVICE =========================================== */

__global__ void clear_forces_dev (system_par sys, double3 *force, double3 *torque)
{
  long    ii, nn;
  double3 zero;
  
  /* Fetches system parameters. */
  nn = sys.npart;

  /* Sets forces and torques to zero. */
  zero.x = zero.y = zero.z = 0.0;

  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (ii < nn)
  {
    /* Saves forces and torques. */
    force[ii] = zero;
    torque[ii] = zero;
  }

  return;
}
