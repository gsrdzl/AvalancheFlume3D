/**************************************************************************************************
 This routine is the final part of the velocity-Verlet integration algorithm for the equations of 
 motion. It corrects the velocities and updates them.
**************************************************************************************************/

#include "cabecera.hu"

/* =========================================== HOST ============================================ */

void finish_Verlet_hst (double3 *accel_i, double3 *ang_acc_i, double3 *accel_f, double3 *ang_acc_f,
                        system_par sys, particle *pp)
{
  long    ii, nn;
  double  hdt, time;
  double3 rr, vv, ww, aai, aaf, angai, angaf, side;

  /* Fetches system parameters */
  nn   = sys.npart;
  time = sys.time;
  side = sys.side;
  hdt  = sys.hdt;
  
  for (ii = 0; ii < nn; ii++)
  {
    /* Fetches particle position. */
    rr = pp[ii].rr;

    if (time > 0.0 && rr.x > side.x)
    {
      pp[ii].vv.x = pp[ii].vv.y = pp[ii].vv.z = 0.0;
      continue;
    }

    /* Fetches particle vectors */
    vv    = pp[ii].vv;
    ww    = pp[ii].ww;
    aai   = accel_i[ii];
    angai = ang_acc_i[ii];
    aaf   = accel_f[ii];
    angaf = ang_acc_f[ii];
    
    /* Corrects velocities. */
    vv.x += hdt*(aaf.x - aai.x);
    vv.y += hdt*(aaf.y - aai.y);
    vv.z += hdt*(aaf.z - aai.z);
    ww.x += hdt*(angaf.x - angai.x);
    ww.y += hdt*(angaf.y - angai.y);
    ww.z += hdt*(angaf.z - angai.z);
    
    /* Saves velocities */
    pp[ii].vv = vv;
    pp[ii].ww = ww;
  }

  return;
}

/* ========================================== DEVICE =========================================== */

__global__ void finish_Verlet_dev (double3 *accel_i, double3 *ang_acc_i, double3 *accel_f,
                                   double3 *ang_acc_f, system_par sys, particle *pp)
{
  long    ii, nn;
  double  hdt, time;
  double3 rr, vv, ww, aai, aaf, angai, angaf, side;

  /* Fetches system parameters */
  nn   = sys.npart;
  hdt  = sys.hdt;
  time = sys.time;
  side = sys.side;

  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;

  if (ii < nn)
  {
    /* Fetches particle position. */
    rr = pp[ii].rr;

    if (time > 0.0 && rr.x > side.x)
    {
      pp[ii].vv.x = pp[ii].vv.y = pp[ii].vv.z = 0.0;
      goto SKIP;
    }
    /* Fetches particle vectors */
    vv    = pp[ii].vv;
    ww    = pp[ii].ww;
    aai   = accel_i[ii];
    angai = ang_acc_i[ii];
    aaf   = accel_f[ii];
    angaf = ang_acc_f[ii];
    
    /* Corrects velocities. */
    vv.x += hdt*(aaf.x - aai.x);
    vv.y += hdt*(aaf.y - aai.y);
    vv.z += hdt*(aaf.z - aai.z);
    ww.x += hdt*(angaf.x - angai.x);
    ww.y += hdt*(angaf.y - angai.y);
    ww.z += hdt*(angaf.z - angai.z);
    
    /* Saves velocities */
    pp[ii].vv = vv;
    pp[ii].ww = ww;

    SKIP:;
  }

  return;
}
