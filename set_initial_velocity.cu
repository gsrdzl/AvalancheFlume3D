/**************************************************************************************************
 This routine sets the initial translational and angular velocities for the particles. This
 routine runs only in HOST.
**************************************************************************************************/
#include "cabecera.hu"

void set_initial_velocity (system_par sys, particle *pp)
{
  int     nn, ii;
  double  v0, rmax, ran[2], inv_rm;
  double3 vv, ww;
  long    idum = 26752658;
  
  /* Fetches system parameters */
  v0   = sys.v0;
  nn   = sys.npart;
  rmax = sys.rmax;
  
  inv_rm = 1.0/(4.0*rmax);
  
  for (ii = 0; ii < nn; ii++)
  {
    ran[0] = 2.0*PI*RAN (&idum) - PI;
    ran[1] = PI*RAN (&idum);
    
    /* Sets translational velocity components. */
    vv.x = v0*cos (ran[0])*cos (ran[1]);
    vv.y = v0*sin (ran[0])*cos (ran[1]);
    vv.z = v0*sin (ran[1]);
    
    /* Sets angular velocity components. */
    ww.x = inv_rm*v0*cos (ran[0])*cos (ran[1]);
    ww.y = inv_rm*v0*sin (ran[0])*cos (ran[1]);
    ww.z = inv_rm*v0*sin (ran[1]);
    
    /* Saves vectors. */
    pp[ii].vv = vv;
    pp[ii].ww = ww;
  }
  
  return;
}
