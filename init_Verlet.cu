/**************************************************************************************************
 This routine is the initial part of the velocity-Verlet integration algorithm for the equations of
 motion. It updates the positions and orientations of particles, given the previous step velocities
 and accelerations and makes the velocity-prediction step.
**************************************************************************************************/

#include "cabecera.hu"

/* =========================================== HOST ============================================ */

void init_Verlet_hst (double3 *accel, double3 *ang_acc, system_par sys, particle *pp)
{
  long    ii, nn;
  double  dt, hdt, mag_dOmega, cdO, sdO, omcdO, norm, rad, pos, time;
  double3 rr, vv, aa, Rx, Ry, Rz, ww, anga, dOmega, ndOmega, dRx, dRy, dRz, Rtmpx, Rtmpy, Rtmpz,
          ra, rb, rc, rd, side;

  /* Fetches system parameters. */
  nn   = sys.npart;
  dt   = sys.dt;
  hdt  = sys.hdt;
  time = sys.time;
  side = sys.side;

  for (ii = 0; ii < nn; ii++)
  {
    /* Fetches particle position. */
    rr   = pp[ii].rr;
    if (time > 0.0 && rr.x > side.x)
    {
      pp[ii].radius = 0.0;
      continue;
    }

    /* Fetches particle vectors. */
    vv   = pp[ii].vv;
    aa   = accel[ii];
    Rx   = pp[ii].Rmatx;
    Ry   = pp[ii].Rmaty;
    Rz   = pp[ii].Rmatz;
    ww   = pp[ii].ww;
    rad  = pp[ii].radius;
    anga = ang_acc[ii];
    
    /* Does traslational part. First calculates new positions. */
    rr.x += dt*(vv.x + hdt*aa.x);
    rr.y += dt*(vv.y + hdt*aa.y);
    rr.z += dt*(vv.z + hdt*aa.z);
    
    /* Then, predicts new velocities. */
    vv.x += dt*aa.x;
    vv.y += dt*aa.y;
    vv.z += dt*aa.z;
    
    /* Now rotational part. First calculates differential rotation vector: components and
       magnitude. */
    dOmega.x   = dt*(ww.x + hdt*anga.x);
    dOmega.y   = dt*(ww.y + hdt*anga.y);
    dOmega.z   = dt*(ww.z + hdt*anga.z);
    mag_dOmega = sqrt (dOmega.x*dOmega.x + dOmega.y*dOmega.y + dOmega.z*dOmega.z);
    
    /* Normalizes dOmega */
    if (mag_dOmega > 0.0)
    {
      ndOmega.x = dOmega.x/mag_dOmega;
      ndOmega.y = dOmega.y/mag_dOmega;
      ndOmega.z = dOmega.z/mag_dOmega;
            
      cdO   = cos (mag_dOmega);
      sdO   = sin (mag_dOmega);
      omcdO = 1.0 - cdO;
      
      /* Differential rotation matrix. */
      dRx.x = omcdO*ndOmega.x*ndOmega.x + cdO;
      dRx.y = omcdO*ndOmega.x*ndOmega.y - sdO*ndOmega.z;
      dRx.z = omcdO*ndOmega.x*ndOmega.z + sdO*ndOmega.y;
      dRy.x = omcdO*ndOmega.y*ndOmega.x + sdO*ndOmega.z;
      dRy.y = omcdO*ndOmega.y*ndOmega.y + cdO;
      dRy.z = omcdO*ndOmega.y*ndOmega.z - sdO*ndOmega.x;
      dRz.x = omcdO*ndOmega.z*ndOmega.x - sdO*ndOmega.y;
      dRz.y = omcdO*ndOmega.z*ndOmega.y + sdO*ndOmega.x;
      dRz.z = omcdO*ndOmega.z*ndOmega.z + cdO;
    }
    else
    {
      dRx.x = 1.0; dRx.y = 0.0; dRx.z = 0.0;
      dRy.x = 0.0; dRy.y = 1.0; dRy.z = 0.0;
      dRz.x = 0.0; dRz.y = 0.0; dRz.z = 1.0;
    }
    
    /* Rotates the particle. */
    Rtmpx.x = dRx.x*Rx.x + dRx.y*Ry.x + dRx.z*Rz.x;
    Rtmpx.y = dRx.x*Rx.y + dRx.y*Ry.y + dRx.z*Rz.y;
    Rtmpx.z = dRx.x*Rx.z + dRx.y*Ry.z + dRx.z*Rz.z;
    Rtmpy.x = dRy.x*Rx.x + dRy.y*Ry.x + dRy.z*Rz.x;
    Rtmpy.y = dRy.x*Rx.y + dRy.y*Ry.y + dRy.z*Rz.y;
    Rtmpy.z = dRy.x*Rx.z + dRy.y*Ry.z + dRy.z*Rz.z;
    Rtmpz.x = dRz.x*Rx.x + dRz.y*Ry.x + dRz.z*Rz.x;
    Rtmpz.y = dRz.x*Rx.y + dRz.y*Ry.y + dRz.z*Rz.y;
    Rtmpz.z = dRz.x*Rx.z + dRz.y*Ry.z + dRz.z*Rz.z;

    /* Normalization correction. */
    norm = sqrt (Rtmpx.x*Rtmpx.x + Rtmpy.x*Rtmpy.x + Rtmpz.x*Rtmpz.x);
    Rx.x = Rtmpx.x/norm;
    Ry.x = Rtmpy.x/norm;
    Rz.x = Rtmpz.x/norm;
    norm = sqrt (Rtmpx.y*Rtmpx.y + Rtmpy.y*Rtmpy.y + Rtmpz.y*Rtmpz.y);
    Rx.y = Rtmpx.y/norm;
    Ry.y = Rtmpy.y/norm;
    Rz.y = Rtmpz.y/norm;
    norm = sqrt (Rtmpx.z*Rtmpx.z + Rtmpy.z*Rtmpy.z + Rtmpz.z*Rtmpz.z);
    Rx.z = Rtmpx.z/norm;
    Ry.z = Rtmpy.z/norm;
    Rz.z = Rtmpz.z/norm;
    
    /* Predicts new angular velocity */
    ww.x += dt*anga.x;
    ww.y += dt*anga.y;
    ww.z += dt*anga.z;
    
    /* Saves vectors. */
    pp[ii].rr    = rr;
    pp[ii].vv    = vv;
    pp[ii].Rmatx = Rx;
    pp[ii].Rmaty = Ry;
    pp[ii].Rmatz = Rz;
    pp[ii].ww    = ww;

    /* Calculates positions of spheres. */
    pos = rad/ROOT2;
    ra.x = rr.x + pos*(Rx.x + Rx.y + Rx.z);
    ra.y = rr.y + pos*(Ry.x + Ry.y + Ry.z);
    ra.z = rr.z + pos*(Rz.x + Rz.y + Rz.z);
    rb.x = rr.x + pos*(Rx.x - Rx.y - Rx.z);
    rb.y = rr.y + pos*(Ry.x - Ry.y - Ry.z);
    rb.z = rr.z + pos*(Rz.x - Rz.y - Rz.z);
    rc.x = rr.x + pos*(-Rx.x - Rx.y + Rx.z);
    rc.y = rr.y + pos*(-Ry.x - Ry.y + Ry.z);
    rc.z = rr.z + pos*(-Rz.x - Rz.y + Rz.z);
    rd.x = rr.x + pos*(-Rx.x + Rx.y - Rx.z);
    rd.y = rr.y + pos*(-Ry.x + Ry.y - Ry.z);
    rd.z = rr.z + pos*(-Rz.x + Rz.y - Rz.z);
    
    /* Saves positions. */
    pp[ii].rr_a = ra;
    pp[ii].rr_b = rb;
    pp[ii].rr_c = rc;
    pp[ii].rr_d = rd;
  }

  return;
}

/* ========================================== DEVICE =========================================== */

__global__ void init_Verlet_dev (double3 *accel, double3 *ang_acc, system_par sys, particle *pp)
{
  long    ii, nn;
  double  dt, hdt, mag_dOmega, cdO, sdO, omcdO, norm, rad, pos, time;
  double3 rr, vv, aa, Rx, Ry, Rz, ww, anga, dOmega, ndOmega, dRx, dRy, dRz, Rtmpx, Rtmpy, Rtmpz,
          ra, rb, rc, rd, side;
 
  /* Fetches system parameters. */
  nn   = sys.npart;
  dt   = sys.dt;
  hdt  = sys.hdt;
  time = sys.time;
  side = sys.side;
  
  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (ii < nn)
  {
    /* Fetches particle position. */
    rr   = pp[ii].rr;
    if (time > 0.0 && rr.x > side.x)
    {
      pp[ii].radius = 0.0;
      goto SKIP;
    }

    /* Fetches particle vectors. */
    vv   = pp[ii].vv;
    aa   = accel[ii];
    Rx   = pp[ii].Rmatx;
    Ry   = pp[ii].Rmaty;
    Rz   = pp[ii].Rmatz;
    ww   = pp[ii].ww;
    rad  = pp[ii].radius;
    anga = ang_acc[ii];
    
    /* Does traslational part. First calculates new positions. */
    rr.x += dt*(vv.x + hdt*aa.x);
    rr.y += dt*(vv.y + hdt*aa.y);
    rr.z += dt*(vv.z + hdt*aa.z);
    
    /* Then, predicts new velocities. */
    vv.x += dt*aa.x;
    vv.y += dt*aa.y;
    vv.z += dt*aa.z;
    
    /* Now rotational part. First calculates differential rotation vector: components and
       magnitude. */
    dOmega.x   = dt*(ww.x + hdt*anga.x);
    dOmega.y   = dt*(ww.y + hdt*anga.y);
    dOmega.z   = dt*(ww.z + hdt*anga.z);
    mag_dOmega = norm3d (dOmega.x, dOmega.y, dOmega.z);
    
    /* Normalizes dOmega */
    if (mag_dOmega > 0.0)
    {
      ndOmega.x = dOmega.x/mag_dOmega;
      ndOmega.y = dOmega.y/mag_dOmega;
      ndOmega.z = dOmega.z/mag_dOmega;
            
      cdO   = cos (mag_dOmega);
      sdO   = sin (mag_dOmega);
      omcdO = 1.0 - cdO;
      
      /* Differential rotation matrix. */
      dRx.x = omcdO*ndOmega.x*ndOmega.x + cdO;
      dRx.y = omcdO*ndOmega.x*ndOmega.y - sdO*ndOmega.z;
      dRx.z = omcdO*ndOmega.x*ndOmega.z + sdO*ndOmega.y;
      dRy.x = omcdO*ndOmega.y*ndOmega.x + sdO*ndOmega.z;
      dRy.y = omcdO*ndOmega.y*ndOmega.y + cdO;
      dRy.z = omcdO*ndOmega.y*ndOmega.z - sdO*ndOmega.x;
      dRz.x = omcdO*ndOmega.z*ndOmega.x - sdO*ndOmega.y;
      dRz.y = omcdO*ndOmega.z*ndOmega.y + sdO*ndOmega.x;
      dRz.z = omcdO*ndOmega.z*ndOmega.z + cdO;
    }
    else
    {
      dRx.x = 1.0; dRx.y = 0.0; dRx.z = 0.0;
      dRy.x = 0.0; dRy.y = 1.0; dRy.z = 0.0;
      dRz.x = 0.0; dRz.y = 0.0; dRz.z = 1.0;
    }
    
    /* Rotates the particle. */
    Rtmpx.x = dRx.x*Rx.x + dRx.y*Ry.x + dRx.z*Rz.x;
    Rtmpx.y = dRx.x*Rx.y + dRx.y*Ry.y + dRx.z*Rz.y;
    Rtmpx.z = dRx.x*Rx.z + dRx.y*Ry.z + dRx.z*Rz.z;
    Rtmpy.x = dRy.x*Rx.x + dRy.y*Ry.x + dRy.z*Rz.x;
    Rtmpy.y = dRy.x*Rx.y + dRy.y*Ry.y + dRy.z*Rz.y;
    Rtmpy.z = dRy.x*Rx.z + dRy.y*Ry.z + dRy.z*Rz.z;
    Rtmpz.x = dRz.x*Rx.x + dRz.y*Ry.x + dRz.z*Rz.x;
    Rtmpz.y = dRz.x*Rx.y + dRz.y*Ry.y + dRz.z*Rz.y;
    Rtmpz.z = dRz.x*Rx.z + dRz.y*Ry.z + dRz.z*Rz.z;

    /* Normalization correction. */
    norm = norm3d (Rtmpx.x, Rtmpy.x, Rtmpz.x);
    Rx.x = Rtmpx.x/norm;
    Ry.x = Rtmpy.x/norm;
    Rz.x = Rtmpz.x/norm;
    norm = norm3d (Rtmpx.y, Rtmpy.y, Rtmpz.y);
    Rx.y = Rtmpx.y/norm;
    Ry.y = Rtmpy.y/norm;
    Rz.y = Rtmpz.y/norm;
    norm = norm3d (Rtmpx.z, Rtmpy.z, Rtmpz.z);
    Rx.z = Rtmpx.z/norm;
    Ry.z = Rtmpy.z/norm;
    Rz.z = Rtmpz.z/norm;
    
    /* Predicts new angular velocity */
    ww.x += dt*anga.x;
    ww.y += dt*anga.y;
    ww.z += dt*anga.z;
    
    /* Saves vectors. */
    pp[ii].rr    = rr;
    pp[ii].vv    = vv;
    pp[ii].Rmatx = Rx;
    pp[ii].Rmaty = Ry;
    pp[ii].Rmatz = Rz;
    pp[ii].ww    = ww;

    /* Calculates positions of spheres. */
    pos = rad/ROOT2;
    ra.x = rr.x + pos*(Rx.x + Rx.y + Rx.z);
    ra.y = rr.y + pos*(Ry.x + Ry.y + Ry.z);
    ra.z = rr.z + pos*(Rz.x + Rz.y + Rz.z);
    rb.x = rr.x + pos*(Rx.x - Rx.y - Rx.z);
    rb.y = rr.y + pos*(Ry.x - Ry.y - Ry.z);
    rb.z = rr.z + pos*(Rz.x - Rz.y - Rz.z);
    rc.x = rr.x + pos*(-Rx.x - Rx.y + Rx.z);
    rc.y = rr.y + pos*(-Ry.x - Ry.y + Ry.z);
    rc.z = rr.z + pos*(-Rz.x - Rz.y + Rz.z);
    rd.x = rr.x + pos*(-Rx.x + Rx.y - Rx.z);
    rd.y = rr.y + pos*(-Ry.x + Ry.y - Ry.z);
    rd.z = rr.z + pos*(-Rz.x + Rz.y - Rz.z);
    
    /* Saves positions. */
    pp[ii].rr_a = ra;
    pp[ii].rr_b = rb;
    pp[ii].rr_c = rc;
    pp[ii].rr_d = rd;

    SKIP:;
  }

  return;
}
