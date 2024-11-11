/**************************************************************************************************
 This routine calculates the acceleration of particles given the total force acting on them.
**************************************************************************************************/

/* ========================================== HOST ============================================= */

#include "cabecera.hu"

void get_accel_hst (double3 *force, double3 *torque, double3 *accel, double3 *ang_acc,
                    system_par sys, particle *pp)
{
  long    ii, nn;
  double  mass, Cxy, Cyz, Czx, rad, adrag, vmag, time;
  double3 ff, aa, ww_lab, ww_body, tq_lab, tq_body, anga_lab, anga_body, Rx, Ry, Rz, Imom, vv, side;

  /* Fetches system parameters */
  nn   = sys.npart;
  time = sys.time;
  side = sys.side;

  for (ii = 0; ii < nn; ii++)
  {
    if (time > 0.0 && pp[ii].rr.x > side.x)
      continue;

    /* Fetches particle parameters. */
    mass   = pp[ii].mass;
    rad    = pp[ii].radius;
    Imom   = pp[ii].moment_of_inertia;
    vv     = pp[ii].vv;
    ww_lab = pp[ii].ww;
    Rx     = pp[ii].Rmatx;
    Ry     = pp[ii].Rmaty;
    Rz     = pp[ii].Rmatz;

    /* Calculates needed parameters. */
    Cxy = (Imom.x - Imom.y)/Imom.z;
    Cyz = (Imom.y - Imom.z)/Imom.x;
    Czx = (Imom.z - Imom.x)/Imom.y;

    /* Fetches forces. */
    ff     = force[ii];
    tq_lab = torque[ii];

    /* Obtains translational accelerations. */
    /* Drag acceleration. */
    adrag  = DRAG*rad*rad/mass;
    vmag   = sqrt (vv.x*vv.x + vv.y*vv.y + vv.z*vv.z);
    adrag *= vmag;

    aa.x = ff.x/mass - adrag*vv.x;
    aa.y = ff.y/mass - adrag*vv.y;
    aa.z = ff.z/mass - adrag*vv.z - GRAV;

    /* Rotates angular velocities and torque vectors to principal axis system */
    ww_body.x = ww_lab.x*Rx.x + ww_lab.y*Ry.x + ww_lab.z*Rz.x;
    ww_body.y = ww_lab.x*Rx.y + ww_lab.y*Ry.y + ww_lab.z*Rz.y;
    ww_body.z = ww_lab.x*Rx.z + ww_lab.y*Ry.z + ww_lab.z*Rz.z;
    tq_body.x = tq_lab.x*Rx.x + tq_lab.y*Ry.x + tq_lab.z*Rz.x;
    tq_body.y = tq_lab.x*Rx.y + tq_lab.y*Ry.y + tq_lab.z*Rz.y;
    tq_body.z = tq_lab.x*Rx.z + tq_lab.y*Ry.z + tq_lab.z*Rz.z;
    
    /* Obtains angular acceleration. */
    anga_body.x = tq_body.x/Imom.x + Cyz*ww_body.y*ww_body.z;
    anga_body.y = tq_body.y/Imom.y + Czx*ww_body.z*ww_body.x;
    anga_body.z = tq_body.z/Imom.z + Cxy*ww_body.x*ww_body.y;

    /* Rotates angular acceleration back to laboratory system. */
    anga_lab.x = anga_body.x*Rx.x + anga_body.y*Rx.y + anga_body.z*Rx.z;
    anga_lab.y = anga_body.x*Ry.x + anga_body.y*Ry.y + anga_body.z*Ry.z;
    anga_lab.z = anga_body.x*Rz.x + anga_body.y*Rz.y + anga_body.z*Rz.z;
    
    /* Saves accelerations. */
    accel[ii]   = aa;
    ang_acc[ii] = anga_lab;
  }
  return;
}

/* ========================================= DEVICE ============================================ */

__global__ void get_accel_dev (double3 *force, double3 *torque, double3 *accel, double3 *ang_acc,
                               system_par sys, particle *pp)
{
  long    ii, nn;
  double  mass, Cxy, Cyz, Czx, rad, adrag, vmag, time;
  double3 ff, aa, ww_lab, ww_body, tq_lab, tq_body, anga_lab, anga_body, Rx, Ry, Rz, Imom, vv, side;

  /* Fetches system parameters */
  nn   = sys.npart;
  time = sys.time;
  side = sys.side;

  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;

  if (ii < nn)
  {
    if (time > 0.0 && pp[ii].rr.x > side.x)
      goto SKIP;

    /* Fetches particle parameters. */
    mass   = pp[ii].mass;
    rad    = pp[ii].radius;
    Imom   = pp[ii].moment_of_inertia;
    vv     = pp[ii].vv;
    ww_lab = pp[ii].ww;
    Rx     = pp[ii].Rmatx;
    Ry     = pp[ii].Rmaty;
    Rz     = pp[ii].Rmatz;

    /* Calculates needed parameters. */
    Cxy = (Imom.x - Imom.y)/Imom.z;
    Cyz = (Imom.y - Imom.z)/Imom.x;
    Czx = (Imom.z - Imom.x)/Imom.y;

    /* Fetches forces. */
    ff     = force[ii];
    tq_lab = torque[ii];

    /* Obtains translational accelerations. */
    /* Drag acceleration. */
    adrag  = DRAG*rad*rad/mass;
    vmag   = norm3d (vv.x, vv.y, vv.z);
    adrag *= vmag;

    aa.x = ff.x/mass - adrag*vv.x;
    aa.y = ff.y/mass - adrag*vv.y;
    aa.z = ff.z/mass - adrag*vv.z - GRAV;

    /* Rotates angular velocities and torque vectors to principal axis system */
    ww_body.x = ww_lab.x*Rx.x + ww_lab.y*Ry.x + ww_lab.z*Rz.x;
    ww_body.y = ww_lab.x*Rx.y + ww_lab.y*Ry.y + ww_lab.z*Rz.y;
    ww_body.z = ww_lab.x*Rx.z + ww_lab.y*Ry.z + ww_lab.z*Rz.z;
    tq_body.x = tq_lab.x*Rx.x + tq_lab.y*Ry.x + tq_lab.z*Rz.x;
    tq_body.y = tq_lab.x*Rx.y + tq_lab.y*Ry.y + tq_lab.z*Rz.y;
    tq_body.z = tq_lab.x*Rx.z + tq_lab.y*Ry.z + tq_lab.z*Rz.z;
    
    /* Obtains angular acceleration. */
    anga_body.x = tq_body.x/Imom.x + Cyz*ww_body.y*ww_body.z;
    anga_body.y = tq_body.y/Imom.y + Czx*ww_body.z*ww_body.x;
    anga_body.z = tq_body.z/Imom.z + Cxy*ww_body.x*ww_body.y;

    /* Rotates angular acceleration back to laboratory system. */
    anga_lab.x = anga_body.x*Rx.x + anga_body.y*Rx.y + anga_body.z*Rx.z;
    anga_lab.y = anga_body.x*Ry.x + anga_body.y*Ry.y + anga_body.z*Ry.z;
    anga_lab.z = anga_body.x*Rz.x + anga_body.y*Rz.y + anga_body.z*Rz.z;
    
    /* Saves accelerations. */
    accel[ii]   = aa;
    ang_acc[ii] = anga_lab;

    SKIP:;
  }
  return;
}
