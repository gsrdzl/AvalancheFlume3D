/***************************************************************************************************
 This routine calculates the dimensions of cell system, based on the maximum and minimum x, y and z
 coordinates of particles.
***************************************************************************************************/
#include "cabecera.hu"

void cell_dim (double3 *rr_max, double3 *rr_min, particle *pp, system_par sys)
{
  int     ii, nn;
  double  rgyr, dt, ddt, vxmax, vzmin;
  double3 rmax, rmin, rr, vv;

  nn   = sys.npart;
  rgyr = 0.5*sys.gyr_max;
  dt   = sys.dt;
  ddt  = 2000.0*dt;

  /* Initializes max and min positions with those of 1st particle. */
  rmax.x = rmin.x = pp[0].rr.x;
  rmax.y = rmin.y = pp[0].rr.y;
  rmax.z = rmin.z = pp[0].rr.z;

  /* Initializes max vx and min vz with those of 1st particle. */
  vxmax = pp[0].vv.x;
  vzmin = pp[0].vv.z;

  for (ii = 1; ii < nn; ii++)
  {
    /* Fetches particle positions and velocities. */
    rr = pp[ii].rr;
    vv = pp[ii].vv;

    /* Compares particle position with maximum and minimum positions. */
    if (rr.x < rmin.x)
      rmin.x = rr.x;
    if (rr.x > rmax.x)
      rmax.x = rr.x;
    if (rr.y < rmin.y)
      rmin.y = rr.y;
    if (rr.y > rmax.y)
      rmax.y = rr.y;
    if (rr.z < rmin.z)
      rmin.z = rr.z;
    if (rr.z > rmax.z)
      rmax.z = rr.z;

    /* Compares particle velocity x and z with max vx and min vy, respectively. */
    if (vv.x > vxmax)
      vxmax = vv.x;
    if (vv.z < vzmin)
      vzmin = vv.z;
  }

  /* Adds particle dimensions. */
  rmax.x += 1.05*rgyr + vxmax*ddt; 
  rmin.x -= 1.05*rgyr;
  rmax.y += 1.05*rgyr; 
  rmin.y -= 1.05*rgyr;
  rmax.z += 1.05*rgyr; 
  rmin.z -= 1.05*rgyr - vzmin*ddt;

  /* Saves. */
  *rr_max = rmax;
  *rr_min = rmin;
  return;
}
