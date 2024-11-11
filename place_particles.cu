#include "cabecera.hu"

/**************************************************************************************************
 This routine places the particles created by 'generate_particles' into a cubic (squared) lattice
 and sets the initial orientations. This routine can only run in HOST.
**************************************************************************************************/

void place_particles (system_par sys, particle *pp)
{
  int     nn, aa, ii, jj, kk;
  double3 side, rr, Rx, Ry, Rz;
  double  rmax, dist, hdist, pos;
  
  /* Fetches system parameters. */
  nn   = sys.npart;
  side = sys.side;
  rmax = sys.rmax;
  
  /* Calculates some needed parameters. */
  dist  = 2.0 + ROOT2 + 0.05*rmax;
  hdist = dist/2.0;

  /* Initializes counters. */
  ii = jj = kk = 0;

  for (aa = 0; aa < nn; aa++)
  {
    /* Sets orientations. Initializes initial rotation matrix. */
    /* NOTE: Rotation matrices in this code are defined so that 'double3' vectors are matrix rows.
       This means that matrix A is represented as follows:
        -           -     -         -
       | Axx Axy Axz |   | double3 Ax |
       | Ayx Ayy Ayz | = | double3 Ay |
       | Azx Azy Azz |   | double3 Az |
        -           -     -         -
       where Axx = Ax.x, Axy = Ax.y and so on. */
    Rx.x = 1.0; Rx.y = 0.0; Rx.z = 0.0;
    Ry.x = 0.0; Ry.y = 1.0; Ry.z = 0.0;
    Rz.x = 0.0; Rz.y = 0.0; Rz.z = 1.0;
    
    /* Puts particles in a cubic lattice */
    rr.x = hdist*rmax + (double) ii*dist*rmax;
    rr.y = hdist*rmax + (double) jj*dist*rmax + 0.5*rmax + YHOP_L;
    rr.z = hdist*rmax + (double) kk*dist*rmax + 0.5*rmax;
    ii++;
    if (rr.x >= side.x - hdist*rmax)
    {
      rr.x  = hdist*rmax;
      rr.y += dist*rmax;
      jj++;
      ii = 1;
      if (rr.y >= YHOP_R - hdist*rmax)
      {
        rr.y  = hdist*rmax + 0.5*rmax + YHOP_L;
        rr.z += dist*rmax;
        kk++;
        jj = 0;
        if (rr.z >= side.z - hdist*rmax)
        {
          printf ("Particles don't fit!\n");
          exit (1);
        }
      }
    }
    rr.z += HOP_HEIGHT;
    pos = pp[aa].radius/ROOT2;

    /* Saves particle positions. */
    pp[aa].rr = rr;
    pp[aa].rr_a.x = rr.x + pos;
    pp[aa].rr_a.y = rr.y + pos;
    pp[aa].rr_a.z = rr.z + pos;
    pp[aa].rr_b.x = rr.x + pos;
    pp[aa].rr_b.y = rr.y - pos;
    pp[aa].rr_b.z = rr.z - pos;
    pp[aa].rr_c.x = rr.x - pos;
    pp[aa].rr_c.y = rr.y - pos;
    pp[aa].rr_c.z = rr.z + pos;
    pp[aa].rr_d.x = rr.x - pos;
    pp[aa].rr_d.y = rr.y + pos;
    pp[aa].rr_d.z = rr.z - pos;

    /* Saves orientations. */
    pp[aa].Rmatx = Rx;
    pp[aa].Rmaty = Ry;
    pp[aa].Rmatz = Rz;
  }

  return;
}
