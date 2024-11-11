/***************************************************************************************************
 This routine rotates position vectors from the laboratory frame to the slope frame.
***************************************************************************************************/
#include "cabecera.hu"

/* =========================================== HOST ============================================= */
void rotate_vectors_hst (particle *pp, system_par sys)
{
  long    nn, ii;
  double  sint, cost;
  double3 rr, rrot, rr_a, rr_b, rr_c, rr_d, rrot_a, rrot_b, rrot_c, rrot_d;

  /* Fetches system parameters. */
  nn   = sys.npart;
  sint = sys.sint;
  cost = sys.cost;

  for (ii = 0; ii < nn; ii++)
  {
    /* Fetches particle vectors. */
    rr   = pp[ii].rr;
    rr_a = pp[ii].rr_a;
    rr_b = pp[ii].rr_b;
    rr_c = pp[ii].rr_c;
    rr_d = pp[ii].rr_d;

    /* Applies rotation matrix. */
    rrot.x = rr.x*cost - rr.z*sint;
    rrot.y = rr.y;
    rrot.z = rr.x*sint + rr.z*cost;

    rrot_a.x = rr_a.x*cost - rr_a.z*sint;
    rrot_a.y = rr_a.y;
    rrot_a.z = rr_a.x*sint + rr_a.z*cost;

    rrot_b.x = rr_b.x*cost - rr_b.z*sint;
    rrot_b.y = rr_b.y;
    rrot_b.z = rr_b.x*sint + rr_b.z*cost;

    rrot_c.x = rr_c.x*cost - rr_c.z*sint;
    rrot_c.y = rr_c.y;
    rrot_c.z = rr_c.x*sint + rr_c.z*cost;

    rrot_d.x = rr_d.x*cost - rr_d.z*sint;
    rrot_d.y = rr_d.y;
    rrot_d.z = rr_d.x*sint + rr_d.z*cost;

   /* Saves. */
    pp[ii].rrot   = rrot;
    pp[ii].rrot_a = rrot_a;
    pp[ii].rrot_b = rrot_b;
    pp[ii].rrot_c = rrot_c;
    pp[ii].rrot_d = rrot_d;
  }
  return;
}

/* ========================================== DEVICE ============================================ */
__global__ void rotate_vectors_dev (particle *pp, system_par sys)
{
  long    nn, ii;
  double  sint, cost;
  double3 rr, rrot, rr_a, rr_b, rr_c, rr_d, rrot_a, rrot_b, rrot_c, rrot_d;

  /* Fetches system parameters. */
  nn   = sys.npart;
  sint = sys.sint;
  cost = sys.cost;

  /* Thread index. */
  ii = threadIdx.x + blockDim.x*blockIdx.x;

  if (ii < nn)
  {
    /* Fetches particle vectors. */
    rr   = pp[ii].rr;
    rr_a = pp[ii].rr_a;
    rr_b = pp[ii].rr_b;
    rr_c = pp[ii].rr_c;
    rr_d = pp[ii].rr_d;

    /* Applies rotation matrix. */
    rrot.x = rr.x*cost - rr.z*sint;
    rrot.y = rr.y;
    rrot.z = rr.x*sint + rr.z*cost;

    rrot_a.x = rr_a.x*cost - rr_a.z*sint;
    rrot_a.y = rr_a.y;
    rrot_a.z = rr_a.x*sint + rr_a.z*cost;

    rrot_b.x = rr_b.x*cost - rr_b.z*sint;
    rrot_b.y = rr_b.y;
    rrot_b.z = rr_b.x*sint + rr_b.z*cost;

    rrot_c.x = rr_c.x*cost - rr_c.z*sint;
    rrot_c.y = rr_c.y;
    rrot_c.z = rr_c.x*sint + rr_c.z*cost;

    rrot_d.x = rr_d.x*cost - rr_d.z*sint;
    rrot_d.y = rr_d.y;
    rrot_d.z = rr_d.x*sint + rr_d.z*cost;

   /* Saves. */
    pp[ii].rrot   = rrot;
    pp[ii].rrot_a = rrot_a;
    pp[ii].rrot_b = rrot_b;
    pp[ii].rrot_c = rrot_c;
    pp[ii].rrot_d = rrot_d;
  }
  return;
}

