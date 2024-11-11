/**************************************************************************************************
 The routines within this module handle the locations and particle ocupations within the cell
 system. Routine 'clear_cell' empties the cells, routine 'clear_ocup' sets occupation number to
 zero, and routine 'cell_locate' places the particle into their corresponding cell.
**************************************************************************************************/

#include "cabecera.hu"

/* ============================================ HOST =========================================== */

void clear_cell_hst (cell_par cc, long *cell_vec)
{
  long ntot, ii;
  
  /* Fetches system parameters. */
  ntot = cc.ntot;
  
  for (ii = 0; ii < ntot; ii++)
    cell_vec[ii] = -1;
  
  return;
}

/*************************************************************************************************/

void clear_ocup_hst (cell_par cc, int *nocup_vec)
{
  long ii, nn;
  
  nn = cc.ncell3;
  
  for (ii = 0; ii < nn; ii++)
    nocup_vec[ii] = 0;
  
  return;
}

/*************************************************************************************************/

void cell_locate_hst (cell_par cc, system_par sys, particle *pp, int *nocup_vec, long *cell_vec)
{
  long    mm, ii, jj, kk, cell_index, nn, ntags;
  int     shift;
  long3   ncell;
  double  time;
  double3 rr, cell_side_inv, rr_min, side;

  /* Fetches system parameters. */
  nn   = sys.npart;
  time = sys.time;
  side = sys.side;
  
  /* Fetches cell parameters. */
  ntags         = cc.ntags;
  ncell         = cc.ncell;
  cell_side_inv = cc.cs_inv;
  rr_min        = cc.rr_min;
   
  /* Places particles into cells. */
  for (mm = 0; mm < nn; mm++)
  {
    /* Fetches particle positions. */
    rr = pp[mm].rr;
    
    if (time > 0.0 && rr.x > side.x)
      continue;

    /* Gets cell location. */
    ii = (long) (cell_side_inv.x*(rr.x - rr_min.x));
    if (ii == ncell.x)
      ii--;
    jj = (long) (cell_side_inv.y*(rr.y - rr_min.y));
    if (jj == ncell.y)
      jj--;
    kk = (long) (cell_side_inv.z*(rr.z - rr_min.z));
    if (kk == ncell.z)
      kk--;
     
    if (ii < 0 || jj < 0 || kk < 0 )
    {
      printf("error in cell_locate: mm = %ld\n", mm);
      printf("i, j, k, ncell_x, ncell_y, ncell_z %ld %ld %ld %ld %ld %ld\n", ii, jj, kk, ncell.x,
             ncell.y, ncell.z);
      printf ("r = (%lf,%lf,%lf); v = (%lf,%lf,%lf)\n", rr.x, rr.y, rr.z, pp[mm].vv.x, pp[mm].vv.y,
              pp[mm].vv.z);
      printf ("time = %lf\n", sys.time);
      exit (1);
    }

    /* Linearized cell vector */
    cell_index = ii + ncell.x*(jj + kk*ncell.y);
    
    shift = nocup_vec[cell_index];

    /* Checks errors and particle misplacements. */
    if (shift == ntags - 1)
    {
      printf ("Error in 'cell_locate': cell_index %ld  shift %d  ntags %ld\n", cell_index, shift,
              ntags);
      printf ("nocup %d\n", nocup_vec[cell_index]);
      exit (1);
    }

    /* Adds 1 to 'nocup'. */
    nocup_vec[cell_index]++;

    /* Saves particle index in cell vector. */
    cell_vec[cell_index*ntags + (long) shift] = mm;
  }
  return;
}

/* =========================================== DEVICE ========================================== */

__global__ void clear_cell_dev (cell_par cc, long *cell_vec)
{
  long ntot, ii;
  
  /* Fetches system parameters. */
  ntot = cc.ntot;
  
  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;
  
  if (ii < ntot)
    cell_vec[ii] = -1;
  
  return;
}

/*************************************************************************************************/

__global__ void clear_ocup_dev (cell_par cc, int *nocup_vec)
{
  long ii, nn;
  
  nn = cc.ncell3;
  
  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;

  if (ii < nn)
    nocup_vec[ii] = 0;
  
  return;
}

/*************************************************************************************************/

__global__ void cell_locate_dev (cell_par cc, system_par sys, particle *pp, int *nocup_vec,
                                 long *cell_vec)
{
  long    mm, ii, jj, kk, cell_index, nn, ntags;
  int     shift;
  long3   ncell;
  double  time;
  double3 rr, cell_side_inv, rr_min, side;

  /* Fetches system parameters. */
  nn   = sys.npart;
  time = sys.time;
  side = sys.side;
  /* Fetches cell parameters. */
  ntags         = cc.ntags;
  ncell         = cc.ncell;
  cell_side_inv = cc.cs_inv;
  rr_min        = cc.rr_min;
   
  /* Thread index. */
  mm = threadIdx.x + blockIdx.x*blockDim.x;

  /* Places particles into cells. */
  if (mm < nn)
  {
    /* Fetches particle positions. */
    rr = pp[mm].rr;
    if (time > 0.0 && rr.x > side.x)
      goto SKIP;
    
    /* Gets cell location. */
    ii = (long) (cell_side_inv.x*(rr.x - rr_min.x));
    if (ii == ncell.x)
      ii--;
    jj = (long) (cell_side_inv.y*(rr.y - rr_min.y));
    if (jj == ncell.y)
      jj--;
    kk = (long) (cell_side_inv.z*(rr.z - rr_min.z));
    if (kk == ncell.z)
      kk--;
     
    /* Linearized cell vector */
    cell_index = ii + ncell.x*(jj + kk*ncell.y);
    
    /* Adds 1 to 'nocup'. */
    shift = atomicAdd (&(nocup_vec[cell_index]), 1);

    /* Saves particle index in cell vector. */
    cell_vec[cell_index*ntags + (long) shift] = mm;
    SKIP:;
  }
  return;
}
