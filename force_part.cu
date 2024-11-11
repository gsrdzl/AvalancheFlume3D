/**************************************************************************************************
 This routine calculates the interaction forces and torques between particles. In this code, 'll'
 stands for the particle on which the force is been exerted and 'mm' stand for the other particles
 exerting the force on 'll'.
**************************************************************************************************/

#include "cabecera.hu"

/* =========================================== HOST ============================================ */

void force_part_hst (double3 *force, double3 *torque, long *cell_vec, int *nocup_vec, particle *pp,
                     system_par sys, cell_par cell, int sph_id)
{
  int     idel, jdel, kdel, nocup, smod, sdiv;
  long    ii, jj, kk, ll, mm, nn, cell_ind, tag, tag_init, tag_end, ntags, in, jn, kn;
  long3   ncell;
  double  dist, dist2, comp, comp_dot, kappa, gamma, mu, dvvcn_mag, dvvct_mag, normal, dt,
          aux_fric, rad_1, rad_2, sep, sep2, xi, xi_dot, ff_sdt, cns_1, t0, massl, massm, rhol,
          rhom, meff, km, gm, dvvct_sig, time;
  double3 rrl, rrm, rr_1, rr_2, rrc, drr, vvl, vvm, wwl, wwm, vvlc, vvmc, dvvc, dvvcn, dvvct, ff,
          ff_sph, tq, nvec, cell_side, tvec, rr_min, side;

  /* Fetches system parameters. */
  nn     = sys.npart;
  mu     = sys.mu_pp;
  kappa  = sys.kappa_pp;
  gamma  = sys.gamma_pp;
  dt     = sys.dt;
  cns_1  = sys.fr_fact_pp;
  time   = sys.time;
  side   = sys.side;
  
  /* Fetches cell parameters. */
  ncell     = cell.ncell;
  ntags     = cell.ntags;
  cell_side = cell.cell_side;
  rr_min    = cell.rr_min;
  
  /* Calculates some needed constants. */
  smod  = (sph_id - 1)%4;
  sdiv  = (sph_id - 1)/4;

  for (ll = 0; ll < nn; ll++)
  {
    /* Fetches first particle's position. */
    rrl = pp[ll].rr;
    if (time > 0.0 && rrl.x > side.x)
      continue;

    /* Fetches first particle's parameters. */
    vvl    = pp[ll].vv;
    wwl    = pp[ll].ww;
    rhol   = pp[ll].density;
    rad_1  = pp[ll].radius;
    ff     = force[ll];
    tq     = torque[ll];

    massl  = SPH_FACT*rhol*rad_1*rad_1*rad_1;

    switch (smod)
    {
      case 0:
        rr_1  = pp[ll].rr_a;
        break;
      case 1:
        rr_1  = pp[ll].rr_b;
        break;
      case 2:
        rr_1  = pp[ll].rr_c;
        break;
      case 3:
        rr_1  = pp[ll].rr_d;
        break;
    }

    /* Locates cells. */
    ii = (long) ((rrl.x - rr_min.x)/cell_side.x);
    if (ii == ncell.x)
      ii--;
    jj = (long) ((rrl.y - rr_min.y)/cell_side.y);
    if (jj == ncell.y)
      jj--;
    kk = (long) ((rrl.z - rr_min.z)/cell_side.z);
    if (kk == ncell.z)
      kk--;

    /* Checks the neighboring cells. */
    for (idel = -1; idel <= 1; idel++)
    for (jdel = -1; jdel <= 1; jdel++)
    for (kdel = -1; kdel <= 1; kdel++)
    {
      in = ii + idel;
      jn = jj + jdel;
      kn = kk + kdel;
      
      if (in >= 0 && jn >= 0 && kn >= 0 && in < ncell.x && jn < ncell.y && kn < ncell.z)
      {
        /* Index for linearized cell system. */
        cell_ind = in + ncell.x*(jn + ncell.y*kn);
        
        /* Ocupation numbers. */
        nocup    = nocup_vec[cell_ind];
        tag_init = cell_ind*ntags;
        tag_end  = tag_init + nocup;
        
        for (tag = tag_init; tag < tag_end; tag++)
        {
          mm = cell_vec[tag];

          if (mm != ll)
          {
            /* Fetches parameters of second particle. */
            rrm    = pp[mm].rr;
            rad_2  = pp[mm].radius;
            vvm    = pp[mm].vv;
            wwm    = pp[mm].ww;
            rhom   = pp[mm].density;
            
            massm = SPH_FACT*rhom*rad_2*rad_2*rad_2;
            meff  = 0.5*(massm + massl);
            km    = 0.5*kappa*meff;
            gm    = 0.5*gamma*meff;
            switch (sdiv)
            {
              case 0:
                rr_2  = pp[mm].rr_a;
                break;
              case 1:
                rr_2  = pp[mm].rr_b;
                break;
              case 2:
                rr_2  = pp[mm].rr_c;
                break;
              case 3:
                rr_2  = pp[mm].rr_d;
                break;
            }

            /* Separation between sphere centers at contact. */
            sep = rad_1 + rad_2;
            sep2 = sep*sep;
            /* Distance between spheres. */
            drr.x = rr_2.x - rr_1.x;
            dist2 = drr.x*drr.x;
            
            if (dist2 < sep2)
            {
              drr.y = rr_2.y - rr_1.y;
              dist2 += drr.y*drr.y;
                
              if (dist2 < sep2)
              {
                drr.z = rr_2.z - rr_1.z;
                dist2 += drr.z*drr.z;
                  
                /* Collision detected. */
                if (dist2 < sep2)
                {
                  dist = sqrt (dist2);

                  /* Compression. */
                  comp = sep - dist;
                    
                  /* Normal vector. */
                  nvec.x = drr.x/dist;
                  nvec.y = drr.y/dist;
                  nvec.z = drr.z/dist;

                  /* Contact point position. */
                  rrc.x = rr_1.x + rad_1*nvec.x;
                  rrc.y = rr_1.y + rad_1*nvec.y;
                  rrc.z = rr_1.z + rad_1*nvec.z;

                  /* Contact point velocities. */
                  vvlc.x = vvl.x + wwl.y*(rrc.z - rrl.z) - wwl.z*(rrc.y - rrl.y);
                  vvlc.y = vvl.y + wwl.z*(rrc.x - rrl.x) - wwl.x*(rrc.z - rrl.z);
                  vvlc.z = vvl.z + wwl.x*(rrc.y - rrl.y) - wwl.y*(rrc.x - rrl.x);

                  vvmc.x = vvm.x + wwm.y*(rrc.z - rrm.z) - wwm.z*(rrc.y - rrm.y);
                  vvmc.y = vvm.y + wwm.z*(rrc.x - rrm.x) - wwm.x*(rrc.z - rrm.z);
                  vvmc.z = vvm.z + wwm.x*(rrc.y - rrm.y) - wwm.y*(rrc.x - rrm.x);

                  /* Relative contact point velocity. */
                  dvvc.x = vvmc.x - vvlc.x;
                  dvvc.y = vvmc.y - vvlc.y;
                  dvvc.z = vvmc.z - vvlc.z;

                  /* Componentes of relative velocity. */
                  dvvcn_mag = dvvc.x*nvec.x + dvvc.y*nvec.y + dvvc.z*nvec.z;
                  dvvcn.x = dvvcn_mag*nvec.x;
                  dvvcn.y = dvvcn_mag*nvec.y;
                  dvvcn.z = dvvcn_mag*nvec.z;
                  dvvct.x = dvvc.x - dvvcn.x;
                  dvvct.y = dvvc.y - dvvcn.y;
                  dvvct.z = dvvc.z - dvvcn.z;
                  dvvct_mag = sqrt (dvvct.x*dvvct.x + dvvct.y*dvvct.y + dvvct.z*dvvct.z);

                  /* Tangential vectors. */
                  tvec.x    = dvvct.x/dvvct_mag;
                  tvec.y    = dvvct.y/dvvct_mag;
                  tvec.z    = dvvct.z/dvvct_mag;
                  dvvct_sig = dvvc.x*tvec.x + dvvc.y*tvec.y + dvvc.z*tvec.z;

                  /* Compression time derivative. */
                  comp_dot = -dvvcn_mag;
                  
                  /* Normal force. Procceeds only if normal force is repulsive. */
                  normal = km*comp + gm*comp_dot;
                  if (normal > 0.0)
                  {
                    ff_sph.x = -normal*nvec.x;
                    ff_sph.y = -normal*nvec.y;
                    ff_sph.z = -normal*nvec.z;

                    /* Friction force. This is the Cundall-Strack model. */
                    xi_dot = dvvct_sig;
                    xi     = dt*xi_dot;
                    if (cns_1*normal < fabs (xi))
                      xi = 0.0;
                    xi_dot = xi/dt;
                    ff_sdt = (1.0e12*massl*xi + 2.0e6*massl*xi_dot);
                    if (dvvct_mag > 0.0)
                    {
                      if (fabs (ff_sdt) > mu*normal)
                      {
                        aux_fric = 0.8*mu*normal/dvvct_mag;
                        t0 = massl/aux_fric;
                        aux_fric /= (1.0 + (dt/t0)*(dt/t0));
                        ff_sph.x += aux_fric*dvvct.x;
                        ff_sph.y += aux_fric*dvvct.y;
                        ff_sph.z += aux_fric*dvvct.z;
                      }
                      else 
                      {
                        if (dvvct.x != 0.0)
                          ff_sph.x += ff_sdt*dvvct.x/dvvct_mag;
                        if (dvvct.y != 0.0)
                          ff_sph.y += ff_sdt*dvvct.y/dvvct_mag;
                        if (dvvct.z != 0.0)
                          ff_sph.z += ff_sdt*dvvct.z/dvvct_mag;
                      }
                    }

                    /* Calculates torques. */
                    tq.x += (rrc.y - rrl.y)*ff_sph.z - (rrc.z - rrl.z)*ff_sph.y;
                    tq.y += (rrc.z - rrl.z)*ff_sph.x - (rrc.x - rrl.x)*ff_sph.z;
                    tq.z += (rrc.x - rrl.x)*ff_sph.y - (rrc.y - rrl.y)*ff_sph.x;

                    ff.x += ff_sph.x;
                    ff.y += ff_sph.y;
                    ff.z += ff_sph.z;
                  }
                }
              }
            }
          }
        }
      }
    }

    /* Saves forces and torques. */
    force[ll]  = ff;
    torque[ll] = tq;
  }
  return;
}

/* ========================================== DEVICE =========================================== */

__global__ void force_part_dev (double3 *force, double3 *torque, long *cell_vec, int *nocup_vec,
                                particle *pp, system_par sys, cell_par cell, int sph_id)
{
  int     idel, jdel, kdel, nocup, smod, sdiv;
  long    ii, jj, kk, ll, mm, nn, cell_ind, tag, tag_init, tag_end, ntags, in, jn, kn;
  long3   ncell;
  double  dist, dist2, comp, comp_dot, kappa, gamma, mu, dvvcn_mag, dvvct_mag, normal, dt,
          aux_fric, rad_1, rad_2, sep, sep2, xi, xi_dot, ff_sdt, cns_1, t0, massl, massm, rhol,
          rhom, meff, km, gm, dvvct_sig, time;
  double3 rrl, rrm, rr_1, rr_2, rrc, drr, vvl, vvm, wwl, wwm, vvlc, vvmc, dvvc, dvvcn, dvvct, ff,
          ff_sph, tq, nvec, cell_side, tvec, rr_min, side;

  /* Fetches system parameters. */
  nn     = sys.npart;
  mu     = sys.mu_pp;
  kappa  = sys.kappa_pp;
  gamma  = sys.gamma_pp;
  dt     = sys.dt;
  cns_1  = sys.fr_fact_pp;
  time   = sys.time;
  side   = sys.side;

  /* Fetches cell parameters. */
  ncell     = cell.ncell;
  ntags     = cell.ntags;
  cell_side = cell.cell_side;
  rr_min    = cell.rr_min;

  /* Calculates some needed constants. */
  smod  = (sph_id - 1)%4;
  sdiv  = (sph_id - 1)/4;

  /* Thread index. */
  ll = threadIdx.x + blockIdx.x*blockDim.x;

  if (ll < nn)
  {
    /* Fetches first particle's position. */
    rrl = pp[ll].rr;
    if (time > 0.0 && rrl.x > side.x)
      goto SKIP;

    /* Fetches first particle's parameters. */
    vvl    = pp[ll].vv;
    wwl    = pp[ll].ww;
    rhol   = pp[ll].density;
    rad_1  = pp[ll].radius;
    ff     = force[ll];
    tq     = torque[ll];

    massl  = SPH_FACT*rhol*rad_1*rad_1*rad_1;

    switch (smod)
    {
      case 0:
        rr_1  = pp[ll].rr_a;
        break;
      case 1:
        rr_1  = pp[ll].rr_b;
        break;
      case 2:
        rr_1  = pp[ll].rr_c;
        break;
      case 3:
        rr_1  = pp[ll].rr_d;
        break;
    }

    /* Locates cells. */
    ii = (long) ((rrl.x - rr_min.x)/cell_side.x);
    if (ii == ncell.x)
      ii--;
    jj = (long) ((rrl.y - rr_min.y)/cell_side.y);
    if (jj == ncell.y)
      jj--;
    kk = (long) ((rrl.z - rr_min.z)/cell_side.z);
    if (kk == ncell.z)
      kk--;

    /* Checks the neighboring cells. */
    for (idel = -1; idel <= 1; idel++)
    for (jdel = -1; jdel <= 1; jdel++)
    for (kdel = -1; kdel <= 1; kdel++)
    {
      in = ii + idel;
      jn = jj + jdel;
      kn = kk + kdel;
      
      if (in >= 0 && jn >= 0 && kn >= 0 && in < ncell.x && jn < ncell.y && kn < ncell.z)
      {
        /* Index for linearized cell system. */
        cell_ind = in + ncell.x*(jn + ncell.y*kn);
        
        /* Ocupation numbers. */
        nocup    = nocup_vec[cell_ind];
        tag_init = cell_ind*ntags;
        tag_end  = tag_init + nocup;
        
        for (tag = tag_init; tag < tag_end; tag++)
        {
          mm = cell_vec[tag];

          if (mm == ll)
            continue;

            /* Fetches parameters of second particle. */
            rrm    = pp[mm].rr;
            rad_2  = pp[mm].radius;
            vvm    = pp[mm].vv;
            wwm    = pp[mm].ww;
            rhom   = pp[mm].density;
            
            massm = SPH_FACT*rhom*rad_2*rad_2*rad_2;
            meff  = 0.5*(massm + massl);
            km    = 0.5*kappa*meff;
            gm    = 0.5*gamma*meff;
            switch (sdiv)
            {
              case 0:
                rr_2  = pp[mm].rr_a;
                break;
              case 1:
                rr_2  = pp[mm].rr_b;
                break;
              case 2:
                rr_2  = pp[mm].rr_c;
                break;
              case 3:
                rr_2  = pp[mm].rr_d;
                break;
            }

            /* Separation between sphere centers at contact. */
            sep = rad_1 + rad_2;
            sep2 = sep*sep;
            /* Distance between spheres. */
            drr.x = rr_2.x - rr_1.x;
            dist2 = drr.x*drr.x;
            
            if (dist2 < sep2)
            {
              drr.y = rr_2.y - rr_1.y;
              dist2 += drr.y*drr.y;
                
              if (dist2 < sep2)
              {
                drr.z = rr_2.z - rr_1.z;
                dist2 += drr.z*drr.z;
                  
                /* Collision detected. */
                if (dist2 < sep2)
                {
                  dist = sqrt (dist2);

                  /* Compression. */
                  comp = sep - dist;
                    
                  /* Normal vector. */
                  nvec.x = drr.x/dist;
                  nvec.y = drr.y/dist;
                  nvec.z = drr.z/dist;

                  /* Contact point position. */
                  rrc.x = rr_1.x + rad_1*nvec.x;
                  rrc.y = rr_1.y + rad_1*nvec.y;
                  rrc.z = rr_1.z + rad_1*nvec.z;

                  /* Contact point velocities. */
                  vvlc.x = vvl.x + wwl.y*(rrc.z - rrl.z) - wwl.z*(rrc.y - rrl.y);
                  vvlc.y = vvl.y + wwl.z*(rrc.x - rrl.x) - wwl.x*(rrc.z - rrl.z);
                  vvlc.z = vvl.z + wwl.x*(rrc.y - rrl.y) - wwl.y*(rrc.x - rrl.x);

                  vvmc.x = vvm.x + wwm.y*(rrc.z - rrm.z) - wwm.z*(rrc.y - rrm.y);
                  vvmc.y = vvm.y + wwm.z*(rrc.x - rrm.x) - wwm.x*(rrc.z - rrm.z);
                  vvmc.z = vvm.z + wwm.x*(rrc.y - rrm.y) - wwm.y*(rrc.x - rrm.x);

                  /* Relative contact point velocity. */
                  dvvc.x = vvmc.x - vvlc.x;
                  dvvc.y = vvmc.y - vvlc.y;
                  dvvc.z = vvmc.z - vvlc.z;

                  /* Componentes of relative velocity. */
                  dvvcn_mag = dvvc.x*nvec.x + dvvc.y*nvec.y + dvvc.z*nvec.z;
                  dvvcn.x = dvvcn_mag*nvec.x;
                  dvvcn.y = dvvcn_mag*nvec.y;
                  dvvcn.z = dvvcn_mag*nvec.z;
                  dvvct.x = dvvc.x - dvvcn.x;
                  dvvct.y = dvvc.y - dvvcn.y;
                  dvvct.z = dvvc.z - dvvcn.z;
                  dvvct_mag = norm3d (dvvct.x, dvvct.y, dvvct.z);

                  /* Tangential vectors. */
                  tvec.x    = dvvct.x/dvvct_mag;
                  tvec.y    = dvvct.y/dvvct_mag;
                  tvec.z    = dvvct.z/dvvct_mag;
                  dvvct_sig = dvvc.x*tvec.x + dvvc.y*tvec.y + dvvc.z*tvec.z;

                  /* Compression time derivative. */
                  comp_dot = -dvvcn_mag;
                  
                  /* Normal force. Procceeds only if normal force is repulsive. */
                  normal = km*comp + gm*comp_dot;
                  if (normal > 0.0)
                  {
                    ff_sph.x = -normal*nvec.x;
                    ff_sph.y = -normal*nvec.y;
                    ff_sph.z = -normal*nvec.z;

                    /* Friction force. This is the Cundall-Strack model. */
                    xi_dot = dvvct_sig;
                    xi     = dt*xi_dot;
                    if (cns_1*normal < fabs (xi))
                      xi = 0.0;
                    xi_dot = xi/dt;
                    ff_sdt = (1.0e12*massl*xi + 2.0e6*massl*xi_dot);
                    if (dvvct_mag > 0.0)
                    {
                      if (fabs (ff_sdt) > mu*normal || xi == 0.0)
                      {
                        aux_fric = 0.8*mu*normal/dvvct_mag;
                        t0 = massl/aux_fric;
                        aux_fric /= (1.0 + (dt/t0)*(dt/t0));
                        ff_sph.x += aux_fric*dvvct.x;
                        ff_sph.y += aux_fric*dvvct.y;
                        ff_sph.z += aux_fric*dvvct.z;
                      }
                      else 
                      {
                        if (dvvct.x != 0.0)
                          ff_sph.x += ff_sdt*dvvct.x/dvvct_mag;
                        if (dvvct.y != 0.0)
                          ff_sph.y += ff_sdt*dvvct.y/dvvct_mag;
                        if (dvvct.z != 0.0)
                          ff_sph.z += ff_sdt*dvvct.z/dvvct_mag;
                      }
                    }

                    /* Calculates torques. */
                    tq.x += (rrc.y - rrl.y)*ff_sph.z - (rrc.z - rrl.z)*ff_sph.y;
                    tq.y += (rrc.z - rrl.z)*ff_sph.x - (rrc.x - rrl.x)*ff_sph.z;
                    tq.z += (rrc.x - rrl.x)*ff_sph.y - (rrc.y - rrl.y)*ff_sph.x;

                    ff.x += ff_sph.x;
                    ff.y += ff_sph.y;
                    ff.z += ff_sph.z;
                  }
                }
              }
            }
          
        }
      }
    }

    /* Saves forces and torques. */
    force[ll]  = ff;
    torque[ll] = tq;
    
    SKIP:;
  }
  return;
}
