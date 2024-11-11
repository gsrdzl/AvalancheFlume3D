/**************************************************************************************************
 This routine calculates the interaction force between the particles and the floor of the
 deposition zone.
**************************************************************************************************/

#include "cabecera.hu"

/* ========================================== HOST ============================================= */

void force_dep_hst (double3 *force, double3 *torque, system_par sys, particle *pp)
{
  int    nn, ii, ns;
  double  kappa, gamma, dist, ri, dvvcn_mag, dvvct_mag, comp, comp_dot, normal, aux_fric, mu, dt,
          xi, xi_dot, ff_sdt, sint, cost, tant, rect, Lf, xh, cns_1, mass, km, gm, rho, t0,
          dvvct_sig;
  double3 rr, rr_a, rr_b, rr_c, rr_d, vv, ww, ff, ff_sph, tq, rrc, nvec, rr_i, vvc, dvvc, dvvcn,
          dvvct, tvec;

  /* Fetches system parameters */
  nn     = sys.npart;
  mu     = sys.mu_pw;
  kappa  = sys.kappa_pw;
  gamma  = sys.gamma_pw;
  dt     = sys.dt;
  Lf     = sys.Lflume;
  xh     = sys.xhopper;
  sint   = sys.sint;
  cost   = sys.cost;
  cns_1  = sys.fr_fact_pw;
  tant   = sint/cost;

  for (ii = 0; ii < nn; ii++)
  {
    /* Fetches particle parameters */
    rr  = pp[ii].rr;
    vv  = pp[ii].vv;
    ww  = pp[ii].ww;
    ff  = force[ii];
    tq  = torque[ii];
    ri  = pp[ii].radius;
    rho = pp[ii].density;

    mass = SPH_FACT*ri*ri*ri*rho;
    km   = 0.5*kappa*mass;
    gm   = 0.5*gamma*mass;

    /* Position of dimer's spheres. */
    rr_a = pp[ii].rr_a;
    rr_b = pp[ii].rr_b;
    rr_c = pp[ii].rr_c;
    rr_d = pp[ii].rr_d;

    /* Checks contact between each of the two spheres and the floor. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rr_i = rr_a;
          break;
        case 1:
          rr_i = rr_b;
          break;
        case 2:
          rr_i = rr_c;
          break;
        case 3:
          rr_i = rr_d;
          break;
      }

      /* Distance to the floor. */
      rect = tant*(rr_i.x - Lf - xh) - rr_i.z;
      dist = cost*(-rect);
      
      if (dist <= ri)
      {
        /* Normal vector. */
        nvec.x = sint;
        nvec.y = 0.0;
        nvec.z = -cost;

        /* Compression. */
        comp = ri - dist;

        /* Contact points. */
        rrc.x = rr_i.x + dist*sint;
        rrc.y = rr_i.y;
        rrc.z = tant*(rrc.x - Lf - xh);

        /* Contact point velocities. */
        vvc.x = vv.x + ww.y*(rrc.z - rr.z) - ww.z*(rrc.y - rr.y);
        vvc.y = vv.y + ww.z*(rrc.x - rr.x) - ww.x*(rrc.z - rr.z);
        vvc.z = vv.z + ww.x*(rrc.y - rr.y) - ww.y*(rrc.x - rr.x);

        /* Relative velocity at contact point. */
        dvvc.x = -vvc.x;
        dvvc.y = -vvc.y;
        dvvc.z = -vvc.z;

        /* Components of relative velocity. */
        dvvcn_mag = dvvc.x*nvec.x + dvvc.y*nvec.y + dvvc.z*nvec.z;
        dvvcn.x   = dvvcn_mag*nvec.x;
        dvvcn.y   = dvvcn_mag*nvec.y;
        dvvcn.z   = dvvcn_mag*nvec.z;
        dvvct.x   = dvvc.x - dvvcn.x;
        dvvct.y   = dvvc.y - dvvcn.y;
        dvvct.z   = dvvc.z - dvvcn.z;
        dvvct_mag = sqrt (dvvct.x*dvvct.x + dvvct.y*dvvct.y + dvvct.z*dvvct.z);

        /* Tangential vectors. */
        tvec.x    = dvvct.x/dvvct_mag;
        tvec.y    = dvvct.y/dvvct_mag;
        tvec.z    = dvvct.z/dvvct_mag;
        dvvct_sig = dvvc.x*tvec.x + dvvc.y*tvec.y + dvvc.z*tvec.z;

        /* Time derivative of compression. */
        comp_dot = -dvvcn_mag;

        /* Normal force. */
        normal = km*comp + gm*comp_dot;

        /* Normal force is always repulsive. */
        if (normal > 0.0)
        {
          ff_sph.x = -normal*nvec.x;
          ff_sph.y = -normal*nvec.y;
          ff_sph.z = -normal*nvec.z;

          /* Friction force. */
          xi_dot = dvvct_sig;
          xi     = dt*xi_dot;
          if (cns_1*normal < fabs (xi))
            xi = 0.0;
          xi_dot = xi/dt;
          ff_sdt = (1.0e12*mass*xi + 2.0e6*mass*xi_dot);
          if (dvvct_mag > 0.0)
          {
            if (fabs (ff_sdt) > mu*normal || xi == 0.0)
            {
              aux_fric = 0.8*mu*normal/dvvct_mag;
              t0 = mass/aux_fric;
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

          /* Torque. */
          tq.x += (rrc.y - rr.y)*ff_sph.z - (rrc.z - rr.z)*ff_sph.y;
          tq.y += (rrc.z - rr.z)*ff_sph.x - (rrc.x - rr.x)*ff_sph.z;
          tq.z += (rrc.x - rr.x)*ff_sph.y - (rrc.y - rr.y)*ff_sph.x;

          ff.x += ff_sph.x;
          ff.y += ff_sph.y;
          ff.z += ff_sph.z;

        }
      }

    /* Saves forces and torques. */
    force[ii]  = ff;
    torque[ii] = tq;
    }
  }

  return;
}

/* ========================================= DEVICE ============================================ */

__global__ void force_dep_dev (double3 *force, double3 *torque, system_par sys, particle *pp)
{
  int     nn, ii, ns;
  double  kappa, gamma, dist, ri, dvvcn_mag, dvvct_mag, comp, comp_dot, normal, aux_fric, mu, dt,
          xi, xi_dot, ff_sdt, sint, cost, tant, rect, Lf, xh, cns_1, mass, km, gm, rho, t0,
          dvvct_sig;
  double3 rr, rr_a, rr_b, rr_c, rr_d, vv, ww, ff, ff_sph, tq, rrc, nvec, rr_i, vvc, dvvc, dvvcn,
          dvvct, tvec;

  /* Fetches system parameters */
  nn     = sys.npart;
  mu     = sys.mu_pw;
  kappa  = sys.kappa_pw;
  gamma  = sys.gamma_pw;
  dt     = sys.dt;
  Lf     = sys.Lflume;
  xh     = sys.xhopper;
  sint   = sys.sint;
  cost   = sys.cost;
  cns_1  = sys.fr_fact_pw;
  tant   = sint/cost;

  /* Thread index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;

  if (ii < nn)
  {
    /* Fetches particle parameters */
    rr  = pp[ii].rr;
    vv  = pp[ii].vv;
    ww  = pp[ii].ww;
    ff  = force[ii];
    tq  = torque[ii];
    ri  = pp[ii].radius;
    rho = pp[ii].density;

    mass = SPH_FACT*ri*ri*ri*rho;
    km   = 0.5*kappa*mass;
    gm   = 0.5*gamma*mass;

    /* Position of dimer's spheres. */
    rr_a = pp[ii].rr_a;
    rr_b = pp[ii].rr_b;
    rr_c = pp[ii].rr_c;
    rr_d = pp[ii].rr_d;

    /* Checks contact between each of the two spheres and the floor. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rr_i = rr_a;
          break;
        case 1:
          rr_i = rr_b;
          break;
        case 2:
          rr_i = rr_c;
          break;
        case 3:
          rr_i = rr_d;
          break;
      }

      /* Distance to the floor. */
      rect = tant*(rr_i.x - Lf - xh) - rr_i.z;
      dist = cost*(-rect);
      
      if (dist <= ri)
      {
        /* Normal vector. */
        nvec.x = sint;
        nvec.y = 0.0;
        nvec.z = -cost;

        /* Compression. */
        comp = ri - dist;

        /* Contact points. */
        rrc.x = rr_i.x + dist*sint;
        rrc.y = rr_i.y;
        rrc.z = tant*(rrc.x - Lf - xh);

        /* Contact point velocities. */
        vvc.x = vv.x + ww.y*(rrc.z - rr.z) - ww.z*(rrc.y - rr.y);
        vvc.y = vv.y + ww.z*(rrc.x - rr.x) - ww.x*(rrc.z - rr.z);
        vvc.z = vv.z + ww.x*(rrc.y - rr.y) - ww.y*(rrc.x - rr.x);

        /* Relative velocity at contact point. */
        dvvc.x = -vvc.x;
        dvvc.y = -vvc.y;
        dvvc.z = -vvc.z;

        /* Components of relative velocity. */
        dvvcn_mag = dvvc.x*nvec.x + dvvc.y*nvec.y + dvvc.z*nvec.z;
        dvvcn.x   = dvvcn_mag*nvec.x;
        dvvcn.y   = dvvcn_mag*nvec.y;
        dvvcn.z   = dvvcn_mag*nvec.z;
        dvvct.x   = dvvc.x - dvvcn.x;
        dvvct.y   = dvvc.y - dvvcn.y;
        dvvct.z   = dvvc.z - dvvcn.z;
        dvvct_mag = norm3d (dvvct.x, dvvct.y, dvvct.z);

        /* Tangential vectors. */
        tvec.x    = dvvct.x/dvvct_mag;
        tvec.y    = dvvct.y/dvvct_mag;
        tvec.z    = dvvct.z/dvvct_mag;
        dvvct_sig = dvvc.x*tvec.x + dvvc.y*tvec.y + dvvc.z*tvec.z;

        /* Time derivative of compression. */
        comp_dot = -dvvcn_mag;

        /* Normal force. */
        normal = km*comp + gm*comp_dot;

        /* Normal force is always repulsive. */
        if (normal > 0.0)
        {
          ff_sph.x = -normal*nvec.x;
          ff_sph.y = -normal*nvec.y;
          ff_sph.z = -normal*nvec.z;

          /* Friction force. */
          xi_dot = dvvct_sig;
          xi     = dt*xi_dot;
          if (cns_1*normal < fabs (xi))
            xi = 0.0;
          xi_dot = xi/dt;
          ff_sdt = (1.0e12*mass*xi + 2.0e6*mass*xi_dot);
          if (dvvct_mag > 0.0)
          {
            if (fabs (ff_sdt) > mu*normal || xi == 0.0)
            {
              aux_fric = 0.8*mu*normal/dvvct_mag;
              t0 = mass/aux_fric;
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

          /* Torque. */
          tq.x += (rrc.y - rr.y)*ff_sph.z - (rrc.z - rr.z)*ff_sph.y;
          tq.y += (rrc.z - rr.z)*ff_sph.x - (rrc.x - rr.x)*ff_sph.z;
          tq.z += (rrc.x - rr.x)*ff_sph.y - (rrc.y - rr.y)*ff_sph.x;

          ff.x += ff_sph.x;
          ff.y += ff_sph.y;
          ff.z += ff_sph.z;

        }
      }

    /* Saves forces and torques. */
    force[ii]  = ff;
    torque[ii] = tq;
    }
  }

  return;
}
