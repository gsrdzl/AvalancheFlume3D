/**************************************************************************************************
 This routine calculates the interaction force between the particles and the walls.
***************************************************************************************************/

#include "cabecera.hu"

/* ========================================== HOST ============================================= */

void force_hop_hst (double3 *force, double3 *torque, system_par sys, particle *pp, char nw)
{
  char    ns, touch_flag;
  long    nn, ii;
  double  xfront, xrear, yleft, yright, zbottom, kappa, gamma, dist, ri, dvvcn_mag, dvvct_mag,
          comp, comp_dot, normal, aux_fric, mu, dt, xi, xi_dot, ff_sdt, cns_1, t0, mass, rho, km,
          gm, dvvct_sig, time;
  double3 rr, vv, ww, ff, ff_sph, tq, rrc, nvec, rr_i, vvc, dvvc, dvvcn, dvvct, tvec;

  /* Fetches system parameters */
  nn    = sys.npart;
  mu    = sys.mu_pw;
  kappa = sys.kappa_pw;
  gamma = sys.gamma_pw;
  dt    = sys.dt;
  cns_1 = sys.fr_fact_pw;
  time  = sys.time;

  /* Calculates positions of hopper walls. */
  xfront  = HOP_WIDTH;
  xrear   = 0.0;
  yleft   = YHOP_L;
  yright  = YHOP_R;
  zbottom = HOP_HEIGHT;

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

    /* Checks contact with all six walls. First checks the four spheres. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rr_i = pp[ii].rr_a;
          break;
        case 1:
          rr_i = pp[ii].rr_b;
          break;
        case 2:
          rr_i = pp[ii].rr_c;
          break;
        case 3:
          rr_i = pp[ii].rr_d;
          break;
      }
      if (rr_i.z < zbottom || rr_i.x > HOP_WIDTH)
        continue;

      touch_flag = 0;
      switch (nw)
      {
        case 0:
          /* Checks contact with rear wall. */
          dist = rr_i.x - xrear;
          if (dist <= ri)
          {
            touch_flag = 1;
            rrc.x      = xrear;
            rrc.y      = rr_i.y;
            rrc.z      = rr_i.z;
            nvec.x     = -1.0;
            nvec.y     =  0.0;
            nvec.z     =  0.0;
            /* Compression. */
            comp = ri - dist;
          }
          /* Checks contact with front wall. */
          dist = xfront - rr_i.x;
          if (dist <= ri)
          {
            touch_flag = 2;
            rrc.x      = xfront;
            rrc.y      = rr_i.y;
            rrc.z      = rr_i.z;
            nvec.x     = 1.0;
            nvec.y     = 0.0;
            nvec.z     = 0.0;
            /* Compression. */
            comp = ri - dist;
          }
          break;

        case 1:
          /* Checks contact with left wall. */
          dist = rr_i.y - yleft;
          if (dist <= ri)
          {
            touch_flag = 3;
            rrc.x      = rr_i.x;
            rrc.y      = yleft;
            rrc.z      = rr_i.z;
            nvec.x = nvec.z = 0.0;
            nvec.y = -1.0;

            /* Compression. */
            comp = ri - dist;
          }
          /* Checks contact with right wall. */
          dist = yright - rr_i.y;
          if (dist <= ri)
          {
            touch_flag = 4;
            rrc.x      = rr_i.x;
            rrc.y      = yright;
            rrc.z      = rr_i.z;
            nvec.x = nvec.z = 0.0;
            nvec.y = 1.0;

            /* Compression. */
            comp = ri - dist;
          }
          break;

        case 2:
          /* Checks contact with bottom wall. */
          dist = rr_i.z - zbottom;
          if (dist <= ri)
          {
            if (time < 0.0)
              touch_flag = 5;
            else
              touch_flag = 0; // After transient, hopper's bottom door is opened. 
            rrc.x      = rr_i.x;
            rrc.y      = rr_i.y;
            rrc.z      = zbottom;

            nvec.x = nvec.y = 0.0;
            nvec.z = -1.0;

            /* Compression. */
            comp = ri - dist;
          }
          break;
      }

      if (touch_flag)
      {
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
              aux_fric  = 0.8*mu*normal/dvvct_mag;
              t0        = mass/aux_fric;
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
    }

    /* Saves forces and torques. */
    force[ii]  = ff;
    torque[ii] = tq;
  }
  return;
}

/* ========================================= DEVICE ============================================ */

__global__ void force_hop_dev (double3 *force, double3 *torque, system_par sys, particle *pp,
                               char nw)
{
  char    ns, touch_flag;
  long    nn, ii;
  double  xfront, xrear, yleft, yright, zbottom, kappa, gamma, dist, ri, dvvcn_mag, dvvct_mag,
          comp, comp_dot, normal, aux_fric, mu, dt, xi, xi_dot, ff_sdt, cns_1, t0, mass, rho, km,
          gm, dvvct_sig, time;
  double3 rr, vv, ww, ff, ff_sph, tq, rrc, nvec, rr_i, vvc, dvvc, dvvcn, dvvct, tvec;

  /* Fetches system parameters */
  nn    = sys.npart;
  mu    = sys.mu_pw;
  kappa = sys.kappa_pw;
  gamma = sys.gamma_pw;
  dt    = sys.dt;
  cns_1 = sys.fr_fact_pw;
  time  = sys.time;

  /* Calculates positions of hopper walls. */
  xfront  = HOP_WIDTH;
  xrear   = 0.0;
  yleft   = YHOP_L;
  yright  = YHOP_R;
  zbottom = HOP_HEIGHT;

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

    /* Checks contact with all six walls. First checks the four spheres. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rr_i = pp[ii].rr_a;
          break;
        case 1:
          rr_i = pp[ii].rr_b;
          break;
        case 2:
          rr_i = pp[ii].rr_c;
          break;
        case 3:
          rr_i = pp[ii].rr_d;
          break;
      }
      if (rr_i.z < zbottom || rr_i.x > HOP_WIDTH)
        continue;

      touch_flag = 0;
      switch (nw)
      {
        case 0:
          /* Checks contact with rear wall. */
          dist = rr_i.x - xrear;
          if (dist <= ri)
          {
            touch_flag = 1;
            rrc.x      = xrear;
            rrc.y      = rr_i.y;
            rrc.z      = rr_i.z;
            nvec.x     = -1.0;
            nvec.y     =  0.0;
            nvec.z     =  0.0;
            /* Compression. */
            comp = ri - dist;
          }
          /* Checks contact with front wall. */
          dist = xfront - rr_i.x;
          if (dist <= ri)
          {
            touch_flag = 2;
            rrc.x      = xfront;
            rrc.y      = rr_i.y;
            rrc.z      = rr_i.z;
            nvec.x     = 1.0;
            nvec.y     = 0.0;
            nvec.z     = 0.0;
            /* Compression. */
            comp = ri - dist;
          }
          break;

        case 1:
          /* Checks contact with left wall. */
          dist = rr_i.y - yleft;
          if (dist <= ri)
          {
            touch_flag = 3;
            rrc.x      = rr_i.x;
            rrc.y      = yleft;
            rrc.z      = rr_i.z;
            nvec.x = nvec.z = 0.0;
            nvec.y = -1.0;

            /* Compression. */
            comp = ri - dist;
          }
          /* Checks contact with right wall. */
          dist = yright - rr_i.y;
          if (dist <= ri)
          {
            touch_flag = 4;
            rrc.x      = rr_i.x;
            rrc.y      = yright;
            rrc.z      = rr_i.z;
            nvec.x = nvec.z = 0.0;
            nvec.y = 1.0;

            /* Compression. */
            comp = ri - dist;
          }
          break;

        case 2:
          /* Checks contact with bottom wall. */
          dist = rr_i.z - zbottom;
          if (dist <= ri)
          {
            if (time < 0.0)
              touch_flag = 5;
            else
              touch_flag = 0; // After transient, hopper's bottom door is opened. 
            rrc.x      = rr_i.x;
            rrc.y      = rr_i.y;
            rrc.z      = zbottom;

            nvec.x = nvec.y = 0.0;
            nvec.z = -1.0;

            /* Compression. */
            comp = ri - dist;
          }
          break;
      }

      if (touch_flag)
      {
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
              aux_fric  = 0.8*mu*normal/dvvct_mag;
              t0        = mass/aux_fric;
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
    }

    /* Saves forces and torques. */
    force[ii]  = ff;
    torque[ii] = tq;
  }
  return;
}
