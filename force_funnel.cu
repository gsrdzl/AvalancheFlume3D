/**************************************************************************************************
 This routine calculates the interaction force between the particles and the walls of the funnel.
***************************************************************************************************/

#include "cabecera.hu"

/* ========================================== HOST ============================================= */

void force_funnel_hst (double3 *force, double3 *torque, system_par sys, particle *pp, Funnel fun)
{
  char    ns, touch_flag;
  long    nn, ii;
  double  kappa, gamma, dist, ri, dvvcn_mag, dvvct_mag, comp, comp_dot, normal, aux_fric, mu, dt,
          xi, xi_dot, ff_sdt, cns_1, t0, mass, rho, km,  gm, dvvct_sig, sinfi, cosfi, tanfi, cost,
          sint;
  double3 rr, vv, ww, ff, ff_sph, tq, rrc, rrotc, nvec, nvecr, rrot_i, vvc, dvvc, dvvcn,
          dvvct, tvec;
  
  /* Fetches system parameters. */
  nn    = sys.npart;
  cost  = sys.cost;
  sint  = sys.sint;
  mu    = sys.mu_pw;
  kappa = sys.kappa_pw;
  gamma = sys.gamma_pw;
  dt    = sys.dt;
  cns_1 = sys.fr_fact_pw;

  /* Fetches funnel parameters. */
  sinfi = fun.sinfi;
  cosfi = fun.cosfi;
  tanfi = fun.tanfi;

  for (ii = 0 ; ii < nn; ii++)
  {    
    /* Fetches particle parameters. */
    rr  = pp[ii].rr;
    ri  = pp[ii].radius;
    vv  = pp[ii].vv;
    ww  = pp[ii].ww;
    ff  = force[ii];
    tq  = torque[ii];
    rho = pp[ii].density;

    /* Calculates some needed parameters.  */
    mass = SPH_FACT*ri*ri*ri*rho;
    km   = 0.5*kappa*mass;
    gm   = 0.5*gamma*mass;

    /* Checks contact with all six walls. First checks the four spheres. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rrot_i = pp[ii].rrot_a;
          break;
        case 1:
          rrot_i = pp[ii].rrot_b;
          break;
        case 2:
          rrot_i = pp[ii].rrot_c;
          break;
        case 3:
          rrot_i = pp[ii].rrot_d;
          break;
      }
      /* If particle is beyond the funnel zone, does nothing. */
      if (rrot_i.x > LFUNNEL - ri*sinfi)
        continue;
      touch_flag = 0;

      /* Checks contact with left funnel wall. */
      dist = fabs (rrot_i.y - rrot_i.x*tanfi)*cosfi;
      if (dist <= ri)
      {
        touch_flag = 1;

        /* Contact point and normal vectors in the rotated frame. */
        rrotc.x  = rrot_i.x + dist*sinfi;
        rrotc.y  = rrot_i.y - dist*cosfi;
        rrotc.z  = rrot_i.z;
        nvecr.x  = sinfi;
        nvecr.y  = -cosfi;
        nvecr.z  = 0.0;

        /* Rotates vectors back to lab. system. */
        rrc.x  = rrotc.x*cost + rrotc.z*sint;
        rrc.y  = rrotc.y;
        rrc.z  = rrotc.z*cost - rrotc.x*sint;
        nvec.x = nvecr.x*cost + nvecr.z*sint;
        nvec.y = nvecr.y;
        nvec.z = nvecr.z*cost - nvecr.x*sint;

        /* Compression. */
        comp = ri - dist;
      }
      /* Checks contact with right funnel side. */
      dist = fabs (YFUNNEL - rrot_i.y - rrot_i.x*tanfi)*cosfi;
      if (dist <= ri)
      {
        touch_flag = 2;

        /* Contact point and normal vectors in the rotated frame. */
        rrotc.x  = rrot_i.x + dist*sinfi;
        rrotc.y  = rrot_i.y + dist*cosfi;
        rrotc.z  = rrot_i.z;
        nvecr.x = sinfi;
        nvecr.y = cosfi;
        nvecr.z = 0.0;

        /* Rotates vectors back to lab. system. */
        rrc.x  = rrotc.x*cost + rrotc.z*sint;
        rrc.y  = rrotc.y;
        rrc.z  = rrotc.z*cost - rrotc.x*sint;
        nvec.x = nvecr.x*cost + nvecr.z*sint;
        nvec.y = nvecr.y;
        nvec.z = nvecr.z*cost - nvecr.x*sint;

        /* Compression. */
        comp = ri - dist;
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

    //SKIP:;
  }
  return;
}

/**************************************************************************************************/
void force_fcorners_hst (double3 *force, double3 *torque, system_par sys, Funnel fun, particle *pp)
{
  char    ns;
  long    nn, ii;
  double  dist, dist2, comp, comp_dot, kappa, gamma, mu, dvvcn_mag, dvvct_mag, normal, dt, aux_fric,
          ri, yc, sep, sep2, xi, xi_dot, ff_sdt, cns_1, t0, km, gm, dvvct_sig, sin_fa, hsy, rho,
          yfun_r, yfun_l, mass, sint, cost;
  double3 rr, vv, ww, ff, ff_sph, tq, rrc, rrotc, nvec, nvecr, rrot_i, drr, vvc, dvvc, dvvcn,
          dvvct, tvec;
  
  /* Fetches system parameters. */
  nn     = sys.npart;
  mu     = sys.mu_pw;
  kappa  = sys.kappa_pw;
  gamma  = sys.gamma_pw;
  dt     = sys.dt;
  cns_1  = sys.fr_fact_pw;
  sin_fa = fun.sinfi;
  yfun_l = sys.dy2;
  cost   = sys.cost;
  sint   = sys.sint;

  yfun_r = sys.side.y + sys.dy2;
  hsy    = sys.dy2 + sys.side.y/2.0;

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

    /* Calculates needed parameters.*/
    mass = SPH_FACT*ri*ri*ri*rho;
    km   = 0.5*kappa*mass;
    gm   = 0.5*gamma*mass;

    /* Checks contact with the four spheres. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rrot_i = pp[ii].rrot_a;
          break;
        case 1:
          rrot_i = pp[ii].rrot_b;
          break;
        case 2:
          rrot_i = pp[ii].rrot_c;
          break;
        case 3:
          rrot_i = pp[ii].rrot_d;
          break;
      }

      /* If the sphere is not near the corner, do nothing. */
      if (rrot_i.x < LFUNNEL - ri*sin_fa || rrot_i.x > LFUNNEL + ri)
        continue;
      if (rrot_i.y > yfun_r || rrot_i.y < yfun_l)
        continue;

      /* Chooses interaction with either right or left corner. */
      yc = (rrot_i.y <= hsy) ? yfun_l : yfun_r;

      /* Distance to corner from sphere center. */
      sep   = ri;
      sep2  = sep*sep;
      drr.x = LFUNNEL - rrot_i.x;
      dist2 = drr.x*drr.x;
      if (dist2 > sep2)
        continue;
      drr.y  = yc - rrot_i.y;
      dist2 += drr.y*drr.y;
      if (dist2 > sep2)
        continue;

      /* If we get here, there is interaction. */
      dist = sqrt (dist2);

      /* Compression. */
      comp = sep - dist;

      /* Normal vector. */
      nvecr.x = drr.x/dist;
      nvecr.y = drr.y/dist;
      nvecr.z = 0.0;
      
      /* Contact point position vector. */
      rrotc.x = rrot_i.x + ri*nvecr.x;
      rrotc.y = rrot_i.y + ri*nvecr.y;
      rrotc.z = rrot_i.z;
      
      /* Rotates vectors back to lab. system. */
      rrc.x  = rrotc.x*cost + rrotc.z*sint;
      rrc.y  = rrotc.y;
      rrc.z  = rrotc.z*cost - rrotc.x*sint;
      nvec.x = nvecr.x*cost + nvecr.z*sint;
      nvec.y = nvecr.y;
      nvec.z = nvecr.z*cost - nvecr.x*sint;
      
      /* Contact point velocity on sphere. */
      vvc.x = vv.x + ww.y*(rrc.z - rr.z) - ww.z*(rrc.y - rr.y);
      vvc.y = vv.y + ww.z*(rrc.x - rr.x) - ww.x*(rrc.z - rr.z);
      vvc.z = vv.z + ww.x*(rrc.y - rr.y) - ww.y*(rrc.x - rr.x);
      
      /* Relative velocity. */
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
      if (normal <= 0.0)
        continue;
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

    /* Saves forces and torques. */
    force[ii]  = ff;
    torque[ii] = tq;
  }
  return;
}

/* ========================================= DEVICE ============================================= */

__global__ void force_funnel_dev (double3 *force, double3 *torque, system_par sys, particle *pp,
                                  Funnel fun)
{
  char    ns, touch_flag;
  long    nn, ii;
  double  kappa, gamma, dist, ri, dvvcn_mag, dvvct_mag, comp, comp_dot, normal, aux_fric, mu, dt,
          xi, xi_dot, ff_sdt, cns_1, t0, mass, rho, km,  gm, dvvct_sig, sinfi, cosfi, tanfi, cost,
          sint;
  double3 rr, vv, ww, ff, ff_sph, tq, rrc, rrotc, nvec, nvecr, rrot_i, vvc, dvvc, dvvcn,
          dvvct, tvec;
  
  /* Fetches system parameters. */
  nn    = sys.npart;
  cost  = sys.cost;
  sint  = sys.sint;
  mu    = sys.mu_pw;
  kappa = sys.kappa_pw;
  gamma = sys.gamma_pw;
  dt    = sys.dt;
  cns_1 = sys.fr_fact_pw;

  /* Fetches funnel parameters. */
  sinfi = fun.sinfi;
  cosfi = fun.cosfi;
  tanfi = fun.tanfi;
  
  /* Thread and block index. */
  ii = threadIdx.x + blockIdx.x*blockDim.x;

  if (ii < nn)
  {    
    /* Fetches particle parameters. */
    rr  = pp[ii].rr;
    ri  = pp[ii].radius;
    vv  = pp[ii].vv;
    ww  = pp[ii].ww;
    ff  = force[ii];
    tq  = torque[ii];
    rho = pp[ii].density;

    /* Calculates some needed parameters.  */
    mass = SPH_FACT*ri*ri*ri*rho;
    km   = 0.5*kappa*mass;
    gm   = 0.5*gamma*mass;

    /* Checks contact with all six walls. First checks the four spheres. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rrot_i = pp[ii].rrot_a;
          break;
        case 1:
          rrot_i = pp[ii].rrot_b;
          break;
        case 2:
          rrot_i = pp[ii].rrot_c;
          break;
        case 3:
          rrot_i = pp[ii].rrot_d;
          break;
      }
      /* If particle is beyond the funnel zone, does nothing. */
      if (rrot_i.x > LFUNNEL - ri*sinfi)
        continue;
      touch_flag = 0;

      /* Checks contact with left funnel wall. */
      dist = fabs (rrot_i.y - rrot_i.x*tanfi)*cosfi;
      if (dist <= ri)
      {
        touch_flag = 1;

        /* Contact point and normal vectors in the rotated frame. */
        rrotc.x  = rrot_i.x + dist*sinfi;
        rrotc.y  = rrot_i.y - dist*cosfi;
        rrotc.z  = rrot_i.z;
        nvecr.x  = sinfi;
        nvecr.y  = -cosfi;
        nvecr.z  = 0.0;

        /* Rotates vectors back to lab. system. */
        rrc.x  = rrotc.x*cost + rrotc.z*sint;
        rrc.y  = rrotc.y;
        rrc.z  = rrotc.z*cost - rrotc.x*sint;
        nvec.x = nvecr.x*cost + nvecr.z*sint;
        nvec.y = nvecr.y;
        nvec.z = nvecr.z*cost - nvecr.x*sint;

        /* Compression. */
        comp = ri - dist;
      }
      /* Checks contact with right funnel side. */
      dist = fabs (YFUNNEL - rrot_i.y - rrot_i.x*tanfi)*cosfi;
      if (dist <= ri)
      {
        touch_flag = 2;

        /* Contact point and normal vectors in the rotated frame. */
        rrotc.x  = rrot_i.x + dist*sinfi;
        rrotc.y  = rrot_i.y + dist*cosfi;
        rrotc.z  = rrot_i.z;
        nvecr.x = sinfi;
        nvecr.y = cosfi;
        nvecr.z = 0.0;

        /* Rotates vectors back to lab. system. */
        rrc.x  = rrotc.x*cost + rrotc.z*sint;
        rrc.y  = rrotc.y;
        rrc.z  = rrotc.z*cost - rrotc.x*sint;
        nvec.x = nvecr.x*cost + nvecr.z*sint;
        nvec.y = nvecr.y;
        nvec.z = nvecr.z*cost - nvecr.x*sint;

        /* Compression. */
        comp = ri - dist;
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

/**************************************************************************************************/
__global__ void force_fcorners_dev (double3 *force, double3 *torque, system_par sys, Funnel fun,
                                    particle *pp)
{
  char    ns;
  long    nn, ii;
  double  dist, dist2, comp, comp_dot, kappa, gamma, mu, dvvcn_mag, dvvct_mag, normal, dt, aux_fric,
          ri, yc, sep, sep2, xi, xi_dot, ff_sdt, cns_1, t0, km, gm, dvvct_sig, sin_fa, hsy, rho,
          yfun_r, yfun_l, mass, sint, cost;
  double3 rr, vv, ww, ff, ff_sph, tq, rrc, rrotc, nvec, nvecr, rrot_i, drr, vvc, dvvc, dvvcn,
          dvvct, tvec;
  
  /* Fetches system parameters. */
  nn     = sys.npart;
  mu     = sys.mu_pw;
  kappa  = sys.kappa_pw;
  gamma  = sys.gamma_pw;
  dt     = sys.dt;
  cns_1  = sys.fr_fact_pw;
  sin_fa = fun.sinfi;
  yfun_l = sys.dy2;
  cost   = sys.cost;
  sint   = sys.sint;

  yfun_r = sys.side.y + sys.dy2;
  hsy    = sys.dy2 + sys.side.y/2.0;

  /* Thread index and blocks. */
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

    /* Calculates needed parameters.*/
    mass = SPH_FACT*ri*ri*ri*rho;
    km   = 0.5*kappa*mass;
    gm   = 0.5*gamma*mass;

    /* Checks contact with the four spheres. */
    for (ns = 0; ns < 4; ns++)
    {
      switch (ns)
      {
        case 0:
          rrot_i = pp[ii].rrot_a;
          break;
        case 1:
          rrot_i = pp[ii].rrot_b;
          break;
        case 2:
          rrot_i = pp[ii].rrot_c;
          break;
        case 3:
          rrot_i = pp[ii].rrot_d;
          break;
      }

      /* If the sphere is not near the corner, do nothing. */
      if (rrot_i.x < LFUNNEL - ri*sin_fa || rrot_i.x > LFUNNEL + ri)
        continue;
      if (rrot_i.y > yfun_r || rrot_i.y < yfun_l)
        continue;

      /* Chooses interaction with either right or left corner. */
      yc = (rrot_i.y <= hsy) ? yfun_l : yfun_r;

      /* Distance to corner from sphere center. */
      sep   = ri;
      sep2  = sep*sep;
      drr.x = LFUNNEL - rrot_i.x;
      dist2 = drr.x*drr.x;
      if (dist2 > sep2)
        continue;
      drr.y  = yc - rrot_i.y;
      dist2 += drr.y*drr.y;
      if (dist2 > sep2)
        continue;

      /* If we get here, there is interaction. */
      dist = sqrt (dist2);

      /* Compression. */
      comp = sep - dist;

      /* Normal vector. */
      nvecr.x = drr.x/dist;
      nvecr.y = drr.y/dist;
      nvecr.z = 0.0;
      
      /* Contact point position vector. */
      rrotc.x = rrot_i.x + ri*nvecr.x;
      rrotc.y = rrot_i.y + ri*nvecr.y;
      rrotc.z = rrot_i.z;
      
      /* Rotates vectors back to lab. system. */
      rrc.x  = rrotc.x*cost + rrotc.z*sint;
      rrc.y  = rrotc.y;
      rrc.z  = rrotc.z*cost - rrotc.x*sint;
      nvec.x = nvecr.x*cost + nvecr.z*sint;
      nvec.y = nvecr.y;
      nvec.z = nvecr.z*cost - nvecr.x*sint;
      
      /* Contact point velocity on sphere. */
      vvc.x = vv.x + ww.y*(rrc.z - rr.z) - ww.z*(rrc.y - rr.y);
      vvc.y = vv.y + ww.z*(rrc.x - rr.x) - ww.x*(rrc.z - rr.z);
      vvc.z = vv.z + ww.x*(rrc.y - rr.y) - ww.y*(rrc.x - rr.x);
      
      /* Relative velocity. */
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
      if (normal <= 0.0)
        continue;
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

    /* Saves forces and torques. */
    force[ii]  = ff;
    torque[ii] = tq;
  }
  return;
}
