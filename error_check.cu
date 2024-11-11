#include "cabecera.hu"

int error_check (system_par sys_par, char *err_msg)
{
  int errn = 0;

  if (sys_par.npart <= 0)
  {
    sprintf (err_msg, "ERROR. No particles in box! n = %ld\n", sys_par.npart);
    errn = 1;
  }

  if (sys_par.side.x <= 0.0 || sys_par.side.y <= 0.0 || sys_par.side.z <= 0.0)
  {
    sprintf (err_msg, "ERROR. Wrong box dimensions. x = %lf, y = %lf, z = %lf\n", sys_par.side.x,
            sys_par.side.y, sys_par.side.z);
    errn = 2;
  }

  if (sys_par.rmin <= 0.0 || sys_par.rmax <= 0.0)
  {
    sprintf (err_msg, "ERROR. Wrong sphere radius. rmax = %lf, rmin = %lf\n", sys_par.rmax,
             sys_par.rmin);
    errn = 3;
  }

  if (sys_par.rmax < sys_par.rmin)
  {
    sprintf (err_msg, "ERROR rmax must be larger than rmin. rmax = %lf, rmin = %lf\n", sys_par.rmax,
            sys_par.rmin);
    errn = 4;
  }

  if (sys_par.eps_pw <= 0.0 || sys_par.eps_pw > 1.0)
  {
    sprintf (err_msg, "ERROR. eps_pw must be 0 < eps_pw <= 1\n");
    errn = 5;
  }

  if (sys_par.eps_pp <= 0.0 || sys_par.eps_pp > 1.0)
  {
    sprintf (err_msg, "ERROR. eps_pp must be 0 < eps_pp <= 1\n");
    errn = 6;
  }

  if (sys_par.mu_pw < 0.0 || sys_par.mu_pw > 1.0)
  {
    sprintf (err_msg, "ERROR. mu_pw must be 0 < mu_pw < 1\n");
    errn = 7;
  }

  if (sys_par.mu_pp < 0.0 || sys_par.mu_pp > 1.5)
  {
    sprintf (err_msg, "ERROR. mu_pp must be 0 < mu_pp < 1\n");
    errn = 8;
  }

  if (tan (sys_par.tilt) < 1.25*sys_par.mu_pw)
  {
    sprintf (err_msg,
            "ERROR. Particles won't go down the slope.\n Please, increase tilt or reduce mu_pw.\n");
    errn = 9;
  }
  
  return errn;
}
