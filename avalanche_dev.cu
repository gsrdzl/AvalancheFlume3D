#include "cabecera.hu"
#define FIRST_TIME 1
#define DEVICE     1

int main ()
{
  int        errn, aval_start, print_flag;
  long       nn, ii, jj, cc, ncell3, ntot, ntags, n3max, ntotmax;
  long3      ncell, ncellmax;
  char       err_msg[100], row[100], *cp_in;
  double     vmean, vm2, wmean, wm2, dist, aux_0, aux_1, xfinal, xcm, zcm, xmass, zmass, hipof,
             vproj;
  double3    cell_side, cs_inv;
  FILE       *fp_log, *fp_xyz, *fp_vmean, *fp_vel, *fp_for, *fp_conf_out, *fp_deb, *fp_xcm;
  system_par sys_par;
  Funnel     fun;
  cell_par   cell;

  /* Vectors in HOST. */
  long       *cell_vec_hst;
  int        *nocup_vec_hst;
  double3    *ff_hst, *torque_hst, *aa_i_hst, *anga_i_hst, *aa_f_hst, *anga_f_hst, *swap_hst;
  particle   *pp_par_hst;
#if DEVICE
  /* Vectors in DEVICE. */
  long       *cell_vec_dev;
  int        *nocup_vec_dev;
  double3    *ff_dev, *torque_dev, *aa_i_dev, *anga_i_dev, *aa_f_dev, *anga_f_dev, *swap_dev;
  particle   *pp_par_dev;
#endif
  fp_log = fopen ("log_dev.dat", "w");
  fp_deb = fopen ("debug.csv", "w");

  /* Scans system parameters from input data file. */
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%ld", &(sys_par.npart));     /* Number of particles. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.side.y));    /* Flume and hopper width. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.side.z));    /* Hopper height. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.rmin));      /* Minimum sphere radius. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.rmax));      /* Maximum sphere radius. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.v0));        /* Initial velocity. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.dt));        /* Time step. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.eps_pw));    /* Particle-wall coefficient of restitution. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.mu_pw));     /* Particle-wall coefficient of friction. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.eps_pp));    /* Particle-particle coefficient of restitution. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.mu_pp));     /* Particle-particle coefficient of friction. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.t_coll));    /* Collision time. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.Lflume));    /* Flume length. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%lf", &(sys_par.tilt));      /* Flume tilt in degrees. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%ld", &(cell.ntags));        /* Number of tagged particles in cells. */
  fprintf (fp_log, "%s", cp_in);
  cp_in = fgets (row, sizeof (row), stdin);
  sscanf (row, "%d", &print_flag);          /* Flag to print .xyz file 1 or 0. */
  fprintf (fp_log, "%s\n", cp_in);
  
  /* Other initial parameters. */
  sys_par.side.x  = HOP_WIDTH;
  sys_par.time    = -1.0;
  sys_par.hdt     = 0.5*sys_par.dt;
  sys_par.tilt   *= PI/180.0;
  sys_par.gyr_max = GYR*sys_par.rmax;
  sys_par.sint    = sin (sys_par.tilt);
  sys_par.cost    = cos (sys_par.tilt);
  sys_par.sintg   = GRAV*sys_par.sint;
  sys_par.costg   = GRAV*sys_par.cost;
  sys_par.xfunnel = LFUNNEL*sys_par.cost;
  sys_par.dy2     = 0.5*(YFUNNEL - sys_par.side.y);
  nn    = sys_par.npart;
  ntags = cell.ntags;

  /* Funnel parameters. */
  hipof     = sqrt (sys_par.dy2*sys_par.dy2 + LFUNNEL*LFUNNEL);
  fun.sinfi = sys_par.dy2/hipof;
  fun.cosfi = LFUNNEL/hipof;
  fun.tanfi = sys_par.dy2/LFUNNEL;
  fprintf (fp_log, "sinf: %lf, cosfi %lf, tanfi %lf\n", fun.sinfi, fun.cosfi, fun.tanfi);

  /* Checks errors.*/
  errn = error_check (sys_par, err_msg);
  if (errn)
  {
    fprintf (stderr, "%sError No. %d.\n", err_msg, errn);
    fprintf (fp_log, "%sError No. %d.\n", err_msg, errn);
    exit (1);
  }

  /* Calculates some needed parameters. */
  dist   = GYR + 0.05f;
  xfinal = (LFUNNEL + sys_par.Lflume)*sys_par.cost;

  fprintf (fp_log, "xfinal: %lf\n", xfinal);
  fflush  (fp_log);

  /* Calculates cell parameters and saves them. */
  ncellmax.x = (long) (xfinal/(dist*sys_par.rmax));
  ncellmax.y = (long) (sys_par.side.y/(dist*sys_par.rmax));
  ncellmax.z = (long) (1.5*sys_par.side.z/(dist*sys_par.rmax));
  n3max      = /*50000000;*/ncellmax.x*ncellmax.y*ncellmax.z;
  ntotmax    = /*400000000;*/n3max*ntags;

  /* Prints parameters out to log. */
  fprintf (fp_log, "ncellmax: (%ld,%ld,%ld).\n", ncellmax.x, ncellmax.y, ncellmax.z);
  fprintf (fp_log, "n3max: %ld, ntotmax: %ld.\n", n3max, ntotmax);

  /* Memory allocation for vectors in HOST. */
  ff_hst     = (double3 *)   malloc (nn*sizeof (double3));
  torque_hst = (double3 *)   malloc (nn*sizeof (double3));
  aa_i_hst   = (double3 *)   malloc (nn*sizeof (double3));
  anga_i_hst = (double3 *)   malloc (nn*sizeof (double3));
  aa_f_hst   = (double3 *)   malloc (nn*sizeof (double3));
  anga_f_hst = (double3 *)   malloc (nn*sizeof (double3));
  swap_hst   = (double3 *)   malloc (nn*sizeof (double3));
  pp_par_hst = (particle *) malloc (nn*sizeof (particle));
#if DEVICE
  /* Memory allocation for vectors in DEVICE. */
  cudaMalloc ((void **) &ff_dev,     nn*sizeof (double3));
  cudaMalloc ((void **) &torque_dev, nn*sizeof (double3));
  cudaMalloc ((void **) &aa_i_dev,   nn*sizeof (double3));
  cudaMalloc ((void **) &anga_i_dev, nn*sizeof (double3));
  cudaMalloc ((void **) &aa_f_dev,   nn*sizeof (double3));
  cudaMalloc ((void **) &anga_f_dev, nn*sizeof (double3));
  cudaMalloc ((void **) &swap_dev,   nn*sizeof (double3));
  cudaMalloc ((void **) &pp_par_dev, nn*sizeof (particle));

  fprintf (fp_deb, "Vector, size (B)\n");
  fprintf (fp_deb, "double3, %ld\n", sizeof (double3));
  fprintf (fp_deb, "long, %ld\n", sizeof (long));
  fprintf (fp_deb, "int, %ld\n", sizeof (int));
  fprintf (fp_deb, "particle, %ld\n", sizeof (particle));
  fflush  (fp_deb);

  /* Threads and blocks. */
  int  nthread = 256;
  long RB_part = 1 + (nn - 1)/((long) nthread);
#endif
  /* Randomly generates the particles. */
  printf ("Creating particles...");
  generate_particles (pp_par_hst, sys_par);
  printf ("\n");

  /* Calculates particle masses, centers of mass and moments of inertia. */
  printf ("Calculating particle parameters...");
  calculate_particle_parameters_hst (sys_par, pp_par_hst);
  printf ("\n");

  /* Calculates total system mass. Saves particle data in log. */
  sys_par.tot_mass = 0.0;
  for (ii = 0; ii < nn; ii++)
    sys_par.tot_mass += pp_par_hst[ii].mass;
  fprintf (fp_log, "Total mass: %lf\n", sys_par.tot_mass); 

  /* Calculates other parameters. */
  aux_0              = PI/sys_par.t_coll;
  aux_1              = log (sys_par.eps_pp)/sys_par.t_coll;
  sys_par.gamma_pp   = -2.0*log (sys_par.eps_pp)/sys_par.t_coll;
  sys_par.kappa_pp   = (aux_0*aux_0 + aux_1*aux_1);
  aux_1              = log (sys_par.eps_pw)/sys_par.t_coll;
  sys_par.gamma_pw   = -2.0*log (sys_par.eps_pw)/sys_par.t_coll;
  sys_par.kappa_pw   = (aux_0*aux_0 + aux_1*aux_1);
  sys_par.fr_fact_pp = sys_par.mu_pp/1.0e12;
  sys_par.fr_fact_pw = sys_par.mu_pw/1.0e12;

  fprintf (fp_log, "gamma_pp: %lf\n", sys_par.gamma_pp);
  fprintf (fp_log, "gamma_pw: %lf\n", sys_par.gamma_pw);
  fprintf (fp_log, "kappa_pp: %lf\n", sys_par.kappa_pp);
  fprintf (fp_log, "kappa_pw: %lf\n", sys_par.kappa_pw);
  fflush  (fp_log);
#if FIRST_TIME
  /* Places particles into a cubic lattice. */
  printf ("Placing particles in their initial positions...");
  place_particles (sys_par, pp_par_hst);
  printf ("\n");

  /* Sets initial velocities */
  printf ("Setting initial velocities...");
  set_initial_velocity (sys_par, pp_par_hst);
  printf ("\n");
  aval_start = 0;
  cc = 0;
  jj = 0;

  /* Writes initial configurations to xyz file. */
  fp_xyz = fopen ("snapshots_dev.xyz", "w");
#else
  int   sec_avflag = 0, ninput;
  double pos;
  FILE  *fp_conf_in;
  printf ("Trying to restart simulation after crash...\n");
  printf ("Reading last configuration...\n");
  fp_conf_in = fopen ("last_config.dat", "r");
  ninput = fscanf (fp_conf_in, "%d\n", &aval_start);
  ninput = fscanf (fp_conf_in, "%lf\n", &(sys_par.time));
  ninput = fscanf (fp_conf_in, "%ld\n", &cc);
  ninput = fscanf (fp_conf_in, "%ld\n", &jj);
  fprintf (fp_deb, "aval_start %d\n", aval_start);
  fprintf (fp_deb, "Time %lf\n", sys_par.time);
  fprintf (fp_deb, "cc %ld\n", cc);
  fprintf (fp_deb, "jj %ld\n", jj);

  while (!feof (fp_conf_in))
  {
    ninput = fscanf (fp_conf_in, "%ld ",&ii);
    ninput = fscanf (fp_conf_in,
            "%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
            &(pp_par_hst[ii].rr.x),    &(pp_par_hst[ii].rr.y),    &(pp_par_hst[ii].rr.z),
            &(pp_par_hst[ii].Rmatx.x), &(pp_par_hst[ii].Rmatx.y), &(pp_par_hst[ii].Rmatx.z),
            &(pp_par_hst[ii].Rmaty.x), &(pp_par_hst[ii].Rmaty.y), &(pp_par_hst[ii].Rmaty.z),
            &(pp_par_hst[ii].Rmatz.x), &(pp_par_hst[ii].Rmatz.y), &(pp_par_hst[ii].Rmatz.z),
            &(pp_par_hst[ii].vv.x),    &(pp_par_hst[ii].vv.y),    &(pp_par_hst[ii].vv.z),
            &(pp_par_hst[ii].ww.x),    &(pp_par_hst[ii].ww.y),    &(pp_par_hst[ii].ww.z));
    if (ninput != 18)
    {
      printf ("Error reading input file\n");
      exit (1);
    }
  }

  fclose (fp_conf_in);

  if (aval_start)
  {
    /* Moves x wall to the end of the deposition zone.*/
    sys_par.side.x = xfinal;

    /* Writes cell parameters in log file. */
    fprintf (fp_log, "Hopper opened.\n");
    fflush  (fp_log);
    sec_avflag = 1;
  }

  /* Calculates positions of spheres. */
  vmean = wmean = 0.0;
  for (ii = 0; ii < nn; ii++)
  {
    pos = pp_par_hst[ii].radius/ROOT2;
    pp_par_hst[ii].rr_a.x = pp_par_hst[ii].rr.x + pos*(pp_par_hst[ii].Rmatx.x +
                            pp_par_hst[ii].Rmatx.y   + pp_par_hst[ii].Rmatx.z);
    pp_par_hst[ii].rr_a.y = pp_par_hst[ii].rr.y + pos*(pp_par_hst[ii].Rmaty.x +
                            pp_par_hst[ii].Rmaty.y   + pp_par_hst[ii].Rmaty.z);
    pp_par_hst[ii].rr_a.z = pp_par_hst[ii].rr.z + pos*(pp_par_hst[ii].Rmatz.x +
                            pp_par_hst[ii].Rmatz.y   + pp_par_hst[ii].Rmatz.z);
    pp_par_hst[ii].rr_b.x = pp_par_hst[ii].rr.x + pos*(pp_par_hst[ii].Rmatx.x -
                            pp_par_hst[ii].Rmatx.y   - pp_par_hst[ii].Rmatx.z);
    pp_par_hst[ii].rr_b.y = pp_par_hst[ii].rr.y + pos*(pp_par_hst[ii].Rmaty.x -
                            pp_par_hst[ii].Rmaty.y   - pp_par_hst[ii].Rmaty.z);
    pp_par_hst[ii].rr_b.z = pp_par_hst[ii].rr.z + pos*(pp_par_hst[ii].Rmatz.x -
                            pp_par_hst[ii].Rmatz.y   - pp_par_hst[ii].Rmatz.z);
    pp_par_hst[ii].rr_c.x = pp_par_hst[ii].rr.x + pos*(-pp_par_hst[ii].Rmatx.x -
                            pp_par_hst[ii].Rmatx.y    + pp_par_hst[ii].Rmatx.z);
    pp_par_hst[ii].rr_c.y = pp_par_hst[ii].rr.y + pos*(-pp_par_hst[ii].Rmaty.x -
                            pp_par_hst[ii].Rmaty.y    + pp_par_hst[ii].Rmaty.z);
    pp_par_hst[ii].rr_c.z = pp_par_hst[ii].rr.z + pos*(-pp_par_hst[ii].Rmatz.x -
                            pp_par_hst[ii].Rmatz.y    + pp_par_hst[ii].Rmatz.z);
    pp_par_hst[ii].rr_d.x = pp_par_hst[ii].rr.x + pos*(-pp_par_hst[ii].Rmatx.x +
                            pp_par_hst[ii].Rmatx.y    - pp_par_hst[ii].Rmatx.z);
    pp_par_hst[ii].rr_d.y = pp_par_hst[ii].rr.y + pos*(-pp_par_hst[ii].Rmaty.x +
                            pp_par_hst[ii].Rmaty.y    - pp_par_hst[ii].Rmaty.z);
    pp_par_hst[ii].rr_d.z = pp_par_hst[ii].rr.z + pos*(-pp_par_hst[ii].Rmatz.x +
                            pp_par_hst[ii].Rmatz.y    - pp_par_hst[ii].Rmatz.z);

    /* Calculates mean velocity. */
    vm2    = pp_par_hst[ii].vv.x*pp_par_hst[ii].vv.x;
    vm2   += pp_par_hst[ii].vv.y*pp_par_hst[ii].vv.y;
    vm2   += pp_par_hst[ii].vv.z*pp_par_hst[ii].vv.z;
    vmean += sqrt (vm2);
    wm2    = pp_par_hst[ii].ww.x*pp_par_hst[ii].ww.x;
    wm2   += pp_par_hst[ii].ww.y*pp_par_hst[ii].ww.y;
    wm2   += pp_par_hst[ii].ww.z*pp_par_hst[ii].ww.z;
    wmean += sqrt (wm2);
  }


  /* Writes initial configurations to xyz file. */
  fp_xyz = fopen ("snapshots_dev.xyz", "a");
#endif
  /* Calculates cell system dimensions. */
  cell_dim (&(cell.rr_max), &(cell.rr_min), pp_par_hst, sys_par);

  /* Calculates cell parameters and saves them. */
  cell.length.x = cell.rr_max.x - cell.rr_min.x;
  cell.length.y = cell.rr_max.y - cell.rr_min.y;
  cell.length.z = cell.rr_max.z - cell.rr_min.z;
  ncell.x       = (long) (cell.length.x/(dist*sys_par.rmax));
  ncell.y       = (long) (cell.length.y/(dist*sys_par.rmax));
  ncell.z       = (long) (cell.length.z/(dist*sys_par.rmax));
  cell_side.x   = cell.length.x/((double) ncell.x);
  cell_side.y   = cell.length.y/((double) ncell.y);
  cell_side.z   = cell.length.z/((double) ncell.z);
  cs_inv.x      = 1.0/cell_side.x;
  cs_inv.y      = 1.0/cell_side.y;
  cs_inv.z      = 1.0/cell_side.z;

  cell.ncell     = ncell;
  cell.cell_side = cell_side;
  cell.cs_inv    = cs_inv;
  cell.ncell3    = ncell.x*ncell.y*ncell.z;
  ncell3         = cell.ncell3;
  cell.ntot      = ncell3*ntags;
  ntot           = cell.ntot;
#if DEVICE
  long RB_cell = 1 + (ntot - 1)/((long) nthread);
  long RB_ocup = 1 + (ncell3 - 1)/((long) nthread);

  /* Memory allocation for cell vectors in DEVICE. */
  cudaMalloc ((void**) &cell_vec_dev,  ntot*sizeof (long));
  cudaMalloc ((void**) &nocup_vec_dev, ncell3*sizeof (int));

  /* Sets forces and torques to zero. */
  clear_forces_dev <<<RB_part, nthread>>> (sys_par, ff_dev, torque_dev);
  cudaCheckErrors ("Routine 'clear_forces' failed! Error 0.");
  
  /* Clears cells. */
  clear_cell_dev <<<RB_cell, nthread>>> (cell, cell_vec_dev);
  clear_ocup_dev <<<RB_ocup, nthread>>> (cell, nocup_vec_dev);

  /* Places particles into cell system. */
  cudaMemcpy (pp_par_dev, pp_par_hst, nn*sizeof (particle), cudaMemcpyHostToDevice);
  cudaCheckErrors ("Copying pp_par HOST to DEVICE failed! Error 0.");
  cell_locate_dev <<<RB_part, nthread>>> (cell, sys_par, pp_par_dev, nocup_vec_dev, cell_vec_dev);
  cudaCheckErrors ("Routine 'cell_locate_dev' failed! Error 0.");

  /* Calculates inter-particle forces and torques. */
  for (ii = 1; ii <= 16; ii++)
    force_part_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, cell_vec_dev, nocup_vec_dev,
                                           pp_par_dev, sys_par, cell, ii);

  /* Calculates hopper wall-particle forces and torques. */
  for (ii = 0; ii < 3; ii++)
    force_hop_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, pp_par_dev, ii);

  /* Calculates accelerations. */
  get_accel_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, aa_i_dev, anga_i_dev, sys_par,
                                        pp_par_dev);
#else
  /* Memory allocation for cell vectors in HOST. */
  cell_vec_hst  = (long *) malloc (ntotmax*sizeof (long));
  nocup_vec_hst = (int *)  malloc (n3max*sizeof (int));

  /* Sets forces and torques to zero. */
  clear_forces_hst (sys_par, ff_hst, torque_hst);
  
  /* Clears cells. */
  clear_cell_hst (cell, cell_vec_hst);
  clear_ocup_hst (cell, nocup_vec_hst);

  /* Places particles into cell system. */
  cell_locate_hst (cell, sys_par, pp_par_hst, nocup_vec_hst, cell_vec_hst);

  /* Calculates inter-particle forces and torques. */
  for (ii = 1; ii <= 16; ii++)
    force_part_hst (ff_hst, torque_hst, cell_vec_hst, nocup_vec_hst, pp_par_hst, sys_par, cell, ii);

  /* Calculates hopper wall-particle forces and torques. */
  for (ii = 0; ii < 3; ii++)
    force_hop_hst (ff_hst, torque_hst, sys_par, pp_par_hst, ii);

  /* Calculates accelerations. */
  get_accel_hst (ff_hst, torque_hst, aa_i_hst, anga_i_hst, sys_par, pp_par_hst);
#endif
  /* ========================== RUN ======================== */
#if FIRST_TIME
  fp_vmean = fopen ("mean_velocity_dev.csv", "w");
  fp_vel   = fopen ("velocities_dev.xyz",    "w");
  fp_for   = fopen ("forces_dev.dat",        "w");
  fp_xcm   = fopen ("xcm_dev.csv",           "w");
#else
  fp_vmean = fopen ("mean_velocity_dev.csv", "a");
  fp_vel   = fopen ("velocities_dev.xyz",    "a");
  fp_for   = fopen ("forces_dev.dat",        "a");
  fp_xcm   = fopen ("xcm_dev.csv",           "a");
#endif  
  if (!aval_start)
    printf ("Filling hopper...\n");

  while (sys_par.time < 0.0 && !aval_start)  /***************** Hopper filling ********************/
  {
#if DEVICE
    /* Initial part of velocity-Verlet integration. */
    init_Verlet_dev <<<RB_part, nthread>>> (aa_i_dev, anga_i_dev, sys_par, pp_par_dev);

    /* Time increment.*/
    sys_par.time += sys_par.dt;

    /* Sets forces and torques to zero. */
    clear_forces_dev <<<RB_part, nthread>>> (sys_par, ff_dev, torque_dev);

    /* Clears cells. */
    clear_cell_dev <<<RB_cell, nthread>>> (cell, cell_vec_dev);
    clear_ocup_dev <<<RB_ocup, nthread>>> (cell, nocup_vec_dev);

    /* Places particles into cell system. */
    cell_locate_dev <<<RB_part, nthread>>> (cell, sys_par, pp_par_dev, nocup_vec_dev,
                                            cell_vec_dev);

    /* Calculates inter-particle forces and torques. */
    for (ii = 1; ii <= 16; ii++)
      force_part_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, cell_vec_dev, nocup_vec_dev,
                                             pp_par_dev, sys_par, cell, ii);

    /* Calculates hopper wall-particle forces and torques. */
    for (ii = 0; ii < 3; ii++)
      force_hop_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, pp_par_dev, ii);

    /* Calculates accelerations. */
    get_accel_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, aa_f_dev, anga_f_dev, sys_par,
                                          pp_par_dev);

    /* Final part of velocity-Verlet integration. */
    finish_Verlet_dev <<<RB_part, nthread>>> (aa_i_dev, anga_i_dev, aa_f_dev, anga_f_dev, sys_par,
                                              pp_par_dev);
#else
    /* Initial part of velocity-Verlet integration. */
    init_Verlet_hst (aa_i_hst, anga_i_hst, sys_par, pp_par_hst);

    /* Time increment.*/
    sys_par.time += sys_par.dt;

    /* Sets forces and torques to zero. */
    clear_forces_hst (sys_par, ff_hst, torque_hst);

    /* Clears cells. */
    clear_cell_hst (cell, cell_vec_hst);
    clear_ocup_hst (cell, nocup_vec_hst);

    /* Places particles into cell system. */
    cell_locate_hst (cell, sys_par, pp_par_hst, nocup_vec_hst, cell_vec_hst);

    /* Calculates inter-particle forces and torques. */
    for (ii = 1; ii <= 16; ii++)
      force_part_hst (ff_hst, torque_hst, cell_vec_hst, nocup_vec_hst, pp_par_hst, sys_par, cell,
                      ii);

    /* Calculates hopper wall-particle forces and torques. */
    for (ii = 0; ii < 3; ii++)
      force_hop_hst (ff_hst, torque_hst, sys_par, pp_par_hst, ii);

    /* Calculates accelerations. */
    get_accel_hst (ff_hst, torque_hst, aa_f_hst, anga_f_hst, sys_par, pp_par_hst);

    /* Final part of velocity-Verlet integration. */
    finish_Verlet_hst (aa_i_hst, anga_i_hst, aa_f_hst, anga_f_hst, sys_par, pp_par_hst);
#endif
    /* Writes snapshots in .xyz format. */
    if (cc%2000 == 0)
    {
#if DEVICE
      cudaMemcpy (pp_par_hst, pp_par_dev, nn*sizeof (particle), cudaMemcpyDeviceToHost);
      cudaCheckErrors ("Copying pp_par DEVICE to HOST failed! Error 1.");
#endif
      /* Recalculates cell system dimensions. */
      cell_dim (&(cell.rr_max), &(cell.rr_min), pp_par_hst, sys_par);

      /* Recalculates cell parameters and saves them. */
      cell.length.x = cell.rr_max.x - cell.rr_min.x;
      cell.length.y = cell.rr_max.y - cell.rr_min.y;
      cell.length.z = cell.rr_max.z - cell.rr_min.z;
      ncell.x       = (long) (cell.length.x/(dist*sys_par.rmax));
      ncell.y       = (long) (cell.length.y/(dist*sys_par.rmax));
      ncell.z       = (long) (cell.length.z/(dist*sys_par.rmax));
      cell_side.x   = cell.length.x/((double) ncell.x);
      cell_side.y   = cell.length.y/((double) ncell.y);
      cell_side.z   = cell.length.z/((double) ncell.z);
      cs_inv.x      = 1.0/cell_side.x;
      cs_inv.y      = 1.0/cell_side.y;
      cs_inv.z      = 1.0/cell_side.z;

      cell.ncell     = ncell;
      cell.cell_side = cell_side;
      cell.cs_inv    = cs_inv;
      cell.ncell3    = ncell.x*ncell.y*ncell.z;
      ncell3         = cell.ncell3;
      cell.ntot      = ncell3*ntags;
      ntot           = cell.ntot;
#if DEVICE
      /* Frees cell and nocup vectors. */
      cudaFree (cell_vec_dev);
      cudaFree (nocup_vec_dev);

      /* Recalculats blocks for cell vectors. */
      RB_cell = 1 + (ntot - 1)/((long) nthread);
      RB_ocup = 1 + (ncell3 - 1)/((long) nthread);

      /* Reallocates cell vectors in DEVICE. */
      cudaMalloc ((void**) &cell_vec_dev,  ntot*sizeof (long));
      cudaMalloc ((void**) &nocup_vec_dev, ncell3*sizeof (int));
#else
      realloc (cell_vec_dev,  ntot*sizeof (long));
      realloc (nocup_vec_dev, ncell3*sizeof (int));
#endif
      printf ("ncell: (%ld,%ld,%ld); ntot: %ld.\n", ncell.x, ncell.y, ncell.z, ntot);

      if (cc%20000 == 0)
      {
        printf ("Time = %lf s, %ld dt\n", sys_par.time - sys_par.dt, jj);
        if (print_flag)
        {
          fprintf (fp_xyz, "%ld\n\n", 4*nn + 2);
          fprintf (fp_vel, "%ld\n\n", nn + 2);
        }

        vmean = wmean = 0.0;
        xcm   = 0.0;

        fp_conf_out = fopen ("last_config.dat", "w");

        fprintf (fp_conf_out, "%d\n%lf\n%ld\n%ld\n", aval_start, sys_par.time, cc + 1, jj + 1);
        for (ii = 0; ii < nn; ii++)
        {
          if (print_flag)
          {
            /* File for rendering in Ovito. */
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_a.x,
                     pp_par_hst[ii].rr_a.y, pp_par_hst[ii].rr_a.z, ii);
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_b.x,
                     pp_par_hst[ii].rr_b.y, pp_par_hst[ii].rr_b.z, ii);
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_c.x,
                     pp_par_hst[ii].rr_c.y, pp_par_hst[ii].rr_c.z, ii);
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_d.x,
                     pp_par_hst[ii].rr_d.y, pp_par_hst[ii].rr_d.z, ii);

            /* Velocity field file. */
            vproj = pp_par_hst[ii].vv.x*sys_par.cost - pp_par_hst[ii].vv.y*sys_par.sint;
            fprintf (fp_vel, "%lf %lf %lf %lf %lf %ld\n", pp_par_hst[ii].rr.x, pp_par_hst[ii].rr.y,
                     pp_par_hst[ii].rr.z, vproj, pp_par_hst[ii].ww.y, ii);
          }
          fprintf (fp_conf_out,
                   "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                   ii,
                   (pp_par_hst[ii].rr.x),    (pp_par_hst[ii].rr.y),    (pp_par_hst[ii].rr.z),
                   (pp_par_hst[ii].Rmatx.x), (pp_par_hst[ii].Rmatx.y), (pp_par_hst[ii].Rmatx.z),
                   (pp_par_hst[ii].Rmaty.x), (pp_par_hst[ii].Rmaty.y), (pp_par_hst[ii].Rmaty.z),
                   (pp_par_hst[ii].Rmatz.x), (pp_par_hst[ii].Rmatz.y), (pp_par_hst[ii].Rmatz.z),
                   (pp_par_hst[ii].vv.x),    (pp_par_hst[ii].vv.y),    (pp_par_hst[ii].vv.z),
                   (pp_par_hst[ii].ww.x),    (pp_par_hst[ii].ww.y),    (pp_par_hst[ii].ww.z));

          /* Calculates mean velocity. */
          vm2    = pp_par_hst[ii].vv.x*pp_par_hst[ii].vv.x;
          vm2   += pp_par_hst[ii].vv.y*pp_par_hst[ii].vv.y;
          vm2   += pp_par_hst[ii].vv.z*pp_par_hst[ii].vv.z;
          vmean += sqrt (vm2);
          wm2    = pp_par_hst[ii].ww.x*pp_par_hst[ii].ww.x;
          wm2   += pp_par_hst[ii].ww.y*pp_par_hst[ii].ww.y;
          wm2   += pp_par_hst[ii].ww.z*pp_par_hst[ii].ww.z;
          wmean += sqrt (wm2);
        }
        fclose (fp_conf_out);

        /* Ghost particles to keep simulation cell size constant in .xyz format. */
        if (print_flag)
        {
          fprintf (fp_xyz, "0.000001 %lf %lf %lf %ld\n", xfinal, YFUNNEL, sys_par.side.z, nn);
          fprintf (fp_xyz, "0.000001 0.000000 0.000000 0.000000 %ld\n", nn);
          fprintf (fp_vel, "%lf %lf %lf 0.000000 0.000000 %ld\n",
                   xfinal, YFUNNEL, sys_par.side.z, nn);
          fprintf (fp_vel, "0.000000 0.000000 0.000000 0.000000 0.000000 %ld\n", nn);
        }

        /* Writes mean velocities to file. */
        vmean /= nn;
        wmean /= nn;
        fprintf (fp_vmean, "%lf, %lf, %lf\n", sys_par.time - sys_par.dt, vmean, wmean);

        /* Flushes data to output files. */
        fflush (fp_vmean);
        fflush (fp_xcm);
        if (print_flag)
        {
          fflush (fp_xyz);
          fflush (fp_vel);
        }
        jj++;
      }
    }

    /* Swaps accelerations */
    swap_dev   = aa_i_dev;
    aa_i_dev   = aa_f_dev;
    aa_f_dev   = swap_dev;
    swap_dev   = anga_i_dev;
    anga_i_dev = anga_f_dev;
    anga_f_dev = swap_dev;

    if (vmean < 0.1*sys_par.v0)
      break;

    cc++;
  }
#if !FIRST_TIME
if (!sec_avflag)
{
#endif
  /* Moves x wall to the end of the deposition zone.*/
  sys_par.side.x = xfinal;

  if (!aval_start)
    printf ("Opening hopper... Avalanche starts!\n");
  aval_start = 1;
  cc = 0;
  sys_par.time = 0.0;
#if !FIRST_TIME
}
#endif
/* ======================== AVALANCHE TIME ======================== */
  while (sys_par.time < 5.0)
  {
#if DEVICE
    /* Initial part of velocity-Verlet integration. */
    init_Verlet_dev <<<RB_part, nthread>>> (aa_i_dev, anga_i_dev, sys_par, pp_par_dev);
    //cudaCheckErrors ("init_Verlet failed!");

    /* Time increment.*/
    sys_par.time += sys_par.dt;

    /* Calculates rotated position vectors. */
    rotate_vectors_dev <<<RB_part, nthread>>> (pp_par_dev, sys_par);
    //cudaCheckErrors ("Rotate_vectors failed!");

    /* Sets forces and torques to zero. */
    clear_forces_dev <<<RB_part, nthread>>> (sys_par, ff_dev, torque_dev);
    //cudaCheckErrors ("clear_forces failed!");

    /* Clears cells. */
    clear_cell_dev <<<RB_cell, nthread>>> (cell, cell_vec_dev);
    clear_ocup_dev <<<RB_ocup, nthread>>> (cell, nocup_vec_dev);
    //cudaCheckErrors ("clear_cell, clear_ocup failed!");

    /* Places particles into cell system. */
    cell_locate_dev <<<RB_part, nthread>>> (cell, sys_par, pp_par_dev, nocup_vec_dev,
                                            cell_vec_dev);
    //cudaCheckErrors ("cell_locate failed!");

    /* Calculates hopper wall-particle forces and torques. */
    for (ii = 0; ii < 3; ii++)
      force_hop_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, pp_par_dev, ii);
    //cudaCheckErrors ("force_hop failed!");

    /* Calculates funnel-particle forces and torques. */
    force_funnel_dev   <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, pp_par_dev, fun);
    //cudaCheckErrors ("force_funnel failed!");

    force_fcorners_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, fun, pp_par_dev);
    //cudaCheckErrors ("force_corners failed!");

    /* Calculates wall-particle forces and torques. */
    force_walls_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, pp_par_dev);
    //cudaCheckErrors ("force_Walls failed!");

    /* Calculates inter-particle forces and torques. */
    for (ii = 1; ii <= 16; ii++)
      force_part_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, cell_vec_dev, nocup_vec_dev,
                                             pp_par_dev, sys_par, cell, ii);
    //cudaCheckErrors ("force_part failed!");

    /* Calculates flume-particle interactions. */
    force_flume_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, sys_par, pp_par_dev);
    //cudaCheckErrors ("force_flume failed!");

    /* Calculates accelerations. */
    get_accel_dev <<<RB_part, nthread>>> (ff_dev, torque_dev, aa_f_dev, anga_f_dev, sys_par,
                                          pp_par_dev);
    //cudaCheckErrors ("get_accel failed!");

    /* Final part of velocity-Verlet integration. */
    finish_Verlet_dev <<<RB_part, nthread>>> (aa_i_dev, anga_i_dev, aa_f_dev, anga_f_dev, sys_par,
                                              pp_par_dev);
    //cudaCheckErrors ("finish_Verlet failed!");
#else
    /* Initial part of velocity-Verlet integration. */
    init_Verlet_hst (aa_i_hst, anga_i_hst, sys_par, pp_par_hst);

    /* Time increment.*/
    sys_par.time += sys_par.dt;

    /* Calculates rotated position vectors. */
    rotate_vectors_hst (pp_par_hst, sys_par);

    /* Sets forces and torques to zero. */
    clear_forces_hst (sys_par, ff_hst, torque_hst);

    /* Clears cells. */
    clear_cell_hst (cell, cell_vec_hst);
    clear_ocup_hst (cell, nocup_vec_hst);

    /* Places particles into cell system. */
    cell_locate_hst (cell, sys_par, pp_par_hst, nocup_vec_hst, cell_vec_hst);

    /* Calculates hopper wall-particle forces and torques. */
    for (ii = 0; ii < 3; ii++)
      force_hop_hst (ff_hst, torque_hst, sys_par, pp_par_hst, ii);

    /* Calculates funnel-particle forces and torques. */
    force_funnel_hst (ff_hst, torque_hst, sys_par, pp_par_hst, fun);
    force_fcorners_hst (ff_hst, torque_hst, sys_par, fun, pp_par_hst);

    /* Calculates wall-particle forces and torques. */
    force_walls_hst (ff_hst, torque_hst, sys_par, pp_par_hst);

    /* Calculates inter-particle forces and torques. */
    for (ii = 1; ii <= 16; ii++)
      force_part_hst (ff_hst, torque_hst, cell_vec_hst, nocup_vec_hst, pp_par_hst, sys_par, cell,
                      ii);

    /* Calculates flume-particle interactions. */
    force_flume_hst (ff_hst, torque_hst, sys_par, pp_par_hst);

    /* Calculates accelerations. */
    get_accel_hst (ff_hst, torque_hst, aa_f_hst, anga_f_hst, sys_par, pp_par_hst);

    /* Final part of velocity-Verlet integration. */
    finish_Verlet_hst (aa_i_hst, anga_i_hst, aa_f_hst, anga_f_hst, sys_par, pp_par_hst);
#endif
    /* Writes snapshots in .xyz format. */	
    if (cc%2000 == 0)
    {
#if DEVICE
      cudaMemcpy (pp_par_hst, pp_par_dev, nn*sizeof (particle), cudaMemcpyDeviceToHost);
      cudaCheckErrors ("Copying pp_par from DEVICE to HOST failed at avalanche time! No load.");
#endif
      /* Recalculates cell system dimensions. */
      cell_dim (&(cell.rr_max), &(cell.rr_min), pp_par_hst, sys_par);

      /* Recalculates cell parameters and saves them. */
      cell.length.x = cell.rr_max.x - cell.rr_min.x;
      cell.length.y = cell.rr_max.y - cell.rr_min.y;
      cell.length.z = cell.rr_max.z - cell.rr_min.z;
      ncell.x       = (long) (cell.length.x/(dist*sys_par.rmax));
      ncell.y       = (long) (cell.length.y/(dist*sys_par.rmax));
      ncell.z       = (long) (cell.length.z/(dist*sys_par.rmax));
      cell_side.x   = cell.length.x/((double) ncell.x);
      cell_side.y   = cell.length.y/((double) ncell.y);
      cell_side.z   = cell.length.z/((double) ncell.z);
      cs_inv.x      = 1.0/cell_side.x;
      cs_inv.y      = 1.0/cell_side.y;
      cs_inv.z      = 1.0/cell_side.z;

      cell.ncell     = ncell;
      cell.cell_side = cell_side;
      cell.cs_inv    = cs_inv;
      cell.ncell3    = ncell.x*ncell.y*ncell.z;
      ncell3         = cell.ncell3;
      cell.ntot      = ncell3*ntags;
      ntot           = cell.ntot;
#if DEVICE
      /* Frees cell vector memory. */
      cudaFree   (cell_vec_dev);
      cudaFree   (nocup_vec_dev);

      /* Recalculates blocks for cell vectors. */
      RB_cell = 1 + (ntot - 1)/((long) nthread);
      RB_ocup = 1 + (ncell3 - 1)/((long) nthread);

      /* Reallocs memory for cell vectors in DEVICE. */
      cudaMalloc ((void**) &cell_vec_dev,  ntot*sizeof (long));
      cudaMalloc ((void**) &nocup_vec_dev, ncell3*sizeof (int));
#endif
      printf ("ncell: (%ld,%ld,%ld); ntot: %ld.\n", ncell.x, ncell.y, ncell.z, ntot);

      if (cc%20000 == 0)
      {
        printf ("Avalanche time = %lf s, %ld dt\n", sys_par.time - sys_par.dt, jj);

        if (print_flag)
        {
          fprintf (fp_xyz, "%ld\n\n", 4*nn + 2);
          fprintf (fp_vel, "%ld\n\n", nn + 2);
        }

        vmean = wmean = 0.0;
        xcm   = zcm   = 0.0;

        fp_conf_out = fopen ("last_config.dat", "w");
        fprintf (fp_conf_out, "%d\n%lf\n%ld\n%ld\n", aval_start, sys_par.time, cc + 1, jj + 1);

        for (ii = 0; ii < nn; ii++)
        {
          if (print_flag)
          {
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_a.x,
                     pp_par_hst[ii].rr_a.y, pp_par_hst[ii].rr_a.z, ii);
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_b.x,
                     pp_par_hst[ii].rr_b.y, pp_par_hst[ii].rr_b.z, ii);
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_c.x,
                     pp_par_hst[ii].rr_c.y, pp_par_hst[ii].rr_c.z, ii);
            fprintf (fp_xyz, "%lf %lf %lf %lf %ld\n", pp_par_hst[ii].radius, pp_par_hst[ii].rr_d.x,
                     pp_par_hst[ii].rr_d.y, pp_par_hst[ii].rr_d.z, ii);

            /* Velocity field file. */
            vproj = pp_par_hst[ii].vv.x*sys_par.cost - pp_par_hst[ii].vv.y*sys_par.sint;
            fprintf (fp_vel, "%lf %lf %lf %lf %lf %ld\n", pp_par_hst[ii].rr.x, pp_par_hst[ii].rr.y,
                     pp_par_hst[ii].rr.z, vproj, pp_par_hst[ii].ww.y, ii);
          }

          fprintf (fp_conf_out,
                   "%ld %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",
                   ii,
                   (pp_par_hst[ii].rr.x),    (pp_par_hst[ii].rr.y),    (pp_par_hst[ii].rr.z),
                   (pp_par_hst[ii].Rmatx.x), (pp_par_hst[ii].Rmatx.y), (pp_par_hst[ii].Rmatx.z),
                   (pp_par_hst[ii].Rmaty.x), (pp_par_hst[ii].Rmaty.y), (pp_par_hst[ii].Rmaty.z),
                   (pp_par_hst[ii].Rmatz.x), (pp_par_hst[ii].Rmatz.y), (pp_par_hst[ii].Rmatz.z),
                   (pp_par_hst[ii].vv.x),    (pp_par_hst[ii].vv.y),    (pp_par_hst[ii].vv.z),
                   (pp_par_hst[ii].ww.x),    (pp_par_hst[ii].ww.y),    (pp_par_hst[ii].ww.z));

          /* Calculates mean velocity. */
          vm2    = pp_par_hst[ii].vv.x*pp_par_hst[ii].vv.x;
          vm2   += pp_par_hst[ii].vv.y*pp_par_hst[ii].vv.y;
          vm2   += pp_par_hst[ii].vv.z*pp_par_hst[ii].vv.z;
          vmean += sqrt (vm2);
          wm2    = pp_par_hst[ii].ww.x*pp_par_hst[ii].ww.x;
          wm2   += pp_par_hst[ii].ww.y*pp_par_hst[ii].ww.y;
          wm2   += pp_par_hst[ii].ww.z*pp_par_hst[ii].ww.z;
          wmean += sqrt (wm2);
        }
        fclose (fp_conf_out);

        vmean /= nn;
        wmean /= nn;

        if (print_flag)
        {
          fprintf (fp_xyz, "0.0001 %lf %lf %lf %ld\n", sys_par.side.x, YFUNNEL, sys_par.side.z, nn);
          fprintf (fp_xyz, "0.000001 0.000000 0.000000 0.000000 %ld\n", nn);
          fprintf (fp_vel, "%lf %lf %lf 0.000000 0.000000 %ld\n", xfinal, YFUNNEL, sys_par.side.z,
                   nn);
          fprintf (fp_vel, "0.000000 0.000000 0.000000 0.000000 0.000000 %ld\n", nn);
        }
        fprintf (fp_vmean, "%lf, %lf, %lf\n", sys_par.time - sys_par.dt, vmean, wmean);

        /* Calculates center of mass position. */
        for (ii = 0; ii < nn; ii ++)
        {
          xmass = pp_par_hst[ii].mass*pp_par_hst[ii].rr.x;
          zmass = pp_par_hst[ii].mass*pp_par_hst[ii].rr.z;
          xcm  += xmass;
          zcm  += zmass;
        }
        xcm /= sys_par.tot_mass;
        zcm /= sys_par.tot_mass;
        fprintf (fp_xcm, "%lf, %lf, %lf\n", sys_par.time - sys_par.dt, xcm, zcm);

        fflush (fp_vmean);
        fflush (fp_xcm);

        if (print_flag)
        {
          fflush (fp_xyz);
          fflush (fp_vel);
        }

        jj++;
      }
    }

    /* Swaps accelerations */
    swap_dev   = aa_i_dev;
    aa_i_dev   = aa_f_dev;
    aa_f_dev   = swap_dev;
    swap_dev   = anga_i_dev;
    anga_i_dev = anga_f_dev;
    anga_f_dev = swap_dev;

    cc++;
    cudaCheckErrors ("General error.");
    
    if (vmean <= 1.0e-3*sys_par.v0)
      break;
  }

  /* Frees memory vectors. */
  cudaFree (aa_i_dev);
  cudaFree (aa_f_dev);
  cudaFree (swap_dev);
  cudaFree (ff_dev);
  cudaFree (torque_dev);
  cudaFree (pp_par_dev);
  cudaFree (cell_vec_dev);
  cudaFree (nocup_vec_dev);

  fclose (fp_vmean);
  fclose (fp_vel);
  fclose (fp_for);
  fclose (fp_log);
  fclose (fp_xyz);
  fclose (fp_deb);
  fclose (fp_xcm);
  return 0;
}
