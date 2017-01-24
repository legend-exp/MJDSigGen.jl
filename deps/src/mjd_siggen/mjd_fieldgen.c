/* program to calculate electric fields and weighting potentials
              of PPC and BEGe Ge detectors by relaxation
   author:           D.C. Radford
   first written:    Nov 2007
   this modified version for MJD:  Oct 2014
      - uses the same (single) config file as modified MJD siggen
      - added intelligent coarse grid / refinement of grid
      - added interpolation of RC and LC positions on the grid
   June 2016: added optional bulletization of point contact

   TO DO:
      - add other bulletizations
      - add dead layer / Li thickness
      - on coarse grids, interpolate the position of L, R, LT, and the ditch
            (as is done now already for RC and LC)
*/

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mjd_siggen.h"

#define MAX_ITS 50000     // default max number of iterations for relaxation
#define MAX_ITS_FACTOR 2  // factor by which max iterations is reduced as grid is refined

int report_config(FILE *fp_out, char *config_file_name);


int main(int argc, char **argv)
{

  MJD_Siggen_Setup setup;

  /* --- default values, normally over-ridden by values in a *.conf file --- */
  int   R = 0;   // radius of detector, in grid lengths
  int   L = 0;   // length of detector, in grid lengths
  int   RC = 0;  // radius of central contact, in grid lengths
  int   LC = 0;  // length of central contact, in grid lengths
  int   LT = 0;  // length of taper, in grid lengths
  int   RO = 0;  // radius of wrap-around outer (Li) contact, in grid lengths
  int   LO = 0;  // length of ditch next to wrap-around outer (Li) contact, in grid lengths
  int   WO = 0;  // width of ditch next to wrap-around outer (Li) contact, in grid lengths
  float BV = 0;  // bias voltage
  float N = 1;   // charge density at z=0 in units of e+10/cm3
  float M = 0;   // charge density gradient, in units of e+10/cm4

  int   WV = 0;  // 0: do not write the V and E values to ppc_ev.dat
                 // 1: write the V and E values to ppc_ev.dat
                 // 2: write the V and E values for both +r, -r (for gnuplot, NOT for siggen)
  int   WP = 0;  // 0: do not calculate the weighting potential
                 // 1: calculate the WP and write the values to ppc_wp.dat
  /* ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  ---  --- */

  double **v[2], **eps, **eps_dr, **eps_dz, **vfraction, *s1, *s2;
  char   **undepleted, config_file_name[256];
  int    **bulk, *rrc;
  float  *drrc, *frrc;
  double eps_sum, v_sum, mean, min, f, f1z, f2z, f1r, f2r;
  double e_over_E = 11.31; // e/epsilon
                           // for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0
  float  dif, sum_dif=0, max_dif, a, b, c, grid = 0.5, dRC, dLC, fLC=0;
  float  E_r, E_z, bubble_volts=0, cs, gridstep[3];
  int    i, j, r, z, iter, old, new=0, zz, rr, istep, max_its;
  FILE   *file;
  time_t t0=0, t1, t2=0;
  double esum, esum2, pi=3.14159, Epsilon=(8.85*16.0/1000.0);  // permittivity of Ge in pF/mm
  double pinched_sum1, pinched_sum2, *imp_ra, *imp_rm, *imp_z, S=0;
  int    gridfact, fully_depleted=0, LL=L, RR=R, zmax, rmax;

  if (argc%2 != 1) {
    printf("Possible options:\n"
	   "      -c config_file_name\n"
	   "      -b bias_volts\n"
	   "      -w {0,1}    (do_not/do write the field file)\n"
	   "      -p {0,1}    (do_not/do write the WP file)\n");
    return 1;
  }

  for (i=1; i<argc-1; i+=2) {
    if (strstr(argv[i], "-c")) {
      if (read_config(argv[i+1], &setup)) return 1;
      strncpy(config_file_name, argv[i+1], sizeof(config_file_name));

      if (setup.xtal_grid < 0.001) setup.xtal_grid = 0.5;
      grid = setup.xtal_grid;

      L  = LL = lrint(setup.xtal_length/grid);
      R  = RR = lrint(setup.xtal_radius/grid);
      // BRT = lrint(setup.top_bullet_radius/grid);
      // BRB = lrint(setup.bottom_bullet_radius/grid);
      LC = lrint(setup.pc_length/grid);
      RC = lrint(setup.pc_radius/grid);
      LT = lrint(setup.taper_length/grid);
      RO = lrint(setup.wrap_around_radius/grid);
      LO = lrint(setup.ditch_depth/grid);
      WO = lrint(setup.ditch_thickness/grid);
      // LiT = lrint(setup.Li_thickness/grid);
      N  = setup.impurity_z0;
      M  = setup.impurity_gradient;
      BV = setup.xtal_HV;
      WV = setup.write_field;
      WP = setup.write_WP;

    } else if (strstr(argv[i], "-b")) {
      BV = atof(argv[i+1]);   // bias volts
    } else if (strstr(argv[i], "-w")) {
      WV = atoi(argv[i+1]);   // write-out options
    } else if (strstr(argv[i], "-p")) {
      WP = atoi(argv[i+1]);   // weighting-potential options
    } else {
      printf("Possible options:\n"
	     "      -c config_file_name\n"
	     "      -b bias_volts\n"
	     "      -w {0,1,2}    (for WV options)\n"
	     "      -p {0,1}      (for WP options)\n");
      return 1;
    }
  }

  if (L <= 1 || R <= 1) {
    printf("ERROR: No configuration file specified.\n"
	   "Possible options:\n"
	   "      -c config_file_name\n"
	   "      -b bias_volts\n"
	   "      -w {0,1,2}    (for WV options)\n"
	   "      -p {0,1}      (for WP options)\n");
    return 1;
  }
  if (L*R > 2500*2500) {
    printf("Error: Crystal size divided by grid size is too large!\n");
    return 1;
  }
  if (WV < 0 || WV > 2) WV = 0;

  if (RO <= 0.0 || RO >= R) {
    RO = R - LT;    // inner radius of taper, in grid lengths
    printf("\n\n"
	   " Crystal: Radius x Length: %.1f x %.1f mm\n"
	   "   Taper: %.1f mm\n"
	   "No wrap-around contact or ditch...\n"
	   "Bias: %.0f V\n"
	   "Impurities: (%.3f + %.3fz) e10/cm3\n\n",
	   grid * (float) R, grid * (float) L, grid * (float) LT,
	   BV, N, M);
  } else {
    printf("\n\n"
	   "    Crystal: Radius x length: %.1f x %.1f mm\n"
	   "      Taper: %.1f mm\n"
	   "Wrap-around: Radius x ditch x gap:  %.1f x %.1f x %.1f mm\n"
	   "       Bias: %.0f V\n"
	   " Impurities: (%.3f + %.3fz) e10/cm3\n\n",
	   grid * (float) R, grid * (float) L, grid * (float) LT,
	   grid * (float) RO, grid * (float) LO, grid * (float) WO, BV, N, M);
  }
  if (setup.bulletize_PC)
    printf("   Contact: Radius x length: %.1f x %.1f mm, bulletized\n\n",
	   grid * (float) RC, grid * (float) LC);
  else
    printf("   Contact: Radius x length: %.1f x %.1f mm, not bulletized\n\n",
	   grid * (float) RC, grid * (float) LC);

    if ((BV < 0 && N < 0) || (BV > 0 && N > 0)) {
    printf("ERROR: Expect bias and impurity to be opposite sign!\n");
    return 1;
  }
  if (N > 0) {
    // swap polarity for n-type material; this lets me assume all voltages are positive
    BV = -BV;
    M = -M;
    N = -N;
  }

  /* malloc arrays
     float v[2][L+5][R+5];
     float eps[L+1][R+1], eps_dr[L+1][R+1], eps_dz[L+1][R+1];
     float vfraction[L+1][R+1], s1[R], s2[R], drrc[LC+2], drrc[LC+2];
     char  undepleted[R+1][L+1];
     int   bulk[L+1][R+1], rrc[LC+2];
  */
  if ((v[0]   = malloc((L+5)*sizeof(*v[0]))) == NULL ||
      (v[1]   = malloc((L+5)*sizeof(*v[1]))) == NULL ||
      (eps    = malloc((L+1)*sizeof(*eps)))  == NULL ||
      (eps_dr = malloc((L+1)*sizeof(*eps_dr))) == NULL ||
      (eps_dz = malloc((L+1)*sizeof(*eps_dz))) == NULL ||
      (bulk   = malloc((L+1)*sizeof(*bulk)))   == NULL ||
      (vfraction  = malloc((L+1)*sizeof(*vfraction)))  == NULL ||
      (undepleted = malloc((R+1)*sizeof(*undepleted))) == NULL ||
      (imp_ra = malloc((R+1)*sizeof(*imp_ra))) == NULL ||
      (imp_rm = malloc((R+1)*sizeof(*imp_rm))) == NULL ||
      (imp_z = malloc((L+1)*sizeof(*imp_z))) == NULL ||
      (rrc   = malloc((LC+2)*sizeof(*rrc)))  == NULL ||
      (drrc  = malloc((LC+2)*sizeof(*drrc))) == NULL ||
      (frrc  = malloc((LC+2)*sizeof(*frrc))) == NULL ||
      (s1 = malloc((R+1)*sizeof(*s1))) == NULL ||
      (s2 = malloc((R+1)*sizeof(*s2))) == NULL) {
    printf("Malloc failed\n");
    return 1;
  }
#define ERR { printf("Malloc failed; j = %d\n", j); return 1; }
  for (j=0; j<L+1; j++) if ((v[0][j] = malloc((R+5)*sizeof(**v[0]))) == NULL) ERR;
  for (j=0; j<L+1; j++) if ((v[1][j] = malloc((R+5)*sizeof(**v[1]))) == NULL) ERR;
  for (j=0; j<L+1; j++) if ((eps[j]  = malloc((R+1)*sizeof(**eps)))  == NULL) ERR;
  for (j=0; j<L+1; j++) if ((eps_dr[j] = malloc((R+1)*sizeof(**eps_dr))) == NULL) ERR;
  for (j=0; j<L+1; j++) if ((eps_dz[j] = malloc((R+1)*sizeof(**eps_dz))) == NULL) ERR;
  for (j=0; j<L+1; j++) if ((bulk[j] = malloc((R+1)*sizeof(**bulk))) == NULL) ERR;
  for (j=0; j<L+1; j++) if ((vfraction[j] = malloc((R+1)*sizeof(**vfraction))) == NULL) ERR;
  for (j=0; j<R+1; j++) {
    if ((undepleted[j] = malloc((L+1)*sizeof(**undepleted))) == NULL) ERR;
    memset(undepleted[j], ' ', (L+1)*sizeof(**undepleted));
  }
  for (r=0; r<R+1; r++) {
    imp_ra[r] = 0.0;
    imp_rm[r] = 1.0;
  }

  /* In the following we divide areas and volumes by pi
    r_bin   rmax  A_top A_outside A_inside  volume  total_surf  out/top  tot/vol
      0     1/2    1/4      1         0       1/4      1.5         4        6  << special case
      1     3/2      2      3         1        2        8        3/2        4
      2     5/2      4      5         3        4       16        5/4        4
      3     7/2      6      7         5        6       24        7/6        4
      r   r+0.5     2r    2r+1      2r-1      2r       8r     (2r+1)/2r     4
                                                              = 1+0.5/r
  */
  // weighting values for the relaxation alg. as a function of r
  s1[0] = 4.0;
  s2[0] = 0.0;
  for (r=1; r<R+1; r++) {
    s1[r] = 1.0 + 0.5 / (double) r;   //  for r+1
    s2[r] = 1.0 - 0.5 / (double) r;   //  for r-1
  }

  /*
    If grid is too small compared to the crystal size, then it will take too
    long for the relaxation to converge. In that case, we use an adaptive
    grid, where we start out coarse and then refine the grid.
  */
  cs = sqrt(setup.xtal_length * setup.xtal_radius);
  i = 1 + ((int) (cs/grid)) / 100;
  if (i < 2) {
    gridstep[0] = grid;
    gridstep[1] = gridstep[2] = 0;
    printf("Single grid size: %.4f\n", grid);
  } else if (i < 6) {
    gridstep[0] = (float) i * grid;
    gridstep[1] = grid;
    gridstep[2] = 0;
    printf("Two grid sizes: %.4f %.4f\n", gridstep[0], grid);
  } else {  // i > 5
    j = (i+4)/5;
    i = (i+j-1)/j;
    gridstep[0] = (float) (i*j) * grid;
    gridstep[1] = (float) j * grid;
    gridstep[2] = grid;
    printf("Three grid sizes: %.4f %.4f %.4f (%d %d)\n",
	   gridstep[0], gridstep[1], grid, i, j);
  }

  /* to be safe, initialize overall potential to bias voltage */
  for (z=0; z<LL+1; z++) {
    for (r=0; r<RR+1; r++) {
      v[0][z][r] = v[1][z][r] = BV;
    }
  }
  if (setup.verbosity >= CHATTY)
    t0 = t2 = time(NULL);  // for calculating elapsed time later...
  max_its = MAX_ITS;
  if (setup.max_iterations > 0) max_its = setup.max_iterations;
  /* now set up and perform the relaxation for each of the grid step sizes in turn */
  for (istep=0; istep<3 && gridstep[istep]>0; istep++) {
    grid = gridstep[istep]; // grid size for this go-around
    old = 1;
    new = 0;
    /*  e/espilon * area of pixel in mm2 / 4
	for 1 mm2, charge units 1e10 e/cm3, espilon = 16*epsilon0
	4.0 = surface area / volume of voxel in cylindrical (2D) coords
	(this would be 6.0 in cartesian coords) */
    e_over_E = 11.31 * grid*grid / 4.0;

    if (istep > 0) {
      /* not the first go-around, so the previous calculation was on a coarser grid...
	 now copy/expand the potential to the new finer grid
      */
      i = (int) (gridstep[istep-1] / gridstep[istep] + 0.5);
      f = 1.0 / (float) i;
      printf("\ngrid %.4f -> %.4f; ratio = %d %.3f\n\n",
	     gridstep[istep-1], gridstep[istep], i, f);
      for (z=0; z<L+1; z++) {
	for (r=0; r<R+1; r++) {
	  f1z = 0.0;
	  zmax = i*z+i;
	  if (zmax > LL+1) zmax = LL+1;
	  for (zz=i*z; zz<zmax; zz++) {
	    f2z = 1.0 - f1z;
	    f1r = 0.0;
	    rmax = i*r+i;
	    if (rmax > RR+1) rmax = RR+1;
	    for (rr=i*r; rr<rmax; rr++) {
	      f2r = 1.0 - f1r;
	      v[0][zz][rr] =      // linear interpolation of potential
		f2z*f2r*v[1][z][r  ] + f1z*f2r*v[1][z+1][r  ] +
		f2z*f1r*v[1][z][r+1] + f1z*f1r*v[1][z+1][r+1];
	      f1r += f;
	    }
	    f1z += f;
	  }
	}
      }
    }

    // recalculate geometry dimensions in units of the current grid size
    L  = lrint(setup.xtal_length/grid);
    R  = lrint(setup.xtal_radius/grid);
    // BRT = lrint(setup.top_bullet_radius/grid);
    // BRB = lrint(setup.bottom_bullet_radius/grid);
    LC = lrint(setup.pc_length/grid);
    // distance in grid units from PC length to the middle of the nearest pixel:
    dLC = setup.pc_length/grid - (float) LC;
    if (dLC < 0.01 && dLC > -0.01) dLC = 0;
    RC = lrint(setup.pc_radius/grid);
    // distance in grid units from PC radius to the middle of the nearest pixel:
    dRC = setup.pc_radius/grid - (float) RC;
    if (dRC < 0.05 && dRC > -0.05) dRC = 0;
    /* set up bulletization inside point contact */
    if (setup.bulletize_PC) {
      for (z=0; z<=LC; z++) {
        if (setup.pc_length <= setup.pc_radius) {  // LC <= RC; use LC as bulletization radius
          a = setup.pc_radius - setup.pc_length;
          b = z * grid;
          c = setup.pc_length*setup.pc_length - b*b;
          if (c < 0.0) c = 0;
          c = a + sqrt(c);
        } else {  // LC > RC; use RC as bulletization radius
          if (z > LC-RC) {
            a = setup.pc_length - setup.pc_radius;
            b = z * grid - a;
            c = setup.pc_radius*setup.pc_radius - b*b;
            if (c < 0.0) c = 0;
            c = sqrt(c);
          } else {
            c = setup.pc_radius;
          }
        }
        rrc[z] = lrint(c/grid);
        drrc[z] = c/grid - (float) rrc[z];
        if (drrc[z] < 0.05 && drrc[z] > -0.05) drrc[z] = 0;
        frrc[z] = 0;
        // printf(">> z rrc drrc: %d %d %f\n", z, rrc[z], drrc[z]);
      }
      drrc[LC+1] = drrc[LC];
    } else {  // no bulletization
      for (z=0; z<=LC+1; z++) {
        rrc[z] = RC;
        drrc[z] = dRC;
        frrc[z] = 0;
      }
    }
    
    LT = lrint(setup.taper_length/grid);
    RO = lrint(setup.wrap_around_radius/grid);
    LO = lrint(setup.ditch_depth/grid);
    WO = lrint(setup.ditch_thickness/grid);
    // LiT = lrint(setup.Li_thickness/grid);

    S = setup.impurity_surface * e_over_E / grid;
    for (z=0; z<L+1; z++) {
      imp_z[z] = (N + 0.1 * M * grid * (double) z + 
		  setup.impurity_quadratic * (1.0 - (double) ((z-L/2)*(z-L/2)) /
					      (double) (L*L/4))) * e_over_E;
    }
    if (setup.impurity_rpower > 0.1) {
      for (r=0; r<R+1; r++) {
	imp_ra[r] = setup.impurity_radial_add * e_over_E *
	  pow((double) r / (double) R, setup.impurity_rpower);
	imp_rm[r] = 1.0 + (setup.impurity_radial_mult - 1.0f) *
	  pow((double) r / (double) R, setup.impurity_rpower);
      }
    }
    if (setup.verbosity >= NORMAL)
      printf("grid = %f  RC = %d  dRC = %f  LC = %d  dLC = %f\n\n",
	     grid, RC, dRC, LC, dLC);
    if (RO <= 0.0 || RO >= R) RO = R - LT;    // inner radius of taper, in grid lengths

    if (istep == 0) {
      // no previous coarse relaxation, so make initial wild guess at potential:
      for (z=0; z<L; z++) {
	a = BV * (float) (z) / (float) L;
	for (r=0; r<R; r++) {
	  v[0][z][r] =  a + (BV - a) * (float) (r) / (float) R;
	}
      }
    }

    /* boundary conditions and permittivity
       boundary condition at Ge-vacuum interface:
       epsilon0 * E_vac = espilon_Ge * E_Ge
    */
    for (z=0; z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	eps[z][r] = eps_dz[z][r] = eps_dr[z][r] = 16;   // permittivity inside Ge
	if (z < LO  && r < RO && r > RO-WO-1) eps[z][r] =  1;  // permittivity inside vacuum
	if (r > 0) eps_dr[z][r-1] = (eps[z][r-1]+eps[z][r])/2.0f;
	if (z > 0) eps_dz[z-1][r] = (eps[z-1][r]+eps[z][r])/2.0f;
      }
    }

    for (z=0; z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	vfraction[z][r] = 1.0;
	if (z < LO && r < RO && r > RO-WO-1) {
	  vfraction[z][r] = 0.0;
	}
	// boundary conditions
	bulk[z][r] = 0;  // flag for normal bulk, no complications
	// outside (HV) contact:
	if (z == L ||
	    r == R ||
	    r >= z + R - LT ||       // taper
	    (z == 0 && r >= RO)) {   // wrap-around
	  bulk[z][r] = -1;               // value of v[*][z][r] is fixed...
	  v[0][z][r] = v[1][z][r] = BV;  // at the bias voltage
	}
	// inside (point) contact, with optional bulletization:
	else if (z <= LC && r <= rrc[z]) {
	  bulk[z][r] = -1;                // value of v[*][z][r] is fixed...
	  v[0][z][r] = v[1][z][r] = 0;    // at zero volts
	  /* radial edge of inside contact; if the PC radius is not in the middle
	     of a pixel, we want to modify interpolation of V in surrounding pixels
	  */
	  if (r == rrc[z] && drrc[z] < -0.05) {
	    bulk[z][r] = 1;  // flag for radial edge of PC
	    frrc[z] = -1.0/drrc[z];  // interpolation weight for pixel at (r-1)
	    // only part of the pixel has volume charge density, the rest is contact
	    vfraction[z][r] *= -2.0*drrc[z];
	  }
	  /* z edge of inside contact; if the PC length is not in the middle
	     of a pixel, we want to modify interpolation of V in surrounding pixels
	  */
	  if (z == LC && dLC < -0.05) {
	    bulk[z][r] = 2;  // flag for z edge of PC
	    fLC = -1.0/dLC;  // interpolation weight for pixel at (z-1)
	    // only part of the pixel has volume charge density, the rest is contact
	    vfraction[z][r] *= -2.0*dLC;
	  }
	}
	/* edges of inside contact; if the PC radius and/or legth is not in the middle
	   of a pixel, we want to modify interpolation of V in surrounding pixels...
	   in this case, the radius/length > grid point, so it modifies the
	   interpolation for the next point out
	*/
	// FIXME: Check for adjacent ditch
	else if (z <= LC && r == rrc[z]+1 && drrc[z] > 0.05) {
	  bulk[z][r] = 1;         // flag for radial edge of PC
	  frrc[z] = 1.0/(1.0 - drrc[z]);  // interpolation weight for pixel at (r-1)
	}
	else if (z == LC+1 && r <= rrc[z] && dLC > 0.05) {
	  bulk[z][r] = 2;         // flag for z edge of PC
	  fLC = 1.0/(1.0 - dLC);  // interpolation weight for pixel at (z-1)
	}
      }
    }

    // now do the actual relaxation
    //for (iter=0; iter<max_its/3; iter++) {
    for (iter=0; iter<max_its; iter++) {
      if (old == 0) {
	old = 1;
	new = 0;
      } else {
	old = 0;
	new = 1;
      }
      sum_dif = 0.0f;
      max_dif = 0.0f;
      bubble_volts = 0.0f;

      for (z=0; z<L; z++) {
	for (r=0; r<R; r++) {
	  if (bulk[z][r] < 0) continue;      // outside or inside contact

	  if (bulk[z][r] == 0) {             // normal bulk, no complications
	    v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r];
	    eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r];
	    min = fminf(v[old][z+1][r], v[old][z][r+1]);
	    if (z > 0) {
	      v_sum += v[old][z-1][r]*eps_dz[z-1][r];
	      eps_sum += eps_dz[z-1][r];
	      min = fminf(min, v[old][z-1][r]);
	    } else {
	      v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	      eps_sum += eps_dz[z][r];
	    }
	    if (r > 0) {
	      v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	      eps_sum += eps_dr[z][r-1]*s2[r];
	      min = fminf(min, v[old][z][r-1]);
	    } else {
	      v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
	      eps_sum += eps_dr[z][r]*s1[r];
	    }

	  } else if (bulk[z][r] == 1) {    // interpolated radial edge of point contact
	    /* since the PC radius is not in the middle of a pixel,
	       use a modified weight for the interpolation to (r-1)
	     */
	    v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r] +
	            v[old][z][r-1]*eps_dr[z][r-1]*s2[r]*frrc[z];
	    eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dr[z][r-1]*s2[r]*frrc[z];
	    min = fminf(v[old][z+1][r], v[old][z][r+1]);
	    min = fminf(min, v[old][z][r-1]);
	    if (z > 0) {
	      v_sum += v[old][z-1][r]*eps_dz[z-1][r];
	      eps_sum += eps_dz[z-1][r];
	      min = fminf(min, v[old][z-1][r]);
	    } else {
	      v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	      eps_sum += eps_dz[z][r];
	    }
	  } else if (bulk[z][r] == 2) {    // interpolated z edge of point contact
	    /* since the PC length is not in the middle of a pixel,
	       use a modified weight for the interpolation to (z-1)
	     */
	    v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r] +
	            v[old][z-1][r]*eps_dz[z-1][r]*fLC;
	    eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dz[z-1][r]*fLC;
	    min = fminf(v[old][z+1][r], v[old][z][r+1]);
	    min = fminf(min, v[old][z-1][r]);
	    if (r > 0) {
	      v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	      eps_sum += eps_dr[z][r-1]*s2[r];
	      min = fminf(min, v[old][z][r-1]);
	    } else {
	      v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
	      eps_sum += eps_dr[z][r]*s1[r];
	    }
	    // check for cases where the PC corner needs modification in both r and z
	    if (z == LC && bulk[z-1][r] == 1) {
	      v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r]*(frrc[z]-1.0);
	      eps_sum += eps_dr[z][r-1]*s2[r]*(frrc[z]-1.0);
	      min = fminf(min, v[old][z][r-1]);
	    }

	  } else {
	    printf(" ERROR! bulk = %d undefined for (z,r) = (%d,%d)\n",
		   bulk[z][r], z, r);
	    return 1;
	  }

	  // calculate the interpolated mean potential and the effect of the space charge
	  mean = v_sum / eps_sum;
	  v[new][z][r] = mean + vfraction[z][r] * (imp_z[z]*imp_rm[r] + imp_ra[r]);
	  if (r == 0)  // special case where volume of voxel is 1/6 of area, not 1/4
	    v[new][z][r] = mean + (vfraction[z][r] * (imp_z[z]*imp_rm[r] + imp_ra[r])) / 1.5;
	  if ((z == 0 && r > RC && r < RO-WO) ||        // passivated surface at z = 0
              (z < LO && (r == RO || r == RO-WO-1)) ||  // passivated surface on sides of ditch
              (z == LO && r <= RO && r >= RO-WO-1))     // passivated surface at top of ditch
	    v[new][z][r] += vfraction[z][r] * S;
	  // check to see if the pixel is undepleted
	  if (vfraction[z][r] > 0.45) undepleted[r][z] = '.';
	  if (v[new][z][r] <= 0.0f) {
	    v[new][z][r] = 0.0f;
	    if (vfraction[z][r] > 0.45) undepleted[r][z] = '*';
	  } else if (v[new][z][r] < min) {
	    if (bubble_volts == 0.0f) bubble_volts = min + 0.1f;
	    v[new][z][r] = bubble_volts;
	    if (vfraction[z][r] > 0.45) undepleted[r][z] = '*';
	  }
	  // calculate difference from last iteration, for convergence check
	  dif = v[old][z][r] - v[new][z][r];
	  if (dif < 0.0f) dif = -dif;
	  sum_dif += dif;
	  if (max_dif < dif) max_dif = dif;
	}
      }
      // report results for some iterations
      if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0)
	printf("%5d %d %d %.10f %.10f\n", iter, old, new, max_dif, sum_dif/(float) (L*R));
      if (max_dif < 0.000000001) break;
    }

    printf("\n>> %d %.16f\n\n", iter, sum_dif);

    fully_depleted = 1;
    for (r=0; r<R+1; r++) {
      for (z=0; z<L+1; z++) {
	if (undepleted[r][z] == '*') {
	  fully_depleted = 0;
	  if (v[new][z][r] > 0.001) undepleted[r][z] = 'B';  // identifies pinch-off
	}
      }
    }
    if (fully_depleted) {
      printf("Detector is fully depleted.\n");
    } else {
      printf("Detector is not fully depleted.\n");
      if (bubble_volts > 0.0f) printf("Pinch-off bubble at %.0f V potential\n", bubble_volts);
    }
    if (setup.verbosity >= CHATTY) {
      t1 = time(NULL);
      printf("\n ^^^^^^^^^^^^^ %d (%d) s elapsed ^^^^^^^^^^^^^^\n",
	     (int) (t1 - t0), (int) (t1 - t2));
      t2 = t1;
    }

    if (istep == 0) {
      // can reduce # of iterations after first go-around
      max_its /= MAX_ITS_FACTOR;
      // report V and E along the axes r=0 and z=0
      if (setup.verbosity >= NORMAL) {
	printf("  z(mm)(r=0)      V   E(V/cm) |  r(mm)(z=0)      V   E(V/cm)\n");
	a = b = v[new][0][0];
	for (z=0; z<L+1; z++) {
	  printf("%10.1f %8.1f %8.1f  |",
		 ((float) z)*grid, v[new][z][0], (v[new][z][0] - a)/(0.1*grid));
	  a = v[new][z][0];
	  if (z > R) {
	    printf("\n");
	  } else {
	    r = z;
	    printf("%10.1f %8.1f %8.1f\n",
		   ((float) r)*grid, v[new][0][r], (v[new][0][r] - b)/(0.1*grid));
	    b = v[new][0][r];
	  }
	}
      }
      // write a little file that shows any undepleted voxels in the crystal
      file = fopen("undepleted.txt", "w");
      for (r=R; r>=0; r--) {
	undepleted[r][L] = '\0';
	fprintf(file, "%s\n", undepleted[r]);
      }
      fclose(file);
    }
  }

  if (WV) {
    if (setup.impurity_z0 > 0) {
      // swap voltages back to negative for n-type material
      for (r=0; r<R+1; r++) {
	for (z=0; z<L+1; z++) {
	  v[new][z][r] = -v[new][z][r];
	}
      }
    }
    // write potential and field to output file
    if (!(file = fopen(setup.field_name, "w"))) {
      printf("ERROR: Cannot open file %s for electric field...\n", setup.field_name);
      return 1;
    } else {
      printf("Writing electric field data to file %s\n", setup.field_name);
    }
    /* copy configuration parameters to output file */
    report_config(file, config_file_name);
    fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", BV);
    if (fully_depleted) {
      fprintf(file, "# Detector is fully depleted.\n");
    } else {
      fprintf(file, "# Detector is not fully depleted.\n");
      if (bubble_volts > 0.0f) fprintf(file, "# Pinch-off bubble at %.0f V potential\n", bubble_volts);
    }
    fprintf(file, "#\n## r (mm), z (mm), V (V),  E (V/cm), E_r (V/cm), E_z (V/cm)\n");

    for (r=0; r<R+1; r++) {
      for (z=0; z<L+1; z++) {
	// calc E in r-direction
	if (r==0) {
	  // E_r = (v[new][z][r] - v[new][z][r+1])/(0.1*grid);
	  E_r = 0;
	} else if (r==R) {
	  E_r = (v[new][z][r-1] - v[new][z][r])/(0.1*grid);
	} else {
	  E_r = (v[new][z][r-1] - v[new][z][r+1])/(0.2*grid);
	}
	// calc E in z-direction
	if (z==0) {
	  E_z = (v[new][z][r] - v[new][z+1][r])/(0.1*grid);
	} else if (z==L) {
	  E_z = (v[new][z-1][r] - v[new][z][r])/(0.1*grid);
	} else {
	  E_z = (v[new][z-1][r] - v[new][z+1][r])/(0.2*grid);
	}
	fprintf(file, "%7.2f %7.2f %7.1f %7.1f %7.1f %7.1f\n",
		((float) r)*grid,  ((float) z)*grid, v[new][z][r],
		sqrt(E_r*E_r + E_z*E_z), E_r, E_z);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }


  if (WP == 0) return 0;
  /*
    -------------------------------------------------------------------------
    now calculate the weighting potential for the central contact
    the WP is also needed for calculating the capacitance
    -------------------------------------------------------------------------
  */

  printf("\nCalculating weighting potential...\n\n");
  if (setup.verbosity >= CHATTY) t0 = t2 = time(NULL);
  max_its = MAX_ITS;
  if (setup.max_iterations > 0) max_its = setup.max_iterations;
  // max_its = 2*MAX_ITS;  // use twice as many iterations for WP; accuracy is more important?
  // if (setup.max_iterations > 0) max_its = 2*setup.max_iterations;

  /* to be safe, initialize overall potential to 0 */
  for (z=0; z<LL+1; z++) {
    for (r=0; r<RR+1; r++) {
      v[0][z][r] = v[1][z][r] = 0;
    }
  }
  for (istep=0; istep<3 && gridstep[istep]>0; istep++) {
    grid = gridstep[istep];
    old = 1;
    new = 0;
    // gridfact = integer ratio of current grid step size to final grid step size
    gridfact = lrintf(grid / setup.xtal_grid);

    if (istep > 0) {
      /* the previous calculation was on a coarser grid...
	 now copy/expand the potential to the new finer grid
      */
      i = (int) (gridstep[istep-1] / gridstep[istep] + 0.5);
      f = 1.0 / (float) i;
      printf("\ngrid %.4f -> %.4f; ratio = %d %.3f\n\n",
	     gridstep[istep-1], gridstep[istep], i, f);
      for (z=0; z<L+1; z++) {
	for (r=0; r<R+1; r++) {
	  f1z = 0.0;
	  zmax = i*z+i;
	  if (zmax > LL+1) zmax = LL+1;
	  for (zz=i*z; zz<zmax; zz++) {
	    f2z = 1.0 - f1z;
	    f1r = 0.0;
	    rmax = i*r+i;
	    if (rmax > RR+1) rmax = RR+1;
	    for (rr=i*r; rr<rmax; rr++) {
	      f2r = 1.0 - f1r;
	      v[0][zz][rr] =      // linear interpolation
		f2z*f2r*v[1][z][r  ] + f1z*f2r*v[1][z+1][r  ] +
		f2z*f1r*v[1][z][r+1] + f1z*f1r*v[1][z+1][r+1];
	      f1r += f;
	    }
	    f1z += f;
	  }
	}
      }
    }

    L  = lrint(setup.xtal_length/grid);
    R  = lrint(setup.xtal_radius/grid);
    // BRT = lrint(setup.top_bullet_radius/grid);
    // BRB = lrint(setup.bottom_bullet_radius/grid);
    LC = lrint(setup.pc_length/grid);
    dLC = setup.pc_length/grid - (float) LC;
    if (dLC < 0.05 && dLC > -0.05) dLC = 0;
    RC = lrint(setup.pc_radius/grid);
    dRC = setup.pc_radius/grid - (float) RC;
    if (dRC < 0.05 && dRC > -0.05) dRC = 0;
    printf("grid = %f  RC = %d  dRC = %f  LC = %d  dLC = %f\n\n",
	   grid, RC, dRC, LC, dLC);
    /* set up bulletization inside point contact */
    if (setup.bulletize_PC) {
      for (z=0; z<=LC; z++) {
        if (setup.pc_length <= setup.pc_radius) {  // LC <= RC; use LC as bulletization radius
          a = setup.pc_radius - setup.pc_length;
          b = z * grid;
          c = setup.pc_length*setup.pc_length - b*b;
          if (c < 0.0) c = 0;
          c = a + sqrt(c);
        } else {  // LC > RC; use RC as bulletization radius
          if (z > LC-RC) {
            a = setup.pc_length - setup.pc_radius;
            b = z * grid - a;
            c = setup.pc_radius*setup.pc_radius - b*b;
            if (c < 0.0) c = 0;
            c = sqrt(c);
          } else {
            c = setup.pc_radius;
          }
        }
        rrc[z] = lrint(c/grid);
        drrc[z] = c/grid - (float) rrc[z];
        if (drrc[z] < 0.05 && drrc[z] > -0.05) drrc[z] = 0;
        frrc[z] = 0;
      }
      drrc[LC+1] = drrc[LC];
    } else {  // no bulletization
      for (z=0; z<=LC+1; z++) {
        rrc[z] = RC;
        drrc[z] = dRC;
        frrc[z] = 0;
      }
    }
    
    LT = lrint(setup.taper_length/grid);
    RO = lrint(setup.wrap_around_radius/grid);
    LO = lrint(setup.ditch_depth/grid);
    WO = lrint(setup.ditch_thickness/grid);
    // LiT = lrint(setup.Li_thickness/grid);
    if (RO <= 0.0 || RO >= R) RO = R - LT;    // inner radius of taper, in grid lengths

    if (istep == 0) {
      // no previous coarse relaxation, so set initial potential:
      for (z=0; z<L+1; z++) {
	for (r=0; r<R+1; r++) {
	  v[0][z][r] = v[1][z][r] = 0.0;
	}
      }
      /*  ----- can comment out this next section to test convergence of WP ----- */
      // perhaps this is a better initial guess than just zero everywhere
      a = LC + RC / 2;
      b = 2.0f*a / (float) (L + R);
      for (z=1; z<L; z++) {
	for (r=1; r<R; r++) {
	  c = a / sqrt(z*z + r*r) - b;
	  if (c < 0) c = 0;
	  if (c > 1) c = 1;
	  v[0][z][r] = v[1][z][r] = c;
	}
      }
      /* ----- ----- */
      // inside contact:
      for (z=0; z<LC+1; z++) {
	for (r=0; r<rrc[z]+1; r++) {
	  v[0][z][r] = v[1][z][r] = 1.0;
	}
      }
    }

    /* boundary conditions and permittivity
       boundary condition at Ge-vacuum interface:
       epsilon0 * E_vac = espilon_Ge * E_Ge
    */
    for (z=0; z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	eps[z][r] = eps_dz[z][r] = eps_dr[z][r] = 16;   // permittivity inside Ge
	if (z < LO  && r < RO && r > RO-WO-1) eps[z][r] =  1;  // permittivity inside vacuum
	if (r > 0) eps_dr[z][r-1] = (eps[z][r-1]+eps[z][r])/2.0f;
	if (z > 0) eps_dz[z-1][r] = (eps[z-1][r]+eps[z][r])/2.0f;
      }
    }

    for (z=0; z<L+1; z++) {
      for (r=0; r<R+1; r++) {
	// boundary conditions
	bulk[z][r] = 0;  // normal bulk, no complications
	// outside (HV) contact:
	if (z == L ||
	    r == R ||
	    r >= z + R - LT ||       // taper
	    (z == 0 && r >= RO)) {   // wrap-around
	  bulk[z][r] = -1;                 // value of v[*][z][r] is fixed...
	  v[0][z][r] = v[1][z][r] = 0.0;   // to zero
	}
	// inside (point) contact:
	else if (z <= LC && r <= rrc[z]) {
	  bulk[z][r] = -1;                 // value of v[*][z][r] is fixed...
	  v[0][z][r] = v[1][z][r] = 1.0;   // to 1.0
	  // radial edge of inside contact:
	  if (r == rrc[z] && drrc[z] < -0.05) {
	    bulk[z][r] = 1;
	    frrc[z] = -1.0/drrc[z];
	  }
	  // z edge of inside contact:
	  if (z == LC && dLC < -0.05) {
	    bulk[z][r] = 2;
	    fLC = -1.0/dLC;
	  }
	}
	// edge of inside contact:
	// FIXME: Check for adjacent ditch
	else if (z <= LC && r == rrc[z]+1 && drrc[z] > 0.05) {
	  bulk[z][r] = 1;
	  frrc[z] = 1.0/(1.0 - drrc[z]);
	}
	else if (z == LC+1 && r <= rrc[z] && dLC > 0.05) {
	  bulk[z][r] = 2;
	  fLC = 1.0/(1.0 - dLC);
	}

	/* determine bulk regions where the detector is undepleted */
	if (!fully_depleted) {
	  if (undepleted[r*gridfact][z*gridfact] == '*') {
	    bulk[z][r] = -1;	            // treat like part of point contact
	    v[0][z][r] = v[1][z][r] = 1.0;  // set WP to one
	  } else if (undepleted[r*gridfact][z*gridfact] == 'B') { // pinch-off
	    bulk[z][r] = 3;
	  }
	}
      }
    }

    // now do the actual relaxation
    for (iter=0; iter<max_its; iter++) {
      if (old == 0) {
	old = 1;
	new = 0;
      } else {
	old = 0;
	new = 1;
      }
      sum_dif = 0.0f;
      max_dif = 0.0f;
      pinched_sum1 = pinched_sum2 = 0.0;

      for (z=0; z<L; z++) {
	for (r=0; r<R; r++) {
	  if (bulk[z][r] < 0) continue;      // outside or inside contact

	  if (bulk[z][r] == 0) {            // normal bulk, no complications
	    v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r];
	    eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r];
	    if (z > 0) {
	      v_sum += v[old][z-1][r]*eps_dz[z-1][r];
	      eps_sum += eps_dz[z-1][r];
	    } else {
	      v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	      eps_sum += eps_dz[z][r];
	    }
	    if (r > 0) {
	      v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	      eps_sum += eps_dr[z][r-1]*s2[r];
	    } else {
	      v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
	      eps_sum += eps_dr[z][r]*s1[r];
	    }

	  } else if (bulk[z][r] == 1) {    // interpolated radial edge of point contact
	    v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r] +
	      v[old][z][r-1]*eps_dr[z][r-1]*s2[r]*frrc[z];
	    eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dr[z][r-1]*s2[r]*frrc[z];
	    if (z > 0) {
	      v_sum += v[old][z-1][r]*eps_dz[z-1][r];
	      eps_sum += eps_dz[z-1][r];
	    } else {
	      v_sum += v[old][z+1][r]*eps_dz[z][r];  // reflection symm around z=0
	      eps_sum += eps_dz[z][r];
	    }
	  } else if (bulk[z][r] == 2) {    // interpolated z edge of point contact
	    v_sum = v[old][z+1][r]*eps_dz[z][r] + v[old][z][r+1]*eps_dr[z][r]*s1[r] +
	      v[old][z-1][r]*eps_dz[z-1][r]*fLC;
	    eps_sum = eps_dz[z][r] + eps_dr[z][r]*s1[r] + eps_dz[z-1][r]*fLC;
	    if (r > 0) {
	      v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	      eps_sum += eps_dr[z][r-1]*s2[r];
	    } else {
	      v_sum += v[old][z][r+1]*eps_dr[z][r]*s1[r];  // reflection symm around r=0
	      eps_sum += eps_dr[z][r]*s1[r];
	    }
	    if (z == LC && bulk[z-1][r] == 1) {
	      v_sum += v[old][z][r-1]*eps_dr[z][r-1]*s2[r]*(frrc[z]-1.0);
	      eps_sum += eps_dr[z][r-1]*s2[r]*(frrc[z]-1.0);
	    }

	  } else if (bulk[z][r] == 3) {   // pinched-off
	    if (bulk[z+1][r] == 0) {
	      pinched_sum1 += v[old][z+1][r]*eps_dz[z][r];
	      pinched_sum2 += eps_dz[z][r];
	    }
	    if (bulk[z][r+1] == 0) {
	      pinched_sum1 += v[old][z][r+1]*eps_dr[z][r]*s1[r];
	      pinched_sum2 += eps_dr[z][r]*s1[r];
	    }
	    if (z > 0 && bulk[z-1][r] == 0) {
	      pinched_sum1 += v[old][z-1][r]*eps_dz[z-1][r];
	      pinched_sum2 += eps_dz[z-1][r];
	    }
	    if (r > 0 && bulk[z][r-1] == 0) {
	      pinched_sum1 += v[old][z][r-1]*eps_dr[z][r-1]*s2[r];
	      pinched_sum2 += eps_dr[z][r-1]*s2[r];
	    }
	    v_sum = pinched_sum1;
	    eps_sum = pinched_sum2;

	  } else {
	    printf(" ERROR! bulk = %d undefined for (z,r) = (%d,%d)\n",
		   bulk[z][r], z, r);
	    return 1;
	  }
	  if (bulk[z][r] != 3) {
	    mean = v_sum / eps_sum;
	    v[new][z][r] = mean;
	    dif = v[old][z][r] - v[new][z][r];
	    if (dif < 0.0f) dif = -dif;
	    sum_dif += dif;
	    if (max_dif < dif) max_dif = dif;
	  }
	}
      }

      if (pinched_sum2 > 0.1) {
	mean = pinched_sum1 / pinched_sum2;
	for (z=0; z<L; z++) {
	  for (r=0; r<R; r++) {
	    if (bulk[z][r] == 3) {
	      v[new][z][r] = mean;
	      dif = v[old][z][r] - v[new][z][r];
	      if (dif < 0.0f) dif = -dif;
	      sum_dif += dif;
	      if (max_dif < dif) max_dif = dif;
	    }
	  }
	}
      }

      // report results for some iterations
      if (iter < 10 || (iter < 600 && iter%100 == 0) || iter%1000 == 0)
	printf("%5d %d %d %.10f %.10f ; %.10f %.10f\n",
	       iter, old, new, max_dif, sum_dif/(float) (L*R),
	       v[new][L/2][R/2], v[new][L-5][R-5]);
      if (max_dif < 0.0000000001) break;
    }
    printf(">> %d %.16f\n\n", iter, sum_dif);
    if (setup.verbosity >= CHATTY) {
      t1 = time(NULL);
      printf(" ^^^^^^^^^^^^^ %d (%d) s elapsed ^^^^^^^^^^^^^^\n",
	     (int) (t1 - t0), (int) (t1 - t2));
      t2 = t1;
    }
    if (istep == 0) max_its /= MAX_ITS_FACTOR;
  }

  /* --------------------- calculate capacitance ---------------------
     1/2 * epsilon * integral(E^2) = 1/2 * C * V^2
     so    C = epsilon * integral(E^2) / V^2
     V = 1 volt
  */
  printf("Calculating integrals of weighting field\n");
  esum = esum2 = j = 0;
  for (z=0; z<L; z++) {
    for (r=1; r<R; r++) {
      E_r = eps_dr[z][r]/16.0 * (v[new][z][r] - v[new][z][r+1])/(0.1*grid);
      E_z = eps_dz[z][r]/16.0 * (v[new][z][r] - v[new][z+1][r])/(0.1*grid);
      esum += (E_r*E_r + E_z*E_z) * (double) r;
      /*
      if ((z <= LC && r == rrc[z]) ||
	  (z == LC && r <= rrc[z]) ||
	  (z <= LC+1 && r == rrc[z]+1) || // average over two different surfaces
	  (z == LC+1 && r <= rrc[z]+1)) {
      */
      if ((r == RC && z <= LC) ||
	  (r <= RC && z == LC) ||
	  (r == RC+1 && z <= LC+1) || // average over two different surfaces
	  (r <= RC+1 && z == LC+1)) {
	if (bulk[z+1][r+1] < 0) j = 1;
	esum2 += 0.5 * sqrt(E_r*E_r + E_z*E_z) * (double) r;  // 0.5 since averaging over 2 surfaces
      }
    }
  }
  esum  *= 2.0 * pi * 0.01 * Epsilon * pow(grid, 3.0);
  // Epsilon is in pF/mm
  // 0.01 converts (V/cm)^2 to (V/mm)^2, pow() converts to grid^3 to mm3
  esum2 *= 2.0 * pi * 0.1 * Epsilon * pow(grid, 2.0);
  // 0.1 converts (V/cm) to (V/mm),  grid^2 to  mm2
  printf("\n  >>  Calculated capacitance at %.0f V: %.3lf pF\n", BV, esum);
  if (j==0) {
    printf("  >>  Alternative calculation of capacitance: %.3lf pF\n\n", esum2);
  } else {
    printf("\n");
  }

  if (WP == 1) {
    // write WP values to output file
    if (!(file = fopen(setup.wp_name, "w"))) {
      printf("ERROR: Cannot open file %s for weighting potential...\n", setup.wp_name);
      return 1;
    } else {
      printf("Writing weighting potential to file %s\n", setup.wp_name);
    }
    /* copy configuration parameters to output file */
    report_config(file, config_file_name);
    fprintf(file, "#\n# HV bias in fieldgen: %.1f V\n", BV);
    if (fully_depleted) {
      fprintf(file, "# Detector is fully depleted.\n");
    } else {
      fprintf(file, "# Detector is not fully depleted.\n");
      if (bubble_volts > 0.0f) fprintf(file, "# Pinch-off bubble at %.0f V potential\n", bubble_volts);
    }
    fprintf(file, "#\n## r (mm), z (mm), WP\n");
    for (r=0; r<R+1; r++) {
      for (z=0; z<L+1; z++) {
	fprintf(file, "%7.2f %7.2f %10.6f\n",
		((float) r)*grid,  ((float) z)*grid, v[new][z][r]);
      }
      fprintf(file, "\n");
    }
    fclose(file);
  }

  return 0;
}

int report_config(FILE *fp_out, char *config_file_name) {

  char  *c, line[256];
  FILE  *file;

  fprintf(fp_out, "# Config file: %s\n", config_file_name);
  if (!(file = fopen(config_file_name, "r"))) return 1;

  while (fgets(line, sizeof(line), file)) {
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    if ((c = strchr(line, '#')) || (c = strchr(line, '\n'))) *c = '\0';
    fprintf(fp_out, "# %s\n", line);
  }
  fclose(file);
  return 0;
}
