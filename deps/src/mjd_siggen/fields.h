/* fields.h -- based on m3d2s.f by I-Yang Lee
 * Karin Lagergren
 *
 * This module handles the electric field and weighting potential and 
 * calculates drift velocities
 */

#ifndef _FIELDS_H
#define _FIELDS_H

/* calculate anisotropic drift velocities? (vel. depends on angle between
   el. field and crystal axis; otherwise the velocity will always be 
   in the direction of the el. field 
*/
#define DRIFT_VEL_ANISOTROPY 1

#include "point.h"
#include "mjd_siggen.h"

/* nearest_field_grid_index
   find existing integer field grid index closest to pt
   returns <0 if outside crystal or too far from a valid grid point
            0 if interpolation is okay
            1 if we can find a point but extrapolation is needed
*/
int nearest_field_grid_index(cyl_pt pt, cyl_int_pt *ipt, MJD_Siggen_Setup *setup);

/* field_setup
   given a field directory file, read electic field and weighting
   potential tables from files listed in directory
   returns 0 for success
*/
int field_setup(MJD_Siggen_Setup *setup);

/* free malloc()'ed memory and do other cleanup*/
int fields_finalize(MJD_Siggen_Setup *setup);

/* wpotential
   gives (interpolated or extrapolated ) weighting potential
   at point pt. These values are stored in wp.
   returns 0 for success, 1 on failure.
*/
int wpotential(point pt, float *wp, MJD_Siggen_Setup *setup);

/* drift_velocity
   calculates drift velocity for charge q at point pt
   returns 0 on success, 1 if successful but extrapolation was needed,
   and -1 for failure
*/
int drift_velocity(point pt, float q, vector *velocity, MJD_Siggen_Setup *setup);

/* TC: get drift_velocity with linear superposition of 2 fields:
		1. Detector field, E
		2. Charge cloud field, Eadd
*/
int drift_velocity_w_Eadd(point pt, float q, vector *velo, MJD_Siggen_Setup *setup, cyl_pt Eadd);


int read_fields(MJD_Siggen_Setup *setup);

/*set detector temperature. 77F (no correction) is the default
   MIN_TEMP & MAX_TEMP defines allowed range*/
void set_temp(float temp, MJD_Siggen_Setup *setup);


#endif /*#ifndef _FIELDS_H*/
