#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "mjd_siggen.h"

int read_config(char *config_file_name, MJD_Siggen_Setup *setup) {

  /* reads and parses configuration file of name config_file_name
     fills in values of MJD_Siggen_Setup setup, defined in mjd_siggen.h
     returns 0 on success, 1 otherwise
  */

  char key_word[][32] = {
    "xtal_length",
    "xtal_radius",
    "top_bullet_radius",
    "bottom_bullet_radius",
    "pc_length",
    "pc_radius",
    "bulletize_PC",
    "taper_length",
    "wrap_around_radius",
    "ditch_depth",
    "ditch_thickness",
    "Li_thickness",
    "xtal_grid",
    "impurity_z0",
    "impurity_gradient",
    "impurity_quadratic",
    "impurity_surface",
    "impurity_radial_add",
    "impurity_radial_mult",
    "impurity_rpower",
    "xtal_HV",
    "drift_name",
    "field_name",
    "wp_name",
    "xtal_temp",
    "preamp_tau",
    "time_steps_calc",
    "step_time_calc",
    "step_time_out",
    "charge_cloud_size",
    "use_diffusion",
    "energy",
    "verbosity_level",
    "max_iterations",
    "write_field",
    "write_WP",
    ""
  };

  int   ii, i, l, n=0, ok, iint = 0;
  float fi;
  char  *c, line[256], name[256];
  FILE  *file;


  /* initialize everything to zero... */
  memset(setup, 0, sizeof(*setup));
  /* ...except for impurity_radial_mult */
  setup->impurity_radial_mult = 1.0f;  // 1.0 is neutral (no radial gradient)

  if (!(file = fopen(config_file_name, "r"))) {
    printf("\nERROR: config file %s does not exist?\n", config_file_name);
    return 1;
  }
  /* read config file */
  printf("\nReading values from config file %s\n", config_file_name);
  while (fgets(line, sizeof(line), file)) {
    n++;
    /* ignore comments and blank lines */
    if (strlen(line) < 3 || *line == ' ' || *line == '\t' || *line == '#') continue;
    for (i=0; (l=strlen(key_word[i])) > 0; i++) {
      if (!strncmp(line, key_word[i], l)) {
	/* line starts with key_word[i] */
	if (line[l] != ' ' && line[l] != '\t') {
	  ok = 0;
	} else {
	  /* find next non-white-space char */
	  for (c = line + l; *c == ' ' || *c == '\t'; c++) ;
	  name[0] = 0;
	  ii = iint = 0;
	  fi = 0;
	  if (strstr(key_word[i], "_name")) {
	    /* extract character string for file name */
	    for (ok=0; ok<256 && *c != ' ' && *c != '\t' && *c != '\n' &&  *c != '\r'; ok++) {
	      name[ok] = *c;
	      c++;
	    }
	    name[ok] = '\0';	    
	  } else if (!strncmp("time_steps_calc", key_word[i], l) ||
		     !strncmp("use_diffusion", key_word[i], l) ||
		     !strncmp("verbosity_level", key_word[i], l) ||
		     !strncmp("max_iterations", key_word[i], l) ||
		     !strncmp("write_field", key_word[i], l) ||
		     !strncmp("write_WP", key_word[i], l) ||
		     !strncmp("bulletize_PC", key_word[i], l)) {
	    /* extract integer value */
	    ok = sscanf(c, "%d", &ii);
	    iint = 1;
	  } else {
	    /* extract float value */
	    ok = sscanf(c, "%f", &fi);
	  }
	}
	if (ok < 1) {
	  printf("ERROR reading %s from config file %s\n"
		 "   ...line number %d is: %s",
		 key_word[i], config_file_name, n, line);
	  return 1;
	}
	if (strstr(key_word[i], "verbosity_level")) {
	  setup->verbosity = ii;
	} else if (strstr(key_word[i], "xtal_length")) {
	  setup->xtal_length = fi;
	} else if (strstr(key_word[i], "xtal_radius")) {
	  setup->xtal_radius = fi;
	} else if (strstr(key_word[i], "top_bullet_radius")) {
	  setup->top_bullet_radius = fi;
	} else if (strstr(key_word[i], "bottom_bullet_radius")) {
	  setup->bottom_bullet_radius = fi;
	} else if (strstr(key_word[i], "pc_length")) {
	  setup->pc_length = fi;
	} else if (strstr(key_word[i], "pc_radius")) {
	  setup->pc_radius = fi;
	} else if (strstr(key_word[i], "bulletize_PC")) {
	  setup->bulletize_PC = ii;
	} else if (strstr(key_word[i], "taper_length")) {
	  setup->taper_length = fi;
	} else if (strstr(key_word[i], "wrap_around_radius")) {
	  setup->wrap_around_radius = fi;
	} else if (strstr(key_word[i], "ditch_depth")) {
	  setup->ditch_depth = fi;
	} else if (strstr(key_word[i], "ditch_thickness")) {
	  setup->ditch_thickness = fi;
	} else if (strstr(key_word[i], "Li_thickness")) {
	  setup->Li_thickness = fi;
	} else if (strstr(key_word[i], "xtal_grid")) {
	  setup->xtal_grid = fi;
	} else if (strstr(key_word[i], "impurity_z0")) {
	  setup->impurity_z0 = fi;
	} else if (strstr(key_word[i], "impurity_gradient")) {
	  setup->impurity_gradient = fi;
	} else if (strstr(key_word[i], "impurity_quadratic")) {
	  setup->impurity_quadratic = fi;
	} else if (strstr(key_word[i], "impurity_surface")) {
	  setup->impurity_surface = fi;
	} else if (strstr(key_word[i], "impurity_radial_add")) {
	  setup->impurity_radial_add = fi;
	} else if (strstr(key_word[i], "impurity_radial_mult")) {
	  setup->impurity_radial_mult = fi;
	} else if (strstr(key_word[i], "impurity_rpower")) {
	  setup->impurity_rpower = fi;
	} else if (strstr(key_word[i], "xtal_HV")) {
	  setup->xtal_HV = fi;
	} else if (strstr(key_word[i], "drift_name")) {
	  strncpy(setup->drift_name, name, 256);
	} else if (strstr(key_word[i], "field_name")) {
	  strncpy(setup->field_name, name, 256);
	} else if (strstr(key_word[i], "wp_name")) {
	  strncpy(setup->wp_name, name, 256);
	} else if (strstr(key_word[i], "xtal_temp")) {
	  setup->xtal_temp = fi;
	} else if (strstr(key_word[i], "preamp_tau")) {
	  setup->preamp_tau = fi;
	} else if (strstr(key_word[i], "time_steps_calc")) {
	  setup->time_steps_calc = ii;
	} else if (strstr(key_word[i], "step_time_calc")) {
	  setup->step_time_calc = fi;
	} else if (strstr(key_word[i], "step_time_out")) {
	  setup->step_time_out = fi;
	} else if (strstr(key_word[i], "charge_cloud_size")) {
	  setup->charge_cloud_size = fi;
	} else if (strstr(key_word[i], "use_diffusion")) {
	  setup->use_diffusion = ii;
	} else if (strstr(key_word[i], "energy")) {
	  setup->energy = fi;
	} else if (strstr(key_word[i], "max_iterations")) {
	  setup->max_iterations = ii;
	} else if (strstr(key_word[i], "write_field")) {
	  setup->write_field = ii;
	} else if (strstr(key_word[i], "write_WP")) {
	  setup->write_WP = ii;
	} else {
	  printf("ERROR; unrecognized keyword %s\n", key_word[i]);
	  return 1;
	}

	if (setup->verbosity >= CHATTY) {
	  // printf("%s", line);
	  if (iint) {
	    printf("%s: %d\n", key_word[i], ii);
	  } else if (strlen(name) > 0) {
	    printf("%s: %s\n", key_word[i], name);
	  } else {
	    printf("%s: %f\n", key_word[i], fi);
	  }
	}
	break;
      }
    }
  }
  fclose(file);

  return 0;
}
