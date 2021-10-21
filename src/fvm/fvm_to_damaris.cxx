/*============================================================================
 * Write a nodal representation associated with a mesh and associated
 * variables to CoProcessing stats daemon output
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2020 EDF S.A.

  This program is free software; you can redistribute it and/or modify it under
  the terms of the GNU General Public License as published by the Free Software
  Foundation; either version 2 of the License, or (at your option) any later
  version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
  FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more
  details.

  You should have received a copy of the GNU General Public License along with
  this program; if not, write to the Free Software Foundation, Inc., 51 Franklin
  Street, Fifth Floor, Boston, MA 02110-1301, USA.
*/

/*----------------------------------------------------------------------------*/

#include "cs_defs.h"

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <assert.h>
#include <errno.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <climits>
/*----------------------------------------------------------------------------
 * Asynchronous I/O and in-situ visulisation library header
 *----------------------------------------------------------------------------*/

#if defined(HAVE_DAMARIS)
#include <Damaris.h>
//#include <damaris/paraview/ParaViewHeaders.hpp>
#endif


/*----------------------------------------------------------------------------
 * Catalyst and VTK library headers
 *----------------------------------------------------------------------------*/


/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "bft_error.h"
#include "bft_mem.h"

#include "fvm_defs.h"
#include "fvm_io_num.h"
#include "fvm_nodal.h"
#include "fvm_nodal_priv.h"
#include "fvm_writer_helper.h"
#include "fvm_writer_priv.h"

#include "cs_map.h"
#include "cs_parall.h"

/*----------------------------------------------------------------------------
 *  Header for the current file
 *----------------------------------------------------------------------------*/

#include "fvm_to_damaris.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*============================================================================
 * Local Type Definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Damaris output writer structure - based on Melissa
 *----------------------------------------------------------------------------*/

typedef struct {

  char        *name;           /* Writer name */

  bool         usm_on;         /* Input xml file switch to initialize unstructured
  	                              mesh in fvm_to_damaris_export_nodal() function.
                                  <writer>
                                  <format name="damaris" options="usm_on"/>
                                  </writer>
                                   Default is off if the string is not present */

  bool         rect_grid_on;    /* If true then run original Damaris rectilinear grid
   	   	   	   	   	   	   	       setup code for coord/xmesh etc. */

  FILE        *tracefile;       /* optional file for tracing */

//  int          rank;          /* Rank of current process in communicator */
//  int          n_ranks;       /* Number of processes in communicator */

//  size_t       buffer_size;   /* buffer size required */

  int          time_step;       /* Latest time step */
  double       time_value;      /* Latest time value */

  cs_map_name_to_id_t  *f_map;  /* field names mapping */
  int         *f_ts;            /* last field output time step */

  bool        ensight_names;    /* Use EnSight rules for
                                   field names */

#if defined(HAVE_MPI)
//  int          min_rank_step;     /* Minimum rank step */
//  int          min_block_size;    /* Minimum block buffer size */
//  MPI_Comm     block_comm;        /* Associated MPI block communicator */
    MPI_Comm     damaris_mpi_comm;  /* Associated MPI communicator */
#endif

  bool          modified;       /* Has output been added since
                                   last coprocessing ? */

} fvm_to_damaris_writer_t;

/*----------------------------------------------------------------------------
 * Context structure for fvm_writer_field_helper_output_* functions.
 *----------------------------------------------------------------------------*/

typedef struct {

  fvm_to_damaris_writer_t  *writer;      /* pointer to writer structure */

  const char               *name;        /* current field name */
  int                       time_step;   /* current_time_step */
//  bool                      call_init;   /* call damaris_init ? */

} _damaris_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Only one Damaris writer */

static fvm_to_damaris_writer_t  *_writer = NULL;


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to VTK.
 *
 * Output fields are non interlaced. Input arrays may be interlaced or not.
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated object handle
 *----------------------------------------------------------------------------*/
static void
_export_field_values_e(const fvm_nodal_t         *mesh,
                       const char                *fieldname,
                       int                        dim,
                       cs_interlace_t             interlace,
                       int                        n_parent_lists,
                       const cs_lnum_t            parent_num_shift[],
                       cs_datatype_t              datatype,
                       const void          *const field_values[])
{
	CS_UNUSED(interlace);
	CS_UNUSED(parent_num_shift);
	CS_UNUSED(datatype);
    int  section_id;

    const int dest_dim = (dim == 6) ? 9 : dim;

	cs_lnum_t start_id = 0;
	cs_lnum_t src_shift = 0;

	int damaris_err = DAMARIS_OK ;

	int n_elements_offset ;
	damaris_err = damaris_parameter_get("n_elements_offset",&n_elements_offset,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
	bft_error(__FILE__, __LINE__, damaris_err,
						 _("ERROR: Damaris damaris_parameter_get():\n"
						   "Parameter: \"%s\"."), "n_elements_offset");
	}

  // loop on sections which should be appended
  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

   // fvm_convert_array(dim,
   //                   0,
   //                   dest_dim,
   //                   src_shift,
   //                   section->n_elements + src_shift,
   //                   interlace,
   //                   datatype,
   //                   CS_DOUBLE,
   //                   n_parent_lists,
   //                   parent_num_shift,
   //                   section->parent_element_num,
   //                   field_values,
   //                   values + start_id);

    // copy input fieldname to lower case
    int t1 = 0;
    char fieldname_lwrcase[12];
    while ( (t1<11)  && fieldname[t1])
    {
    	fieldname_lwrcase[t1] = tolower(fieldname[t1]);
    	t1++;
    }


    if ( strncmp( fieldname_lwrcase, "pressure", 8) == 0 )
    {
		int64_t pos[1];

		pos[0] = n_elements_offset + start_id;

		damaris_err = damaris_set_position("fields/pressure" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/pressure");
		}
		damaris_err = damaris_write("fields/pressure" ,field_values[0]);
		//damaris_err = damaris_write("fields/pressure" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/pressure");
		}

    }
    else if ( strncmp( fieldname_lwrcase, "velocity", 8) == 0 )
    {
    	int64_t pos[2];

		pos[0] = 0;
		pos[1] = n_elements_offset + start_id;

		damaris_err = damaris_set_position("fields/velocity" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		damaris_err = damaris_write("fields/velocity" ,field_values[0]);
		//damaris_err = damaris_write("fields/velocity" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
    }
    else if ( strncmp( fieldname_lwrcase, "mpi_rank_id", 11) == 0 )
    {
		int64_t pos[1];

		pos[0] = n_elements_offset + start_id;

		damaris_err = damaris_set_position("fields/mpi_rank_id" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/mpi_rank_id");
		}
		damaris_err = damaris_write("fields/mpi_rank_id" ,field_values[0]);
		//damaris_err = damaris_write("fields/pressure" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/mpi_rank_id");
		}

	}

    start_id += section->n_elements*dest_dim;
    if (n_parent_lists == 0)
      src_shift += section->n_elements;


  }

  // Special case for symmetric tensors - not implemented

}


/*============================================================================
 * Private function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Write field values associated with element values of a nodal mesh to VTK.
 *
 * Output fields are non interlaced. Input arrays may be interlaced or not.
 *
 * parameters:
 *   mesh             <-- pointer to nodal mesh structure
 *   dim              <-- field dimension
 *   interlace        <-- indicates if field in memory is interlaced
 *   n_parent_lists   <-- indicates if field values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent list to common number index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   field_values     <-- array of associated field value arrays
 *   f                <-- associated object handle
 *----------------------------------------------------------------------------*/
/*static void
_export_field_values_e(const fvm_nodal_t         *mesh,
                       const char                *fieldname,
                       int                        dim,
                       cs_interlace_t             interlace,
                       int                        n_parent_lists,
                       const cs_lnum_t            parent_num_shift[],
                       cs_datatype_t              datatype,
                       const void          *const field_values[])
{

  int  section_id;


  const int dest_dim = (dim == 6) ? 9 : dim;


  // Distribute partition to block values

  cs_lnum_t start_id = 0;
  cs_lnum_t src_shift = 0;

  // loop on sections which should be appended
  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

   // fvm_convert_array(dim,
   //                   0,
   //                   dest_dim,
   //                   src_shift,
   //                   section->n_elements + src_shift,
   //                   interlace,
   //                   datatype,
   //                   CS_DOUBLE,
   //                   n_parent_lists,
   //                   parent_num_shift,
   //                   section->parent_element_num,
   //                   field_values,
   //                   values + start_id);

    // copy input fieldname to lower case
    int t1 = 0;
    char fieldname_lwrcase[9];
    while ( (t1<8)  && fieldname[t1])
    {
    	fieldname_lwrcase[t1] = tolower(fieldname[t1]);
    	t1++;
    }


    int damaris_err = DAMARIS_OK ;

    int param_z ;
	damaris_err = damaris_parameter_get("z",&param_z,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		bft_error(__FILE__, __LINE__, damaris_err,
							 _("ERROR: Damaris damaris_parameter_get():\n"
							   "Parameter: \"%s\"."), "z");
	}
	int param_x ;
	damaris_err = damaris_parameter_get("x",&param_x,sizeof(int));
	int param_y ;
	damaris_err = damaris_parameter_get("y",&param_y,sizeof(int));


    if ( strncmp( fieldname_lwrcase, "pressure", 8) == 0 )
    {
		int64_t pos[3];

		//// N.B. x,y,z == 0,1,2 indices
		pos[0] = 0;
		pos[1] = 0;
		pos[2] = cs_glob_rank_id*param_z/cs_glob_n_ranks;

		damaris_err = damaris_set_position("fields/pressure" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/pressure");
		}
		damaris_err = damaris_write("fields/pressure" ,field_values[0]);
		//damaris_err = damaris_write("fields/pressure" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/pressure");
		}

		// BFT_FREE(mypressure);

    }

    //else if (dim == 3)
    else if ( strncmp( fieldname_lwrcase, "velocity", 8) == 0 )
    {
    	int64_t pos[4];


		pos[0] = 0;
		pos[1] = 0;
		pos[2] = 0;
		pos[3] = (cs_glob_rank_id*param_z/cs_glob_n_ranks);

		damaris_err = damaris_set_position("fields/velocity" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		damaris_err = damaris_write("fields/velocity" ,field_values[0]);
		//damaris_err = damaris_write("fields/velocity" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		// damaris_end_iteration();
		// BFT_FREE(mypressure);
    }
    //else if ( strncmp( fieldname_lwrcase, "mpi_rank_id", 8) == 0 )

    start_id += section->n_elements*dest_dim;
    if (n_parent_lists == 0)
      src_shift += section->n_elements;


  }

  // Special case for symmetric tensors - not implemented

} */
/*
static void
_export_field_values_e(const fvm_nodal_t         *mesh,
                       const char                *fieldname,
                       int                        dim,
                       cs_interlace_t             interlace,
                       int                        n_parent_lists,
                       const cs_lnum_t            parent_num_shift[],
                       cs_datatype_t              datatype,
                       const void          *const field_values[])
{

  int  section_id;


  const int dest_dim = (dim == 6) ? 9 : dim;


  // Distribute partition to block values

  cs_lnum_t start_id = 0;
  cs_lnum_t src_shift = 0;

  // loop on sections which should be appended

  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

   // fvm_convert_array(dim,
    //                  0,
   //                   dest_dim,
   //                   src_shift,
   //                   section->n_elements + src_shift,
   //                   interlace,
   //                   datatype,
  //                    CS_DOUBLE,
  //                    n_parent_lists,
   //                   parent_num_shift,
   //                   section->parent_element_num,
   //                   field_values,
   //                   values + start_id);

    // copy input fieldname to lower case
    int t1 = 0;
    char fieldname_lwrcase[9];
    while ( (t1<8)  && fieldname[t1])
    {
    	fieldname_lwrcase[t1] = tolower(fieldname[t1]);
    	t1++;
    }


    int damaris_err = DAMARIS_OK ;

    int param_z ;
	damaris_err = damaris_parameter_get("z",&param_z,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		bft_error(__FILE__, __LINE__, damaris_err,
							 _("ERROR: Damaris damaris_parameter_get():\n"
							   "Parameter: \"%s\"."), "z");
	}
	int param_x ;
	damaris_err = damaris_parameter_get("x",&param_x,sizeof(int));
	int param_y ;
	damaris_err = damaris_parameter_get("y",&param_y,sizeof(int));

    if ( strncmp( fieldname_lwrcase, "pressure", 8) == 0 )
    {
		int64_t pos[3];

		//double * mypressure ;
		//double * mypressure_cpy ;
		//int x = param_x ;
		//int y = param_y ;
		//int z_persectn = param_z / cs_glob_n_ranks;
		//int size_sectn = x * y * z_persectn ;
		//BFT_MALLOC( mypressure, size_sectn, double );
		//mypressure_cpy = mypressure ;
		//for (int zc = cs_glob_rank_id*z_persectn ; zc < (cs_glob_rank_id+1)*z_persectn ; zc++) {
		//	for (int yc = 0; yc < y ; yc++){
		//		for (int xc = 0 ; xc < x ; xc++) {
		//			*mypressure_cpy = (double) ((zc * y * x)+(yc *x)+xc) ;
		//			 printf("%.1f", (*mypressure_cpy));
		//			mypressure_cpy++ ;
		//		}
		//		printf("\n");
		//	}
		//	printf("\n");
		//}
//
		//printf("\n");
//
		//// N.B. x,y,z == 0,1,2 indices
		pos[0] = 0;
		pos[1] = 0;
		pos[2] = cs_glob_rank_id*param_z/cs_glob_n_ranks;

		damaris_err = damaris_set_position("fields/pressure" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/pressure");
		}
		damaris_err = damaris_write("fields/pressure" ,field_values[0]);
		//damaris_err = damaris_write("fields/pressure" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/pressure");
		}

		// BFT_FREE(mypressure);

    }

    //else if (dim == 3)
    else if ( strncmp( fieldname_lwrcase, "velocity", 8) == 0 )
    {
    	int64_t pos[4];

    	//double * mypressure ;
    	//double * mypressure_cpy ;
    	//int v3 = 3 ;
		//int x = param_x ;
		//int y = param_y ;
		//int z_persectn = param_z / cs_glob_n_ranks;  // cs_glob_rank_id*z_persectn
		//int size_sectn = x * y * z_persectn *v3  ;
		//BFT_MALLOC( mypressure, size_sectn, double );
		//mypressure_cpy = mypressure ;
		//for (int zc =cs_glob_rank_id*z_persectn ; zc < (cs_glob_rank_id+1)*z_persectn ;  zc++)
		//	for (int yc = 0; yc < y ; yc++)
		//		for (int xc = 0 ; xc < x ; xc++)
		//			for (int vc = 0 ; vc < v3 ; vc++) {
		//				*mypressure_cpy = (double) ((zc ) * y * x * v3)+(yc *x*v3)+(xc*v3) + vc;
		//				mypressure_cpy++ ;
		//			}
		//

		pos[0] = 0;
		pos[1] = 0;
		pos[2] = 0;
		pos[3] = (cs_glob_rank_id*param_z/cs_glob_n_ranks);

		damaris_err = damaris_set_position("fields/velocity" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		damaris_err = damaris_write("fields/velocity" ,field_values[0]);
		//damaris_err = damaris_write("fields/velocity" ,mypressure);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		// damaris_end_iteration();
		// BFT_FREE(mypressure);
    }
    //else if ( strncmp( fieldname_lwrcase, "mpi_rank_id", 8) == 0 )

    start_id += section->n_elements*dest_dim;
    if (n_parent_lists == 0)
      src_shift += section->n_elements;


  }

  // Special case for symmetric tensors
  //if (dim == 6) {
  //  cs_lnum_t n_elts = f->GetNumberOfCells();
  //  for (cs_lnum_t i = 0; i < n_elts; i++) {
  //    values[9*i + 8] = values[9*i + 2];
  //    values[9*i + 7] = values[9*i + 4];
  ////    values[9*i + 6] = values[9*i + 5];
  //    values[9*i + 4] = values[9*i + 1];
  //    values[9*i + 2] = values[9*i + 5];
  //    values[9*i + 1] = values[9*i + 3];
  //    values[9*i + 5] = values[9*i + 7];
  //  }
  //}

}
*/

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function definitions
 *============================================================================*/

/*----------------------------------------------------------------------------
 * Initialize FVM to Damaris output writer.
 *
 * Options are:
 *   dry_run               trace output to <name>.log file, but do not
 *                         actually communicate with Damaris server
 *   rank_step=<integer>   MPI rank step
 *   trace                 trace output to <name>.log file
 *
 * parameters:
 *   name           <-- base output case name.
 *   options        <-- whitespace separated, lowercase options list
 *   time_dependecy <-- indicates if and how meshes will change with time
 *   comm           <-- associated MPI communicator.
 *
 * returns:
 *   pointer to opaque Damaris output writer structure.
 *----------------------------------------------------------------------------*/

#if defined(HAVE_MPI)
void *
fvm_to_damaris_init_writer(const char             *name,
                           const char             *path,
                           const char             *options,
                           fvm_writer_time_dep_t   time_dependency,
                           MPI_Comm                comm)
#else
void *
fvm_to_damaris_init_writer(const char             *name,
                           const char             *path,
                           const char             *options,
                           fvm_writer_time_dep_t   time_dependency)
#endif
{
  CS_UNUSED(time_dependency);
  CS_UNUSED(path);

  fvm_to_damaris_writer_t  *w = NULL;

  /* Parse options */
  bool usm_on = false;
  bool rect_grid_on = false;  // not used as too difficult to get resulting switch in cs_solver.c

  if (options != NULL) {

    int i1 = 0, i2 = 0;
    int l_tot = strlen(options);

    const char rs[] = "rank_step=";
    const int l_rs = strlen(rs);

    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      int l_opt = i2 - i1;

      if ((l_opt == 6) && (strncmp(options + i1, "usm_on", l_opt) == 0)) {
    	  usm_on = true;
      }
	  else if ((l_opt == 12) && (strncmp(options + i1, "rect_grid_on", l_opt) == 0)) {
		  rect_grid_on = true;
	  }

      for (i1 = i2 + 1; i1 < l_tot && options[i1] == ' '; i1++);

    }
  }

  /* Initialize writer */
  BFT_MALLOC(w, 1, fvm_to_damaris_writer_t);

  BFT_MALLOC(w->name, strlen(name) + 1, char);
  strcpy(w->name, name);

  w->usm_on = usm_on;
  w->rect_grid_on =rect_grid_on ;
  w->tracefile = NULL;
  w->damaris_mpi_comm  = comm ;
//  w->rank = 0;
//  w->n_ranks = 1;

//  w->buffer_size = 0;

  w->f_map = cs_map_name_to_id_create();
  w->f_ts = NULL;

  w->time_value = 0.0 ;
//  if (w->dry_run == false)
    _writer = w;

  /* Return writer */
  return w;
}

/*----------------------------------------------------------------------------
 * Finalize FVM to Damaris output writer.
 *
 * parameters:
 *   this_writer_p <-- pointer to opaque Damaris writer structure.
 *
 * returns:
 *   NULL pointer
 *----------------------------------------------------------------------------*/

void *
fvm_to_damaris_finalize_writer(void  *this_writer_p)
{
  fvm_to_damaris_writer_t *w = (fvm_to_damaris_writer_t *)this_writer_p;

  if (_writer == w)
    _writer = NULL;

  if (w->tracefile != NULL) {
    if (fclose(w->tracefile) != 0)
      bft_error(__FILE__, __LINE__, errno,
                _("Error closing file: \"%s.log\""), w->name);
  }
  BFT_FREE(w->name);

  cs_map_name_to_id_destroy(&(w->f_map));
  BFT_FREE(w->f_ts);

//  if (w->dry_run == false)
//    damaris_finalize();  // called in cs_base.c  _cs_base_exit()

  BFT_FREE(w);

  return NULL;
}

/*----------------------------------------------------------------------------
 * Associate new time step with a Damaris geometry.
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   time_step     <-- time step number
 *   time_value    <-- time_value number
 *----------------------------------------------------------------------------*/

void
fvm_to_damaris_set_mesh_time(void          *this_writer_p,
                             const int      time_step,
                             const double   time_value)
{
  CS_UNUSED(this_writer_p);
  CS_UNUSED(time_value);
  CS_UNUSED(time_step);
}



/*----------------------------------------------------------------------------
 * Define VTK geometric element type according to FVM element type
 * The definitions come from VTK header vtkCellType.h
 *
 * parameters:
 *   fvm_elt_type <-- pointer to fvm element type.
 *
 * return:
 *   Integer representation of VTK geometrical element type.
 *----------------------------------------------------------------------------*/

static int
_get_vtk_element_type(const fvm_element_t fvm_elt_type)
{
  int  norm_elt_type;

  switch (fvm_elt_type) {

  case FVM_EDGE:
    norm_elt_type = 3 ; // VTK_LINE;
    break;

  case FVM_FACE_TRIA:
    norm_elt_type = 5 ; // VTK_TRIANGLE;
    break;

  case FVM_FACE_QUAD:
    norm_elt_type = 9 ; // VTK_QUAD;
    break;

  case FVM_FACE_POLY:
    norm_elt_type = 7 ; // VTK_POLYGON;
    break;

  case FVM_CELL_TETRA:
    norm_elt_type = 10 ; // VTK_TETRA;
    break;

  case FVM_CELL_PYRAM:
    norm_elt_type = 14 ; // VTK_PYRAMID;
    break;

  case FVM_CELL_PRISM:
    norm_elt_type = 13 ; // VTK_WEDGE;
    break;

  case FVM_CELL_HEXA:
    norm_elt_type = 12 ; // VTK_HEXAHEDRON;
    break;

  case FVM_CELL_POLY:
    norm_elt_type = 42 ; // VTK_POLYHEDRON;
    break;

  default:
    norm_elt_type = 0 ; // VTK_EMPTY_CELL;
    bft_error(__FILE__, __LINE__, 0,
              "_get_vtk_element_type(): "
              "No association with VTK element type has been found\n"
              "FVM element type: \"%i\"\n",
              (int)fvm_elt_type);

  } /* End of switch on element type */

  return norm_elt_type;
}


/*----------------------------------------------------------------------------
 * Get vertex order to describe Catalyst element type.
 *
 * parameters:
 *   norm_elt_type  <-- Catalyst element type.
 *   vertex_order  --> Pointer to vertex order array (0 to n-1).
 *----------------------------------------------------------------------------*/

static void
_get_vertex_order(fvm_element_t norm_elt_type,
                  int          *vertex_order)
{
  switch(norm_elt_type) {

  case FVM_EDGE:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    break;

  case FVM_FACE_TRIA:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    break;

  case FVM_FACE_QUAD:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case FVM_CELL_TETRA:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case FVM_CELL_PYRAM:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    vertex_order[4] = 4;
    break;

  case FVM_CELL_PRISM:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    vertex_order[4] = 5;
    vertex_order[5] = 4;
    break;

  case FVM_CELL_HEXA:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    vertex_order[4] = 4;
    vertex_order[5] = 5;
    vertex_order[6] = 6;
    vertex_order[7] = 7;
    break;

  case FVM_FACE_POLY:
    vertex_order[0] = -1;
    break;

  case FVM_CELL_POLY:
    vertex_order[0] = -1;
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "_get_vertex_order(): No associated FVM element type known\n"
              "VTK element type: \"%i\"\n",
              (int)norm_elt_type);
  }

  return;
}

/*----------------------------------------------------------------------------
 * Write nodal mesh to a Damaris output
 *
 * parameters:
 *   this_writer_p <-- pointer to associated writer
 *   mesh          <-- pointer to nodal mesh structure that should be written
 *----------------------------------------------------------------------------*/

void
fvm_to_damaris_export_nodal(void               *this_writer_p,
                            const fvm_nodal_t  *mesh)
{
	int  section_id;

	fvm_to_damaris_writer_t  *w = (fvm_to_damaris_writer_t *)this_writer_p;

	/* Flag w->usm_on  is set in the Code_Saturne input xml file
	 * in the Damaris writer (format name="damaris")
	   <writer id="1" label="damaris_results">
        <directory name="postprocessing"/>
        <format name="damaris" options="usm_on"/>
        <frequency period="time_step">1</frequency>
        <time_dependency choice="fixed_mesh"/>
      </writer>
	 */
	if (w->usm_on == true)
	{
	const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);


	/* Damaris input xml file parameter values */
	/*----------------*/
	int  vtk_type_stride ;
	int  n_sections = 0 ;
	unsigned long n_connectivity = 0 ;
	unsigned long  n_elements_local = 0;

	/* Damaris parameter values:  n_sections, n_vertices, n_connectivity, n_elements
	 *
	 * Count the sections and total number of elements per section and
	 * size of connectivity array (n_connectivity)
	 *
	 * n_elements is the size of the field array, which I think may need to be
	 * set as the largest in the set of all MPI processes via MPIAllreduce()
	 */
	for (section_id = 0; section_id < mesh->n_sections; section_id++) {
		const fvm_nodal_section_t  *section = mesh->sections[section_id];
		if (section->entity_dim < elt_dim)
		  continue;
	   n_sections++ ;
	   n_elements_local += section->n_elements;
	   vtk_type_stride                = fvm_nodal_n_vertices_element[section->type];
	   // This will need to be changed if these elements are present in a typical mesh
	   if ((section->type != FVM_FACE_POLY) &&  (section->type != FVM_CELL_POLY))
	   {
		   n_connectivity += section->n_elements * vtk_type_stride ;

	   }
	} /* End of loop on sections number 1*/


	const unsigned long n_vertices = mesh->n_vertices;

	/*
	 * Allocate arrays to be passed to Damaris for per section details
	 * (sizes and mesh element types)
	 */
	unsigned int  * sect_sizes ;
	BFT_MALLOC(sect_sizes, n_sections , unsigned int);
	unsigned int  * sect_vtk_type ;
	BFT_MALLOC(sect_vtk_type, n_sections , unsigned int);

	// the section_id index can have null values, so use section_id_real as the index to sect_sizes and sect_vtk_type
	int section_id_real = 0;
	for (section_id = 0; section_id < mesh->n_sections; section_id++) {
		const fvm_nodal_section_t  *section = mesh->sections[section_id];
		if (section->entity_dim < elt_dim)
		  continue;

	   sect_sizes[section_id_real]    =  section->n_elements ;
	   // This could be a uchar vector as element types are stored as unsigned char array (from vtkCellType.h)
	   int vtk_type                   = _get_vtk_element_type(section->type);
	   sect_vtk_type[section_id_real] = vtk_type ;
	   section_id_real++;

	} /* End of loop on sections number 2 */

	/* Allocate local vertex storage array to be passed to Damaris */
	double  * vrtx_umesh_xyz ;
	BFT_MALLOC(vrtx_umesh_xyz, n_vertices * mesh->dim, double);
	/* Allocate local section connectivities array to be passed to Damaris */
	unsigned long * sectn_connectivity ;
	BFT_MALLOC(sectn_connectivity, n_connectivity, unsigned long );

	/* Vertex coordinates */
	/*--------------------*/
	const double  *vertex_coords = mesh->vertex_coords;

	// The 3rd dimension will be added by Damaris if needed
	if (mesh->parent_vertex_num != NULL) {
		const cs_lnum_t  *parent_vertex_num = mesh->parent_vertex_num;
		for (size_t i = 0; i < n_vertices; i++) {
			for (int j = 0; j < mesh->dim; j++)
				vrtx_umesh_xyz[i*mesh->dim + j] = vertex_coords[(parent_vertex_num[i]-1)*mesh->dim + j];
			}
	}
	else {
		/*for (size_t i = 0; i < n_vertices; i++) {
			for (size_t j = 0; j < mesh->dim; j++)
				vrtx_umesh_xyz[i*mesh->dim + j] = vertex_coords[i*mesh->dim + j];
			}*/
		memcpy((void*)vrtx_umesh_xyz,(void*)vertex_coords,n_vertices*mesh->dim*sizeof(double)) ;
	}


/*
    struct _fvm_io_num_t {

	  cs_gnum_t          global_count;    /* Global number of entities
	  cs_lnum_t          global_num_size; /* Local size of global numbering array
	  const cs_gnum_t   *global_num;      /* Global (possibly shared) entity
	                                         numbers (1 to n)
	  cs_gnum_t         *_global_num;     /* Global entity numbers if owner,
	                                         NULL otherwise
    };
*/
	const unsigned long  * vrtx_gid ;
	// GID's are the same data as required by damaris, so do not copy
	// Damaris will add/subtract the offset as specified by:   <vertex_global_id   offset="-1" />
	if (mesh->global_vertex_num != NULL) {
	//	mesh->global_vertex_num ;
		vrtx_gid = static_cast<const unsigned long *>(fvm_io_num_get_global_num(mesh->global_vertex_num));
	}

	/* Element connectivity */
	/*----------------------*/

	int vertex_order[8];
	// These are for conversion of FVM_CELL_PRISM to VTK_WEDGE
	vertex_order[0] = 0;
	vertex_order[1] = 2;
	vertex_order[2] = 1;
	vertex_order[3] = 3;
	vertex_order[4] = 5;
	vertex_order[5] = 4;

	// fvm_nodal_n_vertices_element(section->type)
	// sectn_connectivity
	unsigned long rolling_offset = 0;
	section_id_real = 0 ;
	for (section_id = 0; section_id < mesh->n_sections; section_id++) {

		const fvm_nodal_section_t  *section = mesh->sections[section_id];

		if (section->entity_dim < elt_dim)
		  continue;

		const int  stride = fvm_nodal_n_vertices_element[section->type];

		if (section->stride > 0){
			// 	_write_connect_block(section->type,
			// 		                           section->n_elements,
			// 		                           section->vertex_num,
			// 		                           ugrid);
			cs_lnum_t  i;
			int  j;


			// VTKCellType vtk_type = _get_vtk_element_type(section->type);

			if (section->type != FVM_CELL_PRISM ) {  // == VTK_WEDGE
				for (i = 0; i <  section->n_elements; i++) {
						for (j = 0; j < stride; j++)
							sectn_connectivity[rolling_offset+ i*stride+j] = section->vertex_num[i*stride+j];
				}
				//N.B. memcpy can't be used due to possible size difference between cs_lnum_t and long int
				//memcpy((void*)sectn_connectivity,(void*)section->vertex_num,section->n_elements*stride*sizeof(long int)) ;
			} else {
				 _get_vertex_order(section->type, vertex_order);
				for (i = 0; i <  section->n_elements; i++) {
					for (j = 0; j < stride; j++)
						sectn_connectivity[rolling_offset+i*stride+j] = section->vertex_num[i*stride + vertex_order[j]] ;
				}
			}
		} else {
			// I will not write these connectivities (yet)
			//else if (section->type == FVM_FACE_POLY)
			//  _export_nodal_polygons(section, ugrid);
			//else if (section->type == FVM_CELL_POLY)
			//  _export_nodal_polyhedra(mesh->n_vertices, section, ugrid);
			bft_error(__FILE__, __LINE__, 0,
			   "fvm_to_damaris_export_nodal(): VTK_POLYGON and VTK_POLYHEDRON have not been exported!!! \n"
			   "Needs to be fixed for this mesh\n") ;
		}

		rolling_offset += sect_sizes[section_id_real] * stride ;
		section_id_real++ ;
	} /* End of loop on sections */


   /* Communicate global values which are needed for calls to damaris_set_position() for the various data sets */
   /*---------------------------*/

   /*  n_sections_total == sum of n_sections over all ranks,
	*  n_sections_total is the size of: section_types and section_sizes arrays
	*  sections_offsets == array of size "mpi_n_ranks", each element holding prefix sum of elements
	*  up to its rank position.
	*    i.e. if ranks hold 3,2,4,2 then sections_offsets[] = {0,3,5,9}.
	*  The offsets index into section_types, section_sizes, section_dims and section_vertices arrays
	*
	*  QUESTION: Can MPI_Allreduce and MPI_Allgather be done in Damaris memory using damaris_read()?
	*/
	int n_sections_total ;
	MPI_Allreduce(&n_sections, &n_sections_total, 1, MPI_INT, MPI_SUM, w->damaris_mpi_comm) ;

	int * sectn_offsets ;
	BFT_MALLOC(sectn_offsets, cs_glob_n_ranks, int );
	int * sectn_sizes ;
	BFT_MALLOC(sectn_sizes, cs_glob_n_ranks, int  );
	// This gets the n_sections
	MPI_Allgather(&n_sections, 1, MPI_INT, sectn_sizes, 1, MPI_INT, w->damaris_mpi_comm);
	sectn_offsets[0] = 0ul ;
	for (int t1 = 1 ; t1 < cs_glob_n_ranks; t1++)
	{
		sectn_offsets[t1] = sectn_offsets[t1-1] + sectn_sizes[t1-1];
	}



  /*  The following are used for the vertex coordinate offsets and the
	*  vertex GIDs. Vertex coordinate offsets need to be multiplied by the mesh dimension (2D or 3D)
	*  n_vertices_total is the size of: unstructured_gid and unstructured_mesh_xyz (== n_vertices_totalx3)
	*/
	unsigned long  n_vertices_total ;
	MPI_Allreduce(&n_vertices, &n_vertices_total, 1, MPI_UNSIGNED_LONG, MPI_SUM, w->damaris_mpi_comm) ;

	unsigned long  * vertices_rank_offsets ;
	BFT_MALLOC(vertices_rank_offsets, cs_glob_n_ranks, unsigned long  );
	unsigned long  * vertices_rank_sizes ;
	BFT_MALLOC(vertices_rank_sizes, cs_glob_n_ranks, unsigned long );
	// This gets the n_sections
	MPI_Allgather(&n_vertices, 1, MPI_UNSIGNED_LONG, vertices_rank_sizes, 1, MPI_UNSIGNED_LONG, w->damaris_mpi_comm);
	vertices_rank_offsets[0] = 0ul ;
	for (int t1 = 1 ; t1 < cs_glob_n_ranks; t1++)
	{
		vertices_rank_offsets[t1] = vertices_rank_offsets[t1-1] + vertices_rank_sizes[t1-1];
	}


	/*  The following are used for the vertex connectivity data
	*   The connectivity offsets are the positions of the first element of the ranks subset of conectivities
	*   Each section of the mesh on a rank contributes n_elements * element stride values.
	*   n_connectivity_total is the total size required for the variable to hold the connectivities
	*/
	unsigned long  n_connectivity_total ;
	MPI_Allreduce(&n_connectivity, &n_connectivity_total, 1, MPI_UNSIGNED_LONG, MPI_SUM, w->damaris_mpi_comm) ;

	unsigned long * connectivity_rank_offsets ;
	BFT_MALLOC(connectivity_rank_offsets, cs_glob_n_ranks, unsigned long );
	unsigned long * connectivity_rank_sizes ;
	BFT_MALLOC(connectivity_rank_sizes, cs_glob_n_ranks, unsigned long );
	// This gets the n_connectivity of each rank, which is the size of the connectivity of all sections on the current rank
	MPI_Allgather(&n_connectivity, 1, MPI_UNSIGNED_LONG, connectivity_rank_sizes, 1, MPI_UNSIGNED_LONG, w->damaris_mpi_comm);
	connectivity_rank_offsets[0] = 0ul ;
	for (int t1 = 1 ; t1 < cs_glob_n_ranks; t1++)
	{
		connectivity_rank_offsets[t1] = connectivity_rank_offsets[t1-1] + connectivity_rank_sizes[t1-1];
	}


	// n_elements_local is the total number of local elements in local mesh sections. It is equal to the
	// number of field values for a zonal field on the current client.
	unsigned long  n_elements_total ;
	MPI_Allreduce(&n_elements_local, &n_elements_total, 1, MPI_UNSIGNED_LONG, MPI_SUM, w->damaris_mpi_comm) ;
	// Need to check for overflow of integer type as Damaris paramaters are not long int capable.
	if (n_elements_local > INT_MAX) {
		bft_error(__FILE__, __LINE__, -1, _("ERROR: OVERFLOW Detected. Parameter: n_elements_local \n \
		                                     Damaris Parameters should be in range of C type int"));
	}
	if (n_elements_total > INT_MAX) {
		bft_error(__FILE__, __LINE__, -1, _("ERROR: OVERFLOW Detected. Parameter: n_elements_total \n \
                                         Damaris Parameters should be in range of C type int"));
	}

	unsigned long * elements_rank_offsets ;  // to be computed - Probably could use MPI_Scan()?
	BFT_MALLOC(elements_rank_offsets, cs_glob_n_ranks, unsigned long );
	unsigned long * elements_rank_sizes ;    // to be gathered from each client rank
	BFT_MALLOC(elements_rank_sizes, cs_glob_n_ranks, unsigned long );
  // This gets the vector of elements_rank_sizes of each rank, which is the number of elements of all sections on the current rank
	MPI_Allgather(&n_elements_local, 1, MPI_UNSIGNED_LONG, elements_rank_sizes, 1, MPI_UNSIGNED_LONG, w->damaris_mpi_comm);
	elements_rank_offsets[0] = 0ul ;
  // This is a sum scan of the element sizes to be used as offsets if saving the full element data
	for (int t1 = 1 ; t1 < cs_glob_n_ranks; t1++)
	{
		elements_rank_offsets[t1] = elements_rank_offsets[t1-1] + elements_rank_sizes[t1-1];
	}

	/* Pass everything to Damaris */
	/*----------------------------*/

	int damaris_err = DAMARIS_OK;

	// The offset is used so that output functionality (e.g. HDF5Store) knows global offsets of the data of the rank
	// So, it is not strictly needed, however is here for completeness and is used above, so do not remove.
	int temp_int = elements_rank_offsets[cs_glob_rank_id] ;
	damaris_err = damaris_parameter_set("n_elements_offset",&temp_int, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_elements_offset"));
	}
	temp_int = n_elements_local;
	damaris_err = damaris_parameter_set("n_elements_local",&temp_int, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_elements_local"));
	}
	temp_int = n_elements_total;
	damaris_err = damaris_parameter_set("n_elements_total",&temp_int, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_elements_total"));
	}

	// mesh_dim and n_vertices_local affects the size of the vertices array (unstructured_mesh_xyz) on each rank
	int dim = mesh->dim ;
	damaris_err = damaris_parameter_set("mesh_dim",&dim, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: mesh_dim"));
	}

	// N.B. damaris_parameter_set("cs_glob_n_ranks"...) has been set in cs_base.c in function  _cs_base_mpi_setup()
	/* Global size parameters */
	// n_sections_total, n_vertices_total, n_connectivity_total
	damaris_err = damaris_parameter_set("n_sections_total",&n_sections_total, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_sections_total"));
	}

	damaris_err = damaris_parameter_set("n_vertices_total",&n_vertices_total, sizeof(unsigned long));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_vertices_total"));
	}

	damaris_err = damaris_parameter_set("n_connectivity_total",&n_connectivity_total, sizeof(unsigned long));
	if (damaris_err != DAMARIS_OK ) {
		bft_error(__FILE__, __LINE__, damaris_err,_("ERROR: Damaris damaris_parameter_set():\n  parameter: n_connectivity_total"));
	}

	/* Local (per process) size parameters */
	// n_sections, n_vertices, n_connectivity, they can be different for each MPI process
	damaris_err = damaris_parameter_set("n_sections_local",&n_sections, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_sections_local"));
	}


	damaris_err = damaris_parameter_set("n_vertices_local",&n_vertices, sizeof(unsigned long));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\n  parameter: n_vertices_local"));
	}

	// n_connectivity_local is the total size of all connectivity indices for the sections on a rank
	damaris_err = damaris_parameter_set("n_connectivity_local",&n_connectivity, sizeof(unsigned long));
	if (damaris_err != DAMARIS_OK ) {
		bft_error(__FILE__, __LINE__, damaris_err,_("ERROR: Damaris damaris_parameter_set():\n  parameter: n_connectivity_local"));
	}


	/* Field sizes (per process and global) size parameters */
	/* The field sizes are currently 'hard coded' in the xml file
	int global_size = 262144 ;
	damaris_err = damaris_parameter_set("n_field_array_total",&global_size, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {n_field_array_total}
	int local_size = 131072 ;
	damaris_err = damaris_parameter_set("n_field_array_local",&local_size, sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_parameter_set():\nparamater: n_field_array_local"));
	}
	*/

	/*
	 * <mesh name="fluid_domain_umesh" type="unstructured" topology="3">
	 *
		<coord              name="umesh_vars/unstructured_mesh_xyz" unit="m" />
		<vertex_global_id   name="umesh_vars/unstructured_gid"  offset="-1" />
		<section_types      name="umesh_vars/section_types"  />
		<section_sizes      name="umesh_vars/section_sizes"  />
		<connectivity       name="umesh_vars/section_connectivity"  />
	  </mesh>
	 */

	// damaris_set_position
	int64_t pos;

	pos = vertices_rank_offsets[cs_glob_rank_id] * mesh->dim ;
	damaris_err = damaris_set_position("umesh_vars/unstructured_mesh_xyz" , &pos);
	if (damaris_err != DAMARIS_OK )
		bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_set_position() umesh_vars/unstructured_mesh_xyz"));
	// unstructured_mesh_xyz contains the Mesh coords variable
	damaris_err = damaris_write("umesh_vars/unstructured_mesh_xyz" , vrtx_umesh_xyz);
	if (damaris_err != DAMARIS_OK ) {
		bft_error(__FILE__, __LINE__, damaris_err,_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/unstructured_mesh_xyz" ));
	}

	pos = vertices_rank_offsets[cs_glob_rank_id]  ;
	damaris_err = damaris_set_position("umesh_vars/unstructured_gid" , &pos);
	if (damaris_err != DAMARIS_OK )
		bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_set_position() umesh_vars/unstructured_gid"));
	damaris_err = damaris_write("umesh_vars/unstructured_gid" , vrtx_gid);
	if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/unstructured_gid" ));
		}



	pos = sectn_offsets[cs_glob_rank_id] ;
	damaris_err = damaris_set_position("umesh_vars/section_types" , &pos);
	if (damaris_err != DAMARIS_OK )
		bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_set_position() umesh_vars/section_types"));
	damaris_err = damaris_write("umesh_vars/section_types" , sect_vtk_type);
	if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/section_types" ));
		}

	// The section_sizes are equivalent to the number of elements per section, and thus equals the number of
	// field values for a particular mesh section (for fields that are 'per element' AKA 'zonal'
	damaris_err = damaris_set_position("umesh_vars/section_sizes" , &pos);
	if (damaris_err != DAMARIS_OK )
			bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_set_position() umesh_vars/section_sizes"));
	damaris_err = damaris_write("umesh_vars/section_sizes" , sect_sizes);
	if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/section_sizes" ));
		}


	pos = connectivity_rank_offsets[cs_glob_rank_id] ;
	damaris_err = damaris_set_position("umesh_vars/section_connectivity" , &pos);
		if (damaris_err != DAMARIS_OK )
			bft_error(__FILE__, __LINE__, damaris_err, _("ERROR: Damaris damaris_set_position() umesh_vars/section_connectivity"));
	damaris_err = damaris_write("umesh_vars/section_connectivity" , sectn_connectivity);
	if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/section_connectivity" ));
		}



	BFT_FREE(vrtx_umesh_xyz);
	BFT_FREE(sect_sizes);
	BFT_FREE(sect_vtk_type);
	BFT_FREE(sectn_connectivity);


	BFT_FREE(sectn_offsets);
	BFT_FREE(sectn_sizes);
	BFT_FREE(vertices_rank_offsets);
	BFT_FREE(vertices_rank_sizes);
	BFT_FREE(connectivity_rank_offsets);
	BFT_FREE(connectivity_rank_sizes);

	BFT_FREE(elements_rank_offsets);
	BFT_FREE(elements_rank_sizes);



	w->modified = true;
	} else {
		 printf("DAMARIS INFO: fvm_to_damaris_export_nodal() Damaris has not initialized the VTK Unstructured Mesh data. To enable set <writer><format name=\"damaris\" options=\"usm_on\"/> ");
	}
}


/*----------------------------------------------------------------------------
 * Write field associated with a nodal mesh to a Damaris output.
 *
 * Assigning a negative value to the time step indicates a time-independent
 * field (in which case the time_value argument is unused).
 *
 * parameters:
 *   this_writer_p    <-- pointer to associated writer
 *   mesh             <-- pointer to associated nodal mesh structure
 *   name             <-- variable name
 *   location         <-- variable definition location (nodes or elements)
 *   dimension        <-- variable dimension (0: constant, 1: scalar,
 *                        3: vector, 6: sym. tensor, 9: asym. tensor)
 *   interlace        <-- indicates if variable in memory is interlaced
 *   n_parent_lists   <-- indicates if variable values are to be obtained
 *                        directly through the local entity index (when 0) or
 *                        through the parent entity numbers (when 1 or more)
 *   parent_num_shift <-- parent number to value array index shifts;
 *                        size: n_parent_lists
 *   datatype         <-- indicates the data type of (source) field values
 *   time_step        <-- number of the current time step
 *   time_value       <-- associated time value
 *   field_values     <-- array of associated field value arrays
 *----------------------------------------------------------------------------*/

void
fvm_to_damaris_export_field(void                  *this_writer_p,
                             const fvm_nodal_t     *mesh,
                             const char            *name,
                             fvm_writer_var_loc_t   location,
                             int                    dimension,
                             cs_interlace_t         interlace,
                             int                    n_parent_lists,
                             const cs_lnum_t        parent_num_shift[],
                             cs_datatype_t          datatype,
                             int                    time_step,
                             double                 time_value,
                             const void      *const field_values[])
{

  char _name[128];

  fvm_to_damaris_writer_t *w = (fvm_to_damaris_writer_t *)this_writer_p;

  /* Initialization */
  /*----------------*/

  strncpy(_name, name, 127);
  _name[127] = '\0';
  if (w->ensight_names) {
    for (int i = 0; i < 127 && _name[i] != '\0'; i++) {
      switch (_name[i]) {
      case '(':
      case ')':
      case ']':
      case '[':
      case '+':
      case '-':
      case '@':
      case ' ':
      case '\t':
      case '!':
      case '#':
      case '*':
      case '^':
      case '$':
      case '/':
        _name[i] = '_';
        break;
      default:
        break;
      }
      if (_name[i] == ' ')
        _name[i] = '_';
    }
  }

  /*
  if (mesh_id < 0) {
    mesh_id = _add_catalyst_mesh(w, mesh);
    fvm_to_catalyst_export_nodal(w, mesh);
  }
  */

  int _time_step = (time_step > -1) ? time_step : 0;
  double _time_value = (time_value > 0.0) ? time_value : 0.0;
  if (_time_step > w->time_step) {
    w->time_step = _time_step;
    assert(time_value >= w->time_value);
    w->time_value = _time_value;
  }

  /* Per node variable */
  /*-------------------*/
  if (location == FVM_WRITER_PER_NODE)
    /*_export_field_values_n(mesh,
                           _name,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values,
                           f);*/
	  printf("INFO: Damaris currently is not set up to output per node variables from Code_Saturne");


  /* Per element variable */
  /*----------------------*/
  else if (location == FVM_WRITER_PER_ELEMENT)
    _export_field_values_e(mesh,
                           _name,
                           dimension,
                           interlace,
                           n_parent_lists,
                           parent_num_shift,
                           datatype,
                           field_values);

  /* Update field status */
  /*---------------------*/
/*
  fvm_to_catalyst_set_mesh_time(w, time_step, time_value);
*/
  w->modified = true;
}


/*----------------------------------------------------------------------------*/

END_C_DECLS
