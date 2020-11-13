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

/*----------------------------------------------------------------------------
 * Asynchronous I/O and in-situ visulisation library header
 *----------------------------------------------------------------------------*/

#if defined(HAVE_DAMARIS)
#include <Damaris.h>
#endif


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

  char        *name;               /* Writer name */

//  bool         dry_run;            /* If true, do not connect to Damaris */
  FILE        *tracefile;          /* optional file for tracing */

//  int          rank;               /* Rank of current process in communicator */
//  int          n_ranks;            /* Number of processes in communicator */

//  size_t       buffer_size;        /* buffer size required */

  int          time_step;          /* Latest time step */
  double       time_value;         /* Latest time value */

  cs_map_name_to_id_t  *f_map;     /* field names mapping */
  int         *f_ts;               /* last field output time step */

  bool        ensight_names;   /* Use EnSight rules for
                                                  field names */

#if defined(HAVE_MPI)
//  int          min_rank_step;      /* Minimum rank step */
//  int          min_block_size;     /* Minimum block buffer size */
//  MPI_Comm     block_comm;         /* Associated MPI block communicator */
//  MPI_Comm     comm;               /* Associated MPI communicator */
#endif

  bool          modified;          /* Has output been added since
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
  /*assert(f != NULL);*/

  int  section_id;
/*  double  *values = NULL; */

  const int dest_dim = (dim == 6) ? 9 : dim;

 /* values = vtkDoubleArray::SafeDownCast
    (f->GetCellData()->GetArray(fieldname))
    ->WritePointer(0, dest_dim*f->GetNumberOfCells());

  assert(values != NULL);
  */

  /* Distribute partition to block values */

  cs_lnum_t start_id = 0;
  cs_lnum_t src_shift = 0;

  /* loop on sections which should be appended */

  const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

  for (section_id = 0; section_id < mesh->n_sections; section_id++) {

    const fvm_nodal_section_t  *section = mesh->sections[section_id];

    if (section->entity_dim < elt_dim)
      continue;

   /* fvm_convert_array(dim,
                      0,
                      dest_dim,
                      src_shift,
                      section->n_elements + src_shift,
                      interlace,
                      datatype,
                      CS_DOUBLE,
                      n_parent_lists,
                      parent_num_shift,
                      section->parent_element_num,
                      field_values,
                      values + start_id);
*/
    //if (dim == 1)
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

    if ( strncmp( fieldname_lwrcase, "pressure", 8) == 0 )
    {
		int64_t pos[3];

		pos[0] = cs_glob_rank_id*param_z/cs_glob_n_ranks;
		pos[1] = 0;
		pos[2] = 0;

		damaris_err = damaris_set_position("fields/pressure" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/pressure");
		}
		damaris_err = damaris_write("fields/pressure" ,field_values[0]);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/pressure");
		}

    }
    //else if (dim == 3)
    else if ( strncmp( fieldname_lwrcase, "velocity", 8) == 0 )
    {
    	int64_t pos[4];

		pos[0] = 0;
		pos[1] = (cs_glob_rank_id*param_z/cs_glob_n_ranks);
		pos[2] = 0;
		pos[3] = 0;

		damaris_err = damaris_set_position("fields/velocity" , pos);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_set_position():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		damaris_err = damaris_write("fields/velocity" ,field_values[0]);
		if (damaris_err != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "field: \"%s\"."), "fields/velocity");
		}
		// damaris_end_iteration();
    }
    //else if ( strncmp( fieldname_lwrcase, "mpi_rank_id", 8) == 0 )

    start_id += section->n_elements*dest_dim;
    if (n_parent_lists == 0)
      src_shift += section->n_elements;


  }

  /* Special case for symmetric tensors */
/*
  if (dim == 6) {

    cs_lnum_t n_elts = f->GetNumberOfCells();
    for (cs_lnum_t i = 0; i < n_elts; i++) {
      values[9*i + 8] = values[9*i + 2];
      values[9*i + 7] = values[9*i + 4];
      values[9*i + 6] = values[9*i + 5];
      values[9*i + 4] = values[9*i + 1];
      values[9*i + 2] = values[9*i + 5];
      values[9*i + 1] = values[9*i + 3];
      values[9*i + 5] = values[9*i + 7];
    }
  }
  */
}

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

  fvm_to_damaris_writer_t  *w = NULL;

  /* Parse options */

  int rank_step = 1;
  bool trace = false;
  bool dry_run = false;



	/* Initializing the mesh coordinates for a rectilinear mesh
	*   N.B. This should be done after reading the input mesh so
	*   that the actual mesh dimensions can be used and set with
	*   damaris_parameter_set(). Our problem is that the writers
	*   are initialized before the mesh is read in, so we cannot
	*   use the mesh information here!
	*/
  /*
	int damaris_err = DAMARIS_OK ;
	double x_length ;
	int segments_x, segments_y, segments_z, segments_z_per_rank ;

	damaris_err = damaris_parameter_set("cs_glob_n_ranks",&cs_glob_n_ranks,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
									 _("ERROR: Damaris damaris_parameter_set():\n"
									   "paramater: \"%s\"."), "cs_glob_n_ranks");
	}

	damaris_err = damaris_parameter_get("x",&segments_x,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
							 _("ERROR: Damaris damaris_parameter_get():\n"
							   "Parameter: \"%s\"."), "x");
	}
	damaris_err = damaris_parameter_get("y",&segments_y,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
							 _("ERROR: Damaris damaris_parameter_get():\n"
							   "Parameter: \"%s\"."), "y");
	}
	damaris_err = damaris_parameter_get("z",&segments_z,sizeof(int));
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_parameter_get():\n"
								   "Parameter: \"%s\"."), "z");
	}
	damaris_err = damaris_parameter_get("x_length",&x_length,sizeof(double));
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_parameter_get():\n"
								   "Parameter: \"%s\"."), "x_length");
	}

	double* meshx;
	double* meshy;
	double* meshz;

	/* the z direction is distributed over mpi ranks /
	segments_z_per_rank = segments_z/cs_glob_n_ranks;

	BFT_MALLOC(meshx, segments_x, double);
	BFT_MALLOC(meshy, segments_y, double);
	BFT_MALLOC(meshz, segments_z_per_rank, double);

	/* these dimensions are governed by the mesh creation script mesh_cube_xyz.py /
	double x_step, y_step, z_step ;
	double y_length, z_length ;
	y_length = x_length / 2.0;
	z_length = (segments_z / segments_x) * x_length;

	x_step = x_length / segments_x;
	y_step = y_length / segments_y;
	z_step = z_length / segments_z;

	for(int i=0; i<segments_x+1 ; i++)
		meshx[i] = i*x_step;

	for(int j=0; j<segments_y+1 ; j++)
		meshy[j] = j*y_step;

	double offset = cs_glob_rank_id * segments_z_per_rank * z_step;
	for(int k=0; k<segments_z_per_rank+1 ; k++)
		meshz[k] = offset + (k * z_step);

	damaris_err = damaris_write("coord/meshx" , meshx);
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "Variable: \"%s\"."), "coord/meshx");
	}
	damaris_err = damaris_write("coord/meshy" , meshy);
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "Variable: \"%s\"."), "coord/meshy");
	}
	damaris_err = damaris_write("coord/meshz" , meshz);
	if (damaris_err != DAMARIS_OK ) {
	  bft_error(__FILE__, __LINE__, damaris_err,
								 _("ERROR: Damaris damaris_write():\n"
								   "Variable: \"%s\"."), "coord/meshz");
	}

	BFT_FREE(meshx); ;
	BFT_FREE(meshy);
	BFT_FREE(meshz);
*/


  if (options != NULL) {

    int i1 = 0, i2 = 0;
    int l_tot = strlen(options);

    const char rs[] = "rank_step=";
    const int l_rs = strlen(rs);

    while (i1 < l_tot) {

      for (i2 = i1; i2 < l_tot && options[i2] != ' '; i2++);
      int l_opt = i2 - i1;

      if ((l_opt == 5) && (strncmp(options + i1, "trace", l_opt) == 0))
        trace = true;

      else if ((l_opt == 7) && (strncmp(options + i1, "dry_run", l_opt) == 0)) {
        dry_run = true;
        trace = true;
      }

      else if ((strncmp(options + i1, rs, l_rs) == 0)) {
        if (l_opt < l_rs+32) { /* 32 integers more than enough
                                  for maximum integer string */
          char options_c[32];
          strncpy(options_c, options+i1+l_rs, l_opt-l_rs);
          options_c[l_opt-l_rs] = '\0';
          rank_step = atoi(options_c);

        }
      }

      for (i1 = i2 + 1; i1 < l_tot && options[i1] == ' '; i1++);

    }
  }

  if (dry_run == false && _writer != NULL) {
    bft_error(__FILE__, __LINE__, errno,
              _("Error creating Damaris writer: \"%s\":\n"
                "only one Damaris server may be used, and is already used by\n"
                "writer: \"%s\"."), name, _writer->name);
    return NULL;
  }

  /* Initialize writer */

  BFT_MALLOC(w, 1, fvm_to_damaris_writer_t);

  BFT_MALLOC(w->name, strlen(name) + 1, char);
  strcpy(w->name, name);

//  w->dry_run = dry_run;
  w->tracefile = NULL;

//  w->rank = 0;
//  w->n_ranks = 1;

//  w->buffer_size = 0;

  w->f_map = cs_map_name_to_id_create();
  w->f_ts = NULL;


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
 * Define VTK geometrical element type according to FVM element type
 *
 * parameters:
 *   fvm_elt_type <-- pointer to fvm element type.
 *
 * return:
 *   med geometrical element type.
 *----------------------------------------------------------------------------*/

static VTKCellType
_get_norm_elt_type(const fvm_element_t fvm_elt_type)
{
  VTKCellType  norm_elt_type;

  switch (fvm_elt_type) {

  case FVM_EDGE:
    norm_elt_type = VTK_LINE;
    break;

  case FVM_FACE_TRIA:
    norm_elt_type = VTK_TRIANGLE;
    break;

  case FVM_FACE_QUAD:
    norm_elt_type = VTK_QUAD;
    break;

  case FVM_FACE_POLY:
    norm_elt_type = VTK_POLYGON;
    break;

  case FVM_CELL_TETRA:
    norm_elt_type = VTK_TETRA;
    break;

  case FVM_CELL_PYRAM:
    norm_elt_type = VTK_PYRAMID;
    break;

  case FVM_CELL_PRISM:
    norm_elt_type = VTK_WEDGE;
    break;

  case FVM_CELL_HEXA:
    norm_elt_type = VTK_HEXAHEDRON;
    break;

  case FVM_CELL_POLY:
    norm_elt_type = VTK_POLYHEDRON;
    break;

  default:
    norm_elt_type = VTK_EMPTY_CELL;
    bft_error(__FILE__, __LINE__, 0,
              "_get_norm_elt_type(): "
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
_get_vertex_order(VTKCellType   norm_elt_type,
                  int          *vertex_order)
{
  switch(norm_elt_type) {

  case VTK_LINE:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    break;

  case VTK_TRIANGLE:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    break;

  case VTK_QUAD:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case VTK_TETRA:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    break;

  case VTK_PYRAMID:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    vertex_order[4] = 4;
    break;

  case VTK_WEDGE:
    vertex_order[0] = 0;
    vertex_order[1] = 2;
    vertex_order[2] = 1;
    vertex_order[3] = 3;
    vertex_order[4] = 5;
    vertex_order[5] = 4;
    break;

  case VTK_HEXAHEDRON:
    vertex_order[0] = 0;
    vertex_order[1] = 1;
    vertex_order[2] = 2;
    vertex_order[3] = 3;
    vertex_order[4] = 4;
    vertex_order[5] = 5;
    vertex_order[6] = 6;
    vertex_order[7] = 7;
    break;

  case VTK_POLYGON:
    vertex_order[0] = -1;
    break;

  case VTK_POLYHEDRON:
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
	int  mesh_id, section_id;

	fvm_to_catalyst_t  *w = (fvm_to_catalyst_t *)this_writer_p;

	const int  elt_dim = fvm_nodal_get_max_entity_dim(mesh);

	/* Damaris parameter values */
	/*----------------*/
	int  vtk_type_stride ;
	int  n_sections = 0 ;
	unsigned long n_connectivity = 0 ;
	cs_lnum_t  n_elts = 0;

	/* Damaris parameter values:  n_sections, n_vertices, n_connectivity
	 *
	 * Count the sections and total number of elements per section and
	 * size of connectivity array (n_connectivity)
	 */
	for (section_id = 0; section_id < mesh->n_sections; section_id++) {
		const fvm_nodal_section_t  *section = mesh->sections[section_id];
		if (section->entity_dim < elt_dim)
		  continue;
	   n_sections++ ;
	   n_elts += section->n_elements;
	   vtk_type_stride                = fvm_nodal_n_vertices_element[section->type];
	   // This will need to be changed if these elements are present in a typical mesh
	   if ((section->type != FVM_FACE_POLY) &&  (section->type != FVM_CELL_POLY))
		   n_connectivity += section->n_elements * vtk_type_stride ;
	} /* End of loop on sections number 1*/

	const int n_vertices = mesh->n_vertices;

	/*
	 * Allocate arrays to be passed to Damaris for per section details
	 * (sizes and mesh element types)
	 */
	unsigned int  * sect_sizes ;
	BFT_MALLOC(sect_sizes, n_sections , unsigned long);
	unsigned int  * sect_vtk_type ;
	BFT_MALLOC(sect_vtk_type, n_sections , unsigned int);

	for (section_id = 0; section_id < mesh->n_sections; section_id++) {
		const fvm_nodal_section_t  *section = mesh->sections[section_id];
		if (section->entity_dim < elt_dim)
		  continue;

	   sect_sizes[section_id]         =  section->n_elements ;
	   VTKCellType vtk_type           = _get_norm_elt_type(section->type);
	   sect_vtk_type[section_id]      = vtk_type ;

	} /* End of loop on sections number 2 */

	/* Allocate vertex storage array to be passed to Damaris */
	double  * vrtx_umesh_xyz ;
	BFT_MALLOC(vrtx_umesh_xyz, n_vertices * mesh->dim, double);
	/* Allocate section connectivities array to be passed to Damaris */
	unsigned long * sectn_connectivity ;
	BFT_MALLOC(sectn_connectivity, n_connectivity, unsigned long );

	/* Vertex coordinates */
	/*--------------------*/
	const double  *vertex_coords = mesh->vertex_coords;

	// The 3rd dimension is added by Damaris
	if (mesh->parent_vertex_num != NULL) {
		const cs_lnum_t  *parent_vertex_num = mesh->parent_vertex_num;
		for (i = 0; i < n_vertices; i++) {
			for (j = 0; j < mesh->dim; j++)
				*vrtx_umesh_xyz++ = vertex_coords[(parent_vertex_num[i]-1)*mesh->dim + j];
			}
	}
	else {
		for (i = 0; i < n_vertices; i++) {
			for (j = 0; j < mesh->dim; j++)
				*vrtx_umesh_xyz++ = vertex_coords[i*mesh->dim + j];
			}
	}

	unsigned long  * vrtx_gid ;
	// GID's are the same data in the same order, so do not copy
	// Damaris will add/subtract the offset
	if (mesh->global_vertex_num != NULL) {
		mesh->global_vertex_num ;
		vrtx_gid = (unsigned long *) mesh->global_vertex_num->global_num ;
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

	// sectn_connectivity
	for (section_id = 0; section_id < mesh->n_sections; section_id++) {

		const fvm_nodal_section_t  *section = mesh->sections[section_id];

		if (section->entity_dim < elt_dim)
		  continue;

		if (section->stride > 0){
			// 	_write_connect_block(section->type,
			// 		                           section->n_elements,
			// 		                           section->vertex_num,
			// 		                           ugrid);
			cs_lnum_t  i;
			int  j;

			const int  stride = fvm_nodal_n_vertices_element[section->type];
			// VTKCellType vtk_type = _get_norm_elt_type(section->type);

			if (vtk_type != FVM_CELL_PRISM ) {  // == VTK_WEDGE
				for (i = 0; i <  section->n_elements; i++) {
						for (j = 0; j < stride; j++, sectn_connectivity++)
							sectn_connectivity[j] = section->vertex_num[i*stride];
				}
			} else {
				// _get_vertex_order(vtk_type, vertex_order);
				for (i = 0; i <  section->n_elements; i++) {
					for (j = 0; j < stride; j++, sectn_connectivity++)
						sectn_connectivity[j] = section->vertex_num[i*stride + vertex_order[j]] - 1;
				}
			}
		} else {
			// I will not write these connectivites (yet)
			//else if (section->type == FVM_FACE_POLY)
			//  _export_nodal_polygons(section, ugrid);
			//else if (section->type == FVM_CELL_POLY)
			//  _export_nodal_polyhedra(mesh->n_vertices, section, ugrid);
			bft_error(__FILE__, __LINE__, 0,
			   "fvm_to_damaris_export_nodal(): VTK_POLYGON and VTK_POLYHEDRON have not been exported!!! \n"
			   "Needs to be fixed for this mesh\n") ;
		}
	} /* End of loop on sections */



	/* Pass everything to Damaris */
	/*----------------------------*/

	// n_sections, n_vertices, n_connectivity
	if (damaris_parameter_set("n_vertices",&n_vertices, sizeof(cs_lnum_t)) != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err,
		 _("ERROR: Damaris damaris_parameter_set():\nparamater: \"%s\"."), "n_vertices");
	}
	if (damaris_parameter_set("n_sections",&n_sections, sizeof(int)) != DAMARIS_OK ) {
		  bft_error(__FILE__, __LINE__, damaris_err,
		 _("ERROR: Damaris damaris_parameter_set():\nparamater: \"%s\"."), "n_sections");
	}
	if (damaris_parameter_set("n_connectivity",&n_connectivity, sizeof(unsigned long)) != DAMARIS_OK ) {
			  bft_error(__FILE__, __LINE__, damaris_err,
			 _("ERROR: Damaris damaris_parameter_set():\nparamater: \"%s\"."), "n_connectivity");
	}

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

	if (damaris_write("umesh_vars/unstructured_mesh_xyz" , vrtx_umesh_xyz) != DAMARIS_OK ) {
		bft_error(__FILE__, __LINE__, damaris_err,
			_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/unstructured_mesh_xyz" ));
	}
	if (damaris_write("umesh_vars/unstructured_gid" , vrtx_gid) != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/unstructured_gid" ));
		}
	if (damaris_write("umesh_vars/section_types" , sect_vtk_type) != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/section_types" ));
		}
	if (damaris_write("umesh_vars/section_sizes" , sect_sizes) != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/section_sizes" ));
		}
	if (damaris_write("umesh_vars/section_connectivity" , sectn_connectivity) != DAMARIS_OK ) {
			bft_error(__FILE__, __LINE__, damaris_err,
				_("ERROR: Damaris damaris_write():\nVariable: umesh_vars/section_connectivity" ));
		}



	BFT_FREE(vrtx_umesh_xyz);
	BFT_FREE(sect_sizes);
	BFT_FREE(sect_vtk_type);
	BFT_FREE(sectn_connectivity);


	w->modified = true;
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
  int  mesh_id, field_id;
  char _name[128];

  fvm_to_damaris_writer_t *w = (fvm_to_damaris_writer_t *)this_writer_p;

  /* Initialization */
  /*----------------*/
/*
  mesh_id = _get_catalyst_mesh_id(w, mesh->name);
*/
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
 /*   w->datadesc->SetTimeData(w->time_value, w->time_step);*/
  }

  /* Get field id */
  /*
  field_id = _get_catalyst_field_id(w,
                                    _name,
                                    mesh_id,
                                    dimension,
                                    datatype,
                                    location);

  if (field_id < 0)
    field_id = _add_catalyst_field(w,
                                   _name,
                                   mesh_id,
                                   dimension,
                                   datatype,
                                   location);
*/
  /*vtkUnstructuredGrid  *f = w->fields[field_id]->f;*/

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
