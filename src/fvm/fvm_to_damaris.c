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

/*----------------------------------------------------------------------------
 * Statistics library header
 *----------------------------------------------------------------------------*/


#include "Damaris.h"

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

  bool         dry_run;            /* If true, do not connect to Damaris */
  FILE        *tracefile;          /* optional file for tracing */

  int          rank;               /* Rank of current process in communicator */
  int          n_ranks;            /* Number of processes in communicator */

  size_t       buffer_size;        /* buffer size required */

  int          time_step;          /* Latest time step */
  double       time_value;         /* Latest time value */

  cs_map_name_to_id_t  *f_map;     /* field names mapping */
  int         *f_ts;               /* last field output time step */

  bool        ensight_names;   /* Use EnSight rules for
                                                  field names */

#if defined(HAVE_MPI)
  int          min_rank_step;      /* Minimum rank step */
  int          min_block_size;     /* Minimum block buffer size */
  MPI_Comm     block_comm;         /* Associated MPI block communicator */
  MPI_Comm     comm;               /* Associated MPI communicator */
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
  bool                      call_init;   /* call damaris_init ? */

} _damaris_context_t;

/*============================================================================
 * Static global variables
 *============================================================================*/

/* Since the Damaris API specifies the field name (and communicator)
   only, only one Damaris writer can be active a a given time. */

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
    if ( strncmp( fieldname, "pressure", 8) == 0 )
    {
		int64_t pos[3];
		int our_z ;
		damaris_parameter_get("z",&our_z,sizeof(int));
		pos[0] = cs_glob_rank_id*our_z;
		pos[1] = 0;
		pos[2] = 0;

		damaris_set_position("pressure" , pos);
		damaris_write("pressure" ,field_values[0]);

    }
    //else if (dim == 3)
    else if ( strncmp( fieldname, "velocity", 8) == 0 )
    {
    	int64_t pos[4];
		int our_z ;
		damaris_parameter_get("z",&our_z,sizeof(int));
		pos[0] = cs_glob_rank_id*our_z;
		pos[1] = 0;
		pos[2] = 0;
		pos[3] = 0;

		damaris_set_position("velocity" , pos);
		damaris_write("velocity" ,field_values[0]);
		// damaris_end_iteration();
    }

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

  w->dry_run = dry_run;
  w->tracefile = NULL;

  w->rank = 0;
  w->n_ranks = 1;

  w->buffer_size = 0;

  w->f_map = cs_map_name_to_id_create();
  w->f_ts = NULL;

#if defined(HAVE_MPI)
  {
    int mpi_flag, rank, n_ranks;
    w->min_rank_step = 1;
    w->min_block_size = 1024*1024*8;
    w->block_comm = MPI_COMM_NULL;
    w->comm = MPI_COMM_NULL;
    MPI_Initialized(&mpi_flag);
    if (mpi_flag && comm != MPI_COMM_NULL) {
      w->comm = comm;
      MPI_Comm_rank(w->comm, &rank);
      MPI_Comm_size(w->comm, &n_ranks);
      w->rank = rank;
      w->n_ranks = n_ranks;
      if (rank_step < 1)
        rank_step = 1;
      else if (rank_step > n_ranks)
        rank_step = n_ranks;

#if defined(HAVE_MELISSA_MPI)
      w->min_rank_step = rank_step;
      if (rank_step > 1) {
        w->block_comm = cs_base_get_rank_step_comm(rank_step);
      }
      else
        w->block_comm = comm;
#else
      w->min_rank_step = n_ranks;
      w->block_comm = MPI_COMM_NULL;
#endif

      w->comm = comm;
    }
  }
#endif /* defined(HAVE_MPI) */



  if (w->dry_run == false)
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

  if (w->dry_run == false)
    damaris_finalize();

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
  CS_UNUSED(this_writer_p);
  CS_UNUSED(mesh);
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
