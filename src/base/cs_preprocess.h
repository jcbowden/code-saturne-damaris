#ifndef __CS_PREPROCESS_H__
#define __CS_PREPROCESS_H__

/*============================================================================
 * Handle successive preprocessing operations.
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

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include "cs_defs.h"

#include "cs_halo.h"
#include "cs_mesh.h"
#include "cs_mesh_quantities.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*=============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/*=============================================================================
 * Public function prototypes
 *============================================================================*/
#if defined(HAVE_DAMARIS)
void
cs_preprocess_set_damaris_param_from_mesh(cs_mesh_t *m) ;
#endif

/*----------------------------------------------------------------------------*/
/*!
 * \brief Determine if preprocessing is needed.
 *
 * Preprocessing is ignored when a ./restart/mesh_input file is present but
 * no ./mesh_input file or directory is present. In this case, restart mesh
 * file is read, and all other preprocessing steps are skipped.
 *
 * \return  true if preprocessing is needed, false if only reading is needed.
 */
/*----------------------------------------------------------------------------*/

bool
cs_preprocess_mesh_is_needed(void);

/*----------------------------------------------------------------------------
 * Define all mesh preprocessing operations.
 *----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_define(void);

/*----------------------------------------------------------------------------
 * Apply all mesh preprocessing operations.
 *
 * parameters:
 *   halo_type  <->  type of halo (standard or extended)
 *----------------------------------------------------------------------------*/

void
cs_preprocess_mesh(cs_halo_type_t   halo_type);

/*----------------------------------------------------------------------------*/
/*!
 * \brief Apply numbering changes to ignore selected boundary faces.
 *
 * \param[in, out]  m         pointer to mesh structure
 * \param[in, out]  mq        pointer to mesh quantities structure
 * \param[in]       n_faces   number of selected faces
 * \param[in]       face_ids  ids of selected faces
 */
/*----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_selected_b_faces_ignore(cs_mesh_t             *m,
                                           cs_mesh_quantities_t  *mq,
                                           cs_lnum_t              n_faces,
                                           const cs_lnum_t        face_ids[]);

/*----------------------------------------------------------------------------
 * Update fortran arrays relative to the global mesh.
 *----------------------------------------------------------------------------*/

void
cs_preprocess_mesh_update_fortran(void);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_PREPROCESS_H__ */
