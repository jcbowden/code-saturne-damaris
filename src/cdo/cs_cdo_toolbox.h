#ifndef __CS_CDO_TOOLBOX_H__
#define __CS_CDO_TOOLBOX_H__

/*============================================================================
 * Basic operations: dot product, cross product, sum...
 *============================================================================*/

/*
  This file is part of Code_Saturne, a general-purpose CFD tool.

  Copyright (C) 1998-2017 EDF S.A.

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

/*----------------------------------------------------------------------------
 * Standard C library headers
 *----------------------------------------------------------------------------*/

#include <stdio.h>
#include <math.h>

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*============================================================================
 * Macro definitions
 *============================================================================*/

/*============================================================================
 * Type definitions
 *============================================================================*/

/* Structure enabling the repeated usage of local dense square matrix
   associated  with a local set of entities */
typedef struct {

  int         n_max_ent;  // max number of entities by primal cells
  int         n_ent;      // current number of entities
  cs_lnum_t  *ids;        // list of entity ids (size = n_max_ent)
  double     *val;        // local matrix (size: n_max_ent*n_max_ent)

} cs_locmat_t;

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the euclidean norm 2 of a vector of size len
 *         This algorithm tries to reduce round-off error thanks to
 *         intermediate sums.
 *
 *  \param[in] len     vector dimension
 *  \param[in] v       vector
 *
 * \return  the euclidean norm of a vector
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_euclidean_norm(int               len,
                  const cs_real_t   v[]);

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Compute the weighted sum of square values of an array
 *
 *  \param[in] n      size of arrays v and w
 *  \param[in] x      array of floating-point values
 *  \param[in] w      floating-point values of weights
 *
 * \return the result of this operation
 */
/*----------------------------------------------------------------------------*/

cs_real_t
cs_weighted_sum_squared(cs_lnum_t                   n,
                        const cs_real_t *restrict   x,
                        const cs_real_t *restrict   weight);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Allocate and initialize a cs_locmat_t structure
 *
 * \param[in]  n_max_ent    max number of entities
 *
 * \return  a new allocated cs_locmat_t structure
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_locmat_create(int   n_max_ent);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Free a cs_locmat_t structure
 *
 * \param[in]  lm    pointer to a cs_locmat_t struct. to free
 *
 * \return  a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_locmat_t *
cs_locmat_free(cs_locmat_t  *lm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Copy a cs_locmat_t structure into another cs_locmat_t structure
 *          which has been already allocated
 *
 * \param[in, out]  recv    pointer to a cs_locmat_t struct.
 * \param[in]       send    pointer to a cs_locmat_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_copy(cs_locmat_t        *recv,
               const cs_locmat_t  *send);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Compute a local dense matrix-vector product
 *          matvec has been previously allocated
 *
 * \param[in]      loc    local matrix to use
 * \param[in]      vec    local vector to use
 * \param[in, out] matvec result of the local matrix-vector product
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_matvec(const cs_locmat_t   *loc,
                 const cs_real_t     *vec,
                 cs_real_t           *matvec);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Add two local dense matrices: loc += add
 *
 * \param[in, out] loc   local matrix storing the result
 * \param[in]      add   values to add to loc
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_add(cs_locmat_t        *loc,
              const cs_locmat_t  *add);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Give the result of the following operation: loc = loc + alpha*add
 *
 * \param[in, out] loc    local matrix storing the result
 * \param[in]      alpha  multiplicative coefficient
 * \param[in]      add    values to add to loc
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_mult_add(cs_locmat_t        *loc,
                   cs_real_t           alpha,
                   const cs_locmat_t  *add);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Define a new matrix by adding a local matrix with its transpose.
 *          Keep the transposed matrix for future use.
 *
 * \param[in, out] loc   local matrix to transpose and add
 * \param[in, out] tr    transposed of the local matrix
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_add_transpose(cs_locmat_t  *loc,
                        cs_locmat_t  *tr);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Set the given matrix into its anti-symmetric part
 *
 * \param[in, out] loc   local matrix to transform
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_asymm(cs_locmat_t  *loc);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump a local discrete Hodge operator
 *
 * \param[in]    parent_id  id of the related parent entity
 * \param[in]    lm         pointer to the cs_sla_locmat_t struct.
 */
/*----------------------------------------------------------------------------*/

void
cs_locmat_dump(int                 parent_id,
               const cs_locmat_t  *lm);

/*----------------------------------------------------------------------------*/
/*!
 * \brief   Dump an array in the listing (For DEBUG)
 */
/*----------------------------------------------------------------------------*/

void
cs_dump_array_to_listing(const char        *header,
                         const cs_lnum_t    size,
                         const cs_real_t    array[],
                         int                n_cols);

void
cs_dump_integer_to_listing(const char        *header,
                           const cs_lnum_t    size,
                           const cs_lnum_t    array[],
                           int                n_cols);

/*----------------------------------------------------------------------------*/

END_C_DECLS

#endif /* __CS_TOOLBOX_H__ */
