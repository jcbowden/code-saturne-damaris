/*============================================================================
 * Routines to handle cs_navsto_param_t structure
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
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#if defined(HAVE_MPI)
#include <mpi.h>
#endif

#if defined(HAVE_PETSC)
#include <petscversion.h>
#endif

/*----------------------------------------------------------------------------
 *  Local headers
 *----------------------------------------------------------------------------*/

#include <bft_mem.h>
#include <bft_printf.h>

#include "cs_base.h"
#include "cs_equation.h"
#include "cs_log.h"

/*----------------------------------------------------------------------------
 * Header for the current file
 *----------------------------------------------------------------------------*/

#include "cs_navsto_param.h"

/*----------------------------------------------------------------------------*/

BEGIN_C_DECLS

/*! \cond DOXYGEN_SHOULD_SKIP_THIS */

/*=============================================================================
 * Local Macro definitions
 *============================================================================*/

#define CS_NAVSTO_PARAM_DBG  0

/*============================================================================
 * Type definitions
 *============================================================================*/

/*============================================================================
 * Private variables
 *============================================================================*/

static const char
cs_navsto_param_coupling_name[CS_NAVSTO_N_COUPLINGS][CS_BASE_STRING_LEN] =
  { N_("Artificial compressibility algorithm"),
    N_("Monolithic"),
    N_("Incremental projection algorithm"),
  };

static const char _err_empty_nsp[] =
  N_(" %s: Stop setting an empty cs_navsto_param_t structure.\n"
     " Please check your settings.\n");

static const char
_space_scheme_key[CS_SPACE_N_SCHEMES][CS_BASE_STRING_LEN] =
  { "fv",
    "cdo_vb",
    "cdo_vcb",
    "cdo_eb",
    "cdo_fb",
    "hho_p0",
    "hho_p1",
    "hho_p2"
  };

static const char
_time_scheme_key[CS_TIME_N_SCHEMES][CS_BASE_STRING_LEN] =
  { "steady",
    "euler_implicit",
    "euler_explicit",
    "crank_nicolson",
    "theta_scheme"
  };

static const char
_dof_reduction_key[CS_PARAM_N_REDUCTIONS][CS_BASE_STRING_LEN] =
  { "derham",
    "average"
  };

static const char
_quad_type_key[CS_QUADRATURE_N_TYPES][CS_BASE_STRING_LEN] =
  { "none",
    "bary",
    "bary_subdiv",
    "higher",
    "highest"
  };

static const char
_adv_formulation_key[CS_PARAM_N_ADVECTION_FORMULATIONS][CS_BASE_STRING_LEN] =
  {
    "conservative",
    "non_conservative",
    "skew_symmetric"
  };

static const char
_adv_scheme_key[CS_PARAM_N_ADVECTION_SCHEMES][CS_BASE_STRING_LEN] =
  {
    "centered",
    "cip",
    "cip_cw",
    "mix_centered_upwind",
    "samarskii",
    "sg",
    "upwind"
  };

/*============================================================================
 * Private function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief Check if the prerequisite are fullfilled when a PETSC-related type
 *        of sles strategy is requested
 *
 * \param[in]  val          keyval
 * \param[in]  sles_type    type of SLES strategy
 *
 * \return the same sles_type if ok
 */
/*----------------------------------------------------------------------------*/

static inline cs_navsto_sles_t
_check_petsc_strategy(const char         *val,
                      cs_navsto_sles_t    sles_type)
{
#if defined(HAVE_PETSC)
#if PETSC_VERSION_GE(3,11,0)
  CS_UNUSED(val);
  return sles_type;
#else
  if (sles_type == CS_NAVSTO_SLES_GKB_GMRES ||
      sles_type == CS_NAVSTO_SLES_GKB_PETSC)
    bft_error(__FILE__, __LINE__, 0,
              "%s: PETSc version greater or equal to 3.11 is required"
              " when using the keyval \"%s\"\n", __func__, val);
  return sles_type;
#endif
#else
  bft_error(__FILE__, __LINE__, 0,
            " %s: \"CS_NSKEY_SLES_STRATEGY\" keyval %s requires"
            " an installation with PETSC\n", __func__, val);
  return CS_NAVSTO_SLES_N_TYPES;
#endif
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         momentum equation according to the type of coupling
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the corresponding \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static inline cs_equation_param_t *
_get_momentum_param(cs_navsto_param_t    *nsp)
{
  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    return cs_equation_param_by_name("momentum");

  case CS_NAVSTO_COUPLING_PROJECTION:
    return cs_equation_param_by_name("velocity_prediction");
    break;

  default:
    return NULL;

  }  /* Switch */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         momentum equation according to the type of coupling
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 *
 * \return a pointer to the corresponding \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

static void
_propagate_qtype(cs_navsto_param_t    *nsp)
{
  /* Loop on velocity ICs */
  for (int i = 0; i < nsp->n_velocity_ic_defs; i++)
    cs_xdef_set_quadrature(nsp->velocity_ic_defs[i], nsp->qtype);

  /* Loop on pressure ICs */
  for (int i = 0; i < nsp->n_pressure_ic_defs; i++)
    cs_xdef_set_quadrature(nsp->pressure_ic_defs[i], nsp->qtype);

  /* Loop on velocity BCs */
  for (int i = 0; i < nsp->n_velocity_bc_defs; i++)
    cs_xdef_set_quadrature(nsp->velocity_bc_defs[i], nsp->qtype);

  /* Loop on pressure BCs */
  for (int i = 0; i < nsp->n_pressure_bc_defs; i++)
    cs_xdef_set_quadrature(nsp->pressure_bc_defs[i], nsp->qtype);
}

/*! (DOXYGEN_SHOULD_SKIP_THIS) \endcond */

/*============================================================================
 * Public function prototypes
 *============================================================================*/

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Create a new structure to store all numerical parameters related
 *         to the resolution of the Navier-Stokes (NS) system
 *
 * \param[in]  boundaries       pointer to a cs_boundary_t structure
 * \param[in]  model            model related to the NS system to solve
 * \param[in]  algo_coupling    algorithm used for solving the NS system
 * \param[in]  option_flag      additional high-level numerical options
 * \param[in]  post_flag        predefined post-processings
 *
 * \return a pointer to a new allocated structure
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_create(const cs_boundary_t             *boundaries,
                       cs_navsto_param_model_t          model,
                       cs_navsto_param_coupling_t       algo_coupling,
                       cs_flag_t                        option_flag,
                       cs_flag_t                        post_flag)
{
  cs_navsto_param_t  *param = NULL;
  BFT_MALLOC(param, 1, cs_navsto_param_t);

  /* Flags and indicators */
  param->verbosity = 1;
  param->post_flag = post_flag;

  /* Physical modelling */
  /* ------------------ */

  /* Which equations are solved and which terms are needed */
  param->model = model;
  param->reference_pressure = 0.;
  param->phys_constants = cs_get_glob_physical_constants();

  /* Turbulence modelling (pointer to global structures) */
  param->turbulence = cs_turbulence_param_create();

  /* Main set of properties */
  param->mass_density = cs_property_by_name(CS_PROPERTY_MASS_DENSITY);
  if (param->mass_density == NULL)
    param->mass_density = cs_property_add(CS_PROPERTY_MASS_DENSITY,
                                          CS_PROPERTY_ISO);

  param->lam_viscosity = cs_property_add(CS_NAVSTO_LAM_VISCOSITY,
                                         CS_PROPERTY_ISO);

  if (param->turbulence->model->iturb == CS_TURB_NONE)
    param->tot_viscosity = param->lam_viscosity;
  else
    param->tot_viscosity = cs_property_add(CS_NAVSTO_TOTAL_VISCOSITY,
                                           CS_PROPERTY_ISO);

  /* Default numerical settings */
  /* -------------------------- */

  param->option_flag = option_flag;
  param->coupling = algo_coupling;
  param->space_scheme = CS_SPACE_SCHEME_CDOFB;
  param->dof_reduction_mode = CS_PARAM_REDUCTION_AVERAGE;
  param->qtype = CS_QUADRATURE_BARY;

  param->adv_form   = CS_PARAM_ADVECTION_FORM_NONCONS;
  param->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;

  /* Forcing steady state in order to avoid inconsistencies */
  if (option_flag & CS_NAVSTO_FLAG_STEADY)
    param->time_scheme = CS_TIME_SCHEME_STEADY;
  else
    param->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;
  param->theta = 1.0;

  /* Resolution parameters (inner linear system then the non-linear system )*/
  param->sles_param.strategy = CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK;
  param->sles_param.n_max_il_algo_iter = 100;
  param->sles_param.il_algo_rtol = 1e-08;
  param->sles_param.il_algo_atol = 1e-08;
  param->sles_param.il_algo_dtol = 1e3;
  param->sles_param.il_algo_verbosity = 0;

  param->sles_param.nl_algo = CS_NAVSTO_NL_PICARD_ALGO;
  param->sles_param.n_max_nl_algo_iter = 25;
  param->sles_param.nl_algo_rtol = 1e-5;
  param->sles_param.nl_algo_atol = 1e-5;
  param->sles_param.nl_algo_dtol = 1e3;
  param->sles_param.nl_algo_verbosity = 1;

  /* Management of the outer resolution steps (i.e. the full system including
     the turbulence modelling or the the thermal system) */
  param->n_max_outer_iter = 5;
  param->delta_thermal_tolerance = 1e-2;

  /* Physical boundaries specific to the problem at stake */
  param->boundaries = boundaries; /* shared structure */

  /* Remark: As velocity and pressure may not be associated to an equation
     directly, one stores the definition of initial conditions and boundary
     conditions at this level */

  switch (algo_coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    param->gd_scale_coef = 1.0;    /* Default value if not set by the user */

    param->velocity_ic_is_owner = false;
    param->velocity_bc_is_owner = false;
    param->pressure_ic_is_owner = true;
    param->pressure_bc_is_owner = true;
    break;

  case CS_NAVSTO_COUPLING_MONOLITHIC:
    param->sles_param.strategy = CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK;
    param->gd_scale_coef = 0.0;    /* Default value if not set by the user */

    param->velocity_ic_is_owner = false;
    param->velocity_bc_is_owner = false;
    param->pressure_ic_is_owner = true;
    param->pressure_bc_is_owner = true;
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    param->gd_scale_coef = 0.0;    /* Default value if not set by the user */

    param->velocity_ic_is_owner = false;
    param->velocity_bc_is_owner = false;
    param->pressure_ic_is_owner = false;
    param->pressure_bc_is_owner = false;
    break;

  default:
    /* Nothing done */
    break;
  }

  /* Initial conditions for the pressure field */
  param->n_pressure_ic_defs = 0;
  param->pressure_ic_defs = NULL;

  /* Initial conditions for the velocity field */
  param->n_velocity_ic_defs = 0;
  param->velocity_ic_defs = NULL;

  /* Boundary conditions for the pressure field */
  param->n_pressure_bc_defs = 0;
  param->pressure_bc_defs = NULL;

  /* Boundary conditions for the velocity field */
  param->n_velocity_bc_defs = 0;
  param->velocity_bc_defs = NULL;

  /* Enforcement of a solid zone */
  param->n_solid_cells = 0;
  param->solid_cell_ids = NULL;

  return param;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Free a \ref cs_navsto_param_t structure
 *
 * \param[in, out]  param    pointer to a \ref cs_navsto_param_t structure
 *
 * \return a NULL pointer
 */
/*----------------------------------------------------------------------------*/

cs_navsto_param_t *
cs_navsto_param_free(cs_navsto_param_t    *param)
{
  if (param == NULL)
    return param;

  /* Turbulence modelling */

  BFT_FREE(param->turbulence);

  /* Velocity initial conditions */

  if (param->n_velocity_ic_defs > 0) {

    /* Otherwise this is freed inside the related equation */
    if (param->velocity_ic_is_owner) {
      for (int i = 0; i < param->n_velocity_ic_defs; i++)
        param->velocity_ic_defs[i] = cs_xdef_free(param->velocity_ic_defs[i]);
    }
    BFT_FREE(param->velocity_ic_defs);
    param->velocity_ic_defs = NULL;

  }

  /* Pressure initial conditions */

  if (param->n_pressure_ic_defs > 0) {

    if (param->pressure_ic_is_owner) {
      for (int i = 0; i < param->n_pressure_ic_defs; i++)
        param->pressure_ic_defs[i] = cs_xdef_free(param->pressure_ic_defs[i]);
    }
    BFT_FREE(param->pressure_ic_defs);
    param->pressure_ic_defs = NULL;

  }

  /* Velocity boundary conditions */

  if (param->n_velocity_bc_defs > 0) {
    if (param->velocity_bc_is_owner) {
      for (int i = 0; i < param->n_velocity_bc_defs; i++)
        param->velocity_bc_defs[i] = cs_xdef_free(param->velocity_bc_defs[i]);
    }
    BFT_FREE(param->velocity_bc_defs);
    param->velocity_bc_defs = NULL;
  }

  /* Pressure boundary conditions */

  if (param->n_pressure_bc_defs > 0) {

    if (param->pressure_bc_is_owner) {
      for (int i = 0; i < param->n_pressure_bc_defs; i++)
        param->pressure_bc_defs[i] = cs_xdef_free(param->pressure_bc_defs[i]);
    }
    BFT_FREE(param->pressure_bc_defs);
    param->pressure_bc_defs = NULL;

  }

  BFT_FREE(param->solid_cell_ids);

  /* Free the main structure */
  BFT_FREE(param);

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set a parameter attached to a keyname in a \ref cs_navsto_param_t
 *         structure
 *
 * \param[in, out] nsp      pointer to a \ref cs_navsto_param_t structure to set
 * \param[in]      key      key related to the member of eq to set
 * \param[in]      keyval   accessor to the value to set
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_set(cs_navsto_param_t    *nsp,
                    cs_navsto_key_t       key,
                    const char           *keyval)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  if (keyval == NULL)
    bft_error(__FILE__, __LINE__, 0, "%s: Empty key value.\n", __func__);

  /* Conversion of the string to lower case */
  char val[CS_BASE_STRING_LEN];
  for (size_t i = 0; i < strlen(keyval); i++)
    val[i] = tolower(keyval[i]);
  val[strlen(keyval)] = '\0';

  switch(key) {

  case CS_NSKEY_ADVECTION_FORMULATION:
    if (strcmp(val, "conservative") == 0)
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_CONSERV;
    else if (strcmp(val, "non_conservative") == 0)
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_NONCONS;
    else if (strcmp(val, "skew_symmetric") == 0)
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_SKEWSYM;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key"
                  " CS_NSKEY_ADVECTION_FORMULATION\n"
                  " Choice between conservative, non_conservative"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_ADVECTION_SCHEME:
    if (strcmp(val, "upwind") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_UPWIND;
    else if (strcmp(val, "samarskii") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SAMARSKII;
    else if (strcmp(val, "sg") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_SG;
    else if (strcmp(val, "centered") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CENTERED;
    else if (strcmp(val, "mix_centered_upwind") == 0 ||
             strcmp(val, "hybrid_centered_upwind") == 0)
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_HYBRID_CENTERED_UPWIND;
    else if (strcmp(val, "cip") == 0) {
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP;
      /* Automatically switch to a non-conservative formulation */
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_NONCONS;
    }
    else if (strcmp(val, "cip_cw") == 0) {
      nsp->adv_scheme = CS_PARAM_ADVECTION_SCHEME_CIP_CW;
      /* Automatically switch to a non-conservative formulation */
      nsp->adv_form = CS_PARAM_ADVECTION_FORM_NONCONS;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key"
                  " CS_NSKEY_ADVECTION_SCHEME\n"
                  " Choice between upwind, samarskii, sg, centered, cip, cip_cw,"
                  " mix_centered_upwind"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_DOF_REDUCTION:
    if (strcmp(val, "derham") == 0)
      nsp->dof_reduction_mode = CS_PARAM_REDUCTION_DERHAM;
    else if (strcmp(val, "average") == 0)
      nsp->dof_reduction_mode = CS_PARAM_REDUCTION_AVERAGE;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_DOF_REDUCTION\n"
                  " Choice between \"derham\" or \"average\"."),
                __func__, _val);
    }
    break;

  case CS_NSKEY_GD_SCALE_COEF:
    switch (nsp->coupling) {
    case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
    case CS_NAVSTO_COUPLING_MONOLITHIC:
      nsp->gd_scale_coef = atof(val);
      break;

    case CS_NAVSTO_COUPLING_PROJECTION:
      cs_base_warn(__FILE__, __LINE__);
      bft_printf(" %s: Trying to set the zeta parameter with the %s\n "
                 " although this will not have use in the algorithm.\n",
                 __func__, cs_navsto_param_coupling_name[nsp->coupling]);

    default:
      break;

    } /* End of switch */
    break;

  case CS_NSKEY_IL_ALGO_ATOL:
    nsp->sles_param.il_algo_atol = atof(val);
    if (nsp->sles_param.il_algo_rtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the absolute tolerance\n", __func__);
    break;

  case CS_NSKEY_IL_ALGO_DTOL:
    nsp->sles_param.il_algo_dtol = atof(val);
    if (nsp->sles_param.il_algo_dtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the divergence tolerance\n", __func__);
    break;

  case CS_NSKEY_IL_ALGO_RTOL:
    nsp->sles_param.il_algo_rtol = atof(val);
    if (nsp->sles_param.il_algo_rtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the residual tolerance\n", __func__);
    break;

  case CS_NSKEY_IL_ALGO_VERBOSITY:
    nsp->sles_param.il_algo_verbosity = atoi(val);
    break;

  case CS_NSKEY_MAX_IL_ALGO_ITER:
    nsp->sles_param.n_max_il_algo_iter = atoi(val);
    break;

  case CS_NSKEY_MAX_NL_ALGO_ITER:
    nsp->sles_param.n_max_nl_algo_iter = atoi(val);
    break;

  case CS_NSKEY_MAX_OUTER_ITER:
    nsp->n_max_outer_iter = atoi(val);
    break;

  case CS_NSKEY_NL_ALGO:
    {
      if (strcmp(val, "picard") == 0)
        nsp->sles_param.nl_algo = CS_NAVSTO_NL_PICARD_ALGO;
      else if (strcmp(val, "fixed-point") == 0)
        nsp->sles_param.nl_algo = CS_NAVSTO_NL_PICARD_ALGO;
      else {
        const char *_val = val;
        bft_error(__FILE__, __LINE__, 0,
                  " %s: Invalid value \"%s\" for key CS_NSKEY_NL_ALGO\n"
                  " Valid choices are \"picard\", \"fixed-point\".",
                  __func__, _val);
      }

    }
    break; /* Non-linear algorithm */

  case CS_NSKEY_NL_ALGO_ATOL:
    nsp->sles_param.nl_algo_atol = atof(val);
    if (nsp->sles_param.il_algo_atol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the absolute tolerance of the"
                " non-linear algorithm\n",
                __func__);
    break;

  case CS_NSKEY_NL_ALGO_DTOL:
    nsp->sles_param.nl_algo_dtol = atof(val);
    if (nsp->sles_param.nl_algo_dtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the divergence tolerance of the"
                " non-linear algorithm\n",
                __func__);
    break;

  case CS_NSKEY_NL_ALGO_RTOL:
    nsp->sles_param.nl_algo_rtol = atof(val);
    if (nsp->sles_param.nl_algo_rtol < 0)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for the relative tolerance of the"
                " non-linear algorithm\n",
                __func__);
    break;

  case CS_NSKEY_NL_ALGO_VERBOSITY:
    nsp->sles_param.nl_algo_verbosity = atoi(val);
    break;

  case CS_NSKEY_QUADRATURE:
    {
      nsp->qtype = CS_QUADRATURE_NONE;

      if (strcmp(val, "bary") == 0)
        nsp->qtype = CS_QUADRATURE_BARY;
      else if (strcmp(val, "bary_subdiv") == 0)
        nsp->qtype = CS_QUADRATURE_BARY_SUBDIV;
      else if (strcmp(val, "higher") == 0)
        nsp->qtype = CS_QUADRATURE_HIGHER;
      else if (strcmp(val, "highest") == 0)
        nsp->qtype = CS_QUADRATURE_HIGHEST;
      else {
        const char *_val = val;
        bft_error(__FILE__, __LINE__, 0,
                  _(" %s: Invalid value \"%s\" for key CS_NSKEY_QUADRATURE\n"
                    " Valid choices are \"bary\", \"bary_subdiv\", \"higher\""
                    " and \"highest\"."), __func__, _val);
      }

      _propagate_qtype(nsp);
    }
    break; /* Quadrature */

  case CS_NSKEY_SLES_STRATEGY:
    if (strcmp(val, "no_block") == 0)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK;
    else if (strcmp(val, "by_blocks") == 0)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_BY_BLOCKS;
    else if (strcmp(val, "block_amg_cg") == 0)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG;
    else if (strcmp(val, "gkb_saturne") == 0 ||
             strcmp(val, "gkb") == 0)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_GKB_SATURNE;
    else if (strcmp(val, "uzawa_al") == 0 || strcmp(val, "alu") == 0)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_UZAWA_AL;

    /* All the following options need either PETSC or MUMPS */
    /* ---------------------------------------------------- */

    else if (strcmp(val, "additive_gmres") == 0)
      nsp->sles_param.strategy =
        _check_petsc_strategy(val, CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK);
    else if (strcmp(val, "multiplicative_gmres") == 0)
      nsp->sles_param.strategy =
        _check_petsc_strategy(val,
                              CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK);
    else if (strcmp(val, "diag_schur_gmres") == 0)
      nsp->sles_param.strategy =
        _check_petsc_strategy(val, CS_NAVSTO_SLES_DIAG_SCHUR_GMRES);
    else if (strcmp(val, "upper_schur_gmres") == 0)
      nsp->sles_param.strategy =
        _check_petsc_strategy(val, CS_NAVSTO_SLES_UPPER_SCHUR_GMRES);
    else if (strcmp(val, "gkb_gmres") == 0)
      nsp->sles_param.strategy =
        _check_petsc_strategy(val, CS_NAVSTO_SLES_GKB_GMRES);
    else if (strcmp(val, "gkb_petsc") == 0)
      nsp->sles_param.strategy =
        _check_petsc_strategy(val, CS_NAVSTO_SLES_GKB_PETSC);

    else if (strcmp(val, "mumps") == 0) {
#if defined(HAVE_MUMPS)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_MUMPS;
#else
#if defined(HAVE_PETSC)
#if defined(PETSC_HAVE_MUMPS)
      nsp->sles_param.strategy = CS_NAVSTO_SLES_MUMPS;
#else
      bft_error(__FILE__, __LINE__, 0,
                " %s: Error detected while setting \"%s\" key\n"
                " MUMPS is not available with your installation.\n"
                " Please check your installation settings.\n",
                __func__, "CS_NSKEY_SLES_STRATEGY");
#endif  /* PETSC_HAVE_MUMPS */
#endif  /* HAVE_PETSC */
#endif  /* HAVE_MUMPS */
    }

    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid val %s related to key CS_NSKEY_SLES_STRATEGY\n"
                " Choice between: no_block, by_locks, block_amg_cg,\n"
                " {additive,multiplicative}_gmres, {diag,upper}_schur_gmres,\n"
                " gkb, gkb_petsc, gkb_gmres, gkb_saturne,\n"
                " mumps, uzawa_al or alu", __func__, _val);
    }
    break;

  case CS_NSKEY_SPACE_SCHEME:
    if (strcmp(val, "cdo_fb") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_CDOFB;
    }
    else if (strcmp(val, "hho_p0") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P0;
    }
    else if (strcmp(val, "hho_p1") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P1;
    }
    else if (strcmp(val, "hho_p2") == 0) {
      nsp->space_scheme = CS_SPACE_SCHEME_HHO_P2;
    }
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid val %s related to key CS_NSKEY_SPACE_SCHEME\n"
                  " Choice between hho_{p0, p1, p2} or cdo_fb"),
                __func__, _val);
    }
    break;

  case CS_NSKEY_THERMAL_TOLERANCE:
    nsp->delta_thermal_tolerance = atof(val);
    /* If tolerance is set to a negative value then it stops the outer
       iteration process after the first iteration */
    break;

  case CS_NSKEY_TIME_SCHEME:
    if (strcmp(val, "euler_implicit") == 0) {
      nsp->time_scheme = CS_TIME_SCHEME_EULER_IMPLICIT;
      nsp->theta = 1.;
    }
    else if (strcmp(val, "euler_explicit") == 0) {
      nsp->time_scheme = CS_TIME_SCHEME_EULER_EXPLICIT;
      nsp->theta = 0.;
    }
    else if (strcmp(val, "crank_nicolson") == 0) {
      nsp->time_scheme = CS_TIME_SCHEME_CRANKNICO;
      nsp->theta = 0.5;
    }
    else if (strcmp(val, "theta_scheme") == 0)
      nsp->time_scheme = CS_TIME_SCHEME_THETA;
    else {
      const char *_val = val;
      bft_error(__FILE__, __LINE__, 0,
                _(" %s: Invalid value \"%s\" for CS_EQKEY_TIME_SCHEME\n"
                  " Valid choices are \"euler_implicit\","
                  " \"euler_explicit\"," " \"crank_nicolson\","
                  " and \"theta_scheme\"."), __func__, _val);
    }
    break;

  case CS_NSKEY_TIME_THETA:
    nsp->theta = atof(val);
    if (nsp->theta < 0. - cs_math_zero_threshold ||
        nsp->theta > 1.0 + cs_math_zero_threshold)
      bft_error(__FILE__, __LINE__, 0,
                " %s: Invalid value for theta\n", __func__);
    break;

  case CS_NSKEY_VERBOSITY:
    nsp->verbosity = atoi(val);
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              _(" %s: Invalid key for setting the Navier-Stokes system."),
              __func__);

  }

}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Apply the numerical settings defined for the Navier-Stokes system
 *         to an equation related to this system.
 *
 * \param[in]       nsp    pointer to a \ref cs_navsto_param_t structure
 * \param[in, out]  eqp    pointer to a \ref cs_equation_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_transfer(const cs_navsto_param_t    *nsp,
                         cs_equation_param_t        *eqp)
{
  assert(nsp != NULL && eqp != NULL);

  /*  Set the space discretization scheme */
  const char  *ss_key = _space_scheme_key[nsp->space_scheme];

  cs_equation_set_param(eqp, CS_EQKEY_SPACE_SCHEME, ss_key);

  /*  Set the time discretization scheme */
  const char  *ts_key = _time_scheme_key[nsp->time_scheme];

  cs_equation_set_param(eqp, CS_EQKEY_TIME_SCHEME, ts_key);
  if (nsp->time_scheme == CS_TIME_SCHEME_THETA) {
    char  cvalue[36]; /* include '\0' */
    snprintf(cvalue, 35*sizeof(char), "%g", nsp->theta);
    cs_equation_set_param(eqp, CS_EQKEY_TIME_THETA, cvalue);
  }

  /*  Set the way DoFs are defined */
  const char  *dof_key = _dof_reduction_key[nsp->dof_reduction_mode];

  cs_equation_set_param(eqp, CS_EQKEY_DOF_REDUCTION, dof_key);

  /*  Set quadratures type */
  const char  *quad_key = _quad_type_key[nsp->qtype];

  /* If requested, add advection */
  if ((nsp->model & (CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES |
                     CS_NAVSTO_MODEL_OSEEN)) > 0) {

    /* If different from default value */
    const char *form_key = _adv_formulation_key[nsp->adv_form];
    cs_equation_set_param(eqp, CS_EQKEY_ADV_FORMULATION, form_key);

    const char *scheme_key = _adv_scheme_key[nsp->adv_scheme];
    cs_equation_set_param(eqp, CS_EQKEY_ADV_SCHEME, scheme_key);

  }

  cs_equation_set_param(eqp, CS_EQKEY_BC_QUADRATURE, quad_key);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Summary of the main cs_navsto_param_t structure
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_param_log(const cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    return;

  /* Sanity checks */
  if (nsp->model < 1)
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid model for Navier-Stokes.\n",
              __func__);
  if (nsp->coupling == CS_NAVSTO_N_COUPLINGS)
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid way of coupling the Navier-Stokes equations.\n",
              __func__);

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Verbosity: %d\n", nsp->verbosity);
  if (cs_navsto_param_is_steady(nsp))
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Time status: Steady\n");
  else
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Time status: Unsteady\n");

  /* Describe the physical modelling */

  if (nsp->model & CS_NAVSTO_MODEL_STOKES) {
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  "Stokes velocity-pressure system");
    assert(!(nsp->model & CS_NAVSTO_MODEL_OSEEN) &&
           !(nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES));
  }
  else if (nsp->model & CS_NAVSTO_MODEL_OSEEN) {
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  "Oseen velocity-pressure system");
    assert(!(nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES));
  }
  else if (nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  "Incompressible Navier-Stokes velocity-pressure system");

  if (nsp->model & CS_NAVSTO_MODEL_GRAVITY_EFFECTS)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  "Gravity effect activated");

  if (nsp->model & CS_NAVSTO_MODEL_CORIOLIS_EFFECTS)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  "Coriolis effect activated");

  if (nsp->model & CS_NAVSTO_MODEL_BOUSSINESQ)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  " Boussinesq approximation activated");
  if (nsp->model & CS_NAVSTO_MODEL_SOLIDIFICATION_BOUSSINESQ)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Model: %s\n",
                  " Boussinesq approximation for solidification activated");

  /* Describe the coupling algorithm */
  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Coupling: %s\n",
                cs_navsto_param_coupling_name[nsp->coupling]);

  /* Describe if needed the SLES settings for the non-linear algorithm */
  if (nsp->model & CS_NAVSTO_MODEL_INCOMPRESSIBLE_NAVIER_STOKES) {

    const char algo_name[] = "Picard";
    if (nsp->sles_param.nl_algo != CS_NAVSTO_NL_PICARD_ALGO)
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid non-linear algo.", __func__);

    cs_log_printf(CS_LOG_SETUP, "  * NavSto | %s: rtol: %5.3e;"
                  " atol: %5.3e; dtol: %5.3e\n",
                  algo_name, nsp->sles_param.nl_algo_rtol,
                  nsp->sles_param.nl_algo_atol, nsp->sles_param.nl_algo_dtol);
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | %s: Max.Iterations: %d\n",
                  algo_name, nsp->sles_param.n_max_nl_algo_iter);

  }

  /* Describe the strategy to inverse the (inner) linear system */
  const cs_navsto_param_sles_t  nslesp = nsp->sles_param;

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | SLES.Strategy: ");
  switch (nslesp.strategy) {

  case CS_NAVSTO_SLES_EQ_WITHOUT_BLOCK:
    cs_log_printf(CS_LOG_SETUP, "No specific strategy. System as it is.\n");
    break;
  case CS_NAVSTO_SLES_BLOCK_MULTIGRID_CG:
    cs_log_printf(CS_LOG_SETUP, "Block AMG + CG\n");
    break;
  case CS_NAVSTO_SLES_ADDITIVE_GMRES_BY_BLOCK:
    cs_log_printf(CS_LOG_SETUP, "Additive block preconditioner + GMRES\n");
    break;
  case CS_NAVSTO_SLES_MULTIPLICATIVE_GMRES_BY_BLOCK:
    cs_log_printf(CS_LOG_SETUP, "Multiplicative block preconditioner + GMRES\n");
    break;
  case CS_NAVSTO_SLES_DIAG_SCHUR_GMRES:
    cs_log_printf(CS_LOG_SETUP, "Diag. block preconditioner with Schur approx."
                  " + GMRES\n");
    break;
  case CS_NAVSTO_SLES_UPPER_SCHUR_GMRES:
    cs_log_printf(CS_LOG_SETUP, "Upper block preconditioner with Schur approx."
                  " + GMRES\n");
    break;
  case CS_NAVSTO_SLES_GKB_PETSC:
    cs_log_printf(CS_LOG_SETUP, "GKB algorithm (through PETSc)\n");
    break;
  case CS_NAVSTO_SLES_GKB_GMRES:
    cs_log_printf(CS_LOG_SETUP, "GMRES with a GKB preconditioner\n");
    break;
  case CS_NAVSTO_SLES_GKB_SATURNE:
    cs_log_printf(CS_LOG_SETUP, "GKB algorithm (In-House)\n");
    break;
  case CS_NAVSTO_SLES_MUMPS:
    cs_log_printf(CS_LOG_SETUP, "LU factorization with MUMPS\n");
    break;
  case CS_NAVSTO_SLES_UZAWA_AL:
    cs_log_printf(CS_LOG_SETUP, "Augmented Lagrangian-Uzawa\n");
    break;

  default:
    cs_log_printf(CS_LOG_SETUP, "Not set\n");
    break;
  }

  if (nsp->gd_scale_coef > 0)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Grad-div scaling %e\n",
                  nsp->gd_scale_coef);

  cs_log_printf(CS_LOG_SETUP, "  * NavSto | InnerLinear.Algo.Tolerances:"
                " rtol: %5.3e; atol: %5.3e; dtol: %5.3e\n",
                nslesp.il_algo_rtol, nslesp.il_algo_atol, nslesp.il_algo_dtol);

  const char *space_scheme = cs_param_get_space_scheme_name(nsp->space_scheme);
  if (nsp->space_scheme < CS_SPACE_N_SCHEMES)
    cs_log_printf(CS_LOG_SETUP, "  * NavSto | Space scheme: %s\n",
                  space_scheme);
  else
    bft_error(__FILE__, __LINE__, 0,
              " %s: Undefined space scheme.", __func__);

  if (!cs_navsto_param_is_steady(nsp)) {

    const char  *time_scheme = cs_param_get_time_scheme_name(nsp->time_scheme);
    if (time_scheme != NULL) {
      cs_log_printf(CS_LOG_SETUP, "  * NavSto | Time scheme: %s", time_scheme);
      if (nsp->time_scheme == CS_TIME_SCHEME_THETA)
        cs_log_printf(CS_LOG_SETUP, " with value %f\n", nsp->theta);
      else
        cs_log_printf(CS_LOG_SETUP, "\n");
    }
    else
      bft_error(__FILE__, __LINE__, 0, "%s: Invalid time scheme.", __func__);

  }

  /* Default quadrature type */
  cs_log_printf(CS_LOG_SETUP, "  * NavSto | Default quadrature: %s\n",
                cs_quadrature_get_type_name(nsp->qtype));

  /* Initial conditions for the velocity */
  char  prefix[256];

  cs_log_printf(CS_LOG_SETUP,
                "  * NavSto | Velocity.Init.Cond | Number of definitions %2d\n",
                nsp->n_velocity_ic_defs);

  for (int i = 0; i < nsp->n_velocity_ic_defs; i++) {
    sprintf(prefix, "  * NavSto | Velocity.Init.Cond | Definition %2d", i);
    cs_xdef_log(prefix, nsp->velocity_ic_defs[i]);
  }

  /* Initial conditions for the pressure */
  cs_log_printf(CS_LOG_SETUP,
                "  * NavSto | Pressure.Init.Cond | Number of definitions: %d\n",
                nsp->n_pressure_ic_defs);
  for (int i = 0; i < nsp->n_pressure_ic_defs; i++) {
    sprintf(prefix, "  * NavSto | Pressure.Init.Cond | Definition %2d", i);
    cs_xdef_log(prefix, nsp->pressure_ic_defs[i]);
  }
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the \ref cs_equation_param_t structure related to the
 *         velocity equation (momentum equation in most of the cases)
 *
 * \param[in]  nsp    pointer to a cs_navsto_param_t structure
 *
 * \return a pointer to the set of parameters related to the momentum equation
 */
/*----------------------------------------------------------------------------*/

cs_equation_param_t *
cs_navsto_param_get_velocity_param(const cs_navsto_param_t    *nsp)
{
  cs_equation_param_t  *eqp = NULL;

  switch (nsp->coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
    eqp = cs_equation_param_by_name("momentum");
    break;

  case CS_NAVSTO_COUPLING_PROJECTION:
    eqp = cs_equation_param_by_name("velocity_prediction");
    break;

  default:
    bft_error(__FILE__, __LINE__, 0,
              "%s: Invalid coupling algorithm", __func__);
    break;

  }

  return eqp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Retrieve the name of the coupling algorithm
 *
 * \param[in]     coupling    a \ref cs_navsto_param_coupling_t
 *
 * \return the name
 */
/*----------------------------------------------------------------------------*/

const char *
cs_navsto_param_get_coupling_name(cs_navsto_param_coupling_t  coupling)
{
  switch (coupling) {

  case CS_NAVSTO_COUPLING_ARTIFICIAL_COMPRESSIBILITY:
  case CS_NAVSTO_COUPLING_MONOLITHIC:
  case CS_NAVSTO_COUPLING_PROJECTION:
    return cs_navsto_param_coupling_name[coupling];

  default:
    bft_error(__FILE__, __LINE__, 0, "%s: Invalid coupling.", __func__);
    break;

  }

  return NULL;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the value to consider for the reference pressure
 *
 * \param[in]  nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]  pref      value of the reference pressure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_reference_pressure(cs_navsto_param_t    *nsp,
                                 cs_real_t             pref)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  nsp->reference_pressure = pref;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *val)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_xdef_t  *d = NULL;
  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  if (eqp != NULL) { /* An equation related to the velocity is defined */

    d = cs_equation_add_ic_by_value(eqp, z_name, val);

  }
  else { /* No momentum equation available with the choice of velocity-pressure
            coupling */

    nsp->velocity_ic_is_owner = true;

    /* Add a new cs_xdef_t structure */
    int z_id = 0;
    if (z_name != NULL && strlen(z_name) > 0)
      z_id = (cs_volume_zone_by_name(z_name))->id;

    cs_flag_t  meta_flag = 0;
    if (z_id == 0)
      meta_flag |= CS_FLAG_FULL_LOC;

    d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                              3,  /* dim */
                              z_id,
                              CS_FLAG_STATE_UNIFORM,
                              meta_flag,
                              val);
  }

  int  new_id = nsp->n_velocity_ic_defs;
  nsp->n_velocity_ic_defs += 1;
  BFT_REALLOC(nsp->velocity_ic_defs, nsp->n_velocity_ic_defs, cs_xdef_t *);
  nsp->velocity_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the velocity unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_velocity_ic_by_analytic(cs_navsto_param_t      *nsp,
                                      const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_xdef_t  *d = NULL;
  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  if (eqp != NULL) { /* An equation related to the velocity is defined */

    d = cs_equation_add_ic_by_analytic(eqp, z_name, analytic, input);

  }
  else { /* No momentum equation available with the choice of velocity-pressure
            coupling */

    nsp->velocity_ic_is_owner = true;

    /* Add a new cs_xdef_t structure */
    int z_id = cs_get_vol_zone_id(z_name);

    cs_flag_t  meta_flag = 0;
    if (z_id == 0)
      meta_flag |= CS_FLAG_FULL_LOC;

    cs_xdef_analytic_context_t  anai = { .z_id = z_id,
                                         .func = analytic,
                                         .input = input,
                                         .free_input = NULL };

    d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                              3,  /* dim */
                              z_id,
                              0,  /* state flag */
                              meta_flag,
                              &anai);
  }

  /* Assign the default quadrature type of the Navier-Stokes module to this
   * definition (this can be modified by the user if the same call is performed
   * in cs_user_finalize_setup()) */
  cs_xdef_set_quadrature(d, nsp->qtype);

  int  new_id = nsp->n_velocity_ic_defs;
  nsp->n_velocity_ic_defs += 1;
  BFT_REALLOC(nsp->velocity_ic_defs, nsp->n_velocity_ic_defs, cs_xdef_t *);
  nsp->velocity_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set to a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      val       pointer to the value
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_pressure_ic_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *val)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Add a new cs_xdef_t structure */
  int z_id = 0;
  if (z_name != NULL && strlen(z_name) > 0)
    z_id = (cs_volume_zone_by_name(z_name))->id;

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_VALUE,
                                        1,  /* dim */
                                        z_id,
                                        CS_FLAG_STATE_UNIFORM,
                                        meta_flag,
                                        val);

  int  new_id = nsp->n_pressure_ic_defs;
  nsp->n_pressure_ic_defs += 1;
  BFT_REALLOC(nsp->pressure_ic_defs, nsp->n_pressure_ic_defs, cs_xdef_t *);
  nsp->pressure_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the initial condition for the pressure unknowns.
 *         This definition can be done on a specified mesh location.
 *         By default, the unknown is set to zero everywhere.
 *         Here the initial value is set according to an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" if
 *                           all cells are considered)
 * \param[in]      analytic  pointer to an analytic function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_pressure_ic_by_analytic(cs_navsto_param_t      *nsp,
                                      const char             *z_name,
                                      cs_analytic_func_t     *analytic,
                                      void                   *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  /* Add a new cs_xdef_t structure */
  int z_id = cs_get_vol_zone_id(z_name);

  cs_flag_t  meta_flag = 0;
  if (z_id == 0)
    meta_flag |= CS_FLAG_FULL_LOC;

  cs_xdef_analytic_context_t  ac = { .z_id = z_id,
                                     .func = analytic,
                                     .input = input,
                                     .free_input = NULL };

  cs_xdef_t  *d = cs_xdef_volume_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                        1,  /* dim */
                                        z_id,
                                        0,  /* state flag */
                                        meta_flag,
                                        &ac);

  /* Assign the default quadrature type of the Navier-Stokes module to this
   * definition (this can be modified by the user if the same call is
   * performed in cs_user_finalize_setup()) */
  cs_xdef_set_quadrature(d, nsp->qtype);

  int  new_id = nsp->n_pressure_ic_defs;
  nsp->n_pressure_ic_defs += 1;
  BFT_REALLOC(nsp->pressure_ic_defs, nsp->n_pressure_ic_defs, cs_xdef_t *);
  nsp->pressure_ic_defs[new_id] = d;

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to a fixed wall
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_fixed_walls(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != NULL);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_3_t  zero = {0, 0, 0};

  const cs_boundary_t  *bdy = nsp->boundaries;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (    bdy->types[i] & CS_BOUNDARY_WALL
        && !(bdy->types[i] & CS_BOUNDARY_SLIDING_WALL)) {

      /* Homogeneous Dirichlet on the velocity field. Nothing to enforce on the
         pressure field */
      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              3,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_HMG_DIRICHLET,
                                              (void *)zero);
      int  new_id = nsp->n_velocity_bc_defs;

      nsp->n_velocity_bc_defs += 1;
      BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
      nsp->velocity_bc_defs[new_id] = d;

      cs_equation_add_xdef_bc(eqp, d);

      /* Homogeneous Neumann on the pressure field --> default BC */

    } /* Fixed wall */
  } /* Loop on domain boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to a symmetry
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_symmetries(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != NULL);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_t  zero = 0;

  const cs_boundary_t  *bdy = nsp->boundaries;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (bdy->types[i] & CS_BOUNDARY_SYMMETRY) {

      /* Homogeneous Dirichlet on the normal component of the velocity field
         and homogeneous Neumann on the normal stress (balance between the
         pressure gradient and the viscous term) */
      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              1,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_SLIDING,
                                              (void *)&zero);

      cs_equation_add_xdef_bc(eqp, d);

      int  new_id = nsp->n_velocity_bc_defs;

      nsp->n_velocity_bc_defs += 1;
      BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
      nsp->velocity_bc_defs[new_id] = d;

      /* Homogeneous Neumann on the pressure field --> default BC (Nothing to
         do) */

    } /* Symmetry */
  } /* Loop on domain boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add the definition of boundary conditions related to outlets
 *         into the set of parameters for the management of the Navier-Stokes
 *         system of equations
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_set_outlets(cs_navsto_param_t    *nsp)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);
  assert(nsp->boundaries != NULL);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_real_33_t  zero = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  const cs_boundary_t  *bdy = nsp->boundaries;

  int exclude_filter = CS_BOUNDARY_IMPOSED_P | CS_BOUNDARY_IMPOSED_VEL;

  for (int i = 0; i < bdy->n_boundaries; i++) {
    if (   bdy->types[i] & CS_BOUNDARY_OUTLET
        && ! (bdy->types[i] & exclude_filter)) {

      /* Add the homogeneous Neumann on the normal component */
      cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                              9,    /* dim */
                                              bdy->zone_ids[i],
                                              CS_FLAG_STATE_UNIFORM, /* state */
                                              CS_CDO_BC_HMG_NEUMANN,
                                              (void *)&zero);
      cs_equation_add_xdef_bc(eqp, d);

      int  new_id = nsp->n_velocity_bc_defs;

      nsp->n_velocity_bc_defs += 1;
      BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
      nsp->velocity_bc_defs[new_id] = d;

      /* Homogeneous Neumann on the pressure field --> default BC.
         At the end of the day, we end up with a balance between the pressure
         gradient and the viscous term (and advection term in Navier-Stokes) */

    } /* Symmetry */
  } /* Loop on domain boundaries */
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Set the pressure field on a boundary using a uniform value.
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      value     value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_pressure_bc_by_value(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_real_t            *values)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (!(nsp->boundaries->types[bdy_id] & CS_BOUNDARY_IMPOSED_P))
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" is not related to a pressure boundary.\n"
              " Please check your settings.", __func__, z_name);

  /* Set the boundary condition for the pressure field */
  cs_xdef_t  *dp = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                           1, /* dim */
                                           z_id,
                                           CS_FLAG_STATE_UNIFORM, /* state */
                                           CS_CDO_BC_DIRICHLET,
                                           (void *)values);

  int  pnew_id = nsp->n_pressure_bc_defs;

  nsp->n_pressure_bc_defs += 1;
  BFT_REALLOC(nsp->pressure_bc_defs, nsp->n_pressure_bc_defs, cs_xdef_t *);
  nsp->pressure_bc_defs[pnew_id] = dp;

  if (!nsp->pressure_bc_is_owner) {
    bft_error(__FILE__, __LINE__, 0, "%s: Not implemented yet", __func__);
#if 0 /* TODO */
    /* Retrieve the equation related to the pressure */
    cs_equation_param_t  *p_eqp = NULL;
    assert(p_eqp != NULL);
#endif
  }

  /* Add a new cs_xdef_t structure. For the momentum equation, this is a
   * homogeneous Neumann BC for the velocity */
  cs_real_33_t  zero = {{0, 0, 0}, {0, 0, 0}, {0, 0, 0}};

  cs_xdef_t  *du = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                           9, /* dim */
                                           z_id,
                                           CS_FLAG_STATE_UNIFORM, /* state */
                                           CS_CDO_BC_HMG_NEUMANN,
                                           (void *)zero);

  int  unew_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[unew_id] = du;

  cs_equation_param_t *u_eqp = _get_momentum_param(nsp);
  assert(u_eqp != NULL);
  cs_equation_add_xdef_bc(u_eqp, du);

  return dp;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for a sliding wall boundary using a
 *         uniform value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_wall_by_value(cs_navsto_param_t    *nsp,
                                     const char           *z_name,
                                     cs_real_t            *values)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (! (nsp->boundaries->types[bdy_id] & CS_BOUNDARY_SLIDING_WALL))
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" is not related to a sliding wall boundary.\n"
              " Please check your settings.", __func__, z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          3,    /* dim */
                                          z_id,
                                          CS_FLAG_STATE_UNIFORM, /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          (void *)values);

  int  new_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  assert(eqp != NULL);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using a uniform
 *         value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      values    array of three real values
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_value(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_real_t            *values)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (!(nsp->boundaries->types[bdy_id] & CS_BOUNDARY_IMPOSED_VEL))
    bft_error
      (__FILE__, __LINE__, 0,
       " %s: Zone \"%s\" is not related to an imposed velocity boundary.\n"
       " Please check your settings.", __func__, z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_VALUE,
                                          3,    /* dim */
                                          z_id,
                                          CS_FLAG_STATE_UNIFORM, /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          (void *)values);

  int  new_id = nsp->n_velocity_bc_defs;

  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  assert(eqp != NULL);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define the velocity field for an inlet boundary using an analytical
 *         function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           boundary faces are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_set_velocity_inlet_by_analytic(cs_navsto_param_t    *nsp,
                                         const char           *z_name,
                                         cs_analytic_func_t   *ana,
                                         void                 *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  int  z_id = cs_get_bdy_zone_id(z_name);
  if (z_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not exist.\n"
              " Please check your settings.", __func__, z_name);

  int  bdy_id = cs_boundary_id_by_zone_id(nsp->boundaries, z_id);
  if (bdy_id < 0)
    bft_error(__FILE__, __LINE__, 0,
              " %s: Zone \"%s\" does not belong to an existing boundary.\n"
              " Please check your settings.", __func__, z_name);

  if (!(nsp->boundaries->types[bdy_id] & CS_BOUNDARY_IMPOSED_VEL))
    bft_error
      (__FILE__, __LINE__, 0,
       " %s: Zone \"%s\" is not related to an imposed velocity boundary.\n"
       " Please check your settings.", __func__, z_name);

  /* Add a new cs_xdef_t structure */
  cs_xdef_analytic_context_t  ac = { .z_id = z_id,
                                     .func = ana,
                                     .input = input,
                                     .free_input = NULL };

  cs_xdef_t  *d = cs_xdef_boundary_create(CS_XDEF_BY_ANALYTIC_FUNCTION,
                                          3,    /* dim */
                                          z_id,
                                          0,    /* state */
                                          CS_CDO_BC_DIRICHLET,
                                          &ac);

  /* Assign the default quadrature type of the Navier-Stokes module to this
   * definition (this can be modified by the user if the same call is
   * performed in cs_user_finalize_setup()) */
  cs_xdef_set_quadrature(d, nsp->qtype);

  int  new_id = nsp->n_velocity_bc_defs;
  nsp->n_velocity_bc_defs += 1;
  BFT_REALLOC(nsp->velocity_bc_defs, nsp->n_velocity_bc_defs, cs_xdef_t *);
  nsp->velocity_bc_defs[new_id] = d;

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_equation_add_xdef_bc(eqp, d);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an analytical function
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      ana       pointer to an analytical function
 * \param[in]      input     NULL or pointer to a structure cast on-the-fly
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_analytic(cs_navsto_param_t    *nsp,
                                      const char           *z_name,
                                      cs_analytic_func_t   *ana,
                                      void                 *input)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);
  cs_xdef_t  *d = cs_equation_add_source_term_by_analytic(eqp,
                                                          z_name, ana, input);

  /* Assign the default quadrature type of the Navier-Stokes module to this
   * definition (this can be modified by the user if the same call is
   * performed in cs_user_finalize_setup()) */
  cs_xdef_set_quadrature(d, nsp->qtype);

  return d;
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by a constant value
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      val       pointer to the value to set
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_val(cs_navsto_param_t    *nsp,
                                 const char           *z_name,
                                 cs_real_t            *val)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  return cs_equation_add_source_term_by_val(eqp, z_name, val);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Define a new source term structure defined by an array
 *
 * \param[in]      nsp       pointer to a \ref cs_navsto_param_t structure
 * \param[in]      z_name    name of the associated zone (if NULL or "" all
 *                           cells are considered)
 * \param[in]      loc       information to know where are located values
 * \param[in]      array     pointer to an array
 * \param[in]      is_owner  transfer the lifecycle to the cs_xdef_t structure
 *                           (true or false)
 * \param[in]      index     optional pointer to the array index
 *
 * \return a pointer to the new \ref cs_xdef_t structure
 */
/*----------------------------------------------------------------------------*/

cs_xdef_t *
cs_navsto_add_source_term_by_array(cs_navsto_param_t    *nsp,
                                   const char           *z_name,
                                   cs_flag_t             loc,
                                   cs_real_t            *array,
                                   bool                  is_owner,
                                   cs_lnum_t            *index)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  return cs_equation_add_source_term_by_array(eqp, z_name, loc,
                                              array, is_owner,index);
}

/*----------------------------------------------------------------------------*/
/*!
 * \brief  Add a advection field for the Oseen problem
 *
 * \param[in, out]    nsp        pointer to a \ref cs_navsto_param_t
 * \param[in, out]    adv_fld    pointer to a \ref cs_adv_field_t
 */
/*----------------------------------------------------------------------------*/

void
cs_navsto_add_oseen_field(cs_navsto_param_t   *nsp,
                          cs_adv_field_t      *adv_fld)
{
  if (nsp == NULL)
    bft_error(__FILE__, __LINE__, 0, _err_empty_nsp, __func__);

  if (nsp->model != CS_NAVSTO_MODEL_OSEEN)
    bft_error(__FILE__, __LINE__, 0, " %s: Trying to set an external advection"
                                     " where there should not be one. Stopping",
                                     __func__);

  cs_equation_param_t *eqp = _get_momentum_param(nsp);

  cs_equation_add_advection(eqp, adv_fld);
}

/*----------------------------------------------------------------------------*/

END_C_DECLS
