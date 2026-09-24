//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_shear_cosmic_ray.cpp
//! \brief functions that apply shearing box BCs for cosmic-ray variables
//========================================================================================

// C headers

// C++ headers
#include <algorithm>  // min
#include <cmath>
#include <cstdlib>
#include <cstring>    // memcpy
#include <iomanip>
#include <iostream>   // endl
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../../../coordinates/coordinates.hpp"
#include "../../../cosmic_ray/cosmic_ray.hpp"
#include "../../../eos/eos.hpp"
#include "../../../globals.hpp"
#include "../../../mesh/mesh.hpp"
#include "../../../orbital_advection/orbital_advection.hpp"
#include "../../../parameter_input.hpp"
#include "../../../utils/buffer_utils.hpp"
#include "../../bvals.hpp"
#include "../../bvals_interfaces.hpp"

// MPI header
#ifdef MPI_PARALLEL
#include <mpi.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn void CosmicRayBoundaryVariable::AddCosmicRayShearForInit()
//! \brief Send shearing box boundary buffers for cosmic-ray variables

void CosmicRayBoundaryVariable::AddCosmicRayShearForInit() {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  AthenaArray<Real> &var = *var_cc;

  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  Real qomL = pbval_->qomL_;

  int sign[2]{1, -1};
  int ib[2]{pmb->is - NGHOST, pmb->ie + 1};

  // could call modified ShearQuantities(src=shear_cc_, dst=var, upper), by first loading
  // shear_cc_=var for rho_id, v2_id so that order of v2_id update to var doesn't matter.
  // Would need to reassign src=shear_cc_ to updated dst=var for v2_id after? Is it used?
  for (int upper=0; upper<2; upper++) {
    if (pbval_->is_shear[upper]) {
      // step 1. -- add shear to the periodic boundary values
      for (int n=0; n<NCRS; ++n) {
				int cr_id = n;
				int rho_id  = 4*cr_id;
				int v2_id   = rho_id + 2;
				for (int k=kl; k<=ku; k++) {
					for (int j=jl; j<=ju; j++) {
						for (int i=0; i<NGHOST; i++) {
							// add shear to conservative
							int ii = ib[upper] + i;
							var(v2_id, k, j, ii) += sign[upper]*qomL*var(rho_id, k, j, ii);
						}
					}
				}
			}
    }  // if boundary is shearing
  }  // loop over inner/outer boundaries
  return;
}
//----------------------------------------------------------------------------------------
//! \fn void CosmicRayBoundaryVariable::ShearQuantities(AthenaArray<Real> &shear_cc_,
//!                                                   bool upper)
//! \brief Apply shear to CosmicRay x2 momentum

void CosmicRayBoundaryVariable::ShearQuantities(AthenaArray<Real> &shear_cc_, bool upper) {
  MeshBlock *pmb = pmy_block_;
  Mesh *pmesh = pmb->pmy_mesh;
  int &xgh = pbval_->xgh_;
  int jl = pmb->js - NGHOST;
  int ju = pmb->je + NGHOST+2*xgh+1;
  int kl = pmb->ks;
  int ku = pmb->ke;
  if (pmesh->mesh_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  Real qomL = pbval_->qomL_;
  int sign[2]{1, -1};
  int ib[2]{pmb->is - NGHOST, pmb->ie + 1};

	for (int n=0; n<NCRS; ++n) {
		int cr_id = n;
		int rho_id  = 4*cr_id;
		int v2_id   = rho_id + 2;
		for (int k=kl; k<=ku; k++) {
			for (int i=0; i<NGHOST; i++) {
				for (int j=jl; j<=ju; j++) {
					shear_cc_(v2_id, k, i, j) += + sign[upper]*qomL*shear_cc_(rho_id, k, i, j);
				}
			}
		}
	}
  return;
}
