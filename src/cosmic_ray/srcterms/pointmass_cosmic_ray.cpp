//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file pointmass_cosmic_ray.cpp
//! \brief Adds source terms due to point mass AT ORIGIN

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../cosmic_ray.hpp"
#include "cosmic_ray_srcterms.hpp"

class CosmicRay;
class ParameterInput;

//----------------------------------------------------------------------------------------
//! \fn void CosmicRaySourceTerms::PointMass_CosmicRay
//! \brief Adds source terms due to point mass AT ORIGIN

void CosmicRaySourceTerms::PointMassCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
                                 const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons) {
  MeshBlock *pmb = cr_pointer->pmy_block;
  for (int n=0; n<NCRS; ++n) {
    int cr_id = n;
    int rho_id  = 4*cr_id;
    int v1_id   = rho_id + 1;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real &coord_x1           = pmb->pcoord->x1v(i);
          Real &coord_src1         = pmb->pcoord->coord_src1_i_(i);
          const Real &cr_energy    = cr_prim(rho_id, k, j, i);
          Real src                 = dt*cr_energy*coord_src1*gm_/coord_x1;
          cr_cons(v1_id, k, j, i) -= src;
        }
      }
    }
  }
  return;
}
