//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file constant_acc_cosmic_ray.cpp
//! \brief source terms due to constant acceleration

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../cosmic_ray.hpp"
#include "cosmic_ray_srcterms.hpp"

//----------------------------------------------------------------------------------------
//! \fn void CosmicRay::ConstantAccelerationCosmicRay
//! \brief Adds source terms for constant acceleration to conserved variables

void CosmicRaySourceTerms::ConstantAccelerationCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
                                            const AthenaArray<Real> &cr_prim,
                                            AthenaArray<Real> &cr_cons) {
  MeshBlock *pmb = cr_pointer->pmy_block;

  // acceleration in 1-direction
  if (g1_!=0.0) {
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v1_id   = rho_id + 1;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real src = dt*cr_prim(rho_id, k, j, i)*g1_;
            cr_cons(v1_id, k, j, i) += src;
          }
        }
      }
    }
  }

  // acceleration in 2-direction
  if (g2_!=0.0) {
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v2_id   = rho_id + 2;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real src = dt*cr_prim(rho_id, k, j, i)*g2_;
            cr_cons(v2_id, k, j, i) += src;
          }
        }
      }
    }
  }

  // acceleration in 3-direction
  if (g3_!=0.0) {
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v3_id   = rho_id + 3;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            Real src = dt*cr_prim(rho_id, k, j, i)*g3_;
            cr_cons(v3_id, k, j, i) += src;
          }
        }
      }
    }
  }

  return;
}
