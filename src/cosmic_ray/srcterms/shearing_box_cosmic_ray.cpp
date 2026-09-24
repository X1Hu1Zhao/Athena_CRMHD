//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file shearing_box_cosmic_ray.cpp
//! \brief Adds source terms due to local shearing box approximation
//======================================================================================

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

//--------------------------------------------------------------------------------------
//! \fn void CosmicRaySourceTerms::ShearingBoxSourceTermsCosmicRay(const Real dt,
//!   const AthenaArray<Real> *flux, const AthenaArray<Real> &cr_prim,
//!   AthenaArray<Real> &cr_cons)
//! \brief Shearing Box source terms
//!
//! We add shearing box source term via operator splitting method. The source terms are
//! added after the fluxes are computed in each step of the integration (in
//! FluxDivergence) to give predictions of the conservative variables for either the
//! next step or the final update.

void CosmicRaySourceTerms::ShearingBoxSourceTermsCosmicRay(const Real dt,
            const AthenaArray<Real> *cr_flux, const AthenaArray<Real> &cr_prim,
                                              AthenaArray<Real> &cr_cons) {

  MeshBlock *pmb = cr_pointer->pmy_block;

  // 1) Tidal force:
  //    dM1/dt = 2q\rho\Omega^2 x
  // 2) Coriolis forces:
  //    dM1/dt = 2\Omega(\rho v_y)
  //    dM2/dt = -2\Omega(\rho v_x)
  if (ShBoxCoord_== 1) {
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v1_id   = rho_id + 1;
      int v2_id   = rho_id + 2;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            const Real &den = cr_prim(rho_id, k, j, i);
            const Real qO2  = qshear_*SQR(Omega_0_);
            const Real mom1 = den*cr_prim(v1_id, k, j, i);
            const Real &xc  = pmb->pcoord->x1v(i);
            cr_cons(v1_id, k, j, i) += 2.0*dt*(Omega_0_*(den*cr_prim(v2_id, k, j, i))+qO2*den*xc);
            cr_cons(v2_id, k, j, i) -= 2.0*dt*Omega_0_*mom1;
          }
        }
      }
    }
  } else { // ShBoxCoord_== 2
    int ks = pmb->ks;
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v1_id   = rho_id + 1;
      int v3_id   = rho_id + 3;
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {
        const Real &den = cr_prim(rho_id, ks, j, i);
        const Real qO2  = qshear_*SQR(Omega_0_);
        const Real mom1 = den*cr_prim(v1_id, ks, j, i);
        const Real &xc  = pmb->pcoord->x1v(i);
        cr_cons(v1_id, ks, j, i) += 2.0*dt*(Omega_0_*(den*cr_prim(v3_id, ks, j, i))+qO2*den*xc);
        cr_cons(v3_id, ks, j, i) -= 2.0*dt*Omega_0_*mom1;
        }
      }
    }
  }
  return;
}
