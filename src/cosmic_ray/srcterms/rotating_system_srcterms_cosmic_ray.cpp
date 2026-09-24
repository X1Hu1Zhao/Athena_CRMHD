//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file rotating_system_srcterms_cosmic_ray.cpp
//! \brief Adds coriolis force and centrifugal force
//======================================================================================

// C++ headers
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../cosmic_ray.hpp"
#include "cosmic_ray_srcterms.hpp"

//--------------------------------------------------------------------------------------
//! \fn void CosmicRaySourceTerms::RotatingSystemSourceTermsCosmicRay
//!             (const Real dt, const AthenaArray<Real> *cr_flux,
//!              const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons)
//! \brief source terms for the rotating system

void CosmicRaySourceTerms::RotatingSystemSourceTermsCosmicRay
                 (const Real dt, const AthenaArray<Real> *cr_flux,
                  const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons) {
  MeshBlock *pmb = cr_pointer->pmy_block;
  if(std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    // dM1/dt = 2 \rho vc vp /r +\rho vc^2/r
    // dM2/dt = -2 \rho vc vr /r
    // dE/dt  = \rho vc^2 vr /r
    // vc     = r \Omega
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v1_id   = rho_id + 1;
      int v2_id   = rho_id + 2;
      int v3_id   = rho_id + 3;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            const Real &den = cr_prim(rho_id, k, j, i);
            const Real mom1 = den*cr_prim(v1_id, k, j, i);
            const Real &ri  = pmb->pcoord->coord_src1_i_(i);
            const Real &rv  = pmb->pcoord->x1v(i);
            const Real vc   = rv*Omega_0_;
            const Real src  = SQR(vc); // (rOmega)^2
            const Real flux_c = 0.5*(cr_flux[X1DIR](rho_id, k, j, i)+cr_flux[X1DIR](rho_id, k, j, i+1));
            cr_cons(v1_id, k, j, i) += dt*ri*(2.0*vc*(den*cr_prim(v2_id, k, j, i))+den*src);
            cr_cons(v2_id, k, j, i) -= dt*ri*vc*(mom1 + flux_c);
          }
        }
      }
    }
  } else if(std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    if (pmb->block_size.nx2 == 1) {
      std::stringstream msg;
      msg << "### FATAL ERROR in "
          << "CosmicRaySourceTerms::RotatingSystemSourceTermsCosmicRay"
          << std::endl
          << "Rotating System does not support spherical polar coordinates in 1D."
          << std::endl;
      ATHENA_ERROR(msg);
    }
    // dM1/dt = 2 \rho vc vp / r
    //          +\rho (vc)^2 / r
    // dM2/dt = 2 \rho vp cot(\theta) vc / r
    //          + \rho cot(\theta) (vc)^2 /r
    // dM3/dt = -\rho vr (2 vc)/r
    //          -\rho vt (2 vc) cot(\theta) /r
    // vc     = r sin(\theta)\Omega
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v1_id   = rho_id + 1;
      int v2_id   = rho_id + 2;
      int v3_id   = rho_id + 3;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          const Real &cv1 = pmb->pcoord->coord_src1_j_(j); // cot(theta)
          const Real &cv3 = pmb->pcoord->coord_src3_j_(j); // cot(\theta)
          const Real &sv  = std::sin(pmb->pcoord->x2v(j)); // sin(\theta)
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            const Real &den  = cr_prim(rho_id, k, j, i);
            const Real &rv   = pmb->pcoord->x1v(i);
            const Real &ri   = pmb->pcoord->coord_src1_i_(i); // 1/r
            const Real vc    = rv*sv*Omega_0_;
            const Real src   = SQR(vc); // vc^2
            const Real force = den*ri*(2.0*vc*cr_prim(v3_id,k,j,i)+src);
            const Real flux_xc = 0.5*(cr_flux[X1DIR](rho_id, k, j,   i+1)+cr_flux[X1DIR](rho_id, k, j, i));
            const Real flux_yc = 0.5*(cr_flux[X2DIR](rho_id, k, j+1, i)+cr_flux[X2DIR](rho_id,   k, j, i));
            cr_cons(v1_id, k, j, i) += dt*force;
            cr_cons(v2_id, k, j, i) += dt*force*cv1;
            cr_cons(v3_id, k, j, i) -= dt*ri*vc*(den*cr_prim(v1_id, k, j, i)+flux_xc
                                           +cv3*(den*cr_prim(v2_id, k, j, i)+flux_yc));
          }
        }
      }
    }
  }
  return;
}
