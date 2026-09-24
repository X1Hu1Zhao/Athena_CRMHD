//======================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//======================================================================================
//! \file orbital_advection_srcterms_cosmic_ray.cpp
//! \brief Add source terms for orbital advection
//======================================================================================

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../cosmic_ray.hpp"

// this class header
#include "cosmic_ray_srcterms.hpp"

//--------------------------------------------------------------------------------------
//! \fn void CosmicRaySourceTerms::OrbitalAdvectionSourceTerms
//                          (const Real dt, const AthenaArray<Real> *cr_flux,
//                           const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons)
//  \brief Add source terms for orbital advection
void CosmicRaySourceTerms::OrbitalAdvectionSourceTermsCosmicRay
                       (const Real dt, const AthenaArray<Real> *cr_flux,
                        const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons) {
  MeshBlock *pmb = cr_pointer->pmy_block;
  OrbitalAdvection  *porb = pmb->porb;
  AthenaArray<Real> &vKc  = porb->vKc;
  AthenaArray<Real> &dvKc1 = porb->dvKc1;
  AthenaArray<Real> &dvKc2 = porb->dvKc2;
  if (std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) { // default orbital velocity
      // dM1/dt = 2\Omega \rho vy
      // dM2/dt = -\rho vx (2-q)\Omega
      // dE/dt  = q\Omega0 (\rho vy vx -b1 b2)
      //           + 2\rho\Omega vx (vk+q\Omega x)
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
              const Real &den       = cr_prim(rho_id,k,j,i);
              const Real mom1       = den*cr_prim(v1_id,k,j,i);
              const Real &vy        = cr_prim(v2_id,k,j,i);
              cr_cons(v1_id,k,j,i) += dt*Omega_0_*2.0*(den*vy);
              cr_cons(v2_id,k,j,i) -= dt*Omega_0_*(2.0-qshear_)*mom1;
            }
          }
        }
      }
    } else { // user-defined orbital velocity
      // dM1/dt = 2\Omega \rho vy + 2\Omega \rho vc
      // dM2/dt = -\rho vx (2\Omega+dvk/dx) - \rho vz dvk/dz
      // dE/dt  = (b1 b2 -\rho vy vx) dvk/dx
      //           + (b2 b3 - \rho vy vz) dvk/dz
      //           + 2\rho\Omega vx vc
      // vc     = vk+q\Omega x
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
              const Real &den       = cr_prim(rho_id,k,j,i);
              const Real &dvk_i     = dvKc1(k,i);
              const Real mom1       = den*cr_prim(v1_id,k,j,i);
              const Real &vy        = cr_prim(v2_id,k,j,i);
              const Real vc         = vKc(k,i)+qshear_*Omega_0_*pmb->pcoord->x1v(i);
              // 2D components
              cr_cons(v1_id,k,j,i) += dt*2.0*Omega_0_*(den*(vy+vc));
              cr_cons(v2_id,k,j,i) -= dt*(2.0*Omega_0_+dvk_i)*mom1;
              // 3D components
              if (pmb->block_size.nx3 > 1) {
                const Real &dvk_k  = dvKc2(k,i);
                if (dvk_k != 0.0) {
                  const Real mom3 = den*cr_prim(v3_id,k,j,i);
                  cr_cons(v2_id,k,j,i) -= dt*dvk_k*mom3;
                }
              }
            }
          }
        }
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0) {
    if (pmb->pmy_mesh->OrbitalVelocity_ == nullptr) { // default orbital velocity
      if (pmb->block_size.nx3 > 1) { //3D
        // dM1/dt = 2\rho vp vc/r
        // dM2/dt = -\rho vr (2 vc /r+(dvk/dr-vk/r))-\rho vz dvk/dz
        // dM3/dt = -\rho GM z /d^3
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        //          +(dvk/dz)(b2 b3 - \rho vp vz) -GM z /d^3 \rho vz
        // vc     = vk + r\Omega
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
                const Real &den       = cr_prim(rho_id,k,j,i);
                const Real &rv        = pmb->pcoord->x1v(i);
                const Real &ri        = pmb->pcoord->coord_src1_i_(i); // 1/r
                const Real &zv        = pmb->pcoord->x3v(k);
                const Real dv         = std::sqrt(SQR(rv)+SQR(zv));    // d = sqrt(r^2+z^2)
                const Real vc         = std::sqrt(gm_/dv)*rv/dv;       // vk + r\Omega
                const Real rdok_i     = -1.5*vc*rv/SQR(dv);            // dvk/dr-vk/r
                const Real dvk_k      = -1.5*vc*zv/SQR(dv);            // dvk/dz
                const Real mom1       = den*cr_prim(v1_id,k,j,i);
                const Real &vp        = cr_prim(v2_id,k,j,i);
                const Real mom3       = den*cr_prim(v3_id,k,j,i);
                const Real flux_xc    = 0.5*(cr_flux[X1DIR](rho_id,k,j,i+1)+cr_flux[X1DIR](rho_id,k,j,i));
                const Real flux_zc    = 0.5*(cr_flux[X3DIR](rho_id,k+1,j,i)+cr_flux[X3DIR](rho_id,k,j,i));
                cr_cons(v1_id,k,j,i) += 2.0*dt*vc*ri*(den*vp);
                cr_cons(v2_id,k,j,i) -= 0.5*dt*((2.0*vc*ri+rdok_i)*(mom1+flux_xc)
                                            +dvk_k*(mom3+flux_zc));
                const Real src_k      = -gm_*zv/(dv*dv*dv);
                cr_cons(v3_id,k,j,i) += dt*src_k*den;
              }
            }
          }
        }
      } else { // 2D
        // dM1/dt = 2\rho vp vc/r
        // dM2/dt = -\rho vr (2 vc /r+(dvk/dr-vk/r))
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        // vc     = vk+r\Omega
        int ks = pmb->ks;
        for (int n=0; n<NCRS; ++n) {
          int cr_id = n;
          int rho_id  = 4*cr_id;
          int v1_id   = rho_id + 1;
          int v2_id   = rho_id + 2;
          int v3_id   = rho_id + 3;
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              const Real &den        = cr_prim(rho_id,ks,j,i);
              const Real &rv         = pmb->pcoord->x1v(i);
              const Real &ri         = pmb->pcoord->coord_src1_i_(i);
              const Real vc          = std::sqrt(gm_/rv);
              const Real rdok_i      = -1.5*vc/rv;
              const Real mom1        = den*cr_prim(v1_id,ks,j,i);
              const Real &vp         = cr_prim(v2_id,ks,j,i);
              const Real flux_c      = 0.5*(cr_flux[X1DIR](rho_id,ks,j,i+1)+cr_flux[X1DIR](rho_id,ks,j,i));
              cr_cons(v1_id,ks,j,i) += 2.0*dt*vc*ri*(den*vp);
              cr_cons(v2_id,ks,j,i) -= dt*(vc*ri+0.5*rdok_i)*(mom1+flux_c);
            }
          }
        }
      }
    } else { // user-defined orbital velocity
      if (pmb->block_size.nx3== 1) { // 2D
        // dM1/dt = 2\rho vp vc/r
        //          + \rho (vc^2 - gm/r)/r
        // dM2/dt = -\rho vr (2 vc /r+(dvk/dr-vk/r))
        // dE/dt  = (dvk/dr-vk/r)(b1 b2 - \rho vr vp)
        //          +\rho vr (vc^2 - gm/r)/r
        // vc     = vk+r\Omega
        int ks = pmb->ks;
        for (int n=0; n<NCRS; ++n) {
          int cr_id = n;
          int rho_id  = 4*cr_id;
          int v1_id   = rho_id + 1;
          int v2_id   = rho_id + 2;
          int v3_id   = rho_id + 3;
          for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
            for (int i=pmb->is; i<=pmb->ie; ++i) {
              const Real &den        = cr_prim(rho_id,ks,j,i);
              const Real &rv         = pmb->pcoord->x1v(i);
              const Real &ri         = pmb->pcoord->coord_src1_i_(i);
              const Real &vk         = vKc(ks,i);
              const Real vc          = vk + Omega_0_*rv;
              const Real rdok_i      = dvKc1(ks,i)-vk*ri;
              const Real mom1        = den*cr_prim(v1_id,ks,j,i);
              const Real &vp         = cr_prim(v2_id,ks,j,i);
              const Real src         = SQR(vc)-gm_/rv; // vc^2-gm/r
              const Real flux_c      = 0.5*(cr_flux[X1DIR](rho_id,ks,j,i+1)+cr_flux[X1DIR](rho_id,ks,j,i));
              cr_cons(v1_id,ks,j,i) += dt*ri*(2.0*vc*(den*vp)+src*den);
              cr_cons(v2_id,ks,j,i) -= dt*(vc*ri+0.5*rdok_i)*(mom1+flux_c);
            }
          }
        }
      } else { // 3D
        std::cout << "[OribtalAdvectionSourceTerm]: not compatible to 3D with "
                  << "user-defined vK in cylindrical coordinates ." << std::endl;
        return;
      }
    }
  } else if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) {
    // dM1/dt = 2\rho vp vc/r
    //          + \rho (vc^2/r - gm/r^2)
    // dM2/dt = 2\rho vp vc/r cot(\theta)
    //          + \rho/r vc^2 cot(\theta)
    // dM3/dt = -\rho vr (2 vc /r + rdok_i)
    //          -\rho vt cot(\theta)(2 vc /r + sdvk_j)
    // dE/dt  = rdok_i(b1b3-\rho vr vp)
    //          +sdvk_j (b2 b3 - \rho vt vp)
    //          + \rho vr (vc^2/r - gm/r^2)
    //          + \rho vt/r cot(\theta) vc^2
    // vc     = vk + r sin(\theta) \Omega
    // rdok_i = dvk/dr-vk/r
    // sdvk_j = (dvk/dtheta - vk cot(\theta))/r
    for (int n=0; n<NCRS; ++n) {
      int cr_id = n;
      int rho_id  = 4*cr_id;
      int v1_id   = rho_id + 1;
      int v2_id   = rho_id + 2;
      int v3_id   = rho_id + 3;
      for (int k=pmb->ks; k<=pmb->ke; ++k) {
        for (int j=pmb->js; j<=pmb->je; ++j) {
          const Real &cv1 = pmb->pcoord->coord_src1_j_(j);   // cot(\theta)
          const Real &cv3 = pmb->pcoord->coord_src3_j_(j);   // cot(\theta)
          const Real sv   = std::sin(pmb->pcoord->x2v(j));   // sin(\theta)
#pragma omp simd
          for (int i=pmb->is; i<=pmb->ie; ++i) {
            const Real &den       = cr_prim(rho_id,k,j,i);
            const Real &rv        = pmb->pcoord->x1v(i);
            const Real &ri        = pmb->pcoord->coord_src1_i_(i); // 1/r
            const Real &vk        = vKc(j,i);                      // vk
            const Real vc         = vk + Omega_0_*rv*sv;
            const Real &dvk1      = dvKc1(j,i);
            const Real &dvk2      = dvKc2(j,i);
            const Real mom1       = den*cr_prim(v1_id,k,j,i);
            const Real mom2       = den*cr_prim(v2_id,k,j,i);
            const Real &vp        = cr_prim(v3_id,k,j,i);
            const Real src        = SQR(vc);
            const Real grav       = gm_/rv;
            const Real flux_xc    = 0.5*(cr_flux[X1DIR](rho_id,k,j,i+1)+cr_flux[X1DIR](rho_id,k,j,i));
            const Real flux_yc    = 0.5*(cr_flux[X2DIR](rho_id,k,j+1,i)+cr_flux[X2DIR](rho_id,k,j,i));
            cr_cons(v1_id,k,j,i) += dt*ri*den*(2.0*vc*vp+(src-grav));
            cr_cons(v2_id,k,j,i) += dt*ri*den*cv1*(2.0*vc*vp+src);
            cr_cons(v3_id,k,j,i) -= 0.5*dt*(((vc+Omega_0_*rv*sv)*ri+dvk1)*(mom1+flux_xc)
                                      +((vc+Omega_0_*rv*sv)*cv3+dvk2)*(mom2+flux_yc)*ri);
          }
        }
      }
    }
  }
  return;
}
