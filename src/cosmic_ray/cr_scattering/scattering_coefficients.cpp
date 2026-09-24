//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file scattering_coefficients.cpp
//! \brief Set built-in CR scattering coefficients.

#include <cmath>

#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"

void CRScattering::SetConstantScatteringCoefficients(
    AthenaArray<Real> &sigma_pL, AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL, AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz, AthenaArray<Real> &deltaPcr,
    int is, int ie, int js, int je, int ks, int ke) {
  for (int CR_id=0; CR_id<NCRS; ++CR_id) {
    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          sigma_Lorentz(CR_id, k, j, i) = const_sigma_Lorentz[CR_id];
          sigma_mL(CR_id, k, j, i) = const_sigma_mL[CR_id];
          sigma_mR(CR_id, k, j, i) = const_sigma_mR[CR_id];
          sigma_pL(CR_id, k, j, i) = const_sigma_pL[CR_id];
          sigma_pR(CR_id, k, j, i) = const_sigma_pR[CR_id];
          deltaPcr(CR_id, k, j, i) = 0.0;
        }
      }
    }
  }
  return;
}



void CRScattering::SetStreamingScatteringCoefficients(
    const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
    const AthenaArray<Real> &bcc,
    AthenaArray<Real> &sigma_pL, AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL, AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz, AthenaArray<Real> &deltaPcr,
    int is, int ie, int js, int je, int ks, int ke) {

  const bool f2 = pmb_->block_size.nx2 > 1;
  const bool f3 = pmb_->block_size.nx3 > 1;

      Real EcrGradb = 0.0;
      Real dEcrdx = 0.0;
      Real dEcrdy = 0.0;
      Real dEcrdz = 0.0;
      Real vA = 0.0;
      Real dx, dy, dz;
      Real Dx3v, Dx2v, Dx32v;

  for (int CR_id=0; CR_id<NCRS; ++CR_id) {
    const int cr_energy_id = 4*CR_id;
    for (int k=ks; k<=ke; ++k) {
      Dx3v = pco_->dx3v(k);
      for (int j=js; j<=je; ++j) {
        Dx2v = pco_->dx2v(j);
        Dx32v = Dx3v/pco_->h32v(j);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          dx = pco_->dx1v(i);
          dy = Dx2v/pco_->h2v(i);
          dz = Dx32v/pco_->h31v(i);
          Real &scattcoeff_pL = sigma_pL(CR_id, k, j, i);
          Real &scattcoeff_pR = sigma_pR(CR_id, k, j, i);
          Real &scattcoeff_mL = sigma_mL(CR_id, k, j, i);
          Real &scattcoeff_mR = sigma_mR(CR_id, k, j, i);
          scattcoeff_pL = 0.0;
          scattcoeff_pR = 0.0;
          scattcoeff_mL = 0.0;
          scattcoeff_mR = 0.0;
          sigma_Lorentz(CR_id, k, j, i) = const_sigma_Lorentz[CR_id];

          const Real bx = bcc(IB1, k, j, i);
          const Real by = bcc(IB2, k, j, i);
          const Real bz = bcc(IB3, k, j, i);
          const Real bmag = std::sqrt(SQR(bx) + SQR(by) + SQR(bz));
          deltaPcr(CR_id, k, j, i) = 0.0;
          const Real gas_density = w(IDN, k, j, i);
          const Real Ecr = cr_prim(cr_energy_id, k, j, i);
          if (bmag <= TINY_NUMBER || gas_density <= TINY_NUMBER
              || Ecr <= TINY_NUMBER) {
            continue;
          }

          dEcrdx = (cr_prim(cr_energy_id,k,j,i+1)
                    - cr_prim(cr_energy_id,k,j,i-1))/(2.0*dx);
          dEcrdy = f2 ? (cr_prim(cr_energy_id,k,j+1,i)
                         - cr_prim(cr_energy_id,k,j-1,i))/(2.0*dy) : 0.0;
          dEcrdz = f3 ? (cr_prim(cr_energy_id,k+1,j,i)
                         - cr_prim(cr_energy_id,k-1,j,i))/(2.0*dz) : 0.0;

          EcrGradb = (bx*dEcrdx + by*dEcrdy + bz*dEcrdz)/bmag;
          vA = bmag/std::sqrt(gas_density);

          if (EcrGradb>0.0){
          scattcoeff_mL      = 0.5*sigma0_streaming[CR_id]*fabs(EcrGradb)/Ecr/vA;
          scattcoeff_mR      = scattcoeff_mL;
          }
          else if (EcrGradb<0.0){
          scattcoeff_pL      = 0.5*sigma0_streaming[CR_id]*fabs(EcrGradb)/Ecr/vA;
          scattcoeff_pR      = scattcoeff_pL;
          }
        }
      }
    }
  }
  return;
}
