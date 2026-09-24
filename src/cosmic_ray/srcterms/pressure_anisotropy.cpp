//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file pressure_anisotropy.cpp
//! \brief Adds source terms due to CR pressure anisotropy

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
//! \fn void CosmicRaySourceTerms::PressureAnisotropySource
//! \brief Adds explicit source terms corresponding to the CR pressure anisotropy

void CosmicRaySourceTerms::PressureAnisotropySource(const Real dt,
    const AthenaArray<Real> *cr_flux, const AthenaArray<Real> &cr_prim,
    AthenaArray<Real> &cr_cons, const AthenaArray<Real> &deltaPcr,
    const AthenaArray<Real> &sigma_pL_array,
    const AthenaArray<Real> &sigma_pR_array,
    const AthenaArray<Real> &sigma_mL_array,
    const AthenaArray<Real> &sigma_mR_array,
    const AthenaArray<Real> &bcc) {
  MeshBlock  *pmb  = cr_pointer->pmy_block;
  Real Vm   = cr_pointer->Vm;
  Real Vm2  = Vm*Vm;
  Real dBxdx, dBxdy, dBxdz;
  Real dBydx, dBydy, dBydz;
  Real dBzdx, dBzdy, dBzdz;    
  Real Bx, By, Bz, B2, B;
  Real Vcrx, Vcry, Vcrz;


  for (int n=0; n<NCRS; ++n) {
    int CR_id = n;
    int cr_energy_id = 4*CR_id;
    int cr_flux1_id = cr_energy_id + 1;
    int cr_flux2_id = cr_energy_id + 2;
    int cr_flux3_id = cr_energy_id + 3;
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
#pragma omp simd
        for (int i=pmb->is; i<=pmb->ie; ++i) {

          Real anisoP   = deltaPcr(CR_id,k,j,i);

          Real sigma_pL = sigma_pL_array(CR_id,k,j,i);
          Real sigma_pR = sigma_pR_array(CR_id,k,j,i);
          Real sigma_mL = sigma_mL_array(CR_id,k,j,i);
          Real sigma_mR = sigma_mR_array(CR_id,k,j,i);
               

          Bx   = bcc(IB1,k,j,i);
          By   = bcc(IB2,k,j,i);
          Bz   = bcc(IB3,k,j,i);                    
          B2   = SQR(Bx) + SQR(By) + SQR(Bz);
          B    = sqrt(B2);
          const Real gas_density = pmb->phydro->w(IDN,k,j,i);
          if (B <= TINY_NUMBER || gas_density <= TINY_NUMBER) continue;
          Real vA = B/std::sqrt(gas_density);


          Real Ecr_src    = vA*Vm*(sigma_pL - sigma_pR - sigma_mL + sigma_mR)*anisoP;
          Real Fcr1_src   = Vm2*(sigma_pL - sigma_pR + sigma_mL - sigma_mR)*anisoP*Bx/B;
          Real Fcr2_src   = Vm2*(sigma_pL - sigma_pR + sigma_mL - sigma_mR)*anisoP*By/B;          
          Real Fcr3_src   = Vm2*(sigma_pL - sigma_pR + sigma_mL - sigma_mR)*anisoP*Bz/B;


          cr_cons(cr_energy_id, k, j, i) += Ecr_src*dt;
          cr_cons(cr_flux1_id, k, j, i)  += Fcr1_src*dt;
          cr_cons(cr_flux2_id, k, j, i)  += Fcr2_src*dt;
          cr_cons(cr_flux3_id, k, j, i)  += Fcr3_src*dt;
        }
      }
    }
  }
  return;
}
