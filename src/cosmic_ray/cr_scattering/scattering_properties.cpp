//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file scattering_properties.cpp
//! \brief Set CR scattering properties.

#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"

void CRScattering::SetProperties(const Real time, const AthenaArray<Real> &w,
    const AthenaArray<Real> &cr_prim, const AthenaArray<Real> &bcc) {
  int is = pmb_->is;
  int ie = pmb_->ie;
  int js = pmb_->js;
  int je = pmb_->je;
  int ks = pmb_->ks;
  int ke = pmb_->ke;
  int il = is - NGHOST;
  int iu = ie + NGHOST;
  int jl = js;
  int ju = je;
  int kl = ks;
  int ku = ke;

  if (pmb_->block_size.nx2 > 1) {
    jl -= NGHOST;
    ju += NGHOST;
  }
  if (pmb_->block_size.nx3 > 1) {
    kl -= NGHOST;
    ku += NGHOST;
  }

  if (Scattering_Flag) {
    if (UserDefinedCRScattering == nullptr) {
      if (coefficient_model_ == "streaming"){
        const int scattering_il = is - 1;
        const int scattering_iu = ie + 1;
        const int scattering_jl = (pmb_->block_size.nx2 > 1) ? js - 1 : js;
        const int scattering_ju = (pmb_->block_size.nx2 > 1) ? je + 1 : je;
        const int scattering_kl = (pmb_->block_size.nx3 > 1) ? ks - 1 : ks;
        const int scattering_ku = (pmb_->block_size.nx3 > 1) ? ke + 1 : ke;
        SetStreamingScatteringCoefficients(
            w, cr_prim, bcc,
            cr_pointer->sigma_pL_array, cr_pointer->sigma_pR_array,
            cr_pointer->sigma_mL_array, cr_pointer->sigma_mR_array,
            cr_pointer->sigma_Lorentz_array, cr_pointer->deltaPcr,
            scattering_il, scattering_iu, scattering_jl, scattering_ju,
            scattering_kl, scattering_ku);
      } else {
        SetConstantScatteringCoefficients(
            cr_pointer->sigma_pL_array, cr_pointer->sigma_pR_array,
            cr_pointer->sigma_mL_array, cr_pointer->sigma_mR_array,
            cr_pointer->sigma_Lorentz_array, cr_pointer->deltaPcr,
            il, iu, jl, ju, kl, ku);
      }
    } else {
      UserDefinedCRScattering(cr_pointer, pmb_, w, cr_prim, bcc,
          cr_pointer->sigma_pL_array, cr_pointer->sigma_pR_array,
          cr_pointer->sigma_mL_array, cr_pointer->sigma_mR_array,
          cr_pointer->sigma_Lorentz_array, cr_pointer->deltaPcr, cs_cr_array,
          il, iu, jl, ju, kl, ku);
    }
  }

  if (!cr_pointer->PressureAnisotropyEnabled()) {
    cr_pointer->deltaPcr.ZeroClear();
  }
}
