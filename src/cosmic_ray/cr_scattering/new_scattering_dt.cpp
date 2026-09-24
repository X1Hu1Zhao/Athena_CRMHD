//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file new_scattering_dt.cpp
//! \brief Calculate the CR scattering timestep constraint.

#include <algorithm>
#include <limits>
#include <sstream>

#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"

Real CRScattering::NewScatteringDt() {
  const Real real_max = std::numeric_limits<Real>::max();
  if (!Scattering_Flag) return real_max;

  Real streaming_factor = 1.0;
  if (coefficient_model_ == "streaming" && UserDefinedCRScattering == nullptr) {
    Real max_sigma0_streaming = 0.0;
    for (int CR_id=0; CR_id<NCRS; ++CR_id) {
      max_sigma0_streaming = std::max(max_sigma0_streaming, sigma0_streaming[CR_id]);
    }
    if (max_sigma0_streaming == 0.0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in CRScattering::NewScatteringDt" << std::endl
          << "scattering_model=streaming requires at least one positive "
          << "sigma0_streaming_N" << std::endl;
      ATHENA_ERROR(msg);
    }
    streaming_factor += 4.0*max_sigma0_streaming;
  }

  const bool f2 = pmb_->pmy_mesh->f2;
  const bool f3 = pmb_->pmy_mesh->f3;
  int il = pmb_->is - NGHOST;
  int jl = pmb_->js;
  int kl = pmb_->ks;
  int iu = pmb_->ie + NGHOST;
  int ju = pmb_->je;
  int ku = pmb_->ke;
  Real dt_scattering = real_max;
  AthenaArray<Real> &len = scattering_dx1_;
  AthenaArray<Real> &dx2 = scattering_dx2_;
  AthenaArray<Real> &dx3 = scattering_dx3_;

  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      pco_->CenterWidth1(k, j, il, iu, len);
      pco_->CenterWidth2(k, j, il, iu, dx2);
      pco_->CenterWidth3(k, j, il, iu, dx3);
#pragma omp simd
      for (int i=il; i<=iu; ++i) {
        len(i) = f2 ? std::min(len(i), dx2(i)) : len(i);
        len(i) = f3 ? std::min(len(i), dx3(i)) : len(i);
      }
      for (int i=il; i<=iu; ++i) {
        dt_scattering = std::min(
            dt_scattering, len(i)/streaming_factor/Vm);
      }
    }
  }
  dt_scattering *= pmb_->pmy_mesh->cfl_number;

  return dt_scattering*SQRT3;
}
