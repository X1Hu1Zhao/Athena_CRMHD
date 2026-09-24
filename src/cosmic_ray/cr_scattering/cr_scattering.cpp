//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_scattering.cpp
//! Contains data and functions that implement physical CR scattering and gas-feedback terms

// C++ headers
#include <algorithm>   // min,max
#include <iostream>
#include <limits>
#include <cstring>    // strcmp

// Athena++ headers
#include "../../defs.hpp"
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../hydro/hydro.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"


CRScattering::CRScattering(CosmicRay *cr_pointer, ParameterInput *pin) :
  cr_pointer(cr_pointer), pmb_(cr_pointer->pmy_block),
  pco_(pmb_->pcoord), pmy_hydro_(pmb_->phydro),
  orbital_advection_order_(pmb_->pmy_mesh->orbital_advection) {

  int nc1 = pmb_->ncells1;
  int nc2 = pmb_->ncells2;
  int nc3 = pmb_->ncells3;

  Vm      = cr_pointer->Vm;
  Vm2     = Vm*Vm;
  inv_Vm  = 1.0/Vm;
  inv_Vm2 = inv_Vm*inv_Vm;
  gamma_cr = cr_pointer->gamma_cr;


  force.NewAthenaArray(NCRVARS, nc1);
  delta_mom.NewAthenaArray(NCRVARS, nc1);

  force_Ecr_n.NewAthenaArray(NSPECIES, nc1);
  force_x1_n.NewAthenaArray(NSPECIES, nc1);
  force_x2_n.NewAthenaArray(NSPECIES, nc1);
  force_x3_n.NewAthenaArray(NSPECIES, nc1);

  sigma_p.NewAthenaArray(NCRS, nc1);
  sigma_m.NewAthenaArray(NCRS, nc1);
  Omega.NewAthenaArray(NCRS, nc1);

  sigma_p_n.NewAthenaArray(NCRS, nc1);
  sigma_m_n.NewAthenaArray(NCRS, nc1);
  Omega_n.NewAthenaArray(NCRS, nc1); 

  cosXY.NewAthenaArray(nc1);
  sinXY.NewAthenaArray(nc1);
  cosZ.NewAthenaArray(nc1);
  sinZ.NewAthenaArray(nc1);
  vA.NewAthenaArray(nc1);
  fwdwvspeed.NewAthenaArray(nc1);
  bwdwvspeed.NewAthenaArray(nc1);  


  inv_gas_rho.NewAthenaArray(nc1);
  dens_prim.NewAthenaArray(NSPECIES,nc1);
  dens_prim_n.NewAthenaArray(NSPECIES,nc1);

  inv_gas_rho_n.NewAthenaArray(nc1);

  delta_Ecr.NewAthenaArray(NSPECIES, nc1);
  delta_mom1.NewAthenaArray(NSPECIES, nc1);
  delta_mom2.NewAthenaArray(NSPECIES, nc1);
  delta_mom3.NewAthenaArray(NSPECIES, nc1);

  delta_Ecr_src.NewAthenaArray(NSPECIES, nc1);
  delta_mom1_src.NewAthenaArray(NSPECIES, nc1);
  delta_mom2_src.NewAthenaArray(NSPECIES, nc1);
  delta_mom3_src.NewAthenaArray(NSPECIES, nc1);

  mom1_prim.NewAthenaArray(NSPECIES, nc1);
  mom2_prim.NewAthenaArray(NSPECIES, nc1);
  mom3_prim.NewAthenaArray(NSPECIES, nc1);

  mom1_prim_n.NewAthenaArray(NSPECIES, nc1);
  mom2_prim_n.NewAthenaArray(NSPECIES, nc1);
  mom3_prim_n.NewAthenaArray(NSPECIES, nc1);

  jacobi.NewAthenaArray(NCRVARS, NCRVARS, nc1);
  jacobi_n.NewAthenaArray(NCRVARS, NCRVARS, nc1);
  product.NewAthenaArray(NCRVARS, NCRVARS, nc1);

  lambda.NewAthenaArray(NCRVARS, NCRVARS, nc1);
  lambda_inv.NewAthenaArray(NCRVARS, NCRVARS, nc1);

  biggest_arr.NewAthenaArray(nc1);
  det_arr.NewAthenaArray(nc1);
  mmax_arr.NewAthenaArray(nc1);
  sum_arr.NewAthenaArray(nc1);
  temp_arr.NewAthenaArray(nc1);
  idx_vector.NewAthenaArray(NCRVARS, nc1);
  scale_arr.NewAthenaArray(NCRVARS, nc1);
  xx_arr.NewAthenaArray(NCRVARS, nc1);
  lu_matrix.NewAthenaArray(NCRVARS, NCRVARS, nc1);
  temp_A.NewAthenaArray(NCRVARS, NCRVARS, nc1);
  temp_B.NewAthenaArray(NCRVARS, NCRVARS, nc1);
  temp_C.NewAthenaArray(NCRVARS, NCRVARS, nc1);

  scattering_dx1_.NewAthenaArray(nc1);
  scattering_dx2_.NewAthenaArray(nc1);
  scattering_dx3_.NewAthenaArray(nc1);

  cs_cr_array.NewAthenaArray(NCRS, nc3, nc2, nc1);
  cs_cr_array_n.NewAthenaArray(NCRS, nc3, nc2, nc1);

  Scattering_Flag =
      pin->GetOrAddBoolean("cosmic_ray", "Scattering_Flag", true);
  if (!Scattering_Flag) {
    std::stringstream msg;
    msg << "### FATAL ERROR in CRScattering constructor" << std::endl
        << "Scattering_Flag must be true when the CR module is enabled."
        << std::endl;
    ATHENA_ERROR(msg);
  }
  UserDefinedCRScattering = pmb_->pmy_mesh->CRScatteringCoeff_;
  sigma0_streaming.fill(0.0);

  coefficient_model_ =
      pin->GetOrAddString("cosmic_ray", "scattering_model", "constant");
  if (Scattering_Flag && UserDefinedCRScattering == nullptr) {
    if (coefficient_model_ == "constant") {
      for (int n=0; n<NCRS; ++n) {
        const_sigma_pL[n] =
            pin->GetReal("cosmic_ray", "sigma_pL_" + std::to_string(n+1));
        const_sigma_pR[n] =
            pin->GetReal("cosmic_ray", "sigma_pR_" + std::to_string(n+1));
        const_sigma_mL[n] =
            pin->GetReal("cosmic_ray", "sigma_mL_" + std::to_string(n+1));
        const_sigma_mR[n] =
            pin->GetReal("cosmic_ray", "sigma_mR_" + std::to_string(n+1));
        const_sigma_Lorentz[n] = pin->GetOrAddReal(
            "cosmic_ray", "sigma_Lorentz_" + std::to_string(n+1),
            1.0e8);
      }
    } else if (coefficient_model_ == "streaming") {
      if (!MAGNETIC_FIELDS_ENABLED) {
        std::stringstream msg;
        msg << "### FATAL ERROR in CRScattering constructor" << std::endl
            << "scattering_model=streaming requires magnetic fields" << std::endl;
        ATHENA_ERROR(msg);
      }
      Real max_sigma0_streaming = 0.0;
      for (int n=0; n<NCRS; ++n) {
        sigma0_streaming[n] = pin->GetOrAddReal(
            "cosmic_ray", "sigma0_streaming_" + std::to_string(n+1),
            0.0);
        const_sigma_Lorentz[n] = pin->GetOrAddReal(
            "cosmic_ray", "sigma_Lorentz_" + std::to_string(n+1), 0.0);
        if (sigma0_streaming[n] < 0.0) {
          std::stringstream msg;
          msg << "### FATAL ERROR in CRScattering constructor" << std::endl
              << "sigma0_streaming_" << n+1 << " must be non-negative" << std::endl;
          ATHENA_ERROR(msg);
        }
        max_sigma0_streaming = std::max(max_sigma0_streaming, sigma0_streaming[n]);
      }
      if (max_sigma0_streaming == 0.0) {
        std::stringstream msg;
        msg << "### FATAL ERROR in CRScattering constructor" << std::endl
            << "scattering_model=streaming requires at least one positive "
            << "sigma0_streaming_N" << std::endl;
        ATHENA_ERROR(msg);
      }
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in CRScattering constructor" << std::endl
          << "scattering_model must be constant or streaming, but is "
          << coefficient_model_ << std::endl;
      ATHENA_ERROR(msg);
    }
  }
  CRFeedback_Flag = pin->GetBoolean("cosmic_ray", "CRFeedback_Flag");
  CR_EnergyFeedback_Flag =
      pin->GetOrAddBoolean("cosmic_ray", "CR_EnergyFeedback_Flag", true);
  if (!CRFeedback_Flag && CR_EnergyFeedback_Flag) {
    if (pmb_->gid == 0) {
      std::cout << "### WARNING in CRScattering constructor" << std::endl
                << "CR_EnergyFeedback_Flag is ignored when "
                << "CRFeedback_Flag=false; selecting no-feedback mode."
                << std::endl;
    }
    CR_EnergyFeedback_Flag = false;
  }
  integrator_ = pin->GetOrAddString("time", "integrator", "vl2");
  scattering_integrator_ =
      pin->GetOrAddString("cosmic_ray", "scattering_integrator", "2nd-implicit");

  if (scattering_integrator_ == "2nd-implicit") {
    method_id_ = 1;
  } else if (scattering_integrator_ == "1st-implicit") {
    method_id_ = 2;
  } else {
    std::stringstream msg;
    msg << "scattering_integrator must be 2nd-implicit or 1st-implicit."
        << std::endl;
    ATHENA_ERROR(msg);
  }
}


void CRScattering::ScatteringIntegrator(const int stage, const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim, const AthenaArray<Real> &bcc,
      AthenaArray<Real> &u, AthenaArray<Real> &cr_cons) {

  if (CRFeedback_Flag && !CR_EnergyFeedback_Flag
      && !((method_id_ == 1 && integrator_ == "vl2")
           || (method_id_ == 2 && integrator_ == "rk1"))) {
    std::stringstream msg;
    msg << "CR_EnergyFeedback_Flag=false is supported only by "
        << "VL2 2nd-implicit and RK1 1st-implicit scattering" << std::endl;
    ATHENA_ERROR(msg);
  }

  switch (method_id_) {
    case 1:
      if (integrator_ == "vl2") {
        if (CRFeedback_Flag) {
          if (CR_EnergyFeedback_Flag) {
            VL2ImplicitFeedback(
                stage, dt, w, cr_prim, bcc, u, cr_cons);
          } else {
            VL2ImplicitNoEnergyFeedback(
                stage, dt, w, cr_prim, bcc, u, cr_cons);
          }
        } else {
          VL2ImplicitNoFeedback(
              stage, dt, w, cr_prim, bcc, u, cr_cons);
        }
      } else {
        std::stringstream msg;
        msg << "The integrator combined with 2nd-implicit must be \"vl2\"."
            << std::endl;
        ATHENA_ERROR(msg);
      }
      break;

    case 2:
      if (integrator_ == "rk1") {
        if (CRFeedback_Flag) {
          if (CR_EnergyFeedback_Flag) {
            BackwardEulerFeedback(
                stage, dt, w, cr_prim, bcc, u, cr_cons);
          } else {
            BackwardEulerNoEnergyFeedback(
                stage, dt, w, cr_prim, bcc, u, cr_cons);
          }
        } else {
          BackwardEulerNoFeedback(
              stage, dt, w, cr_prim, bcc, u, cr_cons);
        }
      } else {
        std::stringstream msg;
        msg << "The integrator combined with 1st-implicit must be \"rk1\"."
            << std::endl;
        ATHENA_ERROR(msg);
      }
      break;

    default:
      std::stringstream msg;
      msg << "scattering_integrator must be 2nd-implicit or 1st-implicit"
          << std::endl;
      ATHENA_ERROR(msg);
      break;
  }
  return;
}
