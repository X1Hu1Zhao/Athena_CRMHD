//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cosmic_ray_HLLE_solver.cpp
//! \brief spatially isothermal HLLE Riemann solver for CR
//!
//! Computes 1D fluxes using the Harten-Lax-van Leer (HLL) Riemann solver.  This flux is
//! very diffusive, especially for contacts, and so it is not recommended for use in
//! applications.  However, as shown by Einfeldt et al.(1991), it is positively
//! conservative (cannot return negative densities or pressure), so it is a useful
//! option when other approximate solvers fail and/or when extra dissipation is needed.
//!
//!REFERENCES:
//!- E.F. Toro, "Riemann Solvers and numerical methods for fluid dynamics", 2nd ed.,
//!  Springer-Verlag, Berlin, (1999) chpt. 10.
//!- Einfeldt et al., "On Godunov-type methods near low densities", JCP, 92, 273 (1991)
//!- A. Harten, P. D. Lax and B. van Leer, "On upstream differencing and Godunov-type
//!  schemes for hyperbolic conservation laws", SIAM Review 25, 35-61 (1983).

// C headers

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../eos/eos.hpp"
#include "cr_frame_transform.hpp"
#include "cosmic_ray.hpp"
Real tanh_kernel(Real x, Real b, Real w);

//----------------------------------------------------------------------------------------
//! \fn void CosmicRay::HLLENoCsRiemannSolver_CosmicRay
//! \brief The HLLE Riemann solver for CR (No Sound Speed)

void CosmicRay::HLLENoCsRiemannSolverCosmicRay(const int k, const int j, const int il, const int iu,
                          const int index, const AthenaArray<Real> &b,
                          AthenaArray<Real> &w, AthenaArray<Real> &cr_prim_l, AthenaArray<Real> &cr_prim_r, 
                          AthenaArray<Real> &b_interface_l, AthenaArray<Real> &b_interface_r,
                          AthenaArray<Real> &deltaPcr_l, AthenaArray<Real> &deltaPcr_r,
                          AthenaArray<Real> &wl_,  AthenaArray<Real> &wr_, AthenaArray<Real> &cr_flux) {
  const int ivx = IVX + (index-IVX)%3;
  const int ivy = IVX + ((index-IVX)+1)%3;
  const int ivz = IVX + ((index-IVX)+2)%3;
  const int left_k = k - (index == IVZ ? 1 : 0);
  const int left_j = j - (index == IVY ? 1 : 0);
  const int left_i_offset = index == IVX ? -1 : 0;
  const Real gamma_gas = pmy_block->peos->GetGamma();
  const Real domain_length = index == IVX
      ? pmy_block->pmy_mesh->mesh_size.x1max - pmy_block->pmy_mesh->mesh_size.x1min
      : (index == IVY
         ? pmy_block->pmy_mesh->mesh_size.x2max - pmy_block->pmy_mesh->mesh_size.x2min
         : pmy_block->pmy_mesh->mesh_size.x3max - pmy_block->pmy_mesh->mesh_size.x3min);

  // Compute gas and magnetic quantities shared by all CR species once per interface.
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    const int left_i = i + left_i_offset;
    const Real rho_l = w(IDN, left_k, left_j, left_i);
    const Real rho_r = w(IDN, k, j, i);
    const Real Bx_l = b_interface_l(IB1, i);
    const Real By_l = b_interface_l(IB2, i);
    const Real Bz_l = b_interface_l(IB3, i);
    const Real Bx_r = b_interface_r(IB1, i);
    const Real By_r = b_interface_r(IB2, i);
    const Real Bz_r = b_interface_r(IB3, i);
    const Real B_l = std::sqrt(SQR(Bx_l) + SQR(By_l) + SQR(Bz_l));
    const Real B_r = std::sqrt(SQR(Bx_r) + SQR(By_r) + SQR(Bz_r));

    Real V_bx_l = 0.0;
    Real V_by_l = 0.0;
    Real V_bz_l = 0.0;
    Real V_bx_r = 0.0;
    Real V_by_r = 0.0;
    Real V_bz_r = 0.0;
    if (B_l > TINY_NUMBER) {
      const cr_frame::MagneticFrame frame_l =
          cr_frame::BuildMagneticFrame(Bx_l, By_l, Bz_l);
      V_bx_l = 1.0;
      cr_frame::RotateToLabFrame(frame_l, V_bx_l, V_by_l, V_bz_l);
    }
    if (B_r > TINY_NUMBER) {
      const cr_frame::MagneticFrame frame_r =
          cr_frame::BuildMagneticFrame(Bx_r, By_r, Bz_r);
      V_bx_r = 1.0;
      cr_frame::RotateToLabFrame(frame_r, V_bx_r, V_by_r, V_bz_r);
    }

    const Real b_ivx_l = index == IVX ? V_bx_l
                       : (index == IVY ? V_by_l : V_bz_l);
    const Real b_ivy_l = index == IVX ? V_by_l
                       : (index == IVY ? V_bz_l : V_bx_l);
    const Real b_ivz_l = index == IVX ? V_bz_l
                       : (index == IVY ? V_bx_l : V_by_l);
    const Real b_ivx_r = index == IVX ? V_bx_r
                       : (index == IVY ? V_by_r : V_bz_r);
    const Real b_ivy_r = index == IVX ? V_by_r
                       : (index == IVY ? V_bz_r : V_bx_r);
    const Real b_ivz_r = index == IVX ? V_bz_r
                       : (index == IVY ? V_bx_r : V_by_r);

    HLLE_aux_(HLLE_rho_l, i) = rho_l;
    HLLE_aux_(HLLE_rho_r, i) = rho_r;

    HLLE_aux_(HLLE_gas_V_l, i) = w(ivx, left_k, left_j, left_i);
    HLLE_aux_(HLLE_gas_V_r, i) = w(ivx, k, j, i);

    HLLE_aux_(HLLE_Cs_l, i) = std::sqrt(gamma_gas*w(IPR, left_k, left_j, left_i)/rho_l);
    HLLE_aux_(HLLE_Cs_r, i) = std::sqrt(gamma_gas*w(IPR, k, j, i)/rho_r);

    HLLE_aux_(HLLE_vA_l, i) = B_l/std::sqrt(rho_l);
    HLLE_aux_(HLLE_vA_r, i) = B_r/std::sqrt(rho_r);

    HLLE_aux_(HLLE_b_ivx_l, i) = b_ivx_l;
    HLLE_aux_(HLLE_b_ivy_l, i) = b_ivy_l;
    HLLE_aux_(HLLE_b_ivz_l, i) = b_ivz_l;
    HLLE_aux_(HLLE_b_ivx_r, i) = b_ivx_r;
    HLLE_aux_(HLLE_b_ivy_r, i) = b_ivy_r;
    HLLE_aux_(HLLE_b_ivz_r, i) = b_ivz_r;

    HLLE_aux_(HLLE_dx, i) = index == IVX ? pco_->dx1f(i)
        : (index == IVY ? pco_->dx2f(j)/pco_->h2f(i)
                        : pco_->dx3f(k)/(pco_->h32f(j)*pco_->h31f(i)));
  }


  for (int CR_id=0; CR_id<NCRS; ++CR_id) {
    const int cr_energy_id = 4*CR_id;
    const int cr_ivx = cr_energy_id + ivx;
    const int cr_ivy = cr_energy_id + ivy;
    const int cr_ivz = cr_energy_id + ivz;

#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      const int left_i = i + left_i_offset;
      const Real rho_l = HLLE_aux_(HLLE_rho_l, i);
      const Real rho_r = HLLE_aux_(HLLE_rho_r, i);
      const Real gas_V_l = HLLE_aux_(HLLE_gas_V_l, i);
      const Real gas_V_r = HLLE_aux_(HLLE_gas_V_r, i);
      const Real Cs_l = HLLE_aux_(HLLE_Cs_l, i);
      const Real Cs_r = HLLE_aux_(HLLE_Cs_r, i);
      const Real vA_l = HLLE_aux_(HLLE_vA_l, i);
      const Real vA_r = HLLE_aux_(HLLE_vA_r, i);
      const Real b_ivx_l = HLLE_aux_(HLLE_b_ivx_l, i);
      const Real b_ivy_l = HLLE_aux_(HLLE_b_ivy_l, i);
      const Real b_ivz_l = HLLE_aux_(HLLE_b_ivz_l, i);
      const Real b_ivx_r = HLLE_aux_(HLLE_b_ivx_r, i);
      const Real b_ivy_r = HLLE_aux_(HLLE_b_ivy_r, i);
      const Real b_ivz_r = HLLE_aux_(HLLE_b_ivz_r, i);
      const Real dx = HLLE_aux_(HLLE_dx, i);
      const Real Ecr_l = cr_prim_l(cr_energy_id, i);
      const Real Ecr_r = cr_prim_r(cr_energy_id, i);
      const Real deltaPcr_left = deltaPcr_l(CR_id, i);
      const Real deltaPcr_right = deltaPcr_r(CR_id, i);
      const Real Pcr_perp_l = (Ecr_l - deltaPcr_left)/3.0;
      const Real Pcr_perp_r = (Ecr_r - deltaPcr_right)/3.0;
      const Real Ccr_l = TWO_3RD*std::sqrt(Ecr_l/rho_l);
      const Real Ccr_r = TWO_3RD*std::sqrt(Ecr_r/rho_r);
      const Real Cmax_l = std::sqrt(SQR(Cs_l) + SQR(Ccr_l) + SQR(vA_l));
      const Real Cmax_r = std::sqrt(SQR(Cs_r) + SQR(Ccr_r) + SQR(vA_r));

      const Real sigma_p_l = sigma_pL_array(CR_id, left_k, left_j, left_i) + sigma_pR_array(CR_id, left_k, left_j, left_i);
      const Real sigma_m_l = sigma_mL_array(CR_id, left_k, left_j, left_i) + sigma_mR_array(CR_id, left_k, left_j, left_i);
      const Real sigma_l = sigma_p_l + sigma_m_l;

      const Real sigma_p_r = sigma_pL_array(CR_id, k, j, i) + sigma_pR_array(CR_id, k, j, i);
      const Real sigma_m_r = sigma_mL_array(CR_id, k, j, i) + sigma_mR_array(CR_id, k, j, i);
      const Real sigma_r = sigma_p_r + sigma_m_r;

      const Real Vb_l_dx = sigma_l > 2.0*PI*SQRT3/dx/Vm ? (sigma_l > 0.2*PI/dx/Cmax_l ? 2.0*Cmax_l : 2.0*std::max(vA_l, Cs_l)) : Vm/SQRT3;
      const Real Vb_l_domain = sigma_l > 2.0*PI*SQRT3/domain_length/Vm ? (sigma_l > 0.2*PI/domain_length/Cmax_l ? 2.0*Cmax_l
                                          : 2.0*std::max(vA_l, Cs_l)) : Vm/SQRT3;

      const Real Vb_r_dx = sigma_r > 2.0*PI*SQRT3/dx/Vm ? (sigma_r > 0.2*PI/dx/Cmax_r ? 2.0*Cmax_r : 2.0*std::max(vA_r, Cs_r)) : Vm/SQRT3;
      const Real Vb_r_domain = sigma_r > 2.0*PI*SQRT3/domain_length/Vm ? (sigma_r > 0.2*PI/domain_length/Cmax_r ? 2.0*Cmax_r
                                          : 2.0*std::max(vA_r, Cs_r)) : Vm/SQRT3;

      const Real Vwave_l = std::abs(b_ivx_l) * std::max(Vb_l_dx, Vb_l_domain);
      const Real Vwave_r = std::abs(b_ivx_r) * std::max(Vb_r_dx, Vb_r_domain);

      const Real meanadv = 0.5*(gas_V_l + gas_V_r);
      const Real meandiffv = 0.5*(Vwave_l + Vwave_r);

      Real al = std::min(meanadv - meandiffv, gas_V_l - Vwave_l);
      Real ar = std::max(meanadv + meandiffv, gas_V_r + Vwave_r);
      ar = std::min(ar, Vm/std::sqrt(3.0));
      al = std::max(al, -Vm/std::sqrt(3.0));

      const Real bp = ar > 0.0 ? ar : 0.0;
      const Real bm = al < 0.0 ? al : 0.0;
      const Real crflux_energy_l  = Vm*cr_prim_l(cr_ivx, i) - bm*Ecr_l;
      const Real crflux_energy_r  = Vm*cr_prim_r(cr_ivx, i) - bp*Ecr_r;
      const Real crflux_ivxflux_l =
          Vm*(Pcr_perp_l + deltaPcr_left*SQR(b_ivx_l))
          - bm*cr_prim_l(cr_ivx, i);
      const Real crflux_ivxflux_r =
          Vm*(Pcr_perp_r + deltaPcr_right*SQR(b_ivx_r))
          - bp*cr_prim_r(cr_ivx, i);
      const Real crflux_ivyflux_l =
          Vm*deltaPcr_left*b_ivx_l*b_ivy_l - bm*cr_prim_l(cr_ivy, i);
      const Real crflux_ivyflux_r =
          Vm*deltaPcr_right*b_ivx_r*b_ivy_r - bp*cr_prim_r(cr_ivy, i);
      const Real crflux_ivzflux_l =
          Vm*deltaPcr_left*b_ivx_l*b_ivz_l - bm*cr_prim_l(cr_ivz, i);
      const Real crflux_ivzflux_r =
          Vm*deltaPcr_right*b_ivx_r*b_ivz_r - bp*cr_prim_r(cr_ivz, i);

      Real tmp = 0.0;
      if (bp != bm) tmp = 0.5*(bp + bm)/(bp - bm);

      cr_flux(cr_energy_id, k, j, i) = 0.5*(crflux_energy_l + crflux_energy_r) + (crflux_energy_l - crflux_energy_r)*tmp;
      cr_flux(cr_ivx, k, j, i) = 0.5*(crflux_ivxflux_l + crflux_ivxflux_r) + (crflux_ivxflux_l - crflux_ivxflux_r)*tmp;
      cr_flux(cr_ivy, k, j, i) = 0.5*(crflux_ivyflux_l + crflux_ivyflux_r) + (crflux_ivyflux_l - crflux_ivyflux_r)*tmp;
      cr_flux(cr_ivz, k, j, i) = 0.5*(crflux_ivzflux_l + crflux_ivzflux_r) + (crflux_ivzflux_l - crflux_ivzflux_r)*tmp;
    }
  }
  return;
}
