#ifndef COSMIC_RAY_HPP_
#define COSMIC_RAY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cosmic_ray.hpp
//! \brief definitions for the CosmicRay class

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/cc/cosmic_ray/bvals_cosmic_ray.hpp"
#include "../hydro/hydro_diffusion/hydro_diffusion.hpp"
#include "cr_scattering/cr_scattering.hpp"
#include "srcterms/cosmic_ray_srcterms.hpp"

class MeshBlock;
class ParameterInput;
class Hydro;
class CosmicRaySourceTerms;
class CRScattering;

//! \class CosmicRay
//! \brief CR data and functions

class CosmicRay {
 friend class EquationOfState;
 friend class Hydro;

 public:
  CosmicRay(MeshBlock *pmb, ParameterInput *pin);

  MeshBlock* pmy_block;
  Real Vm;
  Real Vm2;
  Real gamma_cr;
  bool pressure_anisotropy_flag;

  bool PressureAnisotropyEnabled() const {return pressure_anisotropy_flag;}


  // public data:
  // CR conservative variables: energy density and flux components
  AthenaArray<Real> cr_cons, cr_cons1, cr_cons2; // time-integrator memory register #1
  AthenaArray<Real> cr_cons0, cr_cons_fl_div;    // rkl2 STS memory registers;
  AthenaArray<Real> cr_cons_af_src;              // conservatives after explicit source terms

  // CR primitive variables
  AthenaArray<Real> cr_prim, cr_prim1, cr_prim_n;  // time-integrator memory register #3
  AthenaArray<Real> cr_flux[3];                    // face-averaged flux vector

  // storage for mesh refinement, SMR/AMR
  AthenaArray<Real> coarse_cr_cons_, coarse_cr_prim_; // coarse cr_cons and coarse cr_prim, used in mesh refinement
  int refinement_idx{-1};                             // vector of pointers in MeshRefinement class


  AthenaArray<Real> deltaPcr;          // CR pressure anisotropy for each CR species
  AthenaArray<Real> sigma_pL_array;    // Arrays of scattering coefficients by forward left-polarized waves
  AthenaArray<Real> sigma_pR_array;    // Arrays of scattering coefficients by forward right-polarized waves
  AthenaArray<Real> sigma_mL_array;    // Arrays of scattering coefficients by backward left-polarized waves
  AthenaArray<Real> sigma_mR_array;    // Arrays of scattering coefficients by backward right-polarized waves
  AthenaArray<Real> sigma_Lorentz_array;     // Arrays of scattering coefficients by Lorentz force
//   AthenaArray<Real> cs_cr_array;   // CR effective sound speed


  AthenaArray<Real> deltaPcr_n;        // CR pressure anisotropy at stage n
  AthenaArray<Real> sigma_pL_array_n;    // Arrays of scattering coefficients by forward left-polarized waves, at stage n
  AthenaArray<Real> sigma_pR_array_n;    // Arrays of scattering coefficients by forward right-polarized waves, at stage n
  AthenaArray<Real> sigma_mL_array_n;    // Arrays of scattering coefficients by backward left-polarized waves, at stage n
  AthenaArray<Real> sigma_mR_array_n;    // Arrays of scattering coefficients by backward right-polarized waves, at stage n
  AthenaArray<Real> sigma_Lorentz_array_n;     // Arrays of scattering coefficients by Lorentz force, at stage n
//   AthenaArray<Real> cs_cr_array_n; // CR effective sound speed at stage n

  AthenaArray<Real> Stage_I_delta_mom1, Stage_I_delta_mom2, Stage_I_delta_mom3; // Arrays of temporary delta momenta in Stage I
  AthenaArray<Real> Stage_I_vel1, Stage_I_vel2, Stage_I_vel3;                   // Arrays of temporary velocities

  // fourth-order intermediate quantities
  AthenaArray<Real> cr_cons_cc, cr_prim_cc;   // cell-centered approximations
  // Used to convert between cell-averaged and cell-centered CR states at fourth order.

  CosmicRayBoundaryVariable    crbvar;  // CR boundary variables Object (Cell-Centered)
  CRScattering                 crscat;  // CR-gas scattering coupling
  CosmicRaySourceTerms         crsrc;   // CR source terms

  int solver_id;        // 2 for the adaptive CR HLLE solver
  int cr_xorder;      // The reconstruction order of CR

  // Public functions:
  // Calculate CR flux
  void AddCosmicRayFluxDivergence(const Real wght, AthenaArray<Real> &cr_cons);  // Add flux divergence
  void AddCosmicRayFluxDivergence_STS(const Real wght, int stage,
          AthenaArray<Real> &cr_cons_out, AthenaArray<Real> &cr_cons_fl_div_out); // Add flux divergence

  void CalculateCosmicRayFluxes(AthenaArray<Real> &w, AthenaArray<Real> &cr_prim,
      FaceField &b, const AthenaArray<Real> &bcc, const int order);  // Calculate fluxes of CR

  // Riemann Solvers for CR
  // HLLE solver without CR effective sound speed
  void HLLENoCsRiemannSolverCosmicRay(const int k, const int j, const int il, const int iu,
                          const int index, const AthenaArray<Real> &b,
                          AthenaArray<Real> &w, AthenaArray<Real> &cr_prim_l, AthenaArray<Real> &cr_prim_r, 
                          AthenaArray<Real> &b_interface_l, AthenaArray<Real> &b_interface_r,
                          AthenaArray<Real> &deltaPcr_l, AthenaArray<Real> &deltaPcr_r,
                          AthenaArray<Real> &wl_,  AthenaArray<Real> &wr_, AthenaArray<Real> &cr_flux);

  // Compute the CR advection timestep in a MeshBlock
  Real NewAdvectionDt();

 private:
  enum HLLEAuxIndex {
    HLLE_rho_l,
    HLLE_rho_r,
    HLLE_gas_V_l,
    HLLE_gas_V_r,
    HLLE_Cs_l,
    HLLE_Cs_r,
    HLLE_vA_l,
    HLLE_vA_r,
    HLLE_b_ivx_l,
    HLLE_b_ivy_l,
    HLLE_b_ivz_l,
    HLLE_b_ivx_r,
    HLLE_b_ivy_r,
    HLLE_b_ivz_r,
    HLLE_dx,
    nHLLE_aux
  };

  // ptr to coordinates class
  Coordinates *pco_;
  // scratch space used to compute fluxes
  // 2D scratch arrays
  AthenaArray<Real> dt1_, dt2_, dt3_;                     // scratch arrays used in NewAdvectionDt
  AthenaArray<Real> cr_prim_l_, cr_prim_r_, cr_prim_lb_;  // left and right states in reconstruction
  AthenaArray<Real> b_interface_l_, b_interface_r_, b_interface_lb_;
  AthenaArray<Real> deltaPcr_l_, deltaPcr_r_, deltaPcr_lb_;

  // 1D scratch arrays
  AthenaArray<Real> x1face_area_, x2face_area_, x3face_area_; // face area in x1, x2, x3 directions
  AthenaArray<Real> x2face_area_p1_, x3face_area_p1_;         // face area in x2, x3 directions
  AthenaArray<Real> cell_volume_;                             // the volume of the cells
  AthenaArray<Real> crflx_;
  AthenaArray<Real> HLLE_aux_;

  // fourth-order
  // 4D scratch arrays
  AthenaArray<Real> scr1_nkji_, scr2_nkji_;
  AthenaArray<Real> cr_prim_l3d_, cr_prim_r3d_;
  // 1D scratch arrays
  AthenaArray<Real> laplacian_l_cr_fc_, laplacian_r_cr_fc_;
};
#endif // COSMIC_RAY_HPP_
