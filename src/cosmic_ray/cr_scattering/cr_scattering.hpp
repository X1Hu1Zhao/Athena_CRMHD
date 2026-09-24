#ifndef CR_SCATTERING_HPP_
#define CR_SCATTERING_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_scattering.hpp
//! \brief defines class CRScattering
//! Contains data and functions that implement physical CR scattering and gas-feedback terms

// C headers

// C++ headers
#include <array>
#include <cstring>    // strcmp
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../cosmic_ray.hpp"

// Forward declarations
class CosmicRay;
class ParameterInput;

//! \class CRScattering
//! \brief data and functions for CR-gas scattering coupling
class CRScattering {
 public:
  CRScattering(CosmicRay *cr_pointer, ParameterInput *pin);

  // Flag
  bool Scattering_Flag;
  bool CRFeedback_Flag;  // true or false, the flag of CR feedback term
  bool CR_EnergyFeedback_Flag;
  Real Vm, Vm2, inv_Vm, inv_Vm2, gamma_cr;
  AthenaArray<Real> cs_cr_array, cs_cr_array_n;
  std::array<Real, NCRS> const_sigma_pL{};
  std::array<Real, NCRS> const_sigma_pR{};
  std::array<Real, NCRS> const_sigma_mL{};
  std::array<Real, NCRS> const_sigma_mR{};
  std::array<Real, NCRS> const_sigma_Lorentz{};
  std::array<Real, NCRS> sigma0_streaming{};

  // Select the scattering integrators
  void ScatteringIntegrator(const int stage, const Real dt,
    const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim, const AthenaArray<Real> &bcc,
    AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

  void SetConstantScatteringCoefficients(
    AthenaArray<Real> &sigma_pL, AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL, AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz, AthenaArray<Real> &deltaPcr,
    int is, int ie, int js, int je, int ks, int ke);

  void SetStreamingScatteringCoefficients(
    const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
    const AthenaArray<Real> &bcc,
    AthenaArray<Real> &sigma_pL, AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL, AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz, AthenaArray<Real> &deltaPcr,
    int is, int ie, int js, int je, int ks, int ke);

  Real NewScatteringDt();

  void SetProperties(const Real time, const AthenaArray<Real> &w,
    const AthenaArray<Real> &cr_prim, const AthenaArray<Real> &bcc);
 

  // Matrix Addition
  void Add(const AthenaArray<Real> &a_matrix, const Real b_num,
           const AthenaArray<Real> &b_matrix, AthenaArray<Real> &c_matrix);

  void Add(AthenaArray<Real> &a_matrix, const Real b_num, const AthenaArray<Real> &b_matrix);

  void Add(const Real a_num, const Real b_num,
           const AthenaArray<Real> &b_matrix, AthenaArray<Real> &c_matrix);

  void Add(const Real a_num, const Real b_num, AthenaArray<Real> &b_matrix);

  // Matrix Multiplication
  void Multiply(const AthenaArray<Real> &a_matrix,
                const AthenaArray<Real> &b_matrix, AthenaArray<Real> &c_matrix);

  void MultiplyVector(const AthenaArray<Real> &a_matrix,
                      const AthenaArray<Real> &b_vector, AthenaArray<Real> &c_vector);

  void Multiply(const Real a_num, const AthenaArray<Real> &b_matrix, AthenaArray<Real> &c_matrix);

  void Multiply(const Real a_num, AthenaArray<Real> &b_matrix);

  void RotateArraysToFieldAligned(
      const AthenaArray<Real> &cos_xy, const AthenaArray<Real> &sin_xy,
      const AthenaArray<Real> &cos_z, const AthenaArray<Real> &sin_z,
      AthenaArray<Real> &v1, AthenaArray<Real> &v2, AthenaArray<Real> &v3);

  void RotateArraysToLabFrame(
      const AthenaArray<Real> &cos_xy, const AthenaArray<Real> &sin_xy,
      const AthenaArray<Real> &cos_z, const AthenaArray<Real> &sin_z,
      AthenaArray<Real> &v1, AthenaArray<Real> &v2, AthenaArray<Real> &v3);

  // Matrix Inverse
  void LUdecompose(const AthenaArray<Real> &a_matrix, AthenaArray<Real> &index_vector,
                    AthenaArray<Real> &lu_matrix);

  // Solve A*x = b
  void SolveLinearEquation(const AthenaArray<Real> &index_vector, const AthenaArray<Real> &lu_matrix,
                                 AthenaArray<Real> &b_matrix, AthenaArray<Real> &x_matrix);

  void SolveMultipleLinearEquation(const AthenaArray<Real> &index_vector,
      const AthenaArray<Real> &lu_matrix, AthenaArray<Real> &b_matrix, AthenaArray<Real> &x_matrix);

  // Calculate the inverse of matrix
  void Inverse(const AthenaArray<Real> &index_vector, const AthenaArray<Real> &lu_matrix,
                  AthenaArray<Real> &a_matrix, AthenaArray<Real> &a_inv_matrix);

  // Time Integrators
  // Fully Implicit Integartors
  // Backward Euler methods (Backward Differentiation Formula 1, BDF1), 1st order time convergence
  void BackwardEulerFeedback(const int stage, const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim, const AthenaArray<Real> &bcc,
      AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

  void BackwardEulerNoEnergyFeedback(const int stage, const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

  void BackwardEulerNoFeedback(const int stage, const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

  // Van Leer 2 Implicit methods, 2nd order time convergence
  void VL2ImplicitFeedback(const int stage,
      const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

  void VL2ImplicitNoEnergyFeedback(const int stage,
      const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

  void VL2ImplicitNoFeedback(const int stage,
      const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons);

 private:
  std::string integrator_;        // Time Integrator
  int         method_id_;         // The integrator method id
  std::string scattering_integrator_;
  std::string coefficient_model_;
  CosmicRay  *cr_pointer;   // ptr to CosmicRay containing this CRScattering
  MeshBlock   *pmb_;              // ptr to meshblock containing this CRScattering
  Coordinates *pco_;              // ptr to coordinates class
  Hydro       *pmy_hydro_;        // ptr to hydro class
  int         orbital_advection_order_;
  CRScatteringCoeffFunc UserDefinedCRScattering;


  AthenaArray<Real> force, delta_mom;
  AthenaArray<Real> force_Ecr_n, force_x1_n, force_x2_n, force_x3_n;

  AthenaArray<Real> delta_Ecr, delta_mom1,       delta_mom2,       delta_mom3;
  AthenaArray<Real> delta_Ecr_src, delta_mom1_src,   delta_mom2_src,   delta_mom3_src;

  AthenaArray<Real> mom1_prim,   mom2_prim,   mom3_prim;
  AthenaArray<Real> mom1_prim_n, mom2_prim_n, mom3_prim_n;
  AthenaArray<Real> dens_prim, dens_prim_n;
  

  AthenaArray<Real> sigma_p, sigma_m, Omega;
  AthenaArray<Real> sigma_p_n, sigma_m_n, Omega_n;
  AthenaArray<Real> cosXY, sinXY, cosZ, sinZ, vA, fwdwvspeed, bwdwvspeed;

  AthenaArray<Real> idx_vector, lu_matrix;
  AthenaArray<Real> jacobi, jacobi_n, product, lambda, lambda_inv;

  AthenaArray<Real> biggest_arr, temp_arr;
  AthenaArray<Real> det_arr, scale_arr;
  AthenaArray<Real> sum_arr, xx_arr;
  AthenaArray<int>  mmax_arr;

  AthenaArray<Real> inv_gas_rho, inv_gas_rho_n;
  AthenaArray<Real> temp_A, temp_B, temp_C;
  AthenaArray<Real> scattering_dx1_, scattering_dx2_, scattering_dx3_;
};
#endif // CR_SCATTERING_HPP_
