//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file eos_cosmic_ray.cpp
//  \brief implements cosmic-ray functions in the EquationOfState class

// C headers

// C++ headers
#include <cmath>   // sqrt()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../cosmic_ray/cosmic_ray.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "eos.hpp"

//----------------------------------------------------------------------------------------
//! \fn void EquationOfState::CosmicRayConservedToPrimitive(AthenaArray<Real> &cr_cons,
//!           const AthenaArray<Real> &cr_prim_old,
//!           AthenaArray<Real> &r, Coordinates *pco,
//!           int il, int iu, int jl, int ju, int kl, int ku)
//! \brief Converts conserved into primitive cosmic-ray variables

void EquationOfState::CosmicRayConservedToPrimitive(
    AthenaArray<Real> &cr_cons, const AthenaArray<Real> &cr_prim_old,
    AthenaArray<Real> &cr_prim, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  for (int n=0; n<NCRS; ++n) {
    int cr_id = n;
    int cr_energy_id = 4*cr_id;
    int cr_flux1_id = cr_energy_id + 1;
    int cr_flux2_id = cr_energy_id + 2;
    int cr_flux3_id = cr_energy_id + 3;
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          Real &cons_cr_energy = cr_cons(cr_energy_id, k, j, i);
          Real &cons_cr_flux1 = cr_cons(cr_flux1_id, k, j, i);
          Real &cons_cr_flux2 = cr_cons(cr_flux2_id, k, j, i);
          Real &cons_cr_flux3 = cr_cons(cr_flux3_id, k, j, i);

          Real &prim_cr_energy = cr_prim(cr_energy_id, k, j, i);
          Real &prim_cr_flux1 = cr_prim(cr_flux1_id, k, j, i);
          Real &prim_cr_flux2 = cr_prim(cr_flux2_id, k, j, i);
          Real &prim_cr_flux3 = cr_prim(cr_flux3_id, k, j, i);

          cons_cr_energy = (cons_cr_energy > cr_floor_[cr_id])
                               ? cons_cr_energy : cr_floor_[cr_id];
          prim_cr_energy = cons_cr_energy;
          prim_cr_flux1 = cons_cr_flux1;
          prim_cr_flux2 = cons_cr_flux2;
          prim_cr_flux3 = cons_cr_flux3;
        }
      }
    }
  }
  return;
}



void EquationOfState::CosmicRayConservedToPrimitiveCellAverage(
    AthenaArray<Real> &cr_cons, const AthenaArray<Real> &cr_prim_old, AthenaArray<Real> &cr_prim,
    Coordinates *pco, int il, int iu, int jl, int ju, int kl, int ku) {
  MeshBlock  *pmb = pmy_block_;
  CosmicRay *cr_pointer = pmb->cr_pointer;
  int nl = 0; int nu = NCRVARS - 1;

  Real h = pco->dx1f(il);  // pco->dx1f(i); inside loop
  Real C = (h*h)/24.0;

  // Fourth-order accurate approx to cell-centered conserved and primitive variables
  //AthenaArray<Real> &w_cc = ph->w_cc, &w = ph->w; // &u_cc = ph->u_cc;
  AthenaArray<Real> &cr_prim_cc = cr_pointer->cr_prim_cc, &cr_cons_cc = cr_pointer->cr_cons_cc;
  // Laplacians of cell-averaged conserved and 2nd order accurate primitive variables
  AthenaArray<Real> &laplacian_cc = cr_pointer->scr1_nkji_;

  // Compute and store Laplacian of cell-averaged conserved variables
  pco->Laplacian(cr_cons, laplacian_cc, il, iu, jl, ju, kl, ku, nl, nu);

  // Compute fourth-order approximation to cell-centered conserved variables
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // We do not actually need to store all cell-centered conserved variables,
          // but the ConservedToPrimitive() implementation operates on 4D arrays
          cr_cons_cc(n,k,j,i) = cr_cons(n,k,j,i) - C*laplacian_cc(n,k,j,i);
        }
      }
    }
  }

  // Compute Laplacian of 2nd-order approximation to cell-averaged primitive variables
  pco->Laplacian(cr_prim, laplacian_cc, il, iu, jl, ju, kl, ku, nl, nu);

  // Convert cell-centered conserved values to cell-centered primitive values
  CosmicRayConservedToPrimitive(cr_cons_cc, cr_prim_old, cr_prim_cc, pco,
                                il, iu, jl, ju, kl, ku);

  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          // Compute fourth-order approximation to cell-averaged primitive variables
          cr_prim(n,k,j,i) = cr_prim_cc(n,k,j,i) + C*laplacian_cc(n,k,j,i);
        }
      }
    }
  }

  // Reapply primitive variable floors
  // Cannot fuse w/ above loop since floors are applied to all NHYDRO variables at once
  for (int n=nl; n<=nu; ++n) {
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          ApplyCosmicRayPrimitiveConservedFloors(cr_cons, cr_prim, n, k, j, i);
        }
      }
    }
  }

  return;
}


//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::CosmicRayPrimitiveToConserved(const AthenaArray<Real> &cr_prim
//           AthenaArray<Real> &cr_cons, Coordinates *pco,
//           int il, int iu, int jl, int ju, int kl, int ku);
// \brief Converts primitive variables into conservative variables

void EquationOfState::CosmicRayPrimitiveToConserved(
    const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons, Coordinates *pco,
    int il, int iu, int jl, int ju, int kl, int ku) {
  for (int n=0; n<NCRS; ++n) {
    int cr_id = n;
    int cr_energy_id = 4*cr_id;
    int cr_flux1_id = cr_energy_id + 1;
    int cr_flux2_id = cr_energy_id + 2;
    int cr_flux3_id = cr_energy_id + 3;
    for (int k=kl; k<=ku; ++k) {
      for (int j=jl; j<=ju; ++j) {
#pragma omp simd
        for (int i=il; i<=iu; ++i) {
          cr_cons(cr_energy_id, k, j, i) = cr_prim(cr_energy_id, k, j, i);
          cr_cons(cr_flux1_id, k, j, i)  = cr_prim(cr_flux1_id, k, j, i);
          cr_cons(cr_flux2_id, k, j, i)  = cr_prim(cr_flux2_id, k, j, i);
          cr_cons(cr_flux3_id, k, j, i)  = cr_prim(cr_flux3_id, k, j, i);
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
// \!fn void EquationOfState::ApplyCosmicRayFloors(AthenaArray<Real> &prim, int n,
//                                                     int k, int j, int i)
// \brief Apply the CR energy floor to reconstructed left/right interface states.

void EquationOfState::ApplyCosmicRayFloors(AthenaArray<Real> &cr_prim, int n, int k, int j, int i) {


  int cr_id = n/4;
  int cr_energy_id  = 4*cr_id;

  Real &cr_energy_prim = cr_prim(cr_energy_id, k, j, i);
  // Apply the CR energy floor without modifying the conserved cell average.
  cr_energy_prim = (cr_energy_prim > cr_floor_[cr_id]) ?  cr_energy_prim : cr_floor_[cr_id];
  return;
}

void EquationOfState::ApplyCosmicRayPrimitiveConservedFloors(
    AthenaArray<Real> &cr_cons, AthenaArray<Real> &cr_prim,
    int n, int k, int j, int i) {

  int cr_id = n/4;
  int energy_id  = 4*cr_id;

  Real &cr_energy_cons = cr_cons(energy_id, k, j, i);
  Real &cr_energy_prim = cr_prim(energy_id, k, j, i);

  cr_energy_cons = (cr_energy_cons > cr_floor_[cr_id]) ?  cr_energy_cons : cr_floor_[cr_id];
  cr_energy_prim = cr_energy_cons;

  // this next line, when applied indiscriminately, erases the accuracy gains performed in
  // the 4th order stencils, since <r> != <s>*<1/inv_cr_energy>, in general
  // cr_prim_n = cr_cons_n*inv_cr_energy;
  // however, if r_n is riding the variable floor, it probably should be applied so that
  // s_n = rho*r_n is consistent (more concerned with conservation than order of accuracy
  // when quantities are floored)
  return;
}
