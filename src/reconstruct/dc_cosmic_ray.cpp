//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file dc_cosmic_ray.cpp
//! \brief piecewise constant (donor cell) reconstruction
//! Operates on the entire nx4 range of a single AthenaArray<Real> input (no MHD).
//! No assumptions of hydrodynamic fluid variable input; no characteristic projection.

// C headers

// C++ headers

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "reconstruction.hpp"
#include <iostream>

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX1_CosmicRay(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief reconstruct L/R surfaces of the i-th cells

void Reconstruction::DonorCellX1_CosmicRay(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // compute L/R states for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i+1) =  qr(n,i) = q(n,k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX2_CosmicRay(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

void Reconstruction::DonorCellX2_CosmicRay(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // compute L/R states for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qr(n,i) = q(n,k,j,i);
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn Reconstruction::DonorCellX3_CosmicRay(const int k, const int j,
//!                              const int il, const int iu,
//!                              const AthenaArray<Real> &w, const AthenaArray<Real> &bcc,
//!                              AthenaArray<Real> &wl, AthenaArray<Real> &wr)
//! \brief

void Reconstruction::DonorCellX3_CosmicRay(const int k, const int j, const int il, const int iu,
                                 const AthenaArray<Real> &q,
                                 AthenaArray<Real> &ql, AthenaArray<Real> &qr) {
  const int nu = q.GetDim4() - 1;
  // compute L/R states for each variable
  for (int n=0; n<=nu; ++n) {
#pragma omp simd
    for (int i=il; i<=iu; ++i) {
      ql(n,i) = qr(n,i) = q(n,k,j,i);
    }
  }
  return;
}


void Reconstruction::DonorCellX1_MagneticField(const int k, const int j,
                                 const int il, const int iu,
                                 const AthenaArray<Real> &b1,
                                 const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &bl, AthenaArray<Real> &br) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB1,i) = b1(k,j,i);
    br(IB1,i) = b1(k,j,i);
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB2,i+1) = br(IB2,i) = bcc(IB2,k,j,i);
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB3,i+1) = br(IB3,i) = bcc(IB3,k,j,i);
  }
  return;
}

void Reconstruction::DonorCellX2_MagneticField(const int k, const int j,
                                 const int il, const int iu,
                                 const AthenaArray<Real> &b2,
                                 const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &bl, AthenaArray<Real> &br) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB1,i) = br(IB1,i) = bcc(IB1,k,j,i);
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB2,i) = b2(k,j+1,i);
    br(IB2,i) = b2(k,j,i);
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB3,i) = br(IB3,i) = bcc(IB3,k,j,i);
  }
  return;
}

void Reconstruction::DonorCellX3_MagneticField(const int k, const int j,
                                 const int il, const int iu,
                                 const AthenaArray<Real> &b3,
                                 const AthenaArray<Real> &bcc,
                                 AthenaArray<Real> &bl, AthenaArray<Real> &br) {
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB1,i) = br(IB1,i) = bcc(IB1,k,j,i);
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB2,i) = br(IB2,i) = bcc(IB2,k,j,i);
  }
#pragma omp simd
  for (int i=il; i<=iu; ++i) {
    bl(IB3,i) = b3(k+1,j,i);
    br(IB3,i) = b3(k,j,i);
  }
  return;
}
