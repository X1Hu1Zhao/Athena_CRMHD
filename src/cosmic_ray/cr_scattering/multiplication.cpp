//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file multiplication.cpp
//! Compute the multiplications between matrixes.

// C++ headers
#include <algorithm>   // min,max
#include <cstring>    // strcmp
#include <limits>
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../defs.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"


// Matrix Multiplication
void CRScattering::Multiply(const AthenaArray<Real> &a_matrix,
          const AthenaArray<Real> &b_matrix, AthenaArray<Real> &c_matrix) {

  int is = pmb_->is, ie = pmb_->ie;

  for(int m=0; m<NCRVARS; ++m) {
    for(int n=0; n<NCRVARS; ++n) {
      for(int s=0; s<NCRVARS; ++s) {
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          c_matrix(m, n, i) += a_matrix(m, s, i)*b_matrix(s, n, i);
        }
      }
    }
  }
  return;
}


void CRScattering::MultiplyVector(const AthenaArray<Real> &a_matrix,
                      const AthenaArray<Real> &b_vector, AthenaArray<Real> &c_vector) {

  int is = pmb_->is, ie = pmb_->ie;

  for(int m=0; m<NCRVARS; ++m) {
    for(int n=0; n<NCRVARS; ++n) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        c_vector(n, i) += a_matrix(m, n, i)*b_vector(m, i);
      }
    }
  }
  return;
}


void CRScattering::Multiply(const Real a_num, const AthenaArray<Real> &b_matrix,
              AthenaArray<Real> &c_matrix) {

  int is = pmb_->is, ie = pmb_->ie;

  for(int m=0; m<NCRVARS; ++m) {
    for(int n=0; n<NCRVARS; ++n) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        c_matrix(m, n, i) = a_num*b_matrix(m, n, i);
      }
    }
  }
  return;
}


void CRScattering::Multiply(const Real a_num, AthenaArray<Real> &b_matrix) {

  int is = pmb_->is, ie = pmb_->ie;

  for(int m=0; m<NCRVARS; ++m) {
    for(int n=0; n<NCRVARS; ++n) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_matrix(m, n, i) *= a_num;
      }
    }
  }
  return;
}
