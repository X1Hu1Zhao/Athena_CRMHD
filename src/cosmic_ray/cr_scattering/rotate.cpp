//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file rotate.cpp
//! \brief Rotate vector arrays into and out of the magnetic-field-aligned frame.

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr_frame_transform.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"

void CRScattering::RotateArraysToFieldAligned(
    const AthenaArray<Real> &cos_xy, const AthenaArray<Real> &sin_xy,
    const AthenaArray<Real> &cos_z, const AthenaArray<Real> &sin_z,
    AthenaArray<Real> &v1, AthenaArray<Real> &v2, AthenaArray<Real> &v3) {
  const int is = pmb_->is;
  const int ie = pmb_->ie;

  for (int n=0; n<NSPECIES; ++n) {
#pragma omp simd
    for (int i=is; i<=ie; ++i) {
      const cr_frame::MagneticFrame frame =
          {cos_xy(i), sin_xy(i), cos_z(i), sin_z(i)};
      Real rotated_v1 = v1(n, i);
      Real rotated_v2 = v2(n, i);
      Real rotated_v3 = v3(n, i);
      cr_frame::RotateToFieldAligned(
          frame, rotated_v1, rotated_v2, rotated_v3);

      v1(n, i) = rotated_v1;
      v2(n, i) = rotated_v2;
      v3(n, i) = rotated_v3;
    }
  }
  return;
}

void CRScattering::RotateArraysToLabFrame(
    const AthenaArray<Real> &cos_xy, const AthenaArray<Real> &sin_xy,
    const AthenaArray<Real> &cos_z, const AthenaArray<Real> &sin_z,
    AthenaArray<Real> &v1, AthenaArray<Real> &v2, AthenaArray<Real> &v3) {
  const int is = pmb_->is;
  const int ie = pmb_->ie;

  for (int n=0; n<NSPECIES; ++n) {
#pragma omp simd
    for (int i=is; i<=ie; ++i) {
      const cr_frame::MagneticFrame frame =
          {cos_xy(i), sin_xy(i), cos_z(i), sin_z(i)};
      Real rotated_v1 = v1(n, i);
      Real rotated_v2 = v2(n, i);
      Real rotated_v3 = v3(n, i);
      cr_frame::RotateToLabFrame(
          frame, rotated_v1, rotated_v2, rotated_v3);

      v1(n, i) = rotated_v1;
      v2(n, i) = rotated_v2;
      v3(n, i) = rotated_v3;
    }
  }
  return;
}
