//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cr_frame_transform.hpp
//! \brief Inline transformations between lab and magnetic-field-aligned frames.

#ifndef COSMIC_RAY_CR_FRAME_TRANSFORM_HPP_
#define COSMIC_RAY_CR_FRAME_TRANSFORM_HPP_

#include <cmath>

#include "../athena.hpp"

namespace cr_frame {

struct MagneticFrame {
  Real cos_xy;
  Real sin_xy;
  Real cos_z;
  Real sin_z;
};

inline MagneticFrame BuildMagneticFrame(const Real bx, const Real by, const Real bz,
                                        const Real b, const Real bxy_squared) {
  MagneticFrame frame;
  if (bxy_squared > TINY_NUMBER) {
    frame.cos_xy = bx/std::sqrt(bx*bx + by*by);
    frame.sin_xy = by/std::sqrt(bx*bx + by*by);
  } else {
    frame.cos_xy = 1.0;
    frame.sin_xy = 0.0;
  }

  frame.cos_z = bz/b;
  frame.sin_z = std::sqrt(1.0 - frame.cos_z*frame.cos_z);
  return frame;
}

inline MagneticFrame BuildMagneticFrame(const Real bx, const Real by, const Real bz) {
  const Real bxy = std::sqrt(bx*bx + by*by);
  const Real b = std::sqrt(bxy*bxy + bz*bz);

  MagneticFrame frame;
  if (bxy > TINY_NUMBER) {
    frame.cos_xy = bx/bxy;
    frame.sin_xy = by/bxy;
  } else {
    frame.cos_xy = 1.0;
    frame.sin_xy = 0.0;
  }

  frame.cos_z = bz/b;
  frame.sin_z = std::sqrt(1.0 - frame.cos_z*frame.cos_z);
  return frame;
}

inline void RotateToFieldAligned(const MagneticFrame &frame,
                                 Real &v1, Real &v2, Real &v3) {
  const Real new_v1 = frame.cos_xy*frame.sin_z*v1
                    + frame.sin_xy*frame.sin_z*v2 + frame.cos_z*v3;
  const Real new_v2 = -frame.sin_xy*v1 + frame.cos_xy*v2;
  const Real new_v3 = -frame.cos_xy*frame.cos_z*v1
                    - frame.cos_z*frame.sin_xy*v2 + frame.sin_z*v3;

  v1 = new_v1;
  v2 = new_v2;
  v3 = new_v3;
}

inline void RotateToLabFrame(const MagneticFrame &frame,
                             Real &v1, Real &v2, Real &v3) {
  const Real new_v1 = frame.cos_xy*frame.sin_z*v1
                    - frame.sin_xy*v2 - frame.cos_xy*frame.cos_z*v3;
  const Real new_v2 = frame.sin_xy*frame.sin_z*v1
                    + frame.cos_xy*v2 - frame.sin_xy*frame.cos_z*v3;
  const Real new_v3 = frame.cos_z*v1 + frame.sin_z*v3;

  v1 = new_v1;
  v2 = new_v2;
  v3 = new_v3;
}

}  // namespace cr_frame

#endif  // COSMIC_RAY_CR_FRAME_TRANSFORM_HPP_
