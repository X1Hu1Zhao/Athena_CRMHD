//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file ct.cpp
//! \brief

// C headers

// C++ headers
#include <algorithm>  // max(), min()

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../mesh/mesh.hpp"
#include "field.hpp"

// OpenMP header
#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

//----------------------------------------------------------------------------------------
//! \fn  void Field::CT
//! \brief Constrained Transport implementation of dB/dt = -Curl(E), where E=-(v X B)

void Field::CT(const Real wght, FaceField &b_out) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  bool polar_inner = (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("polar"));
  bool polar_outer = (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("polar"));

  bool polar_wedge_inner = (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("polar_wedge"));
  bool polar_wedge_outer = (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("polar_wedge"));

  AthenaArray<Real> &e1 = e.x1e, &e2 = e.x2e, &e3 = e.x3e;
  AthenaArray<Real> &area = face_area_, &len = edge_length_, &len_p1 = edge_length_p1_;
// pre-processing: polar averaging
  if (pmb->block_size.nx3 > 1) {
    if ((polar_inner || polar_wedge_inner)) {
    PolarAveragingEMF(pmb->pfield->e, js,   4, true);
    PolarAveragingEMF(pmb->pfield->e, js+1, 2, true);
    }
    if ((polar_outer || polar_wedge_outer)) {
    PolarAveragingEMF(pmb->pfield->e, je,   4, false);
    PolarAveragingEMF(pmb->pfield->e, je-1, 2, false);
    }
  }

  //---- update B1
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // add curl(E) in 2D and 3D problem
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,area);
        pmb->pcoord->Edge3Length(k,j  ,is,ie+1,len);
        pmb->pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b_out.x1f(k,j,i) -=
              (wght/area(i))*(len_p1(i)*e3(k,j+1,i) - len(i)*e3(k,j,i));
        }

        if (pmb->block_size.nx3 > 1) {
          pmb->pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pmb->pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) +=
                (wght/area(i))*(len_p1(i)*e2(k+1,j,i) - len(i)*e2(k,j,i));
          }
        }
      }
    }
  }

  //---- update B2 (curl terms in 1D and 3D problems)
  for (int k=ks; k<=ke; ++k) {
    // reset loop limits for polar boundary
    int jl=js; int ju=je+1;
    if (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
        || pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge)
      jl=js+1;
    if (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
        || pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge)
      ju=je;
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->Face2Area(k,j,is,ie,area);
      pmb->pcoord->Edge3Length(k,j,is,ie+1,len);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_out.x2f(k,j,i) += (wght/area(i))*(len(i+1)*e3(k,j,i+1) - len(i)*e3(k,j,i));
      }
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Edge1Length(k  ,j,is,ie,len);
        pmb->pcoord->Edge1Length(k+1,j,is,ie,len_p1);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x2f(k,j,i) -=
              (wght/area(i))*(len_p1(i)*e1(k+1,j,i) - len(i)*e1(k,j,i));
        }
      }
    }
  }

  //---- update B3 (curl terms in 1D and 2D problems)
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face3Area(k,j,is,ie,area);
      pmb->pcoord->Edge2Length(k,j,is,ie+1,len);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_out.x3f(k,j,i) -= (wght/area(i))*(len(i+1)*e2(k,j,i+1) - len(i)*e2(k,j,i));
      }
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Edge1Length(k,j  ,is,ie,len);
        pmb->pcoord->Edge1Length(k,j+1,is,ie,len_p1);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x3f(k,j,i) +=
              (wght/area(i))*(len_p1(i)*e1(k,j+1,i) - len(i)*e1(k,j,i));
        }
      }
    }
  }
  return;
}

//----------------------------------------------------------------------------------------
//! \fn  void Field::CT_STS
//! \brief Constrained Transport implementation of dB/dt = -Curl(E), where E=-(v X B)
//! for STS. RKL2 registers are set to -Curl(E) update if first stage of RKL2 STS.

void Field::CT_STS(const Real wght, int stage,
                   FaceField &b_out, FaceField &ct_update_out) {
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int js = pmb->js; int ks = pmb->ks;
  int ie = pmb->ie; int je = pmb->je; int ke = pmb->ke;

  AthenaArray<Real> &e1 = e.x1e, &e2 = e.x2e, &e3 = e.x3e;
  AthenaArray<Real> &area = face_area_, &len = edge_length_, &len_p1 = edge_length_p1_;

  //---- update B1
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      // add curl(E) in 2D and 3D problem
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Face1Area(k,j,is,ie+1,area);
        pmb->pcoord->Edge3Length(k,j  ,is,ie+1,len);
        pmb->pcoord->Edge3Length(k,j+1,is,ie+1,len_p1);
#pragma omp simd
        for (int i=is; i<=ie+1; ++i) {
          b_out.x1f(k,j,i) -=
              (wght/area(i))*(len_p1(i)*e3(k,j+1,i) - len(i)*e3(k,j,i));
          if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
            ct_update_out.x1f(k,j,i) = -1.*((0.5*pmb->pmy_mesh->dt/area(i))
                                            * (len_p1(i)*e3(k,j+1,i)
                                               - len(i)*e3(k,j,i)));
          }
        }

        if (pmb->block_size.nx3 > 1) {
          pmb->pcoord->Edge2Length(k  ,j,is,ie+1,len);
          pmb->pcoord->Edge2Length(k+1,j,is,ie+1,len_p1);
#pragma omp simd
          for (int i=is; i<=ie+1; ++i) {
            b_out.x1f(k,j,i) +=
                (wght/area(i))*(len_p1(i)*e2(k+1,j,i) - len(i)*e2(k,j,i));
            if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
              ct_update_out.x1f(k,j,i) += ((0.5*pmb->pmy_mesh->dt/area(i))
                                           * (len_p1(i)*e2(k+1,j,i)
                                              - len(i)*e2(k,j,i)));
            }
          }
        }
      }
    }
  }

  //---- update B2 (curl terms in 1D and 3D problems)
  for (int k=ks; k<=ke; ++k) {
    // reset loop limits for polar boundary
    int jl=js; int ju=je+1;
    if (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar
        || pmb->pbval->block_bcs[BoundaryFace::inner_x2] == BoundaryFlag::polar_wedge)
      jl=js+1;
    if (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar
        || pmb->pbval->block_bcs[BoundaryFace::outer_x2] == BoundaryFlag::polar_wedge)
      ju=je;
    for (int j=jl; j<=ju; ++j) {
      pmb->pcoord->Face2Area(k,j,is,ie,area);
      pmb->pcoord->Edge3Length(k,j,is,ie+1,len);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_out.x2f(k,j,i) += (wght/area(i))*(len(i+1)*e3(k,j,i+1) - len(i)*e3(k,j,i));
        if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
          ct_update_out.x2f(k,j,i) = ((0.5*pmb->pmy_mesh->dt/area(i))
                                      * (len(i+1)*e3(k,j,i+1)
                                         - len(i)*e3(k,j,i)));
        }
      }
      if (pmb->block_size.nx3 > 1) {
        pmb->pcoord->Edge1Length(k  ,j,is,ie,len);
        pmb->pcoord->Edge1Length(k+1,j,is,ie,len_p1);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x2f(k,j,i) -=
              (wght/area(i))*(len_p1(i)*e1(k+1,j,i) - len(i)*e1(k,j,i));
          if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
            ct_update_out.x2f(k,j,i) -= ((0.5*pmb->pmy_mesh->dt/area(i))
                                         *(len_p1(i)*e1(k+1,j,i)
                                           - len(i)*e1(k,j,i)));
          }
        }
      }
    }
  }

  //---- update B3 (curl terms in 1D and 2D problems)
  for (int k=ks; k<=ke+1; ++k) {
    for (int j=js; j<=je; ++j) {
      pmb->pcoord->Face3Area(k,j,is,ie,area);
      pmb->pcoord->Edge2Length(k,j,is,ie+1,len);
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        b_out.x3f(k,j,i) -= (wght/area(i))*(len(i+1)*e2(k,j,i+1) - len(i)*e2(k,j,i));
        if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
          ct_update_out.x3f(k,j,i) = -1.*((0.5*pmb->pmy_mesh->dt/area(i))
                                           * (len(i+1)*e2(k,j,i+1)
                                              - len(i)*e2(k,j,i)));
        }
      }
      if (pmb->block_size.nx2 > 1) {
        pmb->pcoord->Edge1Length(k,j  ,is,ie,len);
        pmb->pcoord->Edge1Length(k,j+1,is,ie,len_p1);
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
          b_out.x3f(k,j,i) +=
              (wght/area(i))*(len_p1(i)*e1(k,j+1,i) - len(i)*e1(k,j,i));
          if (stage == 1 && pmb->pmy_mesh->sts_integrator == "rkl2") {
            ct_update_out.x3f(k,j,i) += ((0.5*pmb->pmy_mesh->dt/area(i))
                                         *(len_p1(i)*e1(k,j+1,i)
                                           - len(i)*e1(k,j,i)));
          }
        }
      }
    }
  }
  return;
}


// Polar averaging
void Field::PolarAveragingEMF(EdgeField &e, int j, int nlayer, bool north)
{
  MeshBlock *pmb=pmy_block;
  int is = pmb->is; int ks = pmb->ks;
  int ie = pmb->ie; int ke = pmb->ke;
  int dj = north ? 1 : 0;

  Real fac = 1.0/SQR(nlayer);

  for (int n=0; n<3; ++n)
    for (int k=ks; k<=ke; ++k)
#pragma omp simd
      for (int i=is; i<=ie+1; ++i)
        e_avg_(n,k,i)=0.0;

  for (int k=ks; k<=ke; ++k){
    for (int l=-nlayer+1; l<=nlayer-1; ++l){
      int myk = k+l;
      Real wght = (nlayer-fabs(l))*fac;
      myk = myk <= ke ? myk : myk-pmb->block_size.nx3;
      myk = myk >= ks ? myk : myk+pmb->block_size.nx3;
#pragma omp simd
      for (int i=is; i<=ie+1; ++i){
        e_avg_(0,k,i) += e.x1e(myk, j+dj,i)*wght;
        e_avg_(1,k,i) += e.x2e(myk, j,   i)*wght;
        e_avg_(2,k,i) += e.x3e(myk, j+dj,i)*wght;
      }
    }
  }

  for (int k=ks; k<=ke; ++k){
#pragma omp simd
    for (int i=is; i<=ie+1; ++i){
      e.x1e(k,j+dj,i) = e_avg_(0,k,i);
      e.x2e(k,j   ,i) = e_avg_(1,k,i);
      e.x3e(k,j+dj,i) = e_avg_(2,k,i);
    }
  }
#pragma omp simd
  for (int i=is; i<=ie+1; ++i){
    e.x1e(ke+1,j+dj,i) = e.x1e(ks,j+dj,i);
    e.x2e(ke+1,j,   i) = e.x2e(ks,j,   i);
  }
  return;
}
