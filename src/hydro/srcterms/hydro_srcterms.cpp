//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file hydro_srcterms.cpp
//! \brief Class to implement source terms in the hydro equations

// C headers

// C++ headers
#include <cstring>    // strcmp
#include <iostream>
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>     // c_str()

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../mesh/mesh.hpp"
#include "../../orbital_advection/orbital_advection.hpp"
#include "../../parameter_input.hpp"
#include "../hydro.hpp"
#include "hydro_srcterms.hpp"

//! HydroSourceTerms constructor

HydroSourceTerms::HydroSourceTerms(Hydro *phyd, ParameterInput *pin) {
  pmy_hydro_ = phyd;
  hydro_sourceterms_defined = false;

  // read point mass or constant acceleration parameters from input block

  // set the point source only when the coordinate is spherical or 2D
  // It works even for cylindrical with the orbital advection.
  flag_point_mass_ = false;
  gm_ = pin->GetOrAddReal("problem","GM",0.0);
  bool orbital_advection_defined
         = (pin->GetOrAddInteger("orbital_advection","OAorder",0)!=0)?
           true : false;
  if (gm_ != 0.0) {
    if (std::strcmp(COORDINATE_SYSTEM, "spherical_polar") != 0
        && std::strcmp(COORDINATE_SYSTEM, "cylindrical") != 0) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point mass gravity works only in the cylindrical and "
          << "spherical polar coordinates." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (orbital_advection_defined) {
      hydro_sourceterms_defined = true;
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
               && phyd->pmy_block->block_size.nx3>1) {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "The point mass gravity deos not work in the 3D cylindrical "
          << "coordinates without orbital advection." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    } else {
      flag_point_mass_ = true;
      hydro_sourceterms_defined = true;
    }
  }
  g1_ = pin->GetOrAddReal("hydro","grav_acc1",0.0);
  if (g1_ != 0.0) hydro_sourceterms_defined = true;

  g2_ = pin->GetOrAddReal("hydro","grav_acc2",0.0);
  if (g2_ != 0.0) hydro_sourceterms_defined = true;

  g3_ = pin->GetOrAddReal("hydro","grav_acc3",0.0);
  if (g3_ != 0.0) hydro_sourceterms_defined = true;

  // read shearing box parameters from input block
  Omega_0_ = pin->GetOrAddReal("orbital_advection","Omega0",0.0);
  qshear_  = pin->GetOrAddReal("orbital_advection","qshear",0.0);
  ShBoxCoord_ = pin->GetOrAddInteger("orbital_advection","shboxcoord",1);

  // check flag for shearing source
  flag_shearing_source_ = 0;
  if(orbital_advection_defined) { // orbital advection source terms
    if(ShBoxCoord_ == 1) {
      flag_shearing_source_ = 1;
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in HydroSourceTerms constructor" << std::endl
          << "OrbitalAdvection does NOT work with shboxcoord = 2." << std::endl
          << "Check <orbital_advection> shboxcoord parameter in the input file."
          << std::endl;
      ATHENA_ERROR(msg);
    }
  } else if ((Omega_0_ !=0.0) && (qshear_ != 0.0)
             && std::strcmp(COORDINATE_SYSTEM, "cartesian") == 0) {
    flag_shearing_source_ = 2; // shearing box source terms
  } else if ((Omega_0_ != 0.0) &&
             (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
              || std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0)) {
    flag_shearing_source_ = 3; // rotating system source terms
  }

  if (flag_shearing_source_ != 0)
    hydro_sourceterms_defined = true;

  if (SELF_GRAVITY_ENABLED) hydro_sourceterms_defined = true;

  UserSourceTerm = phyd->pmy_block->pmy_mesh->UserSourceTerm_;
  if (UserSourceTerm != nullptr) hydro_sourceterms_defined = true;
  // scratch array for polar averaging
  if ((std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) && (phyd->pmy_block->block_size.nx3>1)){
    int ncells1 = phyd->pmy_block->block_size.nx1 + 2*NGHOST;
    int ncells3 = phyd->pmy_block->block_size.nx3 + 2*NGHOST;
    hyd_avg_.NewAthenaArray(NHYDRO, ncells3, ncells1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void HydroSourceTerms::AddHydroSourceTerms
//! \brief Adds source terms to conserved variables

void HydroSourceTerms::AddHydroSourceTerms(const Real time, const Real dt,
                           const AthenaArray<Real> *flux,
                           const AthenaArray<Real> &prim, const AthenaArray<Real> &cr_prim,
                           const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
                           AthenaArray<Real> &cons, AthenaArray<Real> &cr_cons,
                           AthenaArray<Real> &cons_scalar) {
  MeshBlock *pmb = pmy_hydro_->pmy_block;

  bool polar_inner = (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("polar"));
  bool polar_outer = (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("polar"));

  bool polar_wedge_inner = (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("polar_wedge"));
  bool polar_wedge_outer = (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("polar_wedge"));


  // accleration due to point mass (MUST BE AT ORIGIN)
  if (flag_point_mass_)
    PointMass(dt, flux, prim, cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0)
    ConstantAcceleration(dt, flux, prim, cons);

  // Add new source terms here
  if (SELF_GRAVITY_ENABLED)
    SelfGravity(dt, flux, prim, cons);

  // Sorce terms for orbital advection, shearing box, or rotating system
  if (flag_shearing_source_ == 1)
    OrbitalAdvectionSourceTerms(dt, flux, prim, cons);
  else if (flag_shearing_source_ == 2)
    ShearingBoxSourceTerms(dt, flux, prim, cons);
  else if (flag_shearing_source_ == 3)
    RotatingSystemSourceTerms(dt, flux, prim, cons);

  // MyNewSourceTerms()

  //  user-defined source terms
  if (UserSourceTerm != nullptr) {
    UserSourceTerm(pmb, time, dt, prim, cr_prim, prim_scalar, bcc, cons, cr_cons, cons_scalar);
  }
  // polar averaging
  if ((std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) && (pmb->block_size.nx3 > 1)) {
    if ((polar_inner || polar_wedge_inner)) {
      PolarAveragingHydro(cons, pmb->js,   4);
      PolarAveragingHydro(cons, pmb->js+1, 2);
    }
    if ((polar_outer || polar_wedge_outer)) {
      PolarAveragingHydro(cons, pmb->je,   4);
      PolarAveragingHydro(cons, pmb->je-1, 2);
    }
  }

  return;
}

void HydroSourceTerms::PolarAveragingHydro(AthenaArray<Real> &cons, int j, int nlayer)
{
  MeshBlock *pmb=pmy_hydro_->pmy_block;
  int is = pmb->is; int ks = pmb->ks;
  int ie = pmb->ie; int ke = pmb->ke;
  Real fac = 1.0/SQR(nlayer);

  for (int n=0; n<NHYDRO; ++n)
    for (int k=ks; k<=ke; ++k)
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        hyd_avg_(n,k,i)=0.0;

  for (int k=ks; k<=ke; ++k){
    for (int l=-nlayer+1; l<=nlayer-1; ++l){
      int myk = k+l;
      Real wght = (nlayer-fabs(l))*fac;
      myk = myk <= ke ? myk : myk-pmb->block_size.nx3;
      myk = myk >= ks ? myk : myk+pmb->block_size.nx3;
      for (int n=0; n<NHYDRO; ++n) {
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          hyd_avg_(n, k, i) += cons(n, myk, j, i)*wght;
      }
    }
  }
  for (int n=0; n<NHYDRO; ++n)
    for (int k=ks; k<=ke; ++k)
#pragma omp simd
      for (int i=is; i<=ie; ++i){
        cons(n, k, j, i) = hyd_avg_(n, k, i);
      }
  return;
}
