//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cosmic_ray_srcterms.cpp
//! \brief Class to implement source terms in the CR equations

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
#include "../cosmic_ray.hpp"
#include "cosmic_ray_srcterms.hpp"

//! CosmicRaySourceTerms constructor

CosmicRaySourceTerms::CosmicRaySourceTerms(CosmicRay *cr_pointer, ParameterInput *pin) {
  this->cr_pointer = cr_pointer;
  cr_sourceterms_defined = false;

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
      msg << "### FATAL ERROR in CosmicRaySourceTerms constructor" << std::endl
          << "The point mass gravity works only in the cylindrical and "
          << "spherical polar coordinates." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    }
    if (orbital_advection_defined) {
      cr_sourceterms_defined = true;
    } else if (std::strcmp(COORDINATE_SYSTEM, "cylindrical") == 0
               && cr_pointer->pmy_block->block_size.nx3>1) {
      std::stringstream msg;
      msg << "### FATAL ERROR in CosmicRaySourceTerms constructor" << std::endl
          << "The point mass gravity does not work in the 3D cylindrical "
          << "coordinates without orbital advection." << std::endl
          << "Check <problem> GM parameter in the input file." << std::endl;
      ATHENA_ERROR(msg);
    } else {
      flag_point_mass_ = true;
      cr_sourceterms_defined = true;
    }
  }
  g1_ = pin->GetOrAddReal("cosmic_ray", "grav_acc1", 0.0);
  if (g1_ != 0.0) cr_sourceterms_defined = true;

  g2_ = pin->GetOrAddReal("cosmic_ray", "grav_acc2", 0.0);
  if (g2_ != 0.0) cr_sourceterms_defined = true;

  g3_ = pin->GetOrAddReal("cosmic_ray", "grav_acc3", 0.0);
  if (g3_ != 0.0) cr_sourceterms_defined = true;


  if (cr_pointer->PressureAnisotropyEnabled()) cr_sourceterms_defined = true;

  // read shearing box parameters from input block
  Omega_0_    = pin->GetOrAddReal("orbital_advection",    "Omega0",     0.0);
  qshear_     = pin->GetOrAddReal("orbital_advection",    "qshear",     0.0);
  ShBoxCoord_ = pin->GetOrAddInteger("orbital_advection", "shboxcoord", 1);

  // check flag for shearing source
  flag_shearing_source_ = 0;
  if(orbital_advection_defined) { // orbital advection source terms
    if(ShBoxCoord_ == 1) {
      flag_shearing_source_ = 1;
    } else {
      std::stringstream msg;
      msg << "### FATAL ERROR in CosmicRaySourceTerms constructor" << std::endl
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
    cr_sourceterms_defined = true;

  if (SELF_GRAVITY_ENABLED) cr_sourceterms_defined = true;

  UserSourceTerm = cr_pointer->pmy_block->pmy_mesh->UserSourceTerm_;
  if (UserSourceTerm != nullptr)
    cr_sourceterms_defined = true;
  // scratch array for polar averaging
  if ((std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) && (cr_pointer->pmy_block->block_size.nx3>1)){
    int ncells1 = cr_pointer->pmy_block->block_size.nx1 + 2*NGHOST;
    int ncells3 = cr_pointer->pmy_block->block_size.nx3 + 2*NGHOST;
    cr_avg_.NewAthenaArray(NCRVARS, ncells3, ncells1);
  }
}

//----------------------------------------------------------------------------------------
//! \fn void CosmicRaySourceTerms::AddCosmicRaySourceTerms
//! \brief Adds source terms to conserved variables

void CosmicRaySourceTerms::AddCosmicRaySourceTerms(const Real time, const Real dt,
                     const AthenaArray<Real> *cr_flux, const AthenaArray<Real> &cr_prim,
                     AthenaArray<Real> &cr_cons, const AthenaArray<Real> &bcc) {
  MeshBlock  *pmb = cr_pointer->pmy_block;
  
  bool polar_inner = (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("polar"));
  bool polar_outer = (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("polar"));

  bool polar_wedge_inner = (pmb->pbval->block_bcs[BoundaryFace::inner_x2] == GetBoundaryFlag("polar_wedge"));
  bool polar_wedge_outer = (pmb->pbval->block_bcs[BoundaryFace::outer_x2] == GetBoundaryFlag("polar_wedge"));


  // accleration due to point mass (MUST BE AT ORIGIN)
  if (flag_point_mass_)
    PointMassCosmicRay(dt, cr_flux, cr_prim, cr_cons);

  // constant acceleration (e.g. for RT instability)
  if (g1_ != 0.0 || g2_ != 0.0 || g3_ != 0.0)
    ConstantAccelerationCosmicRay(dt, cr_flux, cr_prim, cr_cons);

  // Add new source terms here
  if (cr_pointer->PressureAnisotropyEnabled())
    PressureAnisotropySource(dt, cr_flux, cr_prim, cr_cons, cr_pointer->deltaPcr, cr_pointer->sigma_pL_array, cr_pointer->sigma_pR_array,
                        cr_pointer->sigma_mL_array, cr_pointer->sigma_mR_array, bcc);


  //if (SELF_GRAVITY_ENABLED) SelfGravity(dt, cr_flux, cr_prim, cr_cons);

  // Sorce terms for orbital advection, shearing box, or rotating system
  if (flag_shearing_source_ == 1)
    OrbitalAdvectionSourceTermsCosmicRay(dt, cr_flux, cr_prim, cr_cons);
  else if (flag_shearing_source_ == 2)
    ShearingBoxSourceTermsCosmicRay(dt, cr_flux, cr_prim, cr_cons);
  else if (flag_shearing_source_ == 3)
    RotatingSystemSourceTermsCosmicRay(dt, cr_flux, cr_prim, cr_cons);

  // polar averaging
  if ((std::strcmp(COORDINATE_SYSTEM, "spherical_polar") == 0) && (pmb->block_size.nx3 > 1)) {
    if ((polar_inner || polar_wedge_inner)) {
      PolarAveragingCosmicRay(cr_cons, pmb->js,   4);
      PolarAveragingCosmicRay(cr_cons, pmb->js+1, 2);
    }
    if ((polar_outer || polar_wedge_outer)) {
      PolarAveragingCosmicRay(cr_cons, pmb->je,   4);
      PolarAveragingCosmicRay(cr_cons, pmb->je-1, 2);
    }
  }

  return;
}

void CosmicRaySourceTerms::PolarAveragingCosmicRay(AthenaArray<Real> &cr_cons, int j, int nlayer)
{
  MeshBlock *pmb=cr_pointer->pmy_block;
  int is = pmb->is; int ks = pmb->ks;
  int ie = pmb->ie; int ke = pmb->ke;
  Real fac = 1.0/SQR(nlayer);

  for (int n=0; n<NCRVARS; ++n)
    for (int k=ks; k<=ke; ++k)
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        cr_avg_(n, k, i)=0.0;

  for (int k=ks; k<=ke; ++k){
    for (int l=-nlayer+1; l<=nlayer-1; ++l){
      int myk = k+l;
      Real wght = (nlayer-fabs(l))*fac;
      myk = myk <= ke ? myk : myk-pmb->block_size.nx3;
      myk = myk >= ks ? myk : myk+pmb->block_size.nx3;
      for (int n=0; n<NCRVARS; ++n){
#pragma omp simd
        for (int i=is; i<=ie; ++i)
          cr_avg_(n, k, i) += cr_cons(n, myk, j, i)*wght;
      }
    }
  }
  for (int n=0; n<NCRVARS; ++n)
    for (int k=ks; k<=ke; ++k)
#pragma omp simd
      for (int i=is; i<=ie; ++i){
        cr_cons(n, k, j, i) = cr_avg_(n, k, i);
      }
  return;
}
