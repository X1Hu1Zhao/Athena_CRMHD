#ifndef COSMIC_RAY_SRCTERMS_HPP_
#define COSMIC_RAY_SRCTERMS_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cosmic_ray_srcterms.hpp
//! \brief defines class CosmicRaySourceTerms
//! Contains data and functions that implement physical (not coordinate) source terms

// C headers

// C++ headers

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../cosmic_ray.hpp"

// Forward declarations
class CosmicRay;
class ParameterInput;

//! \class CosmicRaySourceTerms
//! \brief data and functions for physical CR source terms
class CosmicRaySourceTerms {
 public:
  CosmicRaySourceTerms(CosmicRay *cr_pointer, ParameterInput *pin);

  // accessors
  Real GetGM() const {return gm_;}

  // data
  bool cr_sourceterms_defined;

  // functions
  void AddCosmicRaySourceTerms(const Real time, const Real dt,
      const AthenaArray<Real> *cr_flux, const AthenaArray<Real> &cr_prim,
      AthenaArray<Real> &cr_cons, const AthenaArray<Real> &bcc);

  // Central stellar gravity source term in disk problem
  void PointMassCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
      const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons);
  void ConstantAccelerationCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
                            const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons);
  void PressureAnisotropySource(const Real dt, const AthenaArray<Real> *cr_flux,
      const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons,
      const AthenaArray<Real> &deltaPcr,
      const AthenaArray<Real> &sigma_pL_array,
      const AthenaArray<Real> &sigma_pR_array,
      const AthenaArray<Real> &sigma_mL_array,
      const AthenaArray<Real> &sigma_mR_array,
      const AthenaArray<Real> &bcc);

  // shearing box src terms
  void ShearingBoxSourceTermsCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
      const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons);

  void OrbitalAdvectionSourceTermsCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
      const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons);

  void RotatingSystemSourceTermsCosmicRay(const Real dt, const AthenaArray<Real> *cr_flux,
      const AthenaArray<Real> &cr_prim, AthenaArray<Real> &cr_cons);

  SrcTermFunc UserSourceTerm;
  void PolarAveragingCosmicRay(AthenaArray<Real> &cr_cons, int j, int nlayer);

 private:
  friend class CosmicRay;
  CosmicRay *cr_pointer;      // ptr to CosmicRay containing this CosmicRaySourceTerms
  Real gm_;                         // GM for point mass MUST BE LOCATED AT ORIGIN
  Real g1_, g2_, g3_;               // constant acc'n in each direction
  Real Omega_0_, qshear_;           // Orbital freq and shear rate in shearing box
  int  ShBoxCoord_;                 // ShearCoordinate type: 1=xy (default), 2=xz
  AthenaArray<Real> cr_avg_;      // storage for polar averaging
  bool flag_point_mass_;            // flag for calling PointMass function
  int  flag_shearing_source_;       // 1=orbital advection, 2=shearing box, 3=rotating system
};
#endif // COSMIC_RAY_SRCTERMS_HPP_
