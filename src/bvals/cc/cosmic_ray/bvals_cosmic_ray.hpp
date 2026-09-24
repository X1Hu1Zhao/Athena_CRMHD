#ifndef BVALS_CC_COSMIC_RAY_BVALS_COSMIC_RAY_HPP_
#define BVALS_CC_COSMIC_RAY_BVALS_COSMIC_RAY_HPP_
//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file bvals_cosmic_ray.hpp
//! \brief

// C headers

// C++ headers

// Athena++ headers
#include "../../../athena.hpp"
#include "../../../athena_arrays.hpp"
#include "../bvals_cc.hpp"

//----------------------------------------------------------------------------------------
//! \class CellCenteredBoundaryVariable
//! \brief

class CosmicRayBoundaryVariable : public CellCenteredBoundaryVariable {
 public:
  CosmicRayBoundaryVariable(MeshBlock *pmb,
                        AthenaArray<Real> *var_cosmic_ray, AthenaArray<Real> *coarse_var,
                        AthenaArray<Real> *var_flux,
                        CosmicRayBoundaryQuantity cr_type);
                                                // AthenaArray<Real> &prim);
  virtual ~CosmicRayBoundaryVariable() = default;

  // switch between CosmicRay class members "cr_cons" and "cr_prim" (or "cr_cons" and "cr_cons1", ...)
  void SwapCosmicRayQuantity(AthenaArray<Real> &var_cosmic_ray, CosmicRayBoundaryQuantity cr_type);
  void SelectCoarseBuffer(CosmicRayBoundaryQuantity cr_type);

  void AddCosmicRayShearForInit();
  void ShearQuantities(AthenaArray<Real> &shear_cc_, bool upper) override;

  //!@{
  //! BoundaryPhysics: need to flip sign of velocity vectors for Reflect*()
  void ReflectInnerX1(Real time, Real dt,
                      int il, int jl, int ju, int kl, int ku, int ngh) override;
  void ReflectOuterX1(Real time, Real dt,
                      int iu, int jl, int ju, int kl, int ku, int ngh) override;
  void ReflectInnerX2(Real time, Real dt,
                      int il, int iu, int jl, int kl, int ku, int ngh) override;
  void ReflectOuterX2(Real time, Real dt,
                      int il, int iu, int ju, int kl, int ku, int ngh) override;
  void ReflectInnerX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int kl, int ngh) override;
  void ReflectOuterX3(Real time, Real dt,
                      int il, int iu, int jl, int ju, int ku, int ngh) override;
  //!@}

  //protected:
 private:
  void SetBoundarySameLevel(Real *buf, const NeighborBlock& nb) override;
  //! CosmicRay is a unique cell-centered variable because of the relationship between
  //! CosmicRayBoundaryQuantity::cr_cons and CosmicRayBoundaryQuantity::cr_prim.
  CosmicRayBoundaryQuantity cr_type_;
  int LoadFluxBoundaryBufferSameLevel(Real *buf, const NeighborBlock& nb) final;
  void SetBoundaryFromCoarser(Real *buf, const NeighborBlock& nb) final;
  void SetBoundaryFromFiner(Real *buf, const NeighborBlock& nb) final;
  void PolarWedgeInnerX2( Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) final;
  void PolarWedgeOuterX2( Real time, Real dt, int il, int iu, int jl, int kl, int ku, int ngh) final;
};

#endif // BVALS_CC_COSMIC_RAY_BVALS_COSMIC_RAY_HPP_
