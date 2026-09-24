//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file cosmic_ray.cpp
//! \brief implementation of functions in class CosmicRay

// C headers

// C++ headers
#include <algorithm>
#include <cstring>    // strcmp
#include <sstream>
#include <stdexcept>  // runtime_error
#include <string>
#include <vector>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../mesh/mesh.hpp"
#include "../reconstruct/reconstruction.hpp"
#include "cosmic_ray.hpp"
#include "cr_scattering/cr_scattering.hpp"
#include "srcterms/cosmic_ray_srcterms.hpp"

class CRScattering;
class CosmicRaySourceTerms;

//! constructor, initializes data structures and parameters
CosmicRay::CosmicRay(MeshBlock *pmb, ParameterInput *pin)  :
  pmy_block(pmb),
  Vm(pin->GetOrAddReal("cosmic_ray", "Vmax", 1e2)),
  Vm2(Vm*Vm),
  gamma_cr(pin->GetOrAddReal("cosmic_ray", "gamma_cr", FOUR_3RD)),
  pressure_anisotropy_flag(
      pin->GetOrAddBoolean("cosmic_ray", "pressure_anisotropy_flag", false)),
  cr_cons(NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  cr_cons1(NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  cr_cons_af_src(NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  cr_prim(NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  cr_prim1(NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  cr_prim_n(NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  cr_flux{{NCRVARS, pmb->ncells3, pmb->ncells2, pmb->ncells1+1},
          {NCRVARS, pmb->ncells3, pmb->ncells2+1, pmb->ncells1,
          (pmb->pmy_mesh->f2 ? AthenaArray<Real>::DataStatus::allocated :
          AthenaArray<Real>::DataStatus::empty)},
          {NCRVARS, pmb->ncells3+1, pmb->ncells2, pmb->ncells1,
          (pmb->pmy_mesh->f3 ? AthenaArray<Real>::DataStatus::allocated :
          AthenaArray<Real>::DataStatus::empty)}},
  coarse_cr_cons_(NCRVARS, pmb->ncc3, pmb->ncc2, pmb->ncc1,
          (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
           AthenaArray<Real>::DataStatus::empty)),
  coarse_cr_prim_(NCRVARS, pmb->ncc3, pmb->ncc2, pmb->ncc1,
          (pmb->pmy_mesh->multilevel ? AthenaArray<Real>::DataStatus::allocated :
           AthenaArray<Real>::DataStatus::empty)),
  deltaPcr(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_pL_array(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_pR_array(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_mL_array(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_mR_array(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_Lorentz_array(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),  
  deltaPcr_n(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_pL_array_n(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_pR_array_n(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_mL_array_n(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_mR_array_n(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  sigma_Lorentz_array_n(NCRS, pmb->ncells3, pmb->ncells2, pmb->ncells1), 
  Stage_I_delta_mom1(NSPECIES, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  Stage_I_delta_mom2(NSPECIES, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  Stage_I_delta_mom3(NSPECIES, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  Stage_I_vel1(NSPECIES, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  Stage_I_vel2(NSPECIES, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  Stage_I_vel3(NSPECIES, pmb->ncells3, pmb->ncells2, pmb->ncells1),
  crbvar(pmb, &cr_cons, &coarse_cr_cons_, cr_flux, CosmicRayBoundaryQuantity::cr_cons),
  crscat(this, pin), crsrc(this, pin),
  pco_(pmb->pcoord) {

  int nc1 = pmb->ncells1, nc2 = pmb->ncells2, nc3 = pmb->ncells3;

  if (!MAGNETIC_FIELDS_ENABLED) {
    std::stringstream msg;
    msg << "### FATAL ERROR in CosmicRay constructor" << std::endl
        << "The CR module requires magnetic fields. Reconfigure Athena++ with -b."
        << std::endl;
    ATHENA_ERROR(msg);
  }

  Mesh *pm = pmy_block->pmy_mesh;
  pmb->RegisterMeshBlockData(cr_cons);

  deltaPcr.ZeroClear();
  deltaPcr_n.ZeroClear();
  std::fill_n(sigma_Lorentz_array.data(), sigma_Lorentz_array.GetSize(), 1.0e8);
  std::fill_n(sigma_Lorentz_array_n.data(), sigma_Lorentz_array_n.GetSize(), 1.0e8);

  int xorder  = pmb->precon->xorder;
  cr_xorder = pin->GetOrAddInteger("cosmic_ray", "cr_xorder", xorder);
  if (cr_xorder > xorder)
    cr_xorder = xorder;

  solver_id = pin->GetOrAddInteger("cosmic_ray", "solver_id", 2);
  if (solver_id != 2) {
    std::stringstream msg;
    msg << "The CR module currently supports only solver_id = 2." << std::endl;
    ATHENA_ERROR(msg);
  }

  // Allocate optional cosmic-ray variable memory registers for the time integrator
  if (cr_xorder == 4) {
    // fourth-order cell-centered approximations
    cr_cons_cc.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
    cr_prim_cc.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
  }

  // If user-requested time integrator is type 3S*, allocate additional memory registers
  std::string integrator = pin->GetOrAddString("time", "integrator", "vl2");

  if (integrator == "ssprk5_4" || STS_ENABLED) {
    // future extension may add "int nregister" to Hydro class
    cr_cons2.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
  }

  // If STS RKL2, allocate additional memory registers
  if (STS_ENABLED) {
    std::string sts_integrator = pin->GetOrAddString("time", "sts_integrator", "rkl2");
    if (sts_integrator == "rkl2") {
      cr_cons0.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
      cr_cons_fl_div.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
    }
  }

  // "Enroll" in SMR/AMR by adding to vector of pointers in MeshRefinement class
  if (pm->multilevel) {
    refinement_idx = pmy_block->pmr->AddToRefinement(&cr_cons, &coarse_cr_cons_);
  }

  // enroll CosmicRayBoundaryVariable object
  crbvar.bvar_index = pmb->pbval->bvars.size();
  pmb->pbval->bvars.push_back(&crbvar);
  pmb->pbval->bvars_main_int.push_back(&crbvar);

  // Allocate memory for scratch arrays
  dt1_.NewAthenaArray(nc1);
  dt2_.NewAthenaArray(nc1);
  dt3_.NewAthenaArray(nc1);
  //dx_cr_prim_.NewAthenaArray(nc1);
  cr_prim_l_.NewAthenaArray(NCRVARS, nc1);
  cr_prim_r_.NewAthenaArray(NCRVARS, nc1);
  cr_prim_lb_.NewAthenaArray(NCRVARS, nc1);
  b_interface_l_.NewAthenaArray(3, nc1);
  b_interface_r_.NewAthenaArray(3, nc1);
  b_interface_lb_.NewAthenaArray(3, nc1);
  deltaPcr_l_.NewAthenaArray(NCRS, nc1);
  deltaPcr_r_.NewAthenaArray(NCRS, nc1);
  deltaPcr_lb_.NewAthenaArray(NCRS, nc1);
  x1face_area_.NewAthenaArray(nc1+1);

  if (pm->f2) {
    x2face_area_.NewAthenaArray(nc1);
    x2face_area_p1_.NewAthenaArray(nc1);
  }

  if (pm->f3) {
    x3face_area_.NewAthenaArray(nc1);
    x3face_area_p1_.NewAthenaArray(nc1);
  }

  cell_volume_.NewAthenaArray(nc1);
  crflx_.NewAthenaArray(NCRVARS, nc1);
  HLLE_aux_.NewAthenaArray(nHLLE_aux, nc1);

  // fourth-order integration scheme
  if (cr_xorder == 4) {
    // 4D scratch arrays
    cr_prim_l3d_.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
    cr_prim_r3d_.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
    scr1_nkji_.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
    scr2_nkji_.NewAthenaArray(NCRVARS, nc3, nc2, nc1);
    // store all face-centered mass fluxes (all 3x coordinate directions) from Hydro:

    // 1D scratch arrays
    laplacian_l_cr_fc_.NewAthenaArray(nc1);
    laplacian_r_cr_fc_.NewAthenaArray(nc1);
  }
}
