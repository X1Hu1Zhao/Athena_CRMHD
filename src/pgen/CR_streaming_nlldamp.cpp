//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file CR_streaming_nlldamp.cpp
//  \brief One-dimensional CR streaming test with nonlinear Landau damping
//

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"


namespace{
  Real beta, theta, phi;
}


static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);
static Real minmod  (const Real A, const Real B);

void CRstInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void CRstOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void MyCRScatteringCoefficient(
    CosmicRay *cr_pointer, MeshBlock *pmb,
    const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
    const AthenaArray<Real> &bcc,
    AthenaArray<Real> &sigma_pL, AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL, AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz, AthenaArray<Real> &deltaPcr,
    AthenaArray<Real> &cs_cr,
    int is, int ie, int js, int je, int ks, int ke);



namespace{

}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
  AllocateUserOutputVariables(11);
  return;
}





void Mesh::InitUserMeshData(ParameterInput *pin) {

  beta   = pin->GetOrAddReal("problem","beta",1.0);

  theta  = pin->GetOrAddReal("problem","theta",90.0)*PI/180.0;
  phi    = pin->GetOrAddReal("problem","phi",0.0)*PI/180.0;

  if (mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::inner_x1, CRstInnerX1);
  }

  if (mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(BoundaryFace::outer_x1, CRstOuterX1);
  }
  EnrollCRScatteringCoefficients(MyCRScatteringCoefficient);

  return;
}


void CRstInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  const Real Vm = pmb->cr_pointer->Vm;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      const Real vspeed = 4.0/3.0+0.5;
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        const int ghost_index = is-i;
        Real x = pco->x1v(ghost_index);
        const int sign = x > 0 ? 1 : -1;
        x -= vspeed*time*sign;
        const Real cr_energy = 1.0/(10.0 + 3.0*fabs(x));
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          const int cr_energy_id = 4*CR_id;
          const int cr_flux1_id = cr_energy_id + 1;
          const int cr_flux2_id = cr_energy_id + 2;
          const int cr_flux3_id = cr_energy_id + 3;
          cr_prim(cr_energy_id,k,j,ghost_index) = cr_energy;
          cr_prim(cr_flux1_id,k,j,ghost_index) = -vspeed*cr_energy/Vm;
          cr_prim(cr_flux2_id,k,j,ghost_index) = 0.0;
          cr_prim(cr_flux3_id,k,j,ghost_index) = 0.0;
        }

        prim(IDN,k,j,ghost_index) = 1.0;
        prim(IVX,k,j,ghost_index) = 0.0;
        prim(IVY,k,j,ghost_index) = 0.0;
        prim(IVZ,k,j,ghost_index) = 0.0;
        prim(IPR,k,j,ghost_index) = prim(IPR,k,j,is);
      }
    }
  }
  return;
}



void CRstOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  const Real Vm = pmb->cr_pointer->Vm;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      const Real vspeed = 4.0/3.0+0.5;
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        const int ghost_index = ie+i;
        Real x = pco->x1v(ghost_index);
        const int sign = x > 0 ? 1 : -1;
        x -= vspeed*time*sign;
        const Real cr_energy = 1.0/(10.0 + 3.0*fabs(x));
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          const int cr_energy_id = 4*CR_id;
          const int cr_flux1_id = cr_energy_id + 1;
          const int cr_flux2_id = cr_energy_id + 2;
          const int cr_flux3_id = cr_energy_id + 3;
          cr_prim(cr_energy_id,k,j,ghost_index) = cr_energy;
          cr_prim(cr_flux1_id,k,j,ghost_index) = vspeed*cr_energy/Vm;
          cr_prim(cr_flux2_id,k,j,ghost_index) = 0.0;
          cr_prim(cr_flux3_id,k,j,ghost_index) = 0.0;
        }

        prim(IDN,k,j,ghost_index) = 1.0;
        prim(IVX,k,j,ghost_index) = 0.0;
        prim(IVY,k,j,ghost_index) = 0.0;
        prim(IVZ,k,j,ghost_index) = 0.0;
        prim(IPR,k,j,ghost_index) = prim(IPR,k,j,ie);
      }
    }
  }
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real x,y,z;
  Real P = 1.0;
  Real P_B, B0;
  Real Lx = (pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real Ly = (pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);

  Real gm1 = peos->GetGamma() - 1.0;

  P_B = P/beta;
  // B0  = std::sqrt(2.0*P_B);
  B0 = 1.0;

  Real Ecr0 = 0.1;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        x = pcoord->x1v(i);
        y = pcoord->x2v(j);
        z = pcoord->x3v(k);


        phydro->u(IDN,k,j,i) = 1.0;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NCRS > 0) {
          for (int CR_id=0; CR_id<NCRS; ++CR_id) {
            int cr_energy_id  = 4*CR_id;
            int cr_flux1_id   = cr_energy_id + 1;
            int cr_flux2_id   = cr_energy_id + 2;
            int cr_flux3_id   = cr_energy_id + 3;

            cr_pointer->cr_cons(cr_energy_id, k, j, i) = 1.0/(1.0/Ecr0 + 3.0*fabs(x));
            cr_pointer->cr_cons(cr_flux1_id,  k, j, i) = 0.0;
            cr_pointer->cr_cons(cr_flux2_id,  k, j, i) = 0.0;
            cr_pointer->cr_cons(cr_flux3_id,  k, j, i) = 0.0;
          }



        }



        if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = P/gm1 +(0.5)*
              (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }

        if (MAGNETIC_FIELDS_ENABLED) {
            phydro->u(IEN,k,j,i) += 0.5*B0*B0;
        }

      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          pfield->b.x1f(k,j,i) = B0;
        }
      }
    }

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x2f(k,j,i) = 0.0;
        }
      }
    }

    for (int k=ks; k<=ke+1; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x3f(k,j,i) = 0.0;
        }
      }
    }
  }
  return;
}


//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================




void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin){
Real P_B;
int il, iu, jl, ju, kl, ku;

il = is - NGHOST;
iu = ie + NGHOST;


AthenaArray<Real> w = phydro->w;
AthenaArray<Real> cr_prim = cr_pointer->cr_prim;
AthenaArray<Real> bcc  = pfield->bcc;
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {

      for (int i=il; i<=iu; i++){


        user_out_var(0,k,j,i) = cr_prim(cr_energy_id,k,j,i);
        user_out_var(1,k,j,i) = cr_prim(cr_flux1_id,k,j,i);
        user_out_var(2,k,j,i) = cr_prim(cr_flux2_id,k,j,i);
        user_out_var(3,k,j,i) = cr_prim(cr_flux3_id,k,j,i);

        user_out_var(4,k,j,i) = w(IDN,k,j,i);
        user_out_var(5,k,j,i) = w(IVX,k,j,i);
        user_out_var(6,k,j,i) = w(IVY,k,j,i);
        user_out_var(7,k,j,i) = w(IVZ,k,j,i);

        user_out_var(8,k,j,i) = bcc(IB1,k,j,i);
        user_out_var(9,k,j,i) = bcc(IB2,k,j,i);
        user_out_var(10,k,j,i) = bcc(IB3,k,j,i);




      }
    }
  }
  return;
}







static Real limiter2(const Real A, const Real B)
{
  // van Leer slope limiter
  return vanleer(A,B);
  //return minmod(A,B);

  // monotonized central (MC) limiter
  // return minmod(2.0*minmod(A,B),0.5*(A+B));
}

static Real limiter4(const Real A, const Real B, const Real C, const Real D)
{
  return limiter2(limiter2(A,B),limiter2(C,D));
}

//----------------------------------------------------------------------------
// vanleer: van Leer slope limiter


static Real vanleer(const Real A, const Real B)
{
  if (A*B > 0) {
    return 2.0*A*B/(A+B);
  } else {
    return 0.0;
  }
}

//----------------------------------------------------------------------------
// minmod: minmod slope limiter


static Real minmod(const Real A, const Real B)
{
  if (A*B > 0) {
    if (A > 0) {
      return (((A)<(B))?(A):(B));
    } else {
      return (((A)>(B))?(A):(B));
    }
  } else {
    return 0.0;
  }
}



void MyCRScatteringCoefficient(
    CosmicRay *cr_pointer, MeshBlock *pmb,
    const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
    const AthenaArray<Real> &bcc,
    AthenaArray<Real> &sigma_pL, AthenaArray<Real> &sigma_pR,
    AthenaArray<Real> &sigma_mL, AthenaArray<Real> &sigma_mR,
    AthenaArray<Real> &sigma_Lorentz, AthenaArray<Real> &deltaPcr,
    AthenaArray<Real> &cs_cr,
    int is, int ie, int js, int je, int ks, int ke){

      const bool has_x2 = pmb->block_size.nx2 > 1;
      const bool has_x3 = pmb->block_size.nx3 > 1;
      const int gradient_il = pmb->is - 1;
      const int gradient_iu = pmb->ie + 1;
      const int gradient_jl = has_x2 ? pmb->js - 1 : pmb->js;
      const int gradient_ju = has_x2 ? pmb->je + 1 : pmb->je;
      const int gradient_kl = has_x3 ? pmb->ks - 1 : pmb->ks;
      const int gradient_ku = has_x3 ? pmb->ke + 1 : pmb->ke;

      for(int ncr = 0; ncr < NCRS; ++ncr){
            int cr_energy_id  = 4*ncr;
            int cr_flux1_id   = cr_energy_id + 1;
            int cr_flux2_id   = cr_energy_id + 2;
            int cr_flux3_id   = cr_energy_id + 3;
        for (int k=gradient_kl; k<=gradient_ku; ++k) {
          for (int j=gradient_jl; j<=gradient_ju; ++j) {
            for (int i=gradient_il; i<=gradient_iu; ++i) {

              Real dEcrdx = 0.5*(cr_prim(cr_energy_id,k,j,i+1) - cr_prim(cr_energy_id,k,j,i-1))/(pmb->pcoord->dx1f(i));
              Real dEcrdy = has_x2
                  ? 0.5*(cr_prim(cr_energy_id,k,j+1,i) - cr_prim(cr_energy_id,k,j-1,i))
                      /pmb->pcoord->dx2f(j) : 0.0;
              Real dEcrdz = has_x3
                  ? 0.5*(cr_prim(cr_energy_id,k+1,j,i) - cr_prim(cr_energy_id,k-1,j,i))
                      /pmb->pcoord->dx3f(k) : 0.0;

              Real bdotgradEcr = (bcc(IB1,k,j,i)*dEcrdx + bcc(IB2,k,j,i)*dEcrdy + bcc(IB3,k,j,i)*dEcrdz);
              Real B = sqrt(bcc(IB1,k,j,i)*bcc(IB1,k,j,i) + bcc(IB2,k,j,i)*bcc(IB2,k,j,i) + bcc(IB3,k,j,i)*bcc(IB3,k,j,i));

              sigma_pL(ncr,k,j,i) = 0.0;
              sigma_pR(ncr,k,j,i) = 0.0;
              sigma_mL(ncr,k,j,i) = 0.0;
              sigma_mR(ncr,k,j,i) = 0.0;
              if (bdotgradEcr < 0.0){
                sigma_pL(ncr,k,j,i) = sqrt(-bdotgradEcr/3.0/B);
                sigma_pR(ncr,k,j,i) = sigma_pL(ncr,k,j,i);
              } else{
                sigma_mL(ncr,k,j,i) = sqrt(bdotgradEcr/3.0/B);
                sigma_mR(ncr,k,j,i) = sigma_mL(ncr,k,j,i);
              }
              sigma_Lorentz(ncr,k,j,i) = 0.0;
              deltaPcr(ncr,k,j,i) = 0.0;
            }
          }
        }
        for (int k=ks; k<=ke; ++k) {
          const int source_k = k < gradient_kl ? gradient_kl
                             : (k > gradient_ku ? gradient_ku : k);
          for (int j=js; j<=je; ++j) {
            const int source_j = j < gradient_jl ? gradient_jl
                               : (j > gradient_ju ? gradient_ju : j);
            for (int i=is; i<=ie; ++i) {
              if (k >= gradient_kl && k <= gradient_ku
                  && j >= gradient_jl && j <= gradient_ju
                  && i >= gradient_il && i <= gradient_iu) continue;
              const int source_i = i < gradient_il ? gradient_il
                                 : (i > gradient_iu ? gradient_iu : i);
              sigma_pL(ncr,k,j,i) = sigma_pL(ncr,source_k,source_j,source_i);
              sigma_pR(ncr,k,j,i) = sigma_pR(ncr,source_k,source_j,source_i);
              sigma_mL(ncr,k,j,i) = sigma_mL(ncr,source_k,source_j,source_i);
              sigma_mR(ncr,k,j,i) = sigma_mR(ncr,source_k,source_j,source_i);
              sigma_Lorentz(ncr,k,j,i) = sigma_Lorentz(ncr,source_k,source_j,source_i);
              deltaPcr(ncr,k,j,i) = deltaPcr(ncr,source_k,source_j,source_i);
            }
          }
        }
      }
      return;
    }
