//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file CR_streaming.cpp
//  \brief One-dimensional CR streaming test with a triangular energy profile
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

  return;
}


void CRstInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
  const Real sigma0 = 1.0;
  const Real Vm = pmb->cr_pointer->Vm;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      const Real vspeed = (1.0+4.0*sigma0)/3.0/sigma0;
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        const int ghost_index = is-i;
        Real x = pco->x1v(ghost_index);
        const Real u = 0.0;
        x -= (4.0/3.0*u*time);
        const int sign = x > 0 ? 1 : -1;
        const Real cr_energy = 2.0+vspeed*(time+dt)-fabs(x);
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          const int cr_energy_id = 4*CR_id;
          const int cr_flux1_id = cr_energy_id + 1;
          const int cr_flux2_id = cr_energy_id + 2;
          const int cr_flux3_id = cr_energy_id + 3;
          cr_prim(cr_energy_id,k,j,ghost_index) = cr_energy;
          cr_prim(cr_flux1_id,k,j,ghost_index) =
              cr_energy*(vspeed*sign+4.0/3.0*u)/Vm;
          cr_prim(cr_flux2_id,k,j,ghost_index) = 0.0;
          cr_prim(cr_flux3_id,k,j,ghost_index) = 0.0;
        }

        prim(IDN,k,j,ghost_index) = 1.0;
        prim(IVX,k,j,ghost_index) = u;
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
  const Real sigma0 = 1.0;
  const Real Vm = pmb->cr_pointer->Vm;
  for (int k=ks; k<=ke; ++k) {
    for (int j=js; j<=je; ++j) {
      const Real vspeed = (1.0+4.0*sigma0)/3.0/sigma0;
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        const int ghost_index = ie+i;
        Real x = pco->x1v(ghost_index);
        const Real u = 0.0;
        x -= (4.0/3.0*u*time);
        const int sign = x > 0 ? 1 : -1;
        const Real cr_energy = 2.0+vspeed*(time+dt)-fabs(x);
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          const int cr_energy_id = 4*CR_id;
          const int cr_flux1_id = cr_energy_id + 1;
          const int cr_flux2_id = cr_energy_id + 2;
          const int cr_flux3_id = cr_energy_id + 3;
          cr_prim(cr_energy_id,k,j,ghost_index) = cr_energy;
          cr_prim(cr_flux1_id,k,j,ghost_index) =
              cr_energy*(vspeed*sign+4.0/3.0*u)/Vm;
          cr_prim(cr_flux2_id,k,j,ghost_index) = 0.0;
          cr_prim(cr_flux3_id,k,j,ghost_index) = 0.0;
        }

        prim(IDN,k,j,ghost_index) = 1.0;
        prim(IVX,k,j,ghost_index) = u;
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

            if (fabs(x) < 2.0) cr_pointer->cr_cons(cr_energy_id, k, j, i) = 2.0 - fabs(x);
            else cr_pointer->cr_cons(cr_energy_id, k, j, i) = 0.0;
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

AthenaArray<Real> &w = phydro->w;
AthenaArray<Real> &cr_prim = cr_pointer->cr_prim;
AthenaArray<Real> &bcc = pfield->bcc;
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {

      for (int i=is; i<=ie; i++){


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
