//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//  \brief Problem generator for turbulence generator
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



//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================


void CRstInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;

  for (int k=ks; k<=ke; ++k) {  
    for (int j=js; j<=je; ++j) {
      Real rhos  = cr_prim(cr_energy_id,k,j,is);
      Real rhos1 = cr_prim(cr_energy_id,k,j,is+1);
      Real v1s  = cr_prim(cr_flux1_id,k,j,is);
      Real v2s  = cr_prim(cr_flux2_id,k,j,is);
      Real v3s  = cr_prim(cr_flux3_id,k,j,is);
      Real xs  = pco->x1v(is);
      Real xs1 = pco->x1v(is+1);
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1v(is-i);
        Real rat = (xs*xs*rhos - xs1*xs1*rhos1)/(xs-xs1);
        rat = -1.0;
        Real xxEcr = xs*xs*rhos + rat*(x-xs);
        // cr_prim(cr_energy_id,k,j,is-i) = xxEcr/x/x;
        cr_prim(cr_energy_id,k,j,is-i) = rhos;
        cr_prim(cr_flux1_id,k,j,is-i) = 0.0;
        cr_prim(cr_flux2_id,k,j,is-i) = 0.0;//v2s;//*(x/x0);
        cr_prim(cr_flux3_id,k,j,is-i) = 0.0;

        prim(IDN,k,j,is-i) = 1.0;
        prim(IVX,k,j,is-i) = 0.0;
        prim(IVY,k,j,is-i) = 0.0;
        prim(IVZ,k,j,is-i) = 0.0;
        prim(IPR,k,j,is-i) = 1.0;
    }
  }     
  }
  return;
}



void CRstOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
                 Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;

  for (int k=ks; k<=ke; ++k) {  
    for (int j=js; j<=je; ++j) {
      Real rhoe  = cr_prim(cr_energy_id,k,j,ie);
      Real rhoe1 = cr_prim(cr_energy_id,k,j,ie-1);
      Real v1e  = cr_prim(cr_flux1_id,k,j,ie);
      Real v2e  = cr_prim(cr_flux2_id,k,j,ie);
      Real v3e  = cr_prim(cr_flux3_id,k,j,ie);
      Real xe  = pco->x1v(ie);
      Real xe1 = pco->x1v(ie-1);
#pragma omp simd
      for (int i=1; i<=ngh; ++i) {
        Real x = pco->x1v(ie+i);
        Real rat = (xe*xe*rhoe - xe1*xe1*rhoe1)/(xe-xe1);
        rat = -1.0;
        Real xxEcr = xe*xe*rhoe + rat*(x-xe);
        // cr_prim(cr_energy_id,k,j,ie+i) = xxEcr/x/x;
        cr_prim(cr_energy_id,k,j,ie+i) = rhoe;
        cr_prim(cr_flux1_id,k,j,ie+i) = 0.0;
        cr_prim(cr_flux2_id,k,j,ie+i) = 0.0; //v2s;//*(x/x0);
        cr_prim(cr_flux3_id,k,j,ie+i) = 0.0;

        prim(IDN,k,j,ie+i) = 1.0;
        prim(IVX,k,j,ie+i) = 0.0;
        prim(IVY,k,j,ie+i) = 0.0;
        prim(IVZ,k,j,ie+i) = 0.0;
        prim(IPR,k,j,ie+i) = 1.0;
    }
  }  
  }   

  return;
}



void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
	 AllocateUserOutputVariables(8);
 	 return;
}





void Mesh::InitUserMeshData(ParameterInput *pin) {

  if(mesh_bcs[BoundaryFace::inner_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(inner_x1, CRstInnerX1);
  }

  if(mesh_bcs[BoundaryFace::outer_x1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(outer_x1, CRstOuterX1);
  }
  beta   = pin->GetOrAddReal("problem","beta",1.0);

  theta  = pin->GetOrAddReal("problem","theta",90.0)*PI/180.0;
  phi    = pin->GetOrAddReal("problem","phi",0.0)*PI/180.0;

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  AthenaArray<Real> ax, az;
  az.NewAthenaArray(ncells3, ncells2, ncells1);  
  ax.NewAthenaArray(ncells3, ncells2, ncells1);


  Real x,y,z;
  Real P = 1.0;
  Real P_B, B0;
  Real Lx = (pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real Ly = (pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  Real Lz = (pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);  

  Real kx = 2.0*PI/Lx;
  Real ky = 2.0*PI/Ly;
  Real kz = 2.0*PI/Lz;

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

        phydro->u(IDN,k,j,i) = 1.0/x/x/x/x;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NCRS > 0) {
          for (int CR_id=0; CR_id<NCRS; ++CR_id) {
            int cr_energy_id  = 4*CR_id;
            int cr_flux1_id   = cr_energy_id + 1;
            int cr_flux2_id   = cr_energy_id + 2;
            int cr_flux3_id   = cr_energy_id + 3;

            // cr_pointer->cr_cons(cr_energy_id, k, j, i) = 10.0/x/x-1.0/x;
            if (x<1.0) cr_pointer->cr_cons(cr_energy_id, k, j, i) = 2.0;
            else cr_pointer->cr_cons(cr_energy_id, k, j, i) = 1e-10;
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
            phydro->u(IEN,k,j,i) += 0.5/x/x/x/x;
        }

      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          Real r = pcoord->x1f(i);
          pfield->b.x1f(k,j,i) = 1.0/r/r;
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

AthenaArray<Real> p = phydro->w;
AthenaArray<Real> cr_prim = cr_pointer->cr_prim;
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;
AthenaArray<Real> bcc  = pfield->bcc;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {

      for (int i=is; i<=ie; i++){


        user_out_var(0,k,j,i) = cr_prim(cr_energy_id,k,j,i);
        user_out_var(1,k,j,i) = cr_prim(cr_flux1_id,k,j,i);
        user_out_var(2,k,j,i) = cr_prim(cr_flux2_id,k,j,i);
        user_out_var(3,k,j,i) = cr_prim(cr_flux3_id,k,j,i);
        user_out_var(4,k,j,i) = bcc(IB1,k,j,i);
        user_out_var(5,k,j,i) = bcc(IB2,k,j,i);
        user_out_var(6,k,j,i) = bcc(IB3,k,j,i);
        user_out_var(7,k,j,i) = p(IDN,k,j,i);
 



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
