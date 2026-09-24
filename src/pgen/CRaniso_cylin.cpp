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
  Real beta;
  Real sigma_aniso, lambda, vA0;
}


static Real limiter2(const Real A, const Real B);
static Real limiter4(const Real A, const Real B, const Real C, const Real D);
static Real vanleer (const Real A, const Real B);
static Real minmod  (const Real A, const Real B);


double factorial(int n) {
  if (n <= 1) {
      return 1;
  } else {
      return n * factorial(n - 1);
  }
}

double besselJ(int n, double x, int terms) {
  double result = 0.0;
  for (int k = 0; k < terms; ++k) {
      result += pow(-1, k) / (factorial(k) * factorial(n + k)) * pow(x / 2, 2 * k + n);
  }
  return result;
}


void CRanisoInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);
void CRanisoOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh);



namespace{

}

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================


void MeshBlock::InitUserMeshBlockData(ParameterInput *pin)
{
 

	 AllocateUserOutputVariables(10);
 	 return;
}


void CRanisoInnerX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;

Real Vm = pmb->cr_pointer->Vm;
Real sqrt3 = std::sqrt(3.0);
for (int k=ks; k<=ke; ++k) {  
for (int j=js; j<=je; ++j) {
Real xs = pco->x1v(is);
Real Ecrs = cr_prim(cr_energy_id,k,j,is);
#pragma omp simd
for (int i=1; i<=ngh; ++i) {
Real x = pco->x1v(is-i);
cr_prim(cr_energy_id,k,j,is-i) = besselJ(0, 0.5*x*sqrt3, 100)*std::exp(-(time));
cr_prim(cr_flux1_id,k,j,is-i) = -4.0/3.0*(Ecrs - cr_prim(cr_energy_id,k,j,is-i))/(xs - x)/Vm;
cr_prim(cr_flux2_id,k,j,is-i) = 0.0;
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


void CRanisoOuterX1(MeshBlock *pmb,Coordinates *pco, AthenaArray<Real> &prim, AthenaArray<Real> &cr_prim, FaceField &b,
  Real time, Real dt, int is, int ie, int js, int je, int ks, int ke, int ngh)
{
const int CR_id = 0;
const int cr_energy_id = 4*CR_id;
const int cr_flux1_id = cr_energy_id + 1;
const int cr_flux2_id = cr_energy_id + 2;
const int cr_flux3_id = cr_energy_id + 3;

Real Vm = pmb->cr_pointer->Vm;
Real sqrt3 = std::sqrt(3.0);

for (int k=ks; k<=ke; ++k) {  
for (int j=js; j<=je; ++j) {
Real Ecre = cr_prim(cr_energy_id,k,j,ie);
Real xe = pco->x1v(ie);
#pragma omp simd
for (int i=1; i<=ngh; ++i) {
Real x = pco->x1v(ie+i);
cr_prim(cr_energy_id,k,j,ie+i) = besselJ(0, 0.5*x*sqrt3, 100)*std::exp(-(time));
cr_prim(cr_flux1_id,k,j,ie+i) = -4.0/3.0*(Ecre - cr_prim(cr_energy_id,k,j,ie+i))/(xe - x)/Vm;
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


void Mesh::InitUserMeshData(ParameterInput *pin) {

  beta   = pin->GetOrAddReal("problem","beta",1.0);
  sigma_aniso   = 0.25;
  lambda   = pin->GetOrAddReal("problem","lambda",1.0);
  vA0   = pin->GetOrAddReal("problem","vA0",1.0);


  if(mesh_bcs[INNER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(INNER_X1, CRanisoInnerX1);
  }

  if(mesh_bcs[OUTER_X1] == GetBoundaryFlag("user")) {
    EnrollUserBoundaryFunction(OUTER_X1, CRanisoOuterX1);
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {

  // AthenaArray<Real> ax, az;
  // az.NewAthenaArray(ncells3, ncells2, ncells1);  
  // ax.NewAthenaArray(ncells3, ncells2, ncells1);

  Real sqrt3 = std::sqrt(3.0);

  Real R, theta, phi, inv_R;
  Real inv_Vm = 1.0/cr_pointer->Vm;
  Real P = 1.0;
  Real P_B, B0;
  Real Lx = (pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real Ly = (pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  Real Lz = (pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);  

  Real kx = 2.0*PI/Lx;
  Real ky = 2.0*PI/Ly;
  Real kz = 2.0*PI/Lz;

  Real gm1 = peos->GetGamma() - 1.0;

  // B0  = std::sqrt(2.0*P_B);


  // for (int k=ks; k<=ke+1; ++k) {
  //   for (int j=js; j<=je+1; ++j) {
  //     for (int i=is; i<=ie+1; ++i) {
  //       az(k,j,i) = B0*(1.0+0.5*std::cos(kx*pcoord->x1f(i)))*pcoord->x2f(j);
  //       ax(k,j,i) = -kx*B0*0.25*pcoord->x2f(j)*pcoord->x3f(k)*std::sin(kx*pcoord->x1f(i));
  //     }
  //   }
  // }


  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is-NGHOST; i<=ie+1+NGHOST; i++) {
          Real R = pcoord->x1f(i);
          pfield->b.x1f(k,j,i) = 1.0/R;
        }
      }
    }

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is-NGHOST; i<=ie+NGHOST; i++) {
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


  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is-NGHOST; i<=ie+NGHOST; i++) {
      const Real& b1_i   = pfield->b.x1f(k,j,i  );
      const Real& b1_ip1 = pfield->b.x1f(k,j,i+1);
      const Real& b2_j   = pfield->b.x2f(k,j  ,i);
      const Real& b2_jp1 = pfield->b.x2f(k,j+1,i);
      const Real& b3_k   = pfield->b.x3f(k  ,j,i);
      const Real& b3_kp1 = pfield->b.x3f(k+1,j,i);

      Real& bcc1 = pfield->bcc(IB1,k,j,i);
      Real& bcc2 = pfield->bcc(IB2,k,j,i);
      Real& bcc3 = pfield->bcc(IB3,k,j,i);

      const Real& x1f_i  = pcoord->x1f(i);
      const Real& x1f_ip = pcoord->x1f(i+1);
      const Real& x1v_i  = pcoord->x1v(i);
      const Real& dx1_i  = pcoord->dx1f(i);
      Real lw=(x1f_ip-x1v_i)/dx1_i;
      Real rw=(x1v_i -x1f_i)/dx1_i;
      bcc1 = lw*b1_i + rw*b1_ip1;
      const Real& x2f_j  = pcoord->x2f(j);
      const Real& x2f_jp = pcoord->x2f(j+1);
      const Real& x2v_j  = pcoord->x2v(j);
      const Real& dx2_j  = pcoord->dx2f(j);
      lw=(x2f_jp-x2v_j)/dx2_j;
      rw=(x2v_j -x2f_j)/dx2_j;
      bcc2 = lw*b2_j + rw*b2_jp1;

      const Real& x3f_k  = pcoord->x3f(k);
      const Real& x3f_kp = pcoord->x3f(k+1);
      const Real& x3v_k  = pcoord->x3v(k);
      const Real& dx3_k  = pcoord->dx3f(k);
      lw=(x3f_kp-x3v_k)/dx3_k;
      rw=(x3v_k -x3f_k)/dx3_k;
      bcc3 = lw*b3_k + rw*b3_kp1;
    }}}


      for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is-NGHOST; i<=ie+NGHOST; i++) {

        R = pcoord->x1v(i);
    
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

            cr_pointer->cr_cons(cr_energy_id, k, j, i) = besselJ(0, 0.5*R*sqrt3, 100);
            cr_pointer->cr_cons(cr_flux1_id,  k, j, i) = 2.0/sqrt3*besselJ(1, 0.5*R*sqrt3, 100)*inv_Vm;
            cr_pointer->cr_cons(cr_flux2_id,  k, j, i) = 0.0;
            cr_pointer->cr_cons(cr_flux3_id,  k, j, i) = 0.0;
          }

        }

       P_B = 0.5*(SQR(pfield->bcc(IB1,k,j,i)) + SQR(pfield->bcc(IB2,k,j,i)) + SQR(pfield->bcc(IB3,k,j,i)));

        if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = P/gm1 +(0.5)*
              (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }

        if (MAGNETIC_FIELDS_ENABLED) {
            phydro->u(IEN,k,j,i) += P_B;
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
AthenaArray<Real> bf  = pfield->b.x1f;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {

      for (int i=is-NGHOST; i<=ie+NGHOST; i++){


        user_out_var(0,k,j,i) = cr_prim(cr_energy_id,k,j,i);
        user_out_var(1,k,j,i) = cr_prim(cr_flux1_id,k,j,i);
        user_out_var(2,k,j,i) = cr_prim(cr_flux2_id,k,j,i);
        user_out_var(3,k,j,i) = cr_prim(cr_flux3_id,k,j,i);
        user_out_var(4,k,j,i) = bcc(IB1,k,j,i);
        user_out_var(5,k,j,i) = bcc(IB2,k,j,i);
        user_out_var(6,k,j,i) = bcc(IB3,k,j,i);
        user_out_var(7,k,j,i) = p(IDN,k,j,i); 
        user_out_var(8,k,j,i) = p(IVX,k,j,i);           
        user_out_var(9,k,j,i) = p(IPR,k,j,i);  



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
