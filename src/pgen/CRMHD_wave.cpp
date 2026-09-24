//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file CRMHD_wave.cpp
//  \brief One-dimensional coupled CR-MHD wave eigenmodes
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
  Real vA, Cs, Ccr, Rew, Imw, c, u0, delta_rho_over_rho;
}



Real complex_AtimeB_Re(Real A_Re, Real A_Im, Real B_Re, Real B_Im){
  return (A_Re*B_Re - A_Im*B_Im);
}

Real complex_AtimeB_Im(Real A_Re, Real A_Im, Real B_Re, Real B_Im){
  return (A_Re*B_Im + A_Im*B_Re);
}

Real complex_AoverB_Re(Real A_Re, Real A_Im, Real B_Re, Real B_Im){
  return (A_Re*B_Re + A_Im*B_Im)/(B_Re*B_Re + B_Im*B_Im);
}

Real complex_AoverB_Im(Real A_Re, Real A_Im, Real B_Re, Real B_Im){
  return (A_Im*B_Re - A_Re*B_Im)/(B_Re*B_Re + B_Im*B_Im);
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
  vA     = pin->GetOrAddReal("problem","vA",0.1);
  Cs     = pin->GetOrAddReal("problem","Cs",0.01);
  Ccr    = pin->GetOrAddReal("problem","Ccr",100.0);
  u0     = pin->GetOrAddReal("problem","u0",1.0);
  Rew    = pin->GetOrAddReal("problem","Rew",1.0);
  Imw    = pin->GetOrAddReal("problem","Imw",1.0);
  delta_rho_over_rho = pin->GetOrAddReal("problem","delta_rho_over_rho",1e-3);
  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {  
  Real x,y,z;
  Real Vm = cr_pointer->Vm;

  Real rho = 10.0;
  Real Pcr   = 0.75*Ccr*Ccr*rho;
  Real Lx = (pmy_mesh->mesh_size.x1max - pmy_mesh->mesh_size.x1min);
  Real Ly = (pmy_mesh->mesh_size.x2max - pmy_mesh->mesh_size.x2min);
  Real Lz = (pmy_mesh->mesh_size.x3max - pmy_mesh->mesh_size.x3min);
  Real knumber  = 2.0*PI/Lx;

  Real gamma = peos->GetGamma();
  Real gm1   = gamma - 1.0;
  Real P   = Cs*Cs*rho/gamma;
  Real B = std::sqrt(rho)*vA;

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        x = pcoord->x1v(i);
        y = pcoord->x2v(j);
        z = pcoord->x3v(k);


        Real Re_delta_rho = delta_rho_over_rho*rho*std::cos(knumber*x);
        Real Im_delta_rho = delta_rho_over_rho*rho*std::sin(knumber*x);

        Real Re_delta_u   = ((Rew - u0/vA)*Re_delta_rho - Imw*Im_delta_rho)*vA/rho;
        Real Im_delta_u   = ((Rew - u0/vA)*Im_delta_rho + Imw*Re_delta_rho)*vA/rho;

        Real w2_Re = complex_AtimeB_Re(Rew - u0/vA, Imw, Rew - u0/vA, Imw);
        Real w2_Im = complex_AtimeB_Im(Rew - u0/vA, Imw, Rew - u0/vA, Imw);

        Real temp_Re = complex_AtimeB_Re(3.0*Rew + 1.0 + u0/vA, 3.0*Imw, w2_Re - Cs*Cs/vA/vA, w2_Im);
        Real temp_Im = complex_AtimeB_Im(3.0*Rew + 1.0 + u0/vA, 3.0*Imw, w2_Re - Cs*Cs/vA/vA, w2_Im);

        Real temp2_Re = complex_AtimeB_Re(temp_Re, temp_Im, Re_delta_rho, Im_delta_rho);
        Real temp2_Im = complex_AtimeB_Im(temp_Re, temp_Im, Re_delta_rho, Im_delta_rho);

        w2_Re = complex_AtimeB_Re(Rew, Imw, Rew, Imw);
        w2_Im = complex_AtimeB_Im(Rew, Imw, Rew, Imw);

        Real Re_delta_Fcr = complex_AoverB_Re(vA*vA*vA*temp2_Re, vA*vA*vA*temp2_Im, 1.0 - 3.0*w2_Re*vA*vA/Vm/Vm, -3.0*w2_Im*vA*vA/Vm/Vm);
        Real Im_delta_Fcr = complex_AoverB_Im(vA*vA*vA*temp2_Re, vA*vA*vA*temp2_Im, 1.0 - 3.0*w2_Re*vA*vA/Vm/Vm, -3.0*w2_Im*vA*vA/Vm/Vm);

        temp_Re = complex_AtimeB_Re(Re_delta_Fcr, Im_delta_Fcr, 1.0 + Rew*vA*vA/Vm/Vm*(u0/vA + 1.0), Imw*vA*vA/Vm/Vm*(u0/vA + 1.0))/vA;
        temp_Im = complex_AtimeB_Im(Re_delta_Fcr, Im_delta_Fcr, 1.0 + Rew*vA*vA/Vm/Vm*(u0/vA + 1.0), Imw*vA*vA/Vm/Vm*(u0/vA + 1.0))/vA;        

        Real Re_delta_Pcr = complex_AoverB_Re(temp_Re, temp_Im, 3.0*Rew + 1.0 + u0/vA, 3.0*Imw);
        Real Im_delta_Pcr = complex_AoverB_Im(temp_Re, temp_Im, 3.0*Rew + 1.0 + u0/vA, 3.0*Imw);        

        phydro->u(IDN,k,j,i) = rho + Re_delta_rho;
        phydro->u(IM1,k,j,i) = phydro->u(IDN,k,j,i)*(u0 + Re_delta_u);
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;

        if (NCRS > 0) {
          for (int CR_id=0; CR_id<NCRS; ++CR_id) {
            int cr_energy_id  = 4*CR_id;
            int cr_flux1_id   = cr_energy_id + 1;
            int cr_flux2_id   = cr_energy_id + 2;
            int cr_flux3_id   = cr_energy_id + 3;

            cr_pointer->cr_cons(cr_energy_id, k, j, i) = 3.0*(Pcr + Re_delta_Pcr);
            cr_pointer->cr_cons(cr_flux1_id,  k, j, i) = (4.0*Pcr*(vA+u0) + Re_delta_Fcr)/Vm;
            cr_pointer->cr_cons(cr_flux2_id,  k, j, i) = 0.0;
            cr_pointer->cr_cons(cr_flux3_id,  k, j, i) = 0.0;
          }

        }



        if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = (P+Cs*Cs*Re_delta_rho)/gm1 + 0.5*
              (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }

        if (MAGNETIC_FIELDS_ENABLED) {
            phydro->u(IEN,k,j,i) += 0.5*B*B;
        }

      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          pfield->b.x1f(k,j,i) = B;
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

AthenaArray<Real> w = phydro->w;
AthenaArray<Real> cr_prim = cr_pointer->cr_prim;
AthenaArray<Real> bcc  = pfield->bcc;

for(int ncr=0; ncr<NCRS; ++ncr){

  int cr_energy_id  = 4*ncr;
  int cr_flux1_id   = cr_energy_id + 1;
  int cr_flux2_id   = cr_energy_id + 2;
  int cr_flux3_id   = cr_energy_id + 3;

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

        user_out_var(8,k,j,i)  = bcc(IB1,k,j,i);
        user_out_var(9,k,j,i)  = bcc(IB2,k,j,i);
        user_out_var(10,k,j,i) = bcc(IB3,k,j,i);
      }
    }
  }
}
  return;
}
