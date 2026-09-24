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


  // for(int i=0;i<piece_number;i++){
  //   Lambda[i] = Lambda0*Cooling[i];
  // }


  // TEFInit(T, Lambda, alpha, TEF, piece_number);


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

  AthenaArray<Real> az;
  az.NewAthenaArray(ncells2, ncells1);  


  P_B =P/beta;
  B0 = 1.0;

    for (int j=js; j<=je+1; ++j) {
      for (int i=is; i<=ie+1; ++i) {
        az(j,i) = -std::sqrt(pcoord->x1f(i)*pcoord->x1f(i) + pcoord->x2f(j)*pcoord->x2f(j));
      }
    }

  Real gm1 = peos->GetGamma() - 1.0;


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
            int cr_energy_id = 4*CR_id;

            Real phi = std::atan2(y,x);

            if (((x*x + y*y)>0.25)&&((x*x + y*y)<0.49)&&(fabs(phi)<PI/12.0)){
              cr_pointer->cr_cons(cr_energy_id, k, j, i) = 12.0;
            } else{
              cr_pointer->cr_cons(cr_energy_id, k, j, i) = 10.0;
            }

          }

        }



        if (NON_BAROTROPIC_EOS) {
            phydro->u(IEN,k,j,i) = P/gm1 +(0.5)*
              (SQR(phydro->u(IM1,k,j,i)) + SQR(phydro->u(IM2,k,j,i))
               + SQR(phydro->u(IM3,k,j,i)))/phydro->u(IDN,k,j,i);
        }

      }
    }
  }

  if (MAGNETIC_FIELDS_ENABLED) {

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie+1; i++) {
          pfield->b.x1f(k,j,i) = (az(j+1,i) - az(j,i))/pcoord->dx2f(j);
          // pfield->b.x1f(k,j,i) = -B0*pcoord->x2f(j)/std::sqrt(SQR(pcoord->x1f(i))+SQR(pcoord->x2f(j)));
        }
      }
    }

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
          // pfield->b.x2f(k,j,i) = B0*pcoord->x1f(i)/std::sqrt(SQR(pcoord->x1f(i))+SQR(pcoord->x2f(j)));          
        }
      }
    }

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je+1; j++) {
        for (int i=is; i<=ie; i++) {
          // pfield->b.x2f(k,j,i) = (az(j,i) - az(j,i+1))/pcoord->dx1f(i);
          pfield->b.x3f(k,j,i) = 0.0;          
        }
      }
    }

    for (int k=ks; k<=ke; k++) {
      for (int j=js; j<=je; j++) {
        for (int i=is; i<=ie; i++) {
          phydro->u(IEN,k,j,i) +=
              0.5*(SQR(0.5*(pfield->b.x1f(k,j,i) + pfield->b.x1f(k,j,i+1))) +
                   SQR(0.5*(pfield->b.x2f(k,j,i) + pfield->b.x2f(k,j+1,i))) +
                   SQR(0.5*(pfield->b.x3f(k,j,i) + pfield->b.x3f(k+1,j,i))));
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
        user_out_var(7,k,j,i)   = w(IDN,k,j,i);       
        user_out_var(8,k,j,i)   = w(IVX,k,j,i);
        user_out_var(9,k,j,i)   = w(IVY,k,j,i);
        user_out_var(10,k,j,i)  = w(IVZ,k,j,i);


      }
    }
  }
  return;
}

