//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file backwardEuler_integrator.cpp
//! \brief Backward Euler CR-scattering time integrators.

// C++ headers
#include <algorithm>   // min,max
#include <cstring>    // strcmp
#include <iostream>   // endl
#include <limits>
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../defs.hpp"
#include "../../mesh/mesh.hpp"
#include "../cr_frame_transform.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"




void CRScattering::BackwardEulerFeedback(const int stage,
      const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons) {

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  const AthenaArray<Real> &cr_cons_af_src = cr_pointer->cr_cons_af_src;
  const int gas_id = 0;


    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
          int cr_energy_id = 4*CR_id;
          int cr_flux1_id = cr_energy_id + 1;
          int cr_flux2_id = cr_energy_id + 2;
          int cr_flux3_id = cr_energy_id + 3;
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            const Real &CR_energy = cr_cons_af_src(cr_energy_id, k, j, i);
            const Real &CR_flux1 = cr_cons_af_src(cr_flux1_id, k, j, i);
            const Real &CR_flux2 = cr_cons_af_src(cr_flux2_id, k, j, i);
            const Real &CR_flux3 = cr_cons_af_src(cr_flux3_id, k, j, i);

            dens_prim(Fluid_id, i) = CR_energy;
            mom1_prim(Fluid_id, i) = CR_flux1;
            mom2_prim(Fluid_id, i) = CR_flux2;
            mom3_prim(Fluid_id, i) = CR_flux3;

            sigma_p(CR_id, i) = cr_pointer->sigma_pL_array(CR_id, k, j, i)
                              + cr_pointer->sigma_pR_array(CR_id, k, j, i);
            sigma_m(CR_id, i) = cr_pointer->sigma_mL_array(CR_id, k, j, i)
                              + cr_pointer->sigma_mR_array(CR_id, k, j, i);
            Omega(CR_id, i) = cr_pointer->sigma_Lorentz_array(CR_id, k, j, i);
          }
        }
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
        // Alias the primitives of gas at current stage
        const Real &gas_rho  = w(IDN, k, j, i);


        dens_prim(gas_id,i) = gas_rho;
        inv_gas_rho(i) = 1.0/gas_rho;

        const Real Bx = bcc(IB1,k,j,i);
        const Real By = bcc(IB2,k,j,i);
        const Real Bz = bcc(IB3,k,j,i);
        const Real B2 = Bx*Bx + By*By + Bz*Bz;
        const Real B = std::sqrt(B2);

        vA(i) = std::sqrt(B2*inv_gas_rho(i));

        if (B2>TINY_NUMBER) {
          const cr_frame::MagneticFrame frame =
              cr_frame::BuildMagneticFrame(Bx, By, Bz, B, SQR(Bx) + SQR(By));
          cosXY(i) = frame.cos_xy;
          sinXY(i) = frame.sin_xy;
          cosZ(i) = frame.cos_z;
          sinZ(i) = frame.sin_z;
        }

        const Real &gas_vel1 = w(IVX, k, j, i);
        const Real &gas_vel2 = w(IVY, k, j, i);
        const Real &gas_vel3 = w(IVZ, k, j, i);

        mom1_prim(gas_id, i) = gas_vel1;
        mom2_prim(gas_id, i) = gas_vel2;
        mom3_prim(gas_id, i) = gas_vel3;
      }

        RotateArraysToFieldAligned(
            cosXY, sinXY, cosZ, sinZ, mom1_prim, mom2_prim, mom3_prim);

        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
#pragma omp simd
        for(int i=is; i<=ie; ++i){
        // <<Most general implicit integrator>>
        Real coeff1 = 1.0 + (sigma_p(CR_id,i) + sigma_m(CR_id,i))*Vm2*dt;
        Real coeff2 = 1.0 + Vm2*Vm2*SQR(Omega(CR_id,i)*dt);
        Real coeff3 = (vA(i)*(sigma_p(CR_id,i) - sigma_m(CR_id,i)) + mom1_prim(gas_id,i)*(sigma_p(CR_id,i) + sigma_m(CR_id,i)))*Vm*dt;
        Real coeff4 = mom2_prim(gas_id,i)*Vm*Omega(CR_id,i)*dt;
        Real coeff5 = mom3_prim(gas_id,i)*Vm*Omega(CR_id,i)*dt;
        Real coeff6;

        if ((sigma_p(CR_id,i)+sigma_m(CR_id,i))>TINY_NUMBER)
                coeff6 = (SQR(coeff3*inv_Vm/dt) + 4.0*vA(i)*vA(i)*sigma_p(CR_id,i)*sigma_m(CR_id,i)*coeff1)/(sigma_p(CR_id,i)+sigma_m(CR_id,i));
        else
                coeff6 = 0.0;

        Real coeff0 = gamma_cr*coeff2*dt*coeff6 - coeff1*(coeff2 - gamma_cr*(SQR(coeff4) + SQR(coeff5)));
        coeff0 = 1.0/coeff0;

        delta_Ecr(Fluid_id,i)   = coeff0*(-coeff1*coeff2*dens_prim(Fluid_id,i) + coeff2*coeff3*mom1_prim(Fluid_id,i) + coeff1*(coeff5 + coeff4*Vm2*Omega(CR_id,i)*dt)*mom2_prim(Fluid_id,i)
                + coeff1*(coeff5*Vm2*dt*Omega(CR_id,i) - coeff4)*mom3_prim(Fluid_id,i)) - dens_prim(Fluid_id,i);

        delta_mom1(Fluid_id,i)  = coeff0*(-gamma_cr*coeff3*coeff2*dens_prim(Fluid_id,i) + (gamma_cr*coeff2*dt*coeff6
                      - (coeff2 - gamma_cr*(SQR(coeff4) + SQR(coeff5))))*mom1_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff5 + coeff4*Vm2*dt*Omega(CR_id,i))*mom2_prim(Fluid_id,i)
                + gamma_cr*coeff3*(coeff5*Vm2*Omega(CR_id,i)*dt - coeff4)*mom3_prim(Fluid_id,i)) - mom1_prim(Fluid_id,i);

        delta_mom2(Fluid_id,i)  = coeff0*(-gamma_cr*coeff1*(coeff5 - coeff4*Vm2*dt*Omega(CR_id,i))*dens_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff4*Vm2*dt*Omega(CR_id,i) - coeff5)*mom1_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6 + coeff1*(gamma_cr*SQR(coeff4)-1.0))*mom2_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6*Vm2*dt*Omega(CR_id,i) + coeff1*(gamma_cr*coeff4*coeff5-Vm2*dt*Omega(CR_id,i)))*mom3_prim(Fluid_id,i))
                - mom2_prim(Fluid_id,i);

        delta_mom3(Fluid_id,i)  = coeff0*(-gamma_cr*coeff1*(coeff4 + coeff5*Vm2*dt*Omega(CR_id,i))*dens_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff4 + coeff5*Vm2*dt*Omega(CR_id,i))*mom1_prim(Fluid_id,i)
                - (gamma_cr*coeff2*dt*coeff6*Vm2*dt*Omega(CR_id,i) - coeff1*(gamma_cr*coeff4*coeff5+Vm2*dt*Omega(CR_id,i)))*mom2_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6 + coeff1*(gamma_cr*SQR(coeff5)-1.0))*mom3_prim(Fluid_id,i))
                - mom3_prim(Fluid_id,i);

        }
        }

        RotateArraysToLabFrame(
            cosXY, sinXY, cosZ, sinZ, delta_mom1, delta_mom2, delta_mom3);

        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
          int cr_energy_id = 4*CR_id;
          int cr_flux1_id = cr_energy_id + 1;
          int cr_flux2_id = cr_energy_id + 2;
          int cr_flux3_id = cr_energy_id + 3;
#pragma omp simd
      for(int i=is; i<=ie; ++i){
        // Alias the conserves of gas
        Real &gas_mom1 = u(IM1, k, j, i);
        Real &gas_mom2 = u(IM2, k, j, i);
        Real &gas_mom3 = u(IM3, k, j, i);

        Real &Ecr  = cr_cons(cr_energy_id, k, j, i);
        Real &Fcr1 = cr_cons(cr_flux1_id, k, j, i);
        Real &Fcr2 = cr_cons(cr_flux2_id, k, j, i);
        Real &Fcr3 = cr_cons(cr_flux3_id, k, j, i);


        Ecr   += delta_Ecr(Fluid_id,i);
        Fcr1  += delta_mom1(Fluid_id,i);
        Fcr2  += delta_mom2(Fluid_id,i);
        Fcr3  += delta_mom3(Fluid_id,i);



        Real &gas_energy  = u(IEN, k, j, i);

        gas_energy  -= delta_Ecr(Fluid_id,i);
        gas_mom1 -= (inv_Vm*delta_mom1(Fluid_id,i));
        gas_mom2 -= (inv_Vm*delta_mom2(Fluid_id,i));
        gas_mom3 -= (inv_Vm*delta_mom3(Fluid_id,i));
           }
        }
      }
    }
  return;
}




void CRScattering::BackwardEulerNoEnergyFeedback(const int stage,
      const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons) {

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  const AthenaArray<Real> &cr_cons_af_src = cr_pointer->cr_cons_af_src;
  const int gas_id = 0;


    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
          int cr_energy_id = 4*CR_id;
          int cr_flux1_id = cr_energy_id + 1;
          int cr_flux2_id = cr_energy_id + 2;
          int cr_flux3_id = cr_energy_id + 3;
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            const Real &CR_energy = cr_cons_af_src(cr_energy_id, k, j, i);
            const Real &CR_flux1 = cr_cons_af_src(cr_flux1_id, k, j, i);
            const Real &CR_flux2 = cr_cons_af_src(cr_flux2_id, k, j, i);
            const Real &CR_flux3 = cr_cons_af_src(cr_flux3_id, k, j, i);

            dens_prim(Fluid_id, i) = CR_energy;
            mom1_prim(Fluid_id, i) = CR_flux1;
            mom2_prim(Fluid_id, i) = CR_flux2;
            mom3_prim(Fluid_id, i) = CR_flux3;

            sigma_p(CR_id, i) = cr_pointer->sigma_pL_array(CR_id, k, j, i)
                              + cr_pointer->sigma_pR_array(CR_id, k, j, i);
            sigma_m(CR_id, i) = cr_pointer->sigma_mL_array(CR_id, k, j, i)
                              + cr_pointer->sigma_mR_array(CR_id, k, j, i);
            Omega(CR_id, i) = cr_pointer->sigma_Lorentz_array(CR_id, k, j, i);
          }
        }
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
        // Alias the primitives of gas at current stage
        const Real &gas_rho  = w(IDN, k, j, i);


        dens_prim(gas_id,i) = gas_rho;
        inv_gas_rho(i) = 1.0/gas_rho;

        const Real Bx = bcc(IB1,k,j,i);
        const Real By = bcc(IB2,k,j,i);
        const Real Bz = bcc(IB3,k,j,i);
        const Real B2 = Bx*Bx + By*By + Bz*Bz;
        const Real B = std::sqrt(B2);

        vA(i) = std::sqrt(B2*inv_gas_rho(i));

        if (B2>TINY_NUMBER) {
          const cr_frame::MagneticFrame frame =
              cr_frame::BuildMagneticFrame(Bx, By, Bz, B, SQR(Bx) + SQR(By));
          cosXY(i) = frame.cos_xy;
          sinXY(i) = frame.sin_xy;
          cosZ(i) = frame.cos_z;
          sinZ(i) = frame.sin_z;
        }

        const Real &gas_vel1 = w(IVX, k, j, i);
        const Real &gas_vel2 = w(IVY, k, j, i);
        const Real &gas_vel3 = w(IVZ, k, j, i);

        mom1_prim(gas_id, i) = gas_vel1;
        mom2_prim(gas_id, i) = gas_vel2;
        mom3_prim(gas_id, i) = gas_vel3;
      }

        RotateArraysToFieldAligned(
            cosXY, sinXY, cosZ, sinZ, mom1_prim, mom2_prim, mom3_prim);

        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
#pragma omp simd
        for(int i=is; i<=ie; ++i){
        // <<Most general implicit integrator>>
        Real coeff1 = 1.0 + (sigma_p(CR_id,i) + sigma_m(CR_id,i))*Vm2*dt;
        Real coeff2 = 1.0 + Vm2*Vm2*SQR(Omega(CR_id,i)*dt);
        Real coeff3 = (vA(i)*(sigma_p(CR_id,i) - sigma_m(CR_id,i)) + mom1_prim(gas_id,i)*(sigma_p(CR_id,i) + sigma_m(CR_id,i)))*Vm*dt;
        Real coeff4 = mom2_prim(gas_id,i)*Vm*Omega(CR_id,i)*dt;
        Real coeff5 = mom3_prim(gas_id,i)*Vm*Omega(CR_id,i)*dt;
        Real coeff6;

        if ((sigma_p(CR_id,i)+sigma_m(CR_id,i))>TINY_NUMBER)
                coeff6 = (SQR(coeff3*inv_Vm/dt) + 4.0*vA(i)*vA(i)*sigma_p(CR_id,i)*sigma_m(CR_id,i)*coeff1)/(sigma_p(CR_id,i)+sigma_m(CR_id,i));
        else
                coeff6 = 0.0;

        Real coeff0 = gamma_cr*coeff2*dt*coeff6 - coeff1*(coeff2 - gamma_cr*(SQR(coeff4) + SQR(coeff5)));
        coeff0 = 1.0/coeff0;

        delta_Ecr(Fluid_id,i)   = coeff0*(-coeff1*coeff2*dens_prim(Fluid_id,i) + coeff2*coeff3*mom1_prim(Fluid_id,i) + coeff1*(coeff5 + coeff4*Vm2*Omega(CR_id,i)*dt)*mom2_prim(Fluid_id,i)
                + coeff1*(coeff5*Vm2*dt*Omega(CR_id,i) - coeff4)*mom3_prim(Fluid_id,i)) - dens_prim(Fluid_id,i);

        delta_mom1(Fluid_id,i)  = coeff0*(-gamma_cr*coeff3*coeff2*dens_prim(Fluid_id,i) + (gamma_cr*coeff2*dt*coeff6
                      - (coeff2 - gamma_cr*(SQR(coeff4) + SQR(coeff5))))*mom1_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff5 + coeff4*Vm2*dt*Omega(CR_id,i))*mom2_prim(Fluid_id,i)
                + gamma_cr*coeff3*(coeff5*Vm2*Omega(CR_id,i)*dt - coeff4)*mom3_prim(Fluid_id,i)) - mom1_prim(Fluid_id,i);

        delta_mom2(Fluid_id,i)  = coeff0*(-gamma_cr*coeff1*(coeff5 - coeff4*Vm2*dt*Omega(CR_id,i))*dens_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff4*Vm2*dt*Omega(CR_id,i) - coeff5)*mom1_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6 + coeff1*(gamma_cr*SQR(coeff4)-1.0))*mom2_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6*Vm2*dt*Omega(CR_id,i) + coeff1*(gamma_cr*coeff4*coeff5-Vm2*dt*Omega(CR_id,i)))*mom3_prim(Fluid_id,i))
                - mom2_prim(Fluid_id,i);

        delta_mom3(Fluid_id,i)  = coeff0*(-gamma_cr*coeff1*(coeff4 + coeff5*Vm2*dt*Omega(CR_id,i))*dens_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff4 + coeff5*Vm2*dt*Omega(CR_id,i))*mom1_prim(Fluid_id,i)
                - (gamma_cr*coeff2*dt*coeff6*Vm2*dt*Omega(CR_id,i) - coeff1*(gamma_cr*coeff4*coeff5+Vm2*dt*Omega(CR_id,i)))*mom2_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6 + coeff1*(gamma_cr*SQR(coeff5)-1.0))*mom3_prim(Fluid_id,i))
                - mom3_prim(Fluid_id,i);

        }
        }

        RotateArraysToLabFrame(
            cosXY, sinXY, cosZ, sinZ, delta_mom1, delta_mom2, delta_mom3);

        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
          int cr_energy_id = 4*CR_id;
          int cr_flux1_id = cr_energy_id + 1;
          int cr_flux2_id = cr_energy_id + 2;
          int cr_flux3_id = cr_energy_id + 3;
#pragma omp simd
      for(int i=is; i<=ie; ++i){
        // Alias the conserves of gas
        Real &gas_mom1 = u(IM1, k, j, i);
        Real &gas_mom2 = u(IM2, k, j, i);
        Real &gas_mom3 = u(IM3, k, j, i);

        Real &Ecr  = cr_cons(cr_energy_id, k, j, i);
        Real &Fcr1 = cr_cons(cr_flux1_id, k, j, i);
        Real &Fcr2 = cr_cons(cr_flux2_id, k, j, i);
        Real &Fcr3 = cr_cons(cr_flux3_id, k, j, i);


        Ecr   += delta_Ecr(Fluid_id,i);
        Fcr1  += delta_mom1(Fluid_id,i);
        Fcr2  += delta_mom2(Fluid_id,i);
        Fcr3  += delta_mom3(Fluid_id,i);



        Real &gas_energy  = u(IEN, k, j, i);

        Real Ek_previous = 0.5*(SQR(gas_mom1) + SQR(gas_mom2) + SQR(gas_mom3))*inv_gas_rho(i);

        gas_mom1 -= (inv_Vm*delta_mom1(Fluid_id,i));
        gas_mom2 -= (inv_Vm*delta_mom2(Fluid_id,i));
        gas_mom3 -= (inv_Vm*delta_mom3(Fluid_id,i));

        Real Ek_current  = 0.5*(SQR(gas_mom1) + SQR(gas_mom2) + SQR(gas_mom3))*inv_gas_rho(i);

        gas_energy += (Ek_current - Ek_previous);
          }
        }
      }
    }
  return;
}





void CRScattering::BackwardEulerNoFeedback(const int stage,
      const Real dt,
      const AthenaArray<Real> &w, const AthenaArray<Real> &cr_prim,
      const AthenaArray<Real> &bcc, AthenaArray<Real> &u, AthenaArray<Real> &cr_cons) {

  int is = pmb_->is; int js = pmb_->js; int ks = pmb_->ks;
  int ie = pmb_->ie; int je = pmb_->je; int ke = pmb_->ke;

  const AthenaArray<Real> &cr_cons_af_src = cr_pointer->cr_cons_af_src;
  const int gas_id = 0;


    for (int k=ks; k<=ke; ++k) {
      for (int j=js; j<=je; ++j) {
        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
          int cr_energy_id = 4*CR_id;
          int cr_flux1_id = cr_energy_id + 1;
          int cr_flux2_id = cr_energy_id + 2;
          int cr_flux3_id = cr_energy_id + 3;
#pragma omp simd
          for (int i=is; i<=ie; ++i) {
            const Real &CR_energy = cr_cons_af_src(cr_energy_id, k, j, i);
            const Real &CR_flux1 = cr_cons_af_src(cr_flux1_id, k, j, i);
            const Real &CR_flux2 = cr_cons_af_src(cr_flux2_id, k, j, i);
            const Real &CR_flux3 = cr_cons_af_src(cr_flux3_id, k, j, i);

            dens_prim(Fluid_id, i) = CR_energy;
            mom1_prim(Fluid_id, i) = CR_flux1;
            mom2_prim(Fluid_id, i) = CR_flux2;
            mom3_prim(Fluid_id, i) = CR_flux3;

            sigma_p(CR_id, i) = cr_pointer->sigma_pL_array(CR_id, k, j, i)
                              + cr_pointer->sigma_pR_array(CR_id, k, j, i);
            sigma_m(CR_id, i) = cr_pointer->sigma_mL_array(CR_id, k, j, i)
                              + cr_pointer->sigma_mR_array(CR_id, k, j, i);
            Omega(CR_id, i) = cr_pointer->sigma_Lorentz_array(CR_id, k, j, i);
          }
        }
#pragma omp simd
        for (int i=is; i<=ie; ++i) {
        // Alias the primitives of gas at current stage
        const Real &gas_rho  = w(IDN, k, j, i);


        dens_prim(gas_id,i) = gas_rho;
        inv_gas_rho(i) = 1.0/gas_rho;

        const Real Bx = bcc(IB1,k,j,i);
        const Real By = bcc(IB2,k,j,i);
        const Real Bz = bcc(IB3,k,j,i);
        const Real B2 = Bx*Bx + By*By + Bz*Bz;
        const Real B = std::sqrt(B2);

        vA(i) = std::sqrt(B2*inv_gas_rho(i));

        if (B2>TINY_NUMBER) {
          const cr_frame::MagneticFrame frame =
              cr_frame::BuildMagneticFrame(Bx, By, Bz, B, SQR(Bx) + SQR(By));
          cosXY(i) = frame.cos_xy;
          sinXY(i) = frame.sin_xy;
          cosZ(i) = frame.cos_z;
          sinZ(i) = frame.sin_z;
        }

        const Real &gas_vel1 = w(IVX, k, j, i);
        const Real &gas_vel2 = w(IVY, k, j, i);
        const Real &gas_vel3 = w(IVZ, k, j, i);

        mom1_prim(gas_id, i) = gas_vel1;
        mom2_prim(gas_id, i) = gas_vel2;
        mom3_prim(gas_id, i) = gas_vel3;
      }

        RotateArraysToFieldAligned(
            cosXY, sinXY, cosZ, sinZ, mom1_prim, mom2_prim, mom3_prim);

        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
#pragma omp simd
        for(int i=is; i<=ie; ++i){
        // <<Most general implicit integrator>>
        Real coeff1 = 1.0 + (sigma_p(CR_id,i) + sigma_m(CR_id,i))*Vm2*dt;
        Real coeff2 = 1.0 + Vm2*Vm2*SQR(Omega(CR_id,i)*dt);
        Real coeff3 = (vA(i)*(sigma_p(CR_id,i) - sigma_m(CR_id,i)) + mom1_prim(gas_id,i)*(sigma_p(CR_id,i) + sigma_m(CR_id,i)))*Vm*dt;
        Real coeff4 = mom2_prim(gas_id,i)*Vm*Omega(CR_id,i)*dt;
        Real coeff5 = mom3_prim(gas_id,i)*Vm*Omega(CR_id,i)*dt;
        Real coeff6;

        if ((sigma_p(CR_id,i)+sigma_m(CR_id,i))>TINY_NUMBER)
                coeff6 = (SQR(coeff3*inv_Vm/dt) + 4.0*vA(i)*vA(i)*sigma_p(CR_id,i)*sigma_m(CR_id,i)*coeff1)/(sigma_p(CR_id,i)+sigma_m(CR_id,i));
        else
                coeff6 = 0.0;

        Real coeff0 = gamma_cr*coeff2*dt*coeff6 - coeff1*(coeff2 - gamma_cr*(SQR(coeff4) + SQR(coeff5)));
        coeff0 = 1.0/coeff0;

        delta_Ecr(Fluid_id,i)   = coeff0*(-coeff1*coeff2*dens_prim(Fluid_id,i) + coeff2*coeff3*mom1_prim(Fluid_id,i) + coeff1*(coeff5 + coeff4*Vm2*Omega(CR_id,i)*dt)*mom2_prim(Fluid_id,i)
                + coeff1*(coeff5*Vm2*dt*Omega(CR_id,i) - coeff4)*mom3_prim(Fluid_id,i)) - dens_prim(Fluid_id,i);

        delta_mom1(Fluid_id,i)  = coeff0*(-gamma_cr*coeff3*coeff2*dens_prim(Fluid_id,i) + (gamma_cr*coeff2*dt*coeff6
                      - (coeff2 - gamma_cr*(SQR(coeff4) + SQR(coeff5))))*mom1_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff5 + coeff4*Vm2*dt*Omega(CR_id,i))*mom2_prim(Fluid_id,i)
                + gamma_cr*coeff3*(coeff5*Vm2*Omega(CR_id,i)*dt - coeff4)*mom3_prim(Fluid_id,i)) - mom1_prim(Fluid_id,i);

        delta_mom2(Fluid_id,i)  = coeff0*(-gamma_cr*coeff1*(coeff5 - coeff4*Vm2*dt*Omega(CR_id,i))*dens_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff4*Vm2*dt*Omega(CR_id,i) - coeff5)*mom1_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6 + coeff1*(gamma_cr*SQR(coeff4)-1.0))*mom2_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6*Vm2*dt*Omega(CR_id,i) + coeff1*(gamma_cr*coeff4*coeff5-Vm2*dt*Omega(CR_id,i)))*mom3_prim(Fluid_id,i))
                - mom2_prim(Fluid_id,i);

        delta_mom3(Fluid_id,i)  = coeff0*(-gamma_cr*coeff1*(coeff4 + coeff5*Vm2*dt*Omega(CR_id,i))*dens_prim(Fluid_id,i) + gamma_cr*coeff3*(coeff4 + coeff5*Vm2*dt*Omega(CR_id,i))*mom1_prim(Fluid_id,i)
                - (gamma_cr*coeff2*dt*coeff6*Vm2*dt*Omega(CR_id,i) - coeff1*(gamma_cr*coeff4*coeff5+Vm2*dt*Omega(CR_id,i)))*mom2_prim(Fluid_id,i)
                + (gamma_cr*coeff2*dt*coeff6 + coeff1*(gamma_cr*SQR(coeff5)-1.0))*mom3_prim(Fluid_id,i))
                - mom3_prim(Fluid_id,i);

        }
        }

        RotateArraysToLabFrame(
            cosXY, sinXY, cosZ, sinZ, delta_mom1, delta_mom2, delta_mom3);

        for (int CR_id=0; CR_id<NCRS; ++CR_id) {
          int Fluid_id = CR_id + 1;
          int cr_energy_id = 4*CR_id;
          int cr_flux1_id = cr_energy_id + 1;
          int cr_flux2_id = cr_energy_id + 2;
          int cr_flux3_id = cr_energy_id + 3;
#pragma omp simd
      for(int i=is; i<=ie; ++i){
        // Alias the conserves of gas
        Real &gas_mom1 = u(IM1, k, j, i);
        Real &gas_mom2 = u(IM2, k, j, i);
        Real &gas_mom3 = u(IM3, k, j, i);

        Real &Ecr  = cr_cons(cr_energy_id, k, j, i);
        Real &Fcr1 = cr_cons(cr_flux1_id, k, j, i);
        Real &Fcr2 = cr_cons(cr_flux2_id, k, j, i);
        Real &Fcr3 = cr_cons(cr_flux3_id, k, j, i);


        Fcr1  += delta_mom1(Fluid_id,i);
        Fcr2  += delta_mom2(Fluid_id,i);
        Fcr3  += delta_mom3(Fluid_id,i);
          }
        }
      }
    }
  return;
}
