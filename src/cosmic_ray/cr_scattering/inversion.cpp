//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file inverse.cpp
//! Compute the ludecompose, inverse of matrixes.

// C++ headers
#include <algorithm>   // min,max
#include <cstring>    // strcmp
#include <limits>
#include <sstream>

// Athena++ headers
#include "../../athena.hpp"
#include "../../athena_arrays.hpp"
#include "../../coordinates/coordinates.hpp"
#include "../../defs.hpp"
#include "../cosmic_ray.hpp"
#include "cr_scattering.hpp"


// Reference: "Nurmerical Recipes, 3ed", Charpter 2.3, William H. Press et al. 2007
// LUdecompose and Matrix Inversion
void CRScattering::LUdecompose(const AthenaArray<Real> &a_matrix,
                AthenaArray<Real> &index_vector, AthenaArray<Real> &lu_matrix) {

  int is = pmb_->is, ie = pmb_->ie;

  //scale_arr.ZeroClear();
  index_vector.ZeroClear();

  for (int m=0; m<NCRVARS; ++m) {
    for (int n=0; n<NCRVARS; ++n) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        lu_matrix(m, n, i) = a_matrix(m, n, i);
      }
    }
  }

  //biggest_arr.ZeroClear();
  //for (int m=0; m<NCRVARS; ++m) {     // Loop over rows to get the implicit scaling information
    //for (int n=0; n<NCRVARS; ++n) {
//#pragma omp simd
      //for (int i=is; i<=ie; ++i) {
        //if ((temp_arr(i) = std::abs(lu_matrix(m, n, i)) > biggest_arr(i)))
          //biggest_arr(i) = temp_arr(i);
      //}
    //}

//#pragma omp simd
    //for (int i=is; i<=ie; ++i) {              // biggest_arr(i) must be larger than 0.0
      //if (biggest_arr(i) == 0.0) {
        //std::stringstream msg;
        //msg << "### FATAL ERROR in Singular matrix in LU decomposition" << std::endl;
        //ATHENA_ERROR(msg);                    // No nonzero largest element.
      //}
      //scale_arr(m, i) = 1.0/biggest_arr(i);   // Save the scaling.
    //}
  //}

  for (int l=0; l<NCRVARS; ++l) {   // This is the outermost lmn loop
    biggest_arr.ZeroClear();         // Initialize for the search for largest pivot element.
    for (int m=l; m<NCRVARS; ++m) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        //temp_arr(i) = scale_arr(m, i)*std::abs(lu_matrix(m, l, i));
        temp_arr(i) = std::abs(lu_matrix(m, l, i));
        if (temp_arr(i) > biggest_arr(i)) {
          biggest_arr(i) = temp_arr(i);
          mmax_arr(i) = m;
        }
      }
    }

//#pragma omp simd
    for (int i=is; i<=ie; ++i) {
      if (l != mmax_arr(i)) {              // Interchange Rows
        for (int n=0; n<NCRVARS; ++n) {
          temp_arr(i)                  = lu_matrix(mmax_arr(i), n, i);
          lu_matrix(mmax_arr(i), n, i) = lu_matrix(l, n, i);
          lu_matrix(l, n, i)           = temp_arr(i);
        }

        det_arr(i) = -det_arr(i);                       // change the parity of det_arr
        //scale_arr(mmax_arr(i), i) = scale_arr(l, i);  // Also interchange the scale factor
      }

      index_vector(l, i) = mmax_arr(i);
      if (lu_matrix(l, l, i) == 0.0)
        lu_matrix(l, l, i) = TINY_NUMBER;
    }

    for (int m=l+1; m<NCRVARS; ++m) {         // Divide by the pivot element
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        temp_arr(i) = lu_matrix(m, l, i) /= lu_matrix(l, l, i);

      for (int n=l+1; n<NCRVARS; ++n) {       // Innermost loop: reduce remaining submatrix.
#pragma omp simd
        for (int i=is; i<=ie; ++i){
          lu_matrix(m, n, i) -= temp_arr(i)*lu_matrix(l, n, i);
        }
      }
    }
  }
  return;
}


void CRScattering::SolveLinearEquation(const AthenaArray<Real> &index_vector,
    const AthenaArray<Real> &lu_matrix, AthenaArray<Real> &b_vector, AthenaArray<Real> &x_matrix) {

  int mm = 0;
  int is = pmb_->is, ie = pmb_->ie;

  for (int m=0; m<NCRVARS; ++m)
#pragma omp simd
    for (int i=is; i<=ie; ++i)
      x_matrix(m, i) = b_vector(m, i);

  for (int m=0; m<NCRVARS; ++m) {       // When mm is set to a positive value,
//#pragma omp simd
    for (int i=is; i<=ie; ++i) {
      int mp = index_vector(m, i);       // it will become the index of the first nonvanishing element of b.
      sum_arr(i) = x_matrix(mp, i);      // We now do the forward substitution
      x_matrix(mp, i) = x_matrix(m, i);  // The only new wrinkle is to unscramble the permutation

      if (mm!=0)
        for (int n=mm-1; n<m; ++n)
          sum_arr(i) -= lu_matrix(m, n, i)*x_matrix(n, i);
      else if (sum_arr(i) != 0.0)        // A nonzero element was encountered, so from now on we
        mm = m+1;                        // will have to do the sums in the loop above.

      x_matrix(m, i) = sum_arr(i);
    }
  }

  for (int m=NCRVARS-1; m>=0; m--) { // Now we do the backsubstitution,
#pragma omp simd
    for (int i=is; i<=ie; ++i)
      sum_arr(i) = x_matrix(m, i);

    for (int n=m+1; n<NCRVARS; ++n)
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        sum_arr(i) -= lu_matrix(m, n, i)*x_matrix(n, i);

#pragma omp simd
    for (int i=is; i<=ie; ++i)
      x_matrix(m, i) = sum_arr(i)/lu_matrix(m, m, i); // Store a component of the solution vector X
  }
  return;
}


void CRScattering::SolveMultipleLinearEquation(const AthenaArray<Real> &index_vector,
      const AthenaArray<Real> &lu_matrix, AthenaArray<Real> &b_vector, AthenaArray<Real> &x_matrix) {

  int is = pmb_->is, ie = pmb_->ie;

  xx_arr.ZeroClear();

  for (int n=0; n<NCRVARS; ++n) {  // Copy and solve each column in turn.
    for (int m=0; m<NCRVARS; ++m) {
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        xx_arr(m, i) = b_vector(m, n, i);
    }

    SolveLinearEquation(index_vector, lu_matrix, xx_arr, xx_arr);

    for (int m=0; m<NCRVARS; ++m)
#pragma omp simd
      for (int i=is; i<=ie; ++i)
        x_matrix(m, n, i) = xx_arr(m, i);
  }
  return;
}


void CRScattering::Inverse(const AthenaArray<Real> &index_vector, const AthenaArray<Real> &lu_matrix,
                  AthenaArray<Real> &a_matrix, AthenaArray<Real> &a_inv_matrix) {

  int is = pmb_->is, ie = pmb_->ie;

  //Using the stored LU decomposition, return in ainv the matrix inverse A^-1.
  for (int m=0; m<NCRVARS; ++m) {
    for (int n=0; n<NCRVARS; ++n) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        a_inv_matrix(m, n, i) = a_matrix(m, n, i);
      }
    }
  }

  a_matrix.ZeroClear();

  for (int m=0; m<NCRVARS; ++m) {
    for (int n=0; n<NCRVARS; ++n) {
#pragma omp simd
      for (int i=is; i<=ie; ++i) {
        a_matrix(m, m, i) = 1.0;
      }
    }
  }

  SolveMultipleLinearEquation(index_vector, lu_matrix, a_matrix, a_inv_matrix);
  return;
}
