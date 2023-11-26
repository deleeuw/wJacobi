
#include "R.h"

int jacobi( const size_t n, double** a, double** vecs, double* vals, const size_t m )
// truncated eigenvalue decomposition of symmetric n by n matrix a
// eigenvectors returned in vecs and eigenvalues returned in vals
// number of components of vecs and proper values in vals equals m 
{
  // convergence constants  
  const size_t MAXITER = 128;                                              // 128 should be enough
  const double EPS = DBL_EPSILON;                                          // 2.2204460492503131e-16
  const double TOL = sqrt( EPS );                                          // 1.4901161193847656e-08
  const double TINY = pow( 10.0, ( log10( EPS ) + log10( TOL ) ) / 2.0 );  // 1.8189894035458617e-12

  // initialize vecs to identity matrix
  for ( size_t i = 1; i <= n; i++ ) {
    for ( size_t j = 1; j <= n; j++ ) vecs[i][j] = 0.0;
    vecs[i][i] = 1.0;
  }

  // administration
  double fold = 0.0;
  for ( size_t i = 1; i <= m; i++ ) fold += a[i][i] * a[i][i];

  // main loop
  size_t iter = 0;
  for ( iter = 1; iter <= MAXITER; iter++ ) {

    // loop over columns / components
    for ( size_t j = 1; j <= m; j++ ) {

      // loop over rows
      for ( size_t i = j + 1; i <= n; i++ ) {

        const double p = a[i][j];
        if ( fabs( p ) < TINY ) continue;

        // set constants
        const double q = a[i][i];
        const double r = a[j][j];
        const double d = 0.5 * ( q - r );
        const double t = -1.0 * d / sqrt( d * d + p * p );
        const double uu = 0.5 * ( 1.0 + t );
        const double u = sqrt( uu );
        const double vv = 0.5 * ( 1.0 - t );
        const double v = ( p < 0.0 ? -1.0 : 1.0 ) * sqrt( vv );
        const double twouvp = 2.0 * u * v * p;

        // rotate, using information from the lower triangle of a only
        // case of rotations 1 <= k <= j
        for ( size_t k = 1; k <= j; k++ ) {
          const double ik = a[i][k];
          const double jk = a[j][k];
          a[i][k] = u * ik - v * jk;
          a[j][k] = v * ik + u * jk;
        }

        // case of rotations j < k < i
        for ( size_t k = j + 1; k < i; k++ ) {
          const double ik = a[i][k];
          const double kj = a[k][j];
          a[i][k] = u * ik - v * kj;
          a[k][j] = v * ik + u * kj;
        }

        // case of rotations i <= k <= n
        for ( size_t k = i; k <= n; k++ ) {
          const double ki = a[k][i];
          const double kj = a[k][j];
          a[k][i] = u * ki - v * kj;
          a[k][j] = v * ki + u * kj;
        }

        // set eigenvectors
        for ( size_t k = 1; k <= n; k++ ) {
          const double ki = vecs[k][i];
          const double kj = vecs[k][j];
          vecs[k][i] = u * ki - v * kj;
          vecs[k][j] = v * ki + u * kj;
        }

        // update a
        a[i][i] = uu * q + vv * r - twouvp;
        a[j][j] = vv * q + uu * r + twouvp;
        a[i][j] = u * v * ( q - r ) + ( uu - vv ) * p;
      }
    }

    // administration
    double fnew = 0.0;
    for ( size_t i = 1; i <= m; i++ ) fnew += a[i][i] * a[i][i];
    if ( ( fnew - fold ) < TINY ) break;
    fold = fnew;
  }

  // set eigenvalues
  for ( size_t i = 1; i <= n; i++ ) vals[i] = a[i][i];

  // proper return value
  return( iter <= MAXITER ? 0 : 1 );
} // jacobi

void Cjacobi( int* rn, double* ra, double* rvecs, double* rvals, int* rm )
// Function Cjacobi() performs truncated eigenvalue decomposition.
// Copyright (C) 2020 Frank M.T.A. Busing (e-mail: busing at fsw dot leidenuniv dot nl)
// This function is free software:
// you can redistribute it and/or modify it under the terms of the GNU Lesser General Public License
// as published by the Free Software Foundation.
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY;
// without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// See the GNU General Public License for more details.
// You should have received a copy of the GNU General Public License along with this function.
// If not, see <https://www.gnu.org/licenses/>.
{
  // transfer to C; set one-based
  size_t n = *rn;
  size_t m = *rm;
  double** restrict a = ( double ** )malloc( n * sizeof( double ) ); a--;
  for ( size_t i = 1, im1 = 0; i <= n; i++, im1++ ) a[i] = &ra[im1 * n - 1];
  double** restrict vecs = ( double ** )malloc( n * sizeof( double ) ); vecs--;
  for ( size_t i = 1, im1 = 0; i <= n; i++, im1++ ) vecs[i] = &rvecs[im1 * n - 1];
  double* restrict vals = &rvals[0]; vals--;

  // run function one-based
  int rvalue = jacobi( n, a, vecs, vals, m );

  // transfer to R; transpose relevant matrices
  if ( rvalue == 0 ) {
    for ( size_t i = 2; i <= n; i++ ) for ( size_t j = 1; j < i; j++ ) {
      const double work = vecs[i][j];
      vecs[i][j] = vecs[j][i];
      vecs[j][i] = work;
    }
  }  
  else m = 0;  // return on error

  // de-allocate memory
  a++; free( a );
  vecs++; free( vecs );

} // Cjacobi
