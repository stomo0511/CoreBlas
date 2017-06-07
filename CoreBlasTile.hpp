/*
 * CoreBlas.hpp
 *
 *  Created on: 2015/09/11
 *      Author: stomo
 */

#ifndef COREBLASTILE_HPP_
#define COREBLASTILE_HPP_

#include <iostream>
#include <plasma.h>
#include <core_blas.h>

#include "BMatrix.hpp"
#include "TMatrix.hpp"

// for QR
void GEQRT( BMatrix *A, BMatrix *T );
void TSQRT( BMatrix *A1, BMatrix *A2, BMatrix *T );
void LARFB( plasma_enum_t side, plasma_enum_t trans,
			BMatrix *A, BMatrix *T, BMatrix *C );
void SSRFB( plasma_enum_t side, plasma_enum_t trans,
		   BMatrix *A, BMatrix *T, BMatrix *C1, BMatrix *C2 );

// for TSQR
void TTQRT( BMatrix *A1, BMatrix *A2, BMatrix *T );
void TTMQR( plasma_enum_t side, plasma_enum_t trans,
		   BMatrix *A, BMatrix *T, BMatrix *C1, BMatrix *C2 );

// for check QR
void dorgqr( const TMatrix A, const TMatrix T, TMatrix& Q );

// for Cholesky
//void POTRF( BMatrix *A );
//void SYRK( BMatrix *A, BMatrix *B );
//void TRSM( BMatrix *A, BMatrix *B );
//void GEMM( BMatrix *A, BMatrix *B, BMatrix *C );

// for LU
//void GETRF( BMatrix *A, int *PIV );
//void TSTRF( BMatrix *U, BMatrix *A, BMatrix *L, int *PIV );
//void GESSM( BMatrix *L, BMatrix *A, const int *PIV );
//void SSSSM( BMatrix *L1, BMatrix *L2, BMatrix *A1, BMatrix *A2, const int *PIV );

#endif /* COREBLASTILE_HPP_ */
