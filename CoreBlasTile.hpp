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

void GEQRT( BMatrix *A, BMatrix *T );
void TSQRT( BMatrix *A1, BMatrix *A2, BMatrix *T );
void LARFB( PLASMA_enum side, PLASMA_enum trans,
			BMatrix *A, BMatrix *T, BMatrix *C );
void SSRFB( PLASMA_enum side, PLASMA_enum trans,
		   BMatrix *A, BMatrix *T, BMatrix *C1, BMatrix *C2 );

void TTQRT( BMatrix *A1, BMatrix *A2, BMatrix *T );
void STRFB( PLASMA_enum side, PLASMA_enum trans,
		   BMatrix *A, BMatrix *T, BMatrix *C1, BMatrix *C2 );

void dorgqr( const TMatrix A, const TMatrix T, TMatrix& Q );


#endif /* COREBLASTILE_HPP_ */
