//============================================================================
// Name        : CoreBlas.cpp
// Author      : T. Suzuki
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <cassert>
#include <omp.h>
#include "CoreBlas.hpp"

#ifndef __DEF__MIN__
#define __DEF__MIN__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __DEF__MIN__

#ifndef  __DEF_MAX__
#define  __DEF_MAX__
#define max(a,b) (((a)>(b)) ? (a) : (b))
#endif

/*
 * GEQRT conputes a QR factorization of a tile A: A = Q * R
 *
 * @param A (M x N) tile matrix
 * @param T (IB x N) upper triangular block reflector
 */
void GEQRT( BMatrix *A, BMatrix *T )
{
	const int M = A->m();
	const int N = A->n();
	const int IB = A->ib();
	const int LDA = A->m();
	const int LDT = T->m();

	const int NB = max(LDA,LDT);

	double* WORK = new double[ IB*NB ];
	double* TAU = new double[ NB ];

	CORE_dgeqrt( M, N, IB,
			A->top(), LDA,
			T->top(), LDT,
			TAU, WORK );

	delete [] WORK;
	delete [] TAU;
}

/*
 * TSQRT conputes a QR factorization of a rectangular matrix formed by cupling (N x N) upper triangular tile A1 on top of (M x N) tile A2
 *
 * @param A1 (N x N) upper triangular tile matrix
 * @param A2 (M x N) tile matrix
 * @param T (IB x N) upper triangular block reflector
 */
void TSQRT( BMatrix *A1, BMatrix *A2, BMatrix *T )
{
	const int M = A2->m();
	const int N = A1->n();

	assert( N == A2->n() );

	const int IB = A1->ib();
	const int LDA1 = A1->m();
	const int LDA2 = A2->m();
	const int LDT = T->m();

	const int NB = max(LDA1,LDT);

	double* WORK = new double[ IB*NB ];
	double* TAU = new double[ NB ];

	CORE_dtsqrt( M, N, IB,
			A1->top(), LDA1,
			A2->top(), LDA2,
			T->top(), LDT,
			TAU, WORK );

	delete [] WORK;
	delete [] TAU;
}

/*
 * LARFB updates (M x N) tile C with the transformation formed with A and T
 *
 * @param[in] side
 *		@arg PlasmaLeft:  apply transformation from the left
 *		@arg PlasmaRight: apply transformation from the right
 *
 * @param[in] trans
 *		@arg PlasmaNoTrans: no transpose the matrix
 *		@arg PlasmaTrans: transpose the matrix
 *
 * @param[in] A (LDA x K) tile matrix
 * @param[in] T (IB x K) upper triangular block reflector
 * @param[in,out] C (M x N) tile matrix
 */
void LARFB( PLASMA_enum side, PLASMA_enum trans,
		BMatrix *A, BMatrix *T, BMatrix *C )
{
	assert( (side==PlasmaLeft) || (side==PlasmaRight) );
	assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );

	const int M = C->m();
	const int N = C->n();
	const int K = A->n();

	if (side == PlasmaLeft)
		assert( M >= K );
	else // (side == PlasmaRight)
		assert( N >= K );

	const int IB = A->ib();
	const int LDA = A->m();
	const int LDT = T->m();
	const int LDC = C->m();

	const int NB = max(LDA,LDT);

	double* WORK = new double[ IB*NB ];

	CORE_dormqr( side, trans,
			M, N, K, IB,
			A->top(), LDA,
			T->top(), LDT,
			C->top(), LDC,
			WORK, NB );

	delete [] WORK;
}

/*
 * SSRFB updates (M1 x N1) tile C1 and (M2 x N2) tile C2 with the transformation formed with A and T
 *
 * @param[in] side
 *		@arg PlasmaLeft:  apply transformation from the left
 *		@arg PlasmaRight: apply transformation from the right
 *
 * @param[in] trans
 *		@arg PlasmaNoTrans: no transpose the matrix
 *		@arg PlasmaTrans: transpose the matrix
 *
 * @param[in] A (LDA x K) tile matrix
 * @param[in] T (IB x K) upper triangular block reflector
 * @param[in,out] C1 (M1 x N1) tile matrix
 * @param[in,out] C2 (M2 x N2) tile matrix
 */
void SSRFB( PLASMA_enum side, PLASMA_enum trans,
		BMatrix *A, BMatrix *T,
	    BMatrix *C1, BMatrix *C2 )
{
	assert( (side==PlasmaLeft) || (side==PlasmaRight) );
	assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );

	const int M1 = C1->m();
	const int M2 = C2->m();

	if (side == PlasmaRight)
		assert( M2 == M1);

	const int N1 = C1->n();
	const int N2 = C2->n();

	if (side == PlasmaLeft)
		assert( N2 == N1);

	const int K = A->n();

	const int IB = C1->ib();
	const int LDA1 = C1->m();
	const int LDA2 = C2->m();
	const int LDV = A->m();
	const int LDT = T->m();

	int LDWORK;
	if (side == PlasmaLeft)
		LDWORK = IB;
	else // side == PlasmaRight
		LDWORK = M1;

	int WSIZE;
	if (side == PlasmaLeft)
		WSIZE = N1;
	else // side == PlasmaRight
		WSIZE = IB;

	double* WORK = new double[ LDWORK * WSIZE ];

	CORE_dtsmqr( side, trans,
			M1, N1, M2, N2, K, IB,
			C1->top(), LDA1,
			C2->top(), LDA2,
			A->top(), LDV,
			T->top(), LDT,
			WORK, LDWORK);

	delete [] WORK;
}

/*
 * dorgqr: genarates (M x N) orthogonal matrix Q: A = Q x R
 *
 * @param[in] A tile matrix
 * @param[in] T tile matrix
 * @param[in] Q tile matirx
 *
 */

/*
 * TTQRT conputes a QR factorization of a rectangular matrix formed by cupling (N x N) upper triangular tile A1 on top of (M x N) upper trapezoidal tile A2
 *
 * @param A1 (N x N) upper triangular tile matrix
 * @param A2 (M x N) upper trapezoidal tile matrix
 * @param T (IB x N) upper triangular block reflector
 */
void TTQRT( BMatrix *A1, BMatrix *A2, BMatrix *T )
{
	const int M = A2->m();
	const int N = A1->n();

	assert( N == A2->n() );

	const int IB = A1->ib();
	const int LDA1 = A1->m();
	const int LDA2 = A2->m();
	const int LDT = T->m();

	const int NB = max(LDA1,LDT);

	double* WORK = new double[ IB*NB ];
	double* TAU = new double[ NB ];

	CORE_dttqrt( M, N, IB,
			A1->top(), LDA1,
			A2->top(), LDA2,
			T->top(), LDT,
			TAU, WORK );

	delete [] WORK;
	delete [] TAU;
}

/*
 * STRFB updates (M1 x N1) tile C1 and (M2 x N2) tile C2 with the transformation formed with A and T
 *
 * @param[in] side
 *		@arg PlasmaLeft:  apply transformation from the left
 *		@arg PlasmaRight: apply transformation from the right
 *
 * @param[in] trans
 *		@arg PlasmaNoTrans: no transpose the matrix
 *		@arg PlasmaTrans: transpose the matrix
 *
 * @param[in] A (LDA x K) tile matrix
 * @param[in] T (IB x K) upper triangular block reflector
 * @param[in,out] C1 (M1 x N1) tile matrix
 * @param[in,out] C2 (M2 x N2) tile matrix
 */
void STRFB( PLASMA_enum side, PLASMA_enum trans,
		BMatrix *A, BMatrix *T,
	    BMatrix *C1, BMatrix *C2 )
{
	assert( (side==PlasmaLeft) || (side==PlasmaRight) );
	assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );

	const int M1 = C1->m();
	const int M2 = C2->m();

	const int N1 = C1->n();
	const int N2 = C2->n();

	const int K = A->n();

	const int IB = C1->ib();
	const int LDA1 = C1->m();
	const int LDA2 = C2->m();
	const int LDV = A->m();
	const int LDT = T->m();

	int LDWORK;
	if (side == PlasmaLeft)
		LDWORK = IB;
	else // side == PlasmaRight
		LDWORK = M1;

	int WSIZE;
	if (side == PlasmaLeft)
		WSIZE = N1;
	else // side == PlasmaRight
		WSIZE = IB;

	double* WORK = new double[ LDWORK * WSIZE ];

	CORE_dttmqr( side, trans,
			M1, N1, M2, N2, K, IB,
			C1->top(), LDA1,
			C2->top(), LDA2,
			A->top(), LDV,
			T->top(), LDT,
			WORK, LDWORK);

	delete [] WORK;
}

/*
 * dorgqr: genarates (M x N) orthogonal matrix Q: A = Q x R
 *
 * @param[in] A tile matrix
 * @param[in] T tile matrix
 * @param[in] Q tile matirx
 *
 */
void dorgqr( const TMatrix A, const TMatrix T, TMatrix& Q )
{
	assert( A.N() == Q.N() );

	const int aMT = A.mt();
	const int aNT = A.nt();
	const int qMT = Q.mt();
	const int qNT = Q.nt();

	for (int tk = min(aMT, aNT)-1; tk+1 >= 1; tk--)
	{
		for (int ti = qMT - 1; ti > tk; ti--)
		{
			#pragma omp parallel for
			for (int tj = tk; tj < qNT; tj++)
			{
				SSRFB( PlasmaLeft, PlasmaNoTrans,
						A(ti,tk), T(ti,tk), Q(tk,tj), Q(ti,tj) );
			}
		}
		#pragma omp parallel for
		for (int tj = tk; tj < qNT; tj++)
		{
			LARFB( PlasmaLeft, PlasmaNoTrans,
					A(tk,tk), T(tk,tk), Q(tk,tj) );
		}
	}
}
