//============================================================================
// Name        : CoreBlas.cpp
// Author      : T. Suzuki
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C, Ansi-style
//============================================================================

#include <cassert>
#include <omp.h>

#include "CoreBlasTile.hpp"

#ifndef __DEF__MIN__
#define __DEF__MIN__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __DEF__MIN__

#ifndef  __DEF_MAX__
#define  __DEF_MAX__
#define max(a,b) (((a)>(b)) ? (a) : (b))
#endif

/*
 * GEQRT computes a QR factorization of a tile A: A = Q * R
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
 * TTMQR updates (M1 x N1) tile C1 and (M2 x N2) tile C2 with the transformation formed with A and T
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
void TTMQR( PLASMA_enum side, PLASMA_enum trans,
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
	assert( A.M() == Q.M() );

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

/*
 * POTRF Computes the Cholesky factorization of a symmetric positive definite matrix A
 *
 * @param[in,out] A
 *
 */
/*void POTRF( BMatrix *A )
{
	const int M = A->m();
	const int N = A->n();
	const int LDA = A->m();

	int info;

	CORE_dpotrf( PlasmaLower,
			M,
			A->top(), LDA,
			&info );
}*/

/*
 * SYRK Performs one of the hermitian rank k operations
 *
 *
 * @param[in] A
 * @param[in,out] C
 *
 *
 *    C = \alpha [ A \times A' ] + \beta C ,
 *
  *
 *  where alpha and beta are real scalars, C is an N-by-N hermitian
 *  matrix and A is an N-by-K matrix.
 *
 * @param[in] uplo
 *          = PlasmaLower: Lower triangle of C is stored.
 *
 * @param[in] trans
 *          = PlasmaNoTrans:   A is not transposed;
 *
 * @param[in] N
 *          N specifies the order of the matrix C. N must be at least zero.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix A.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          A is a LDA-by-K matrix.
  *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA must be at least max( 1, N )
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array uplo part of the matrix is overwritten
 *          by the uplo part of the updated matrix.
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max( 1, N ).
 *
 */
/*void SYRK( BMatrix *A, BMatrix *C )
{

	const int M = A->m();

	assert( M == C->m() );

    CORE_dsyrk(
         PlasmaLower, PlasmaNoTrans,
		 M, M,
         -1.0, A->top(), M,
          1.0, C->top(), M);
}*/

/*
 *  TRSM - Computes triangular solve X*A = B.
 *
 *******************************************************************************
 *
 * @param[in] side
 *          = PlasmaRight: X*A = B
 *
 * @param[in] uplo
 *          = PlasmaLower: Lower triangle of A is stored.
 *
 * @param[in] transA
 *          = PlasmaTrans:     A is not transposed;
 *
 * @param[in] diag
 *          = PlasmaNonUnit: A is non unit;
 *
 * @param[in] M
 *          The order of the matrix A. M >= 0.
 *
 * @param[in] N
 *          The number of right hand sides, i.e., the number of columns of the matrix B. N >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha.
 *
 * @param[in] A
 *          The triangular matrix A. If uplo = PlasmaUpper, the leading M-by-M upper triangular
 *          part of the array A contains the upper triangular matrix, and the strictly lower
 *          triangular part of A is not referenced. If uplo = PlasmaLower, the leading M-by-M
 *          lower triangular part of the array A contains the lower triangular matrix, and the
 *          strictly upper triangular part of A is not referenced. If diag = PlasmaUnit, the
 *          diagonal elements of A are also not referenced and are assumed to be 1.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in,out] B
 *          On entry, the M-by-N right hand side matrix B.
 *          On exit, if return value = 0, the M-by-N solution matrix X.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,M).
  */
/*void TRSM( BMatrix *A, BMatrix *B )
{
	const int M = A->m();
	const int N = B->n();

	assert ( M == B->m() );

	CORE_dtrsm( PlasmaRight, PlasmaLower, PlasmaTrans, PlasmaNonUnit,
			M, N,
			1.0,
			A->top(), M,
			B->top(), M );
}*/

/*
 *  *  GEMM - Performs one of the matrix-matrix operations
 *
 *    C = \alpha [ A \times B' ] + \beta C ,
 *
 *  alpha and beta are scalars, and A, B and C  are matrices, with op( A )
 *  an m by k matrix, op( B ) a k by n matrix and C an m by n matrix.
 *
 *******************************************************************************
 *
 * @param[in] transA
 *          Specifies whether the matrix A is transposed, not transposed or ugate transposed:
 *          = PlasmaNoTrans:   A is not transposed;
 *          = PlasmaTrans:     A is transposed;
 *          = PlasmaTrans: A is ugate transposed.
 *
 * @param[in] transB
 *          Specifies whether the matrix B is transposed, not transposed or ugate transposed:
 *          = PlasmaNoTrans:   B is not transposed;
 *          = PlasmaTrans:     B is transposed;
 *          = PlasmaTrans: B is ugate transposed.
 *
 * @param[in] M
 *          M specifies the number of rows of the matrix op( A ) and of the matrix C. M >= 0.
 *
 * @param[in] N
 *          N specifies the number of columns of the matrix op( B ) and of the matrix C. N >= 0.
 *
 * @param[in] K
 *          K specifies the number of columns of the matrix op( A ) and the number of rows of
 *          the matrix op( B ). K >= 0.
 *
 * @param[in] alpha
 *          alpha specifies the scalar alpha
 *
 * @param[in] A
 *          A is a LDA-by-ka matrix, where ka is K when  transA = PlasmaNoTrans,
 *          and is  M  otherwise.
 *
 * @param[in] LDA
 *          The leading dimension of the array A. LDA >= max(1,M).
 *
 * @param[in] B
 *          B is a LDB-by-kb matrix, where kb is N when  transB = PlasmaNoTrans,
 *          and is  K  otherwise.
 *
 * @param[in] LDB
 *          The leading dimension of the array B. LDB >= max(1,N).
 *
 * @param[in] beta
 *          beta specifies the scalar beta
 *
 * @param[in,out] C
 *          C is a LDC-by-N matrix.
 *          On exit, the array is overwritten by the M by N matrix ( alpha*op( A )*op( B ) + beta*C )
 *
 * @param[in] LDC
 *          The leading dimension of the array C. LDC >= max(1,M).
 *

 */
/*void GEMM( BMatrix *A, BMatrix *B, BMatrix *C )
{

    CORE_dgemm( PlasmaNoTrans, PlasmaTrans,
        tempmn, A.nb, A.nb,
        -1.0,
		A(m, n), ldam,
               A(k, n), ldak,
        1.0, A(m, k), ldam);

}*/

// for LU
/*
 * GETRF computes an LU factorization of a tile A: A = P * L * U
 * using partial tile pivoting with row interchanges (incremental pivoting)
 *
 * @param A (M x N) tile matrix
 * @param PIV (M) pivot indices the define the permutations
 */
void GETRF( BMatrix *A, int *PIV )
{
	const int M = A->m();
	const int N = A->n();
	const int IB = A->ib();
	const int LDA = A->m();

	int info;

	int ret = CORE_dgetrf_incpiv( M, N, IB, A->top(), LDA, PIV, &info );

	assert(info==0);
	//  0 on successful exit
	// <0 if -i, the i-th argument had an illegal value
	// >0 if i, U(i,i) is exactly zero

	assert(ret==PLASMA_SUCCESS);
}

 /*
  * TSTRF computes an LU factorization of a tile formed by an upper triangular (NB x N) tile U
  * on top of a (M x N) tile A using partial tile pivoting with row interchanges
  *
  * @param A (M x N) tile matrix
  * @param PIV (M) pivot indices the define the permutations
  */
void TSTRF( BMatrix *U, BMatrix *A, BMatrix *L, int *PIV )
{
	const int M = A->m();
	const int N = A->n();
	const int IB = A->ib();
	const int LDA = A->m();

	int info;

	// PLASMA tile LU with incpiv では連立方程式を解くために L が必要
	int ret = CORE_dtstrf(int M, int N, int IB, int NB,
            double *U, int LDU,
            double *A, int LDA,
            double *L, int LDL,
            int PIV,
            double *WORK, int LDWORK,
            &info );

	assert(info==0);
	//  0 on successful exit
	// <0 if -i, the i-th argument had an illegal value
	// >0 if i, U(i,i) is exactly zero

	assert(ret==PLASMA_SUCCESS);
}

/*
 * GESSM applies the factors L computed by GETRF to a (M,N) tile A
 */
void GESSM( BMatrix *L, BMatrix *A, const int *PIV )
{
	const int M = A->m();
	const int N = A->n();
	const int K = L->n();
	const int IB = A->ib();
	const int LDA = A->m();
	const int LDL = L->m();

	int ret = CORE_dgessm( M, N, K, IB, PIV, L->top(), LDL, A->top(), LDA );

	assert(ret==PLASMA_SUCCESS);
}

void SSSM( BMatrix *A, BMatrix *B, BMatrix *C )
{

}
