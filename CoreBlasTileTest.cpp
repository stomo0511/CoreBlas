/*
 * CoreBlasTest.cpp
 *
 *  Created on: 2015/10/08
 *      Author: stomo
 */

#include <iostream>
#include <omp.h>
#include <cassert>

#include "TMatrix.hpp"

#include "CoreBlasTile.hpp"

using namespace std;

int main(int argc, char* argv[])
{
	int m = 8;
	int n = 8;
	int b = 4;
	int s = 2;

	cout << "Test for QR kernels\n";
	// for Single tile
	cout << "For Single Tile:\n";
	{
		BMatrix A1(b,b,s), A2(b,b,s);
		A1.Set_Rnd(20151008);
		A2.Set_Rnd(20160213);
		A1.Show_all();
		A2.Show_all();

		BMatrix T1(s,b,s), T2(s,b,s);

		GEQRT( &A1, &T1 );
		A1.Show_all();
		T1.Show_all();

		TSQRT( &A1, &A2, &T2 );
		A1.Show_all();
		A2.Show_all();
		T2.Show_all();
	}

	// for Tile matrix
	cout << endl << "For Tile Matrix:\n";
	{
		TMatrix A(m,n,b,b,s);
		A.Set_Rnd(20160213);

		Matrix Ma(m,n);
		A.Mat_Copy(Ma);
		cout << "Input matrix:\n";
		Ma.Show_all();

		int p = m/b;

		TMatrix T(p*s,n,s,b,s);

		GEQRT( A(0,0), T(0,0) );
		TSQRT( A(0,0), A(1,0), T(1,0) );
		LARFB( PlasmaLeft, PlasmaTrans, A(0,0), T(0,0), A(0,1) );
		SSRFB( PlasmaLeft, PlasmaTrans, A(1,0), T(1,0), A(0,1), A(1,1) );
		GEQRT( A(1,1), T(1,1) );

		A.Mat_Copy(Ma);
		cout << "Output matrix (after k=0 step):\n";
		Ma.Show_all();
	}

	cout << "Test for LU kernels\n";
	// for Single tile
	cout << "For Single Tile:\n";
	{
		BMatrix A1(b,b,s);
		BMatrix A2(b,b,s);
		BMatrix L(s,b,s);
		int* PIV = new int[ b ];

		A1.Set_Rnd(20151008);

		cout << "A1: \n";
		A1.Show_all();

		GETRF( &A1, PIV );

		cout << "U: \n";
		A1.Show_all();

		cout << "PIV: \n";
		for (int i=0; i<b; i++)
			cout << PIV[i] << ", ";
		cout << endl << endl;

		A2.Set_Rnd(20170107);

		cout << "A2: \n";
		A2.Show_all();

		TSTRF( &A1, &A2, &L, PIV );

		cout << "U:\n";
		A1.Show_all();
		cout << "L:\n";
		A2.Show_all();
		cout << "L?:\n";
		L.Show_all();

		cout << "PIV: \n";
		for (int i=0; i<b; i++)
			cout << PIV[i] << ", ";
		cout << endl << endl;

		delete [] PIV;
	}

	// for Tile matrix
	cout << endl << "For Tile Matrix:\n";
	{
		TMatrix A(m,n,b,b,s);
		A.Set_Rnd(20160213);

		BMatrix L(s,b,s);
		int* PIV = new int[ b ];

		cout << "Input matrix:\n";
		A(0,0)->Show_all();
		A(0,1)->Show_all();
		A(1,0)->Show_all();
		A(1,1)->Show_all();

		GETRF( A(0,0), PIV );
		GESSM( A(0,0), A(0,1), PIV );

		cout << "PIV: \n";
				for (int i=0; i<b; i++)
					cout << PIV[i] << ", ";
				cout << endl << endl;

		TSTRF( A(0,0), A(1,0), &L, PIV );
		SSSSM( &L, A(1,0), A(0,1), A(1,1), PIV );

		cout << "PIV: \n";
				for (int i=0; i<b; i++)
					cout << PIV[i] << ", ";
				cout << endl << endl;

		GETRF( A(1,1), PIV );

		cout << "Output matrix:\n";
		A(0,0)->Show_all();
		A(0,1)->Show_all();
		A(1,0)->Show_all();
		A(1,1)->Show_all();


		cout << "PIV: \n";
		for (int i=0; i<b; i++)
			cout << PIV[i] << ", ";
		cout << endl << endl;

		delete [] PIV;

	}
}
