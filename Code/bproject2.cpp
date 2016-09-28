#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include "lib.h"
//#include "jacobi.h"

//Function to find maximum matrix element.

double maxoffdiag ( double **A, int * k, int * l, int n)
{
	double max = 0.0;

	for ( int i = 0; i < n; i++ ) {
		for ( int j = 0; j < n; j++ ) {
			if ( fabs(A[i][j]) > max ) {
				max  = fabs(A[i][j]);
				*l = i;
				*k = j;
			}
		}
	}
	return max;
}

//Function to find the values of cos and sin
void rotate ( double ** A, double ** R, int k, int l, int n)
{
	double s, c;
	if ( A[k][l] != 0.0) {
		double t, tau;
		tau = ( A[l][l] - A[k][k] ) / ( 2*A[k][l] );
		if ( tau > 0 ) {
			t = 1.0 / ( tau + sqrt( 1.0 + tau*tau ) );
		}
		else {
			t = -1.0 / ( -tau + sqrt( 1.0 + tau*tau ) );
		}

		c = 1 / sqrt( 1 + t*t );
		s = c * t;
	}
	else {
		c = 1.0;
		s = 0.0;
	}
	double a_kk, a_ll, a_ik, a_il, r_ik, r_il;
	a_kk = A[k][k];
	a_ll = A[l][l];
	//changing the matrix elements with indices k and l
	A[k][k] = c*c*a_kk - 2.0*c*s*A[k][l] + s*s*a_ll;
	A[l][l] =  s*s*a_kk + 2.0*c*s*A[k][l] + c*c*a_ll;
	A[k][l] = 0.0; //hardcoding of zeros
	A[l][k] = 0.0;
	//change remaining elements
	for ( int i = 0; i < n; i++ ) {
		if ( i != k && i != l ) {
			a_ik = A[i][k];
			a_il = A[i][l];
			A[i][k] = c*a_ik - s*a_il;
			A[k][i] = A[i][k];
			A[i][l] = c*a_il + s*a_ik;
			A[l][i] = A[i][l];
		}
		r_ik = R[i][k];
		r_il = R[i][l];
		R[i][k] = c*r_ik - s*r_il;
		R[i][l] = c*r_il + s*r_ik;
	}
	return;
}

void jacobi_method ( double ** A, double ** R, int n )
{
//Settup of eigenvector
	for ( int i = 0; i < n; i++){
		for ( int j = 0; j < n; j++){
			if ( i == j){
				R[i][j] = 1.0;
			}
			else{
				R[i][j] = 0.0;
			}
		}
	}
	int k, l;

	double epsilon = 1.0e-8;
	double max_number_iterations = (double) n * (double) n * (double) n;
	int iterations = 0;
	double max_offdiag = maxoffdiag ( A, &k, &l, n );

	while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ){
		max_offdiag = maxoffdiag( A, &k, &l, n );
		rotate (A, R, k, l, n);
		iterations++;
	}
	std::cout <<"Number of iterations: " << iterations << "\n";
	return;
}


int main( int argc, char * argv[] )
{
	int n =  atoi(argv[1]);
	double rhoN =  atof(argv[1]);
	double rho0 = 0.0;
	double h = ( rhoN - rho0 ) / (double) n;
	double* rho = new double[n];
	
	//define our tridiagonal matrix
	double **A;
	double **R;
	A = new double*[n];
	R = new double*[n];
	for ( int i = 0; i < n; i++ ) {
		A[i] = new double [n]; 
		R[i] = new double [n];
	}

	for ( int i = 0; i < n; i++) {
		rho[i] = rho0 + i*h;
		double Vi = rho[i]*rho[i];
		for (int j = 0; j < n; j++) {
			if ( i == j ) {
				A[i][j] = 2 / ( h*h ) + Vi;
			}
			if ( fabs( i - j ) == 2 ) {
				A[i][j] = -1 / ( h*h );
				A[j][i] = A[i][j];
			}
		}
	}
	jacobi_method ( A, R, n );
}