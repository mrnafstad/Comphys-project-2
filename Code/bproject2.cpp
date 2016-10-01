#include <math.h>
#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
//#include "lib.h"

//Function to find maximum matrix element.

double maxoffdiag ( double **A, int * k, int * l, int n)
{
	double max = 0.0;

	for ( int i = 0; i < n; i++ ) {
		for ( int j = 0; j < n; j++ ) {
			if ( fabs(A[i][j]) > max && i != j) {
				max  = fabs(A[i][j]);
				*l = i;
				*k = j;
				//cout << i << j << "\n";
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
			t = 1.0 / ( double ) ( tau + sqrt( 1.0 + tau*tau ) );
		}
		else {
			t = -1.0 / ( double ) ( -tau + sqrt( 1.0 + tau*tau ) );
		}

		c = 1 / ( double ) sqrt( 1 + t*t );
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

void jacobi_method ( double ** A, double ** R, int n)
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
	double max_offdiag = maxoffdiag ( A, &k, &l, n);

	while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ){
		max_offdiag = maxoffdiag( A, &k, &l, n);
		rotate (A, R, k, l, n);
		iterations++;
	}
	std::cout <<"Number of iterations: " << iterations << "\n";
	return;
}

void three_low ( double * v, int n)
{
//v is array to find the lowest values of, n is length of array

	double k, l, total_time;
	double * smallest = new double [3];
	clock_t s, f;
	s = clock(); 
	for ( int i = 0; i < 3; i++) smallest[i] = v[i];

	for ( int i = 0; i < n; i++ ) {
		
		if ( fabs(v[i]) < fabs(smallest[0]) ) {
			smallest[0] = v[i];
			k = i;
			//cout << v[i] << " k = " << i << "\n";
		}
	}

	for ( int i = 0; i < n; i++ ) {
		if ( ( fabs(v[i]) < fabs(smallest[1]) ) && i != k ) {
			smallest[1] = v[i];
			l = i;
			//cout << v[i] << " l = " << i << "\n";
		}
	}

	for ( int i = 0; i < n; i++ ) {
		if ( ( fabs(v[i]) < fabs(smallest[2]) ) && ( i !=k ) && i != l) {
			smallest[2] = v[i];
			//cout << v[i] << " m = " << i << "\n";
		}
	}

	for ( int i = 0; i < 3; i++) std::cout << smallest[i] << "\n";
	f = clock();
	total_time = ( ( double ) ( f - s ) / CLOCKS_PER_SEC );
	printf("Time spent on extracting eigenvalues:  %1.3e \n", total_time);
	return;
}

void unit_test(double **R){
	double **Matr, max_test;
	int n_test = 5, k, l;

	Matr = new double*[n_test];

	for(int i = 0; i < n_test; i++){
		Matr[i] = new double [n_test];
	}

	/*Creating a symmetric tridiagonal 5x5 matrix with 1's on the diagonal
	and the numbers 1-4 on the diagonals above and below. 
	Known eigenvalues:
	1: -4.16352, 2: -0.827046 3: 1.0, 4: 2.82705, 5: 6.16352
	https://www.wolframalpha.com/input/?i=eigenvalues+{{1,+1,+0,+0,+0},{1,+1,+2,+0,0},+{0,+2,+1,+3,+0},+{0,+0,+3,+1,+4},+{0,+0,+0,+4,+1}}*/

	for(int i = 0; i < n_test; i++){
		for(int j = 0; j < n_test; j++){
			if(i == j){
				Matr[i][j] = 1;
			}
			else if(j == i +1){
				Matr[i][j] = i + 1;
			}
			else if(j == i -1){
				Matr[i][j] = i;
			}
			else{
				Matr[i][j] = 0;
			}
		}
	}
	// Print Matrix elements to see that it's correct
	/*for(int i = 0; i < n_test; i++){
		for(int j = 0; j < n_test; j++){
			printf("%f\n", Matr[i][j]);
		}
	}*/

	printf("-------------------------------------------\nUNIT TEST\n\n");

	//Known max offdiag element: 4
	max_test = maxoffdiag(Matr, &k, &l, n_test);
	printf("Max off diagonal element of Matr: %f\n\n", max_test);
	
	jacobi_method(Matr, R, n_test);

	for(int i = 0; i < n_test; i++){
		printf("Eigenvalue %i: %f\n", i+1, Matr[i][i]);
	}
	printf("-------------------------------------------\n");

	delete[] Matr;

	return;

}


double potential(double omega, double r, int coloumb )
{
	double V;
	if ( coloumb == 0.0 ) V =  omega*r*r;
	if ( coloumb == 1.0 ) V = omega*r*r + 1.0/(double)r;
	return V;
}


int main( int argc, char * argv[] )
{
	int n =  atoi(argv[1]);
	double rhoN =  atof(argv[2]);
	double rho0 = 0.0;
	double h = ( rhoN - rho0 ) / (double) n;
	double* rho = new double[n];
	double total_time;
	clock_t start, finish;
	
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
		double Vi = potential( 1.0,  rho[i], 1.0 );
		//double Vi = rho[i]*rho[i];
		for (int j = 0; j < n; j++) {
			if ( i == j ) {
				A[i][j] = 2.0 / (double) ( h*h ) + Vi;
			}
			else if ( fabs( i - j ) == 1 ) {
				A[i][j] = -1.0 / (double)( h*h );
				A[j][i] = A[i][j];
			}
			else A[i][j] = 0.0;
			//cout << "A_ "<< i << j << " " << A[i][j] << "\n";
		}
		//cout << "--\n";
	}

	//RUNNING UNIT TESTS
	unit_test(R);


	//MAIN RUN
	//------------------------------------
	start = clock();
	jacobi_method ( A, R, n );
	finish = clock();
	total_time = ( ( double ) ( finish - start ) / CLOCKS_PER_SEC );
	printf("Time spent on algorithm:  %1.3e \n", total_time);
	double * lam = new double[n];
	for ( int i = 0; i < n; i++) {
		lam[i] = A[i][i];
	}

	//------------------------------------


		//cout <<" | " << lam[i];
		/*for ( int j = 0; j < n; j++) {
			cout << A[i][j] << "\n";
		}
	*/

	three_low ( lam, n);

	delete[] A; delete[] R; delete[] rho;

	return 0;
	
}