#include <cmath>
#include <iostream>

using namespace.std;

void jacobi_method ( double ** A, double ** R int n )
{
//Settup of eigenvector
	for ( int i = 0; i < n; i++){
		for ( int j = 0; j < n; j++){
			if ( i == j){
				R[i][j] = 1.0;
			}
			else{
				R[i][j] = 0.0
			}
		}
	}
	int k, l;

	double epsilon = 1.0e-8;
	double max_number_iterations = (double) n * (double) n * (double) n;
	int iterations = 0;
	double max_offdiag = maxoffdiag ( A, &k, &l, n );

	while ( fabs(max_offdiag) > epsilon && (double) iterations < max_number_iterations ){
		max:offdiag = maxoffdiag( A, &k, &l, n );
		rotate (A, R, k, l, n);
		iterations++;
	}
	cout <<"Number of iterations: " << iterations << "\n";
	return;
}


