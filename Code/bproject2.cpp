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

//Function to find the three lowest eigenvalues of a set
double *three_low ( double * v, int n)
{
//v is array to find the lowest values of, n is length of array. Randomly chosen variables xyz here

	int x, y, z;
	double total_time;
	double * smallest = new double [3];
	double *eigpos = new double[3];
	clock_t s, f;
	s = clock(); 
	for ( int i = 0; i < 3; i++) smallest[i] = v[i];

	for ( int i = 0; i < n; i++ ) {
		
		if ( fabs(v[i]) < fabs(smallest[0]) ) {
			smallest[0] = v[i];
			x = i;
			//cout << v[i] << " k = " << i << "\n";
		}
	}
	eigpos[0] = x;

	for ( int i = 0; i < n; i++ ) {
		if ( ( fabs(v[i]) < fabs(smallest[1]) ) && i != x ) {
			smallest[1] = v[i];
			y = i;
			//cout << v[i] << " l = " << i << "\n";
		}
	}
	eigpos[1] = y;

	for ( int i = 0; i < n; i++ ) {
		if ( ( fabs(v[i]) < fabs(smallest[2]) ) && ( i !=x ) && i != y) {
			smallest[2] = v[i];
			z = i;
			//cout << v[i] << " m = " << i << "\n";
		}
	}
	eigpos[2] = z;

	for ( int i = 0; i < 3; i++) std::cout << smallest[i] << "\n";
	f = clock();
	total_time = ( ( double ) ( f - s ) / CLOCKS_PER_SEC );
	printf("Time spent on extracting eigenvalues:  %1.3e \n", total_time);
	return eigpos;
}

//Unit tests for eigenvalues and ortogonality of basis R
void unit_test(double **R){
	double **Matr, **Trans, **Prod;
	int n_test = 5, k, l;

	Matr = new double*[n_test];   //Test matrix
	Trans = new double*[n_test];  //Tranpose of the test matrix
	Prod = new double*[n_test];   //Multiplication of the two above

	for(int i = 0; i < n_test; i++){
		Matr[i] = new double [n_test];
		Trans[i] = new double [n_test];
		Prod[i] = new double [n_test];
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

	//PERFOMING THE TRANSFORMATION
	jacobi_method(Matr, R, n_test);


	printf("-------------------------------------------\nUNIT TEST\n\n");



	//------------------------------------------------
	//CHECK FOR ORTOGONALITY
	//----------------------------------------------
	int c, d, e;
	double eps = 1e-6, sum = 0;



    for( c = 0 ; c < n_test ; c++ )
    {
       for( d = 0 ; d < n_test ; d++ )
       {
          Trans[d][c] = R[c][d];
       }
    }
 
	// Multiplication of R and it's transpose
    for ( c = 0 ; c < n_test ; c++ )
    {
    	for ( d = 0 ; d < n_test ; d++ )
        {
        	for ( e = 0 ; e < n_test ; e++ )
        	{
            	sum = sum + R[c][e]*Trans[e][d];
        	}
 
        	Prod[c][d] = sum;
        	sum = 0;
        }
    }
 
 
    //Checks if the result is the identity matrix with a small margin for error, eps
    for ( c = 0 ; c < n_test ; c++ )
    {
        for ( d = 0 ; d < n_test ; d++ )
        {
        	if ( c == d )
            {
                if ( Prod[c][d] <= (1 - eps) || Prod[c][d] >= (1 + eps))
                	break;
            }
            else
            {
                if ( Prod[c][d] <= (- eps) || Prod[c][d] >= eps)
                	break;
            }
        }
        if ( d != n_test )
        	break;
    }
    if ( c != n_test )
    	printf("Orthogonality of the basis R is not preserved.\n\n");
    else
    	printf("Orthogonality of the basis R is preserved.\n\n");
    //-------------------------------------------------

    printf("The known eigenvalues for the sample matrix Matr:\n");
    printf("1: -4.16352, 2: -0.827046 3: 1.0, 4: 2.82705, 5: 6.16352\n");
    printf("Computed eigenvalues:\n");


	for(int i = 0; i < n_test; i++){
		printf("Eigenvalue %i: %f\n", i+1, Matr[i][i]);
	}
	printf("-------------------------------------------\n");

	delete[] Matr; delete[] Trans; delete[] Prod;

	return;

}

//Function for calculating the potential for the non interacting case (coloumb = 0.0)
// and the interacting case (coloumb = 1.0). Takes omega, r (aka rho) and one of the two coloumbs.
double potential(double omega, double r, double coloumb )
{
	double V;
	if ( coloumb == 0.0 ) V = r*r;
	else if ( coloumb == 1.0 ) V = omega*omega*r*r + 1.0/(double)r;
	return V;
}

//Computes the square of each element of the transformed basis R -> U
double *usquared( double **R, int n, int t)
{
	double * u2 = new double[n];
	for ( int j = 0; j < n; j++ ) {
		u2[j] = R[j][t]*R[j][t];
	
	}
	//printf("U(r) = %1.3e\n", u2[i] );
	return u2;
}

// When running the program you need to enter 4 values in the command line as well:
// n, rhoN, coloumb and file (listed below)
int main( int argc, char * argv[] )
{
	int n =  atoi(argv[1]);          // Dimension of arrays
	double rhoN =  atof(argv[2]);    // Value for rho_max
	double coloumb = atof(argv[3]);  // Value to determine if there is a coloumb potential present, 1.0 for yes, 0.0 for no.
	std::string file = ( argv[4] );  // Give filename for the output file
	double omega;    // Frequency
	//if coloumb = 1 we need a omega value, if not omega is set to 1.
	if(coloumb == 1.0){
		printf("Enter omega value: ");
		scanf("%lf", &omega);
	}
	else if(coloumb == 0.0){
		omega = 1.0;
	}
	else{
		printf("Error: Wrong value for coloumb.\n");
		printf("usage: argv[3] = 1.0 to include Coloumb potential, else argv[3] = 0.0.");
		return 0;
	}
	const char* filename = file.c_str();
	double rho0 = 1e-10;
	double h = ( rhoN - rho0 ) / (double) n;  //stepelength
	double* rho = new double[n];
	double total_time;
	clock_t start, finish;
	
	//define our tridiagonal matrix A and basis R
	double **A;
	double **R;
	A = new double*[n];
	R = new double*[n];
	for ( int i = 0; i < n; i++ ) {
		A[i] = new double [n]; 
		R[i] = new double [n];
	}

	// Computes the elements of A
	for ( int i = 0; i < n; i++) {
		rho[i] = rho0 + i*h;
		double Vi = potential( omega,  rho[i], coloumb );
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

	double * lam = new double[n];  // Array for the eigenvalues, lambda
	
	//Retrieving the eigenvalues from the diagonal
	for ( int i = 0; i < n; i++) {
		lam[i] = A[i][i];
	}



	double *t = three_low ( lam, n);    //Finding the indices of the three lowest eigenvalues
	double *sqru1 = usquared( R, n, t[0] );  //eigenvector for the lowest eigenvalue in the transformed basis R

	
	FILE *fp;

	fp = fopen(filename, "w+");

	//Writes the resulting u^2 to file for plotting in python
	for (int i = 0; i < n; i++){
		fprintf(fp, "%f\n", sqru1[i]);
	}

	fclose(fp);


	

	delete[] A; delete[] R; delete[] rho; delete [] t;
	delete[] sqru1;

	return 0;
	
}