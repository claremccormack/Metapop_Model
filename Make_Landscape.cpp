#include<stdio.h>
#include<math.h>
#include<stdlib.h>
#include<iostream>
#include<iomanip>
#include<cstdlib>
#include<cmath>
#include<random>
#include<algorithm>
#include<time.h>
#include<fstream>
#include<string>
#include<sstream>
#include<direct.h>
#include<windows.h>
#include"randlib_par.h"
#include"Functions.h"
#include"Make_Landscape.h"

using namespace std;

random_device rd2;
mt19937 gen2(rd2());

/*make landscape with mean det_K, std. deviation sigma and covariane corr*/
void Make_Landscape(int capacity, double***K, double***old_K, int r, int n_row, int n_col, int alpha, double det_K, double sigma, double corr) {

	double sig, z1, z2, mu, z3;
	int i, j,k,sites;

	sites = n_row*n_col;

	if (capacity == 2) {

		/*standard deviation for landscape coefficient of variation=sigma*/
		sig = sigma*det_K;

		/*find mean and sd of the normal distribution corresponding to a logNormal with mean det_K and std dev sig*/
		z1 = pow(det_K, 2);
		z2 = sqrt(pow(sig, 2) + z1);
		mu = log(z1 / z2);/*mean of corresponding Normal*/
		z3 = sqrt(log(1 + (pow(sig, 2) / z1)));/*sd of corresponding Normal*/

	}

	else { mu = 0; z3 = 0; }

	normal_distribution<double> normal(mu, z3);

	/*set carrying capacity of each patch*/

	switch (capacity)

	{

	case 1:

		if (alpha == 1) {

			for (i = 0; i < n_row; i++) {

				for (j = 0; j < n_col; j++) {

					K[r][i][j] = det_K;

				}
			}
		}

		else if (alpha > 1) {

			Aggregate_Values(n_row, n_col, old_K[r], K[r]);
		}

		break;

	case 2:

		if (alpha == 1) {

			sig = sigma * det_K;
			double init_b1, init_b2;
			double*params;
			params = new double[2];

			for (i = 0; i < n_row; i++) {

				for (j = 0; j < n_col; j++) {

					if (corr == 0 || sigma == 0) { K[r][i][j] = exp(normal(gen2)); } /*carrying capacity of each patch drawn from logNormal distribution with mean at deterministic equilibrium*/

					else K[r][i][j] = snorm(); /*generate carrying capacity from std normal before adjusting to give correlated values*/

				}
			}

			if (corr > 0 && sigma > 0) {

				double**D;
				double**C;
				double**L;
				double*temp;
				double*list;

				temp = new double[sites];
				list = new double[sites];

				D = new double*[sites];
				C = new double*[sites];
				L = new double*[sites];

				for (i = 0; i < sites; i++) { D[i] = new double[sites]; C[i] = new double[sites]; L[i] = new double[sites]; }

				Create_Distance_Matrix(n_row, n_col, D);/*matrix describing distance between each patch*/
				Create_Correlation_Matrix(n_row, n_col, corr, D, C); /*matrix describing correlation between each patch*/
				Cholesky_Decomposition(n_row, n_col, C, L); /*Cholesky decomposition of correlation matrix*/

				k = 0;

				/*put carrying capacities in a (sites x 1) dim array*/

				for (i = 0; i < n_row; i++) {

					for (j = 0; j < n_col; j++) {

						list[k] = K[r][i][j];

						k = k + 1;
					}
				}

				/*multiply carrying capacities by Cholesky decomposition matrix*/

				for (i = 0; i < sites; i++) {

					temp[i] = 0;

					for (j = 0; j < sites; j++) {

						temp[i] += L[i][j] * list[j];

					}

				}

				k = 0;

				/*replace old carrying capacity values with adjusted values (temp)*/

				for (i = 0; i < n_row; i++) {

					for (j = 0; j < n_col; j++) {

						K[r][i][j] = exp(temp[k]);

						k = k + 1;
					}

				}
	
				for (i = 0; i < sites; i++) { delete D[i]; }
				delete[] D;

				for (i = 0; i < sites; i++) { delete C[i]; }
				delete[] C;

				for (i = 0; i < sites; i++) { delete L[i]; }
				delete[] L;

				delete[] temp;
				delete[] list;
			}
				/*find paramters a and b such that applying the transformation a*pow(K,b) gives a set of values with fixed mean, sd and sum*/

				/*intial guesses b*/
				init_b1 = 1;
				init_b2 = 5;

				Num_Solution_Interpolation(K[r], params, n_row, n_col, det_K, sig, init_b1, init_b2);

				for (i = 0; i < n_row; i++) {

					for (j = 0; j < n_col; j++) {

						K[r][i][j] = params[0] * pow(K[r][i][j], params[1]);
					}
				}

				delete[] params;
		}

		else if (alpha > 1) {

			Aggregate_Values(n_row, n_col, old_K[r], K[r]);
		}

		break;

	}

}

/*aggregate values across 4 grid squares*/
void Aggregate_Values(int n_row, int n_col, double**X, double**Y) {

	int x, y;
	double sum;

	x = 0;
	y = 0;

	for (int i = 0; i < (2 * n_row); i += 2) {

		for (int j = 0; j < (2 * n_col); j += 2) {

			sum = 0;
			sum += X[i][j] + X[i][j + 1] + X[i + 1][j] + X[i + 1][j + 1];

			Y[x][y] = sum;

			y = y + 1;

		}
		x = x + 1;
		y = 0;
	}
}

/*Find values a and b such that expected value of aX^b is mu, and variance of aX^b is sigma using linear interpolation*/
void Num_Solution_Interpolation(double**X, double*result, int nrow, int ncol, double mu, double sigma, double lower, double upper) {

	int sites;
	double sum, value_lower, value_upper, temp, value_temp, temp_product, product, epsilon, b, tolerance;
	sites = nrow * ncol;
	sum = 0;
	epsilon = 1;
	tolerance = .000001;

	value_lower = Func_Interpolation(nrow, ncol, X, mu, sigma, lower); /*value of function we wish to solve for lower guess of b*/
	value_upper = Func_Interpolation(nrow, ncol, X, mu, sigma, upper); /*value of function we wish to solve for upper guess of b*/

	product = value_lower * value_upper; /*product of initial guesses*/

	if (product == 0) {

		/*we have solution*/

		if (value_lower == 0)

			result[1] = lower;

		else
			result[1] = upper;

		for (int i = 0; i < nrow; i++) {

			for (int j = 0; j < ncol; j++) {

				sum += pow(X[i][j], result[1]);
			}
		}

		result[0] = sites * mu / sum;

	}

	else {

		/*if product>0, adjust guesses of solution until product<0*/

		while (product > 0) {

			if (value_lower < 0) {

				temp = upper + .5;

				value_temp = Func_Interpolation(nrow, ncol, X, mu, sigma, temp);

				temp_product = value_lower * value_temp;

				upper = temp;

			}

			else {

				temp = lower - .5;

				value_temp = Func_Interpolation(nrow, ncol, X, mu, sigma, temp);

				temp_product = value_temp * value_upper;

				lower = temp;

			}

			product = temp_product;

		}

		while (fabs(epsilon) > tolerance) {

			/*when product<0, set b as follows and continue until we find value of b such that the value of the function at b lies within specified tolerance range*/

			b = (lower*value_upper - upper * value_lower) / (value_upper - value_lower);
			
			epsilon = Func_Interpolation(nrow, ncol, X, mu, sigma, b); /*value of function for b*/

			if (fabs(epsilon) > tolerance) {

				if (epsilon > 0) {

					upper = b;
					value_upper = Func_Interpolation(nrow, ncol, X, mu, sigma, upper);
				}

				else

					lower = b;
					value_lower = Func_Interpolation(nrow, ncol, X, mu, sigma, lower);
			}

		}
	}

	result[1] = b; /*final value of b*/

	for (int i = 0; i < nrow; i++) {

		for (int j = 0; j < ncol; j++) {

			sum += pow(X[i][j], result[1]); /*calculate value of a for given b*/

		}
	}

	result[0] = sites * mu / sum; /*value of a*/

}

/*function we need to solve to give values with a fixed mean, sd, and sum*/
double Func_Interpolation(int nrow, int ncol, double**X, double mu, double sigma, double b) {

	int i, j, sites;
	double sum1, sum2, sum3, result;

	sites = nrow * ncol;

	sum1 = 0;
	sum2 = 0;
	sum3 = 0;

	for (i = 0; i < nrow; i++) {

		for (j = 0; j < ncol; j++) {

			sum1 += pow(X[i][j], 2 * b);

			sum2 += pow(X[i][j], b);
			sum3 = pow(sum2, 2);

		}
	}

	result = (sum1 / sum3) - (1 / double(sites))*((pow(sigma, 2) / pow(mu, 2)) + 1);

	return result;
}

/*function to calculate distance between each pair of sites*/
void Create_Distance_Matrix(int nrow, int ncol, double**D) {

	int i, j, r, c, sites;
	double c1, c2, r1, r2;
	sites = nrow * ncol;

	double**patches;
	patches = new double*[sites];

	for (i = 0; i < sites; i++) { patches[i] = new double[3]; }

	for (i = 0; i < sites; i++) { for (j = 0; j < sites; j++) D[i][j] = 0; }

	i = 0;

	for (r = 0; r < nrow; r++) {

		for (c = 0; c < ncol; c++) {

			patches[i][0] = r;
			patches[i][1] = c;

			i = i + 1;
		}
	}
	for (i = 0; i < sites; i++) {

		for (j = 0; j < sites; j++) {

			if (i == j) { D[i][j] = 0; }

			else

			{
				if (D[i][j] == 0) {

					r1 = patches[i][0];
					r2 = patches[j][0];

					c1 = patches[i][1];
					c2 = patches[j][1];

					D[i][j] = sqrt(pow(double(r1 - r2), 2) + pow(double(c1 - c2), 2));

					D[j][i] = D[i][j];
				}
			}
		}
	}
}

/*create correlation matrix*/
void Create_Correlation_Matrix(int nrow, int ncol, double corr, double**D, double**C) {

	int i, j, sites;
	sites = nrow * ncol;

	for (i = 0; i < sites; i++) { for (j = 0; j < sites; j++) C[i][j] = 0; }

	for (i = 0; i < sites; i++) {

		for (j = 0; j < sites; j++) {

			C[i][j] = exp(-corr * D[i][j]);

		}
	}
}

/*find Cholesky Decomposition of correlation matrix*/
void Cholesky_Decomposition(int nrow, int ncol, double**X, double**L) {

	int i, j, k, sites;
	double sum;
	sites = nrow * ncol;

	for (i = 0; i < sites; i++) {

		for (j = 0; j < sites; j++) {

			L[i][j] = 0;

		}
	}

	for (i = 0; i < sites; i++) {

		for (j = 0; j < (i + 1); j++) {

			if (i == j) {

				sum = 0;

				for (k = 0; k < i; k++) {

					sum += pow(L[i][k], 2);
				}

				L[i][j] = sqrt(X[i][i] - sum);
			}

			else {

				if (i > j) {

					sum = 0;

					for (k = 0; k < j; k++) {

						sum += L[i][k] * L[j][k];
					}

					L[i][j] = (1 / L[j][j])*(X[i][j] - sum);
				}
			}
		}
	}

}