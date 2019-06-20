#include<math.h>
#include<stdlib.h>
#include<iostream>
#include"randlib_par.h"
#include"Functions.h"

using namespace std;

/*find min values of a matrix*/
double Find_Min_Value(double**X, int n_row, int n_col) {

	int i, j;
	double min;

	min = X[0][0];

	for (i = 0; i < n_row; i++) {

		for (j = 0; j < n_col; j++) {

			if (X[i][j] < min) min = X[i][j];

		}

	}
	return min;
}

/*fine mean value of array of integers*/
double Mean_Value_Int(int*X, int**Y, int steps, int runs, int n) {

	double mean;
	int r, total;

	total = 0;

	for (r = 0; r < runs; r++) {

		if (Y[steps][r] != 0) {

			total += X[r];
		}
	}

	mean = total / double(n);

	return mean;

}

/*fine standard deviation of array of integers*/
double Std_Dev_Int(int*X, int**Y, int steps, double mean, int runs, int n) {

	int r;
	double diff, var, sd, sum, x;

	sum = 0;

	for (r = 0; r < runs; r++) {

		if (Y[steps][r] != 0) {

			diff = X[r] - mean;
			x = pow(diff, 2);
			sum += x;

		}
	}

	var = sum / double(n);
	sd = sqrt(var);

	return sd;
}

/*fine standard deviation of array of doubles*/
double Std_Dev_Doub(double*X, int**Y, int steps, double mean, int runs, int n) {

	int r;
	double diff, var, sd, x, sum;

	sum = 0;

	for (r = 0; r < runs; r++) {

		if (Y[steps][r] != 0) {

			diff = X[r] - mean;
			x = pow(diff, 2);
			sum += x;

		}
	}

	var = sum / double(n);
	sd = sqrt(var);

	return sd;
}

