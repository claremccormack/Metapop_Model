#include<math.h>
#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include"randlib_par.h"
#include"Dispersal_Functions.h"

using namespace std;

/*create matrix describing location of neighbours of each patch within distance 'dist'*/
int **Create_Offset_Matrix(int dist) {

	static int **offset;
	int *max;

	int a, max_total, d, k, r1;

	max = new int[dist + 1];

	/*cout max no. of neighbours a patch can have within each distance*/

	max[0] = 1;

	for (k = 1; k < dist + 1; k++) {

		max[k] = 8 * k + max[k - 1]; /*max eight neighbours within distance 1*/
	}

	max_total = max[dist]; /*max total no of neighbours within distance dist*/

	offset = new int*[max_total];
	for (k = 0; k < max_total; k++) { offset[k] = new int[3]; }

	/*create offset matrix*/

	offset[0][0] = 0;
	offset[0][1] = 0;
	offset[0][2] = 0;

	for (d = 1; d < dist + 1; d++) {

		if (d == 1) r1 = 1;
		else r1 = max[d - 1];

		/*top row neighbours*/

		a = -d;

		int z = 2 * d + 1;

		for (k = r1; k < r1 + z; k++) {

			offset[k][0] = d;
			offset[k][1] = a;
			offset[k][2] = d;

			a = a + 1;

		}

		/*bottom row neighbours*/

		a = -d;

		for (k = r1 + z; k < r1 + 2 * z; k++) {

			offset[k][0] = d;
			offset[k][1] = a;
			offset[k][2] = -d;

			a = a + 1;

		}

		/*right side neighbours*/

		a = -(d - 1);

		for (k = r1 + 2 * z; k < r1 + 3 * z - 2; k++) {

			offset[k][0] = d;
			offset[k][1] = d;
			offset[k][2] = a;

			a = a + 1;

		}

		/*left side neighbours*/

		a = -(d - 1);

		for (k = r1 + 3 * z - 2; k < max[d]; k++) {

			offset[k][0] = d;
			offset[k][1] = -d;
			offset[k][2] = a;

			a = a + 1;

		}

	}

	delete[] max;

	return offset;

}

/*find total number of neighbours of patch (i,j) within distance 'dist'*/
int Find_Number_of_Neighbours(int i, int j, int n_row, int n_col, int max_dist) {

	int *n;

	int d, x1, x2, y1, y2, patches, n_total;

	n = new int[max_dist + 1];

	/*count how many neighbours the patch has at each distance*/

	int sum = 0;

	n[0] = 1;

	for (d = 1; d < max_dist + 1; d++) {

		if (i - d >= 0) x1 = d; else x1 = i;
		if (i + d <= n_row - 1) x2 = d; else x2 = n_row - 1 - i;

		if (j - d >= 0) y1 = d; else y1 = j;
		if (j + d <= n_col - 1) y2 = d; else y2 = n_col - 1 - j;

		patches = (x1 + x2 + 1)*(y1 + y2 + 1);

		if (d == 1) { n[d] = patches - 1; sum += n[d]; }

		else

		{
			n[d] = patches - 1 - sum;
			sum += n[d];
		}

	}

	n_total = 0;

	for (d = 0; d < max_dist + 1; d++) { n_total += n[d]; } /*total no. of neighbours*/

	return n_total;

}

/*create matrix giving neighbours of each patch, distance to neighbours and probability of dispersing to those neighbours*/
void Create_Neighbours_Matrix(int i, int j, int n_row, int n_col, int n_total, int kernel, double mean_dist, int max_dist, double****neighbours, int**offset, double rate) {

	int k, n, r, c, max_total;
	double a, b, eucl_d, x, f, temp1, temp2, temp3, temp4;
	float g;
	double pi = 22 / double(7);

	k = 0;
	int *max;
	max = new int[max_dist + 1];

	/*max no. of neighbours a patch can have within each distance*/

	max[0] = 1;

	for (k = 1; k < max_dist + 1; k++) {

		max[k] = 8 * k + max[k - 1]; /*max eight neighbours within distance 1*/
	}

	max_total = max[max_dist]; /*max total no of neighbours within distance dist*/

	/*Calculate Euclidean distance to each neighbour*/

	k = 0;

	for (n = 0; n < max_total; n++) {

		r = i + offset[n][1];
		c = j + offset[n][2];

		if (r >= 0 && r <= n_row - 1 && c >= 0 && c <= n_col - 1) {

			if (r == i && c == j) { eucl_d = 0.5214; } /*average distance between two random points in a square*/

			else {
				x = pow(double(r - i), 2) + pow(double(c - j), 2);
				eucl_d = sqrt(x);
			}

			switch (kernel) {

			case 0:

				/*negative exponential kernel*/
				a = double(mean_dist) / 2;
				f = (1 / double(2 * pi*a*a))*exp(-eucl_d / a);
				break;

			case 1:

				/*Normal*/
				a = double(2 * mean_dist) / sqrt(pi);
				g = float((1 / double(pi*a*a))*exp(-pow(eucl_d / a, 2)));
				f = double(floorf(g * 100000000000000) / 100000000000000);
				break;


			case 2:

				/*exponenetial power kernel with shape=4*/
				a = mean_dist * sqrt(pi) / tgamma(double(3) / 4);
				g = float((4 / double(2 * pi*sqrt(pi)*pow(a, 2)))*exp(-pow(eucl_d / a, 4)));
				f = double(floorf(g * 100000000000000) / 100000000000000);
				break;

			case 3:

				/*inverse power law*/
				b = 3.5;
				a = double((mean_dist)*(b - 3)) / 2;
				f = ((b - 2)*(b - 1) / double(2 * pi*pow(a, 2)))*pow(1 + (eucl_d / a), -b);
				break;

			case 4:

				/*2d t-distribution with shape=2*/
				b = 2;
				a = (2 * mean_dist*tgamma(b - 1)) / double(sqrt(pi)*tgamma((b - (double(3) / 2))));
				f = ((b - 1) / double(pi*pow(a, 2)))*pow(1 + (pow(eucl_d / a, 2)), -b);
				break;

			}

			neighbours[i][j][k][0] = eucl_d;
			neighbours[i][j][k][1] = f;
			neighbours[i][j][k][2] = r;
			neighbours[i][j][k][3] = c;

			k += 1;

		}

	}

	/*sort in ascending order based on distance value*/

	for (n = 0; n < n_total; n++) {

		for (k = 0; k < n_total - 1; k++) {

			if (neighbours[i][j][k][0] > neighbours[i][j][k + 1][0]) {

				temp1 = neighbours[i][j][k][0];
				temp2 = neighbours[i][j][k][1];
				temp3 = neighbours[i][j][k][2];
				temp4 = neighbours[i][j][k][3];

				neighbours[i][j][k][0] = neighbours[i][j][k + 1][0];
				neighbours[i][j][k][1] = neighbours[i][j][k + 1][1];
				neighbours[i][j][k][2] = neighbours[i][j][k + 1][2];
				neighbours[i][j][k][3] = neighbours[i][j][k + 1][3];

				neighbours[i][j][k + 1][0] = temp1;
				neighbours[i][j][k + 1][1] = temp2;
				neighbours[i][j][k + 1][2] = temp3;
				neighbours[i][j][k + 1][3] = temp4;
			}
		}

	}
}

/*calculate number of adult mosquitoes dispersing and number of adult mosquito deaths*/
void Dispersal(int i, int j, int n_row, int n_col, int runs, int*result, int**population, int current_pop, double****neighbours, int**neighbours_total, double****dispersed, int***no_neighbours, int**no_dist, int max_dist, double mu_m, int t, double dt, double rate) {

	int  d, k, r, c, sites, n_total, num_dist, final_d;
	int out_m, deaths_m, dispersed_m;
	int thread = 0;
	double a, p;

	sites = n_row*n_col;
	n_total = neighbours_total[i][j];/*total no. of neighbours of (i,j) within dispersal range*/

	double n1, n2, n3, n4, n7, n8;
	int n6;

	final_d = 0;

	/*total hazard of leaving site (i,j) */

	double sum_h = rate;
	num_dist = no_dist[i][j];

	double p_m = 1-exp(-(mu_m + sum_h)*dt);

	/*count no. of deaths and total no of mosquitoes dispersing*/

	out_m = ignbin_mt(long(current_pop), p_m, thread);
	deaths_m = ignbin_mt(long(out_m), mu_m / (mu_m + sum_h), thread);
	dispersed_m = out_m - deaths_m ;

	result[0] = deaths_m;
	result[1] = dispersed_m;

	for (d = 0; d < num_dist; d++) { no_neighbours[i][j][d] = 0; }

	/*no of sites at each possible distance*/

	dispersed[i][j][0][0] = 0;
	dispersed[i][j][0][1] = neighbours[i][j][0][1]; /*hazard of dispersing distance zero*/

	d = 1;

	for (k = 1; k < n_total; k++) {

		n1 = neighbours[i][j][k][0];
		n2 = neighbours[i][j][k - 1][0];
		n3 = neighbours[i][j][k][1];

		no_neighbours[i][j][d - 1] += 1;/*count no. of neighbours at distance (d-1)*/

		if (n1 != n2) {

			dispersed[i][j][d][0] = n1; /*value of distance*/
			dispersed[i][j][d][1] = n3; /*hazard of dispersing to that distance*/

			d = d + 1;
		}
	}

	/*count no. of neighbours at final distance*/

	a = dispersed[i][j][num_dist - 1][0]; /*distance value*/
	no_neighbours[i][j][num_dist - 1] = 0;

	for (k = 0; k < n_total; k++) {

		n4 = neighbours[i][j][k][0];

		if (n4 == a) { no_neighbours[i][j][num_dist - 1] += 1; }

	}

	/*no. of mosquitoes dispersing to each of these distances-competing hazards*/

	if (sum_h > 0) { p = dispersed[i][j][0][1] / sum_h; }
	else { p = 0; }

	if (p > 1) { p = 1; }

	dispersed[i][j][0][2] = ignbin_mt(long(dispersed_m), p, thread); /*no. dispersing within same patch*/
	dispersed_m = dispersed_m - int(dispersed[i][j][0][2]);
	sum_h = sum_h - dispersed[i][j][0][1];

	for (d = 1; d < num_dist; d++) {

		if (dispersed_m > 0) {

			n6 = no_neighbours[i][j][d];
			n7 = dispersed[i][j][d][1];
			n8 = 0;

			if (sum_h > 0) { p = (n6 * n7) / sum_h; }

			else { p = 0; }

			if (p > 1) { p = 1; }

			n8 = ignbin_mt(long(dispersed_m), p, thread);/*no. dispersing distance d*/

			dispersed[i][j][d][2] = n8;

			sum_h = sum_h - (n6*n7);

			dispersed_m = dispersed_m - int(n8);
		}

		else {

			if (final_d == 0) { final_d += d; }

			dispersed[i][j][d][2] = 0;
		}
	}

	/*For mosquitoes dispersing distance d, randomly select patches at distance d to move to (including (i,j))*/

	int k1, k2, x, z, m;
	double y;

	for (d = 0; d < num_dist; d++) {

		if (d == 0) { k1 = 0; k2 = 0; }

		else {

			k1 = k1 + no_neighbours[i][j][d - 1];

			k2 = k2 + no_neighbours[i][j][d];
		}

		for (m = 0; m < dispersed[i][j][d][2]; m++) {

			y = floor(ranf()*(k2 + 1));
			z = int(y);
			x = (z % (k2 - k1 + 1)) + k1;

			r = int(neighbours[i][j][x][2]);
			c = int(neighbours[i][j][x][3]);

			population[r][c] += 1;

		}
	}
}