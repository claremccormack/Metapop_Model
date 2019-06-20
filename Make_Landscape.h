#ifndef MAKE_LANDSCAPE_H
#define MAKE_LANDSCAPE_H

void Make_Landscape(int capacity, double***K, double***old_K,int r, int n_row, int n_col, int alpha, double det_K, double sigma, double corr);
void Aggregate_Values(int n_row, int n_col, double**X, double**Y);
void Num_Solution_Interpolation(double**X, double*result, int nrow, int ncol, double mu, double sigma, double lower, double upper);
double Func_Interpolation(int nrow, int ncol, double**X, double mu, double sigma, double b);
void Create_Distance_Matrix(int, int, double**);
void Create_Correlation_Matrix(int, int, double, double**, double**);
void Cholesky_Decomposition(int, int, double**, double**);

#endif
