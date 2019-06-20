#ifndef FUNCTIONS_H
#define FUNCTIONS_H

double Find_Min_Value(double**X, int n_row, int n_col);
double Mean_Value_Int(int*X, int**Y, int steps, int runs, int n);
double Std_Dev_Int(int*X, int**Y, int steps, double mean, int runs, int n);
double Std_Dev_Doub(double*X, int**Y, int steps, double mean, int runs, int n);

#endif

