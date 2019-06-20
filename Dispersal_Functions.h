#ifndef DISPERSAL_FUNCTIONS_H
#define DISPERSAL_FUNCTIONS_H

int **Create_Offset_Matrix(int dist);
int Find_Number_of_Neighbours(int i, int j, int n_row, int n_col, int max_dist);
void Create_Neighbours_Matrix(int i, int j, int n_row, int n_col, int n_total, int kernel, double mean_dist, int max_dist, double****neighbours, int**offset, double rate);
void Dispersal(int i, int j, int n_row, int n_col, int runs, int*result, int**population, int current_pop, double****neighbours, int**neighbours_total, double****dispersed, int***no_neighbours, int**no_dist, int max_dist, double mu_m, int t, double dt, double rate);

#endif