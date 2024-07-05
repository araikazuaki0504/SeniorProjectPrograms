#pragma once

#define _USE_MATH_DEFINES

#include <iostream>
#include <math.h>
#include <iomanip>
#include <limits>
#include <algorithm>
#include <random>
#include <memory>

void showMatrix(double* Matrix, int colunm, int row, bool allShow = false);
void showMatrix(bool* Matrix, int colunm, int row);
double norm(double* vector, int vector_size);
double frobenius_norm(double* Matrix, int Matrix_colunm, int Matrix_row);
double determinant(double* Matrix, int Matrix_colunm_, int Matrix_row_);
double* transpose(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* Matrix_t);
double* add(double* Matrix_A, double* Matrix_B, int Matrix_colunm, int Matrix_row, double* Matrix);
double* sub(double* Matrix_A, double* Matrix_B, int Matrix_colunm, int Matrix_row, double* Matrix);
double dot(double* vector_A, double* vector_B, int vectorSize);
double* product(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_, double* Matrix);
bool QR_Decompose(double* Matrix, int Matrix_colunm, int Matrix_row, double* Q_Matrix,double* R_Matrix);
bool HouseHolderTransform(double* Matrix, int Matrix_colunm, int Matrix_row, double* R_Matrix,double* Q_Matrix);
bool checkDiagonal(double* Matrix, int Matrix_colunm, int Matrix_row, double epsilon = 1e-10);
double* eig(double* Matrix, int Matrix_colunm, int Matrix_row,double* eigValue);
double* inverse_lower_triangular_matrix(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix);
double* inverse_upper_triangular_matrix(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix);
bool LU_Decompose(int mode, double* Matrix, int Matrix_colunm_, int Matrix_row_, double* L, double* U, double* Q);
double* inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix);
double* Pseudo_inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* pinv_Matrix);
double* equal_solve(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* vector, int vector_size_, double* solution_vector);
void SVD(double* Matrix, int Matrix_colunm, int Matrix_row, double* U, double* S, double* V,int maxIter=50);
void NMF(double* Matrix, int Matrix_colunm, int Matrix_row, int n_components, double* W, double* H, double epsilon = 1e-10);
int Matrix_rank(double* Matrix, int Matrix_colunm, int Matrix_row, bool DoSVD = false, double* U = nullptr, double* S = nullptr, double* V = nullptr);
int eig_rank(double* eig_value, int eig_value_size);
std::vector<std::pair<int,int>> get_transposition(std::vector<int> index);