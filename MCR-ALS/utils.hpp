#pragma once
#include <iostream>
#include <vector>
#include <algorithm>
#include <limits>
#include <math.h>

#include "linear_algebra.hpp"

using namespace std;

typedef struct lstsq_result
{
    double* _lstsq_solution;
    double* _residues = nullptr;//-1のときは無効
    int _rank = -1;//-1のときは無効
    double* _single_value;
    int _solution_colunm_size;
    int _solution_row_size;
    int _current_solution_row_index = 0;
    lstsq_result(int solution_colunm_size, int solution_row_size);
    ~lstsq_result();
    void showData();
}lstsq_result;

bool check_P(bool* P, int P_size);
bool check_w(double* w, bool* P, int P_size, double tol);
bool check_s(double* s, bool* P, int P_size);
int maxValueIndex(double* w, bool* P, int P_size);
double get_miniumAlphaWithInds(double* x, double* s, bool* inds, int inds_size);
int lstsq(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_ ,lstsq_result* cal_result);
int nnls(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row, lstsq_result* cal_result);
double meanSquareError(double* calMatrix, double* acutuallyMatrix, int Matrix_colunm, int Matrix_row);