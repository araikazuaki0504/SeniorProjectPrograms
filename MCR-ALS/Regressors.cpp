#include "Regressors.hpp"

LinerRegression::LinerRegression()
{
    
}

LinerRegression::~LinerRegression()
{
    if(result != nullptr)delete[] result;
    if(coef != nullptr)delete[] coef;
}

double* LinerRegression::getCoef(bool DoTranspose)
{
    if(coef == nullptr)coef = new double[result->_solution_colunm_size * result->_solution_row_size];

    if(DoTranspose)transpose(result->_lstsq_solution,result->_solution_colunm_size,result->_solution_row_size,coef);
    else copyMatrix(result->_lstsq_solution,result->_solution_colunm_size,result->_solution_row_size,coef);
    
    return coef;
}

void OLS::fit(double* Matrix_A, int Matrix_A_colunm, int Matrix_A_row,double* Matrix_B, int Matrix_B_colunm, int Matrix_B_row)
{
    if(result == nullptr)result = new lstsq_result(Matrix_B_colunm,Matrix_B_row);
    lstsq(Matrix_A,Matrix_A_colunm,Matrix_A_row,Matrix_B,Matrix_B_colunm,Matrix_B_row,result);
    //result->showData();
}

void NNLS::fit(double* Matrix_A, int Matrix_A_colunm, int Matrix_A_row,double* Matrix_B, int Matrix_B_colunm, int Matrix_B_row)
{
    if(result == nullptr)result = new lstsq_result(Matrix_B_colunm,Matrix_B_row);
    nnls(Matrix_A,Matrix_A_colunm,Matrix_A_row,Matrix_B,Matrix_B_colunm,Matrix_B_row,result);
    //result->showData();
}
