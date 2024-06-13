#pragma once
#include "utils.hpp"

class MCR_ALS;

class LinerRegression
{
    public:
        LinerRegression();
        ~LinerRegression();
        double* getCoef();
    protected:
        lstsq_result* result = nullptr;
        double* coef = nullptr;
        virtual void fit(double* Matrix_A, int Matrix_A_colunm, int Matrix_A_row,double* Matrix_B, int Matrix_B_colunm, int Matrix_B_row){};
    friend MCR_ALS;
};

class OLS : public LinerRegression
{
    void fit(double* Matrix_A, int Matrix_A_colunm, int Matrix_A_row,double* Matrix_B, int Matrix_B_colunm, int Matrix_B_row) override;
};

class NNLS : public LinerRegression
{
    void fit(double* Matrix_A, int Matrix_A_colunm, int Matrix_A_row,double* Matrix_B, int Matrix_B_colunm, int Matrix_B_row) override;
};