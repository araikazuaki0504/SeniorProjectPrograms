#pragma once

#include <iostream>
#include <cstring>

#include "Regressors.hpp"

class MCR_ALS{
    public:
        MCR_ALS(int maxIter = 50);
        ~MCR_ALS();
        MCR_ALS& changeRegressorType(std::string RegressorType);
        MCR_ALS& fit(double* D, int D_colunm, int D_row, double* C = nullptr, double* ST = nullptr, int targetMatrix_colunm = 0, int targetMatrix_row = 0);
        double* get_C();
        double* get_St();
    private:
        const int tol_n_above_min = 10;
        const int tol_n_increase = 10;
        const double tol_increase = 0.0;
        LinerRegression* regr = nullptr; 
        double* _C = nullptr;
        double* _ST = nullptr;
        int _maxIter = 0;
        std::vector<double> errs;
        bool ismin_err(double val);
};