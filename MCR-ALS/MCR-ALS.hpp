#pragma once

#include <iostream>

#include "Regressors.hpp"

class MCR_ALS{
    public:
        MCR_ALS(int maxIter = 50);
        ~MCR_ALS();
        MCR_ALS& changeRegressorType(const char* RegressorType);
        LinerRegression& fit(double* D, int D_colunm, int D_row, double* C = nullptr, double* ST = nullptr, int targetMatrix_colunm = 0, int targetMatrix_row = 0);
    private:
        LinerRegression* regr = nullptr; 
        double* _C = nullptr;
        double* _ST = nullptr;
        int _maxIter = 0;
};