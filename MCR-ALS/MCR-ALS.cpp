#include "MCR-ALS.hpp"

MCR_ALS::MCR_ALS(int maxIter):_maxIter{maxIter}
{
    regr = new OLS();
}

MCR_ALS::~MCR_ALS()
{
    delete regr;
}

MCR_ALS& MCR_ALS::changeRegressorType(const char* RegressorType)
{
    if(RegressorType == "OLS")
    {
        delete regr;
        regr = new OLS();
    }
    else if (RegressorType == "NNLS")
    {
        delete regr;
        regr = new NNLS();
    }
    
    std::cout << "changeTo" << RegressorType << std::endl;

    return *this;
}

LinerRegression& MCR_ALS::fit(double* D, int D_colunm, int D_row, double* C, double* ST, int targetMatrix_colunm, int targetMatrix_row)
{
    if(C == nullptr && ST == nullptr)return *regr;
    double* D_T = new double[D_colunm * D_row];
    double* Matrix_T = new double[targetMatrix_colunm * targetMatrix_row];
    double* cal_D = new double[D_colunm * D_row];

    for(int i = 0; i < _maxIter; i++)
    {
        if(C == nullptr)
        {
            _ST = ST; 
            transpose(D_T,D_colunm,D_row,D_T);
            transpose(ST,targetMatrix_colunm,targetMatrix_row,Matrix_T);
            regr->fit(Matrix_T,targetMatrix_row,targetMatrix_colunm,D_T,D_row,D_colunm);
            _C = regr->getCoef();
        }
        else if (ST == nullptr)
        {
            _C = C;
            regr->fit(C,targetMatrix_colunm,targetMatrix_row,D,D_colunm,D_row);
            _ST = regr->getCoef();
        }

        cal_D = dot(_C,D_colunm,targetMatrix_row,_ST,targetMatrix_row,D_row,cal_D);
        if(meanSquareError(cal_D,D,D_colunm,D_row) > 1e-5)break;
    }

    delete D_T;
    delete Matrix_T;
    delete cal_D;

    return *regr;
}