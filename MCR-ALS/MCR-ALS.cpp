#include "MCR-ALS.hpp"

MCR_ALS::MCR_ALS(int maxIter):_maxIter{maxIter}
{
    regr = new OLS();
}

MCR_ALS::~MCR_ALS()
{
    delete[] _C;
    delete[] _ST;
    delete[] regr;
}

MCR_ALS& MCR_ALS::changeRegressorType(std::string RegressorType)
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

MCR_ALS& MCR_ALS::fit(double* D, int D_colunm, int D_row, double* C, double* ST, int targetMatrix_colunm, int targetMatrix_row)
{
    if(C == nullptr && ST == nullptr)return *this;
    if(_C != nullptr)delete[] _C;
    if(_ST != nullptr)delete[] _ST;

    double* D_T = new double[D_colunm * D_row];
    double* Matrix_T = new double[targetMatrix_colunm * targetMatrix_row];
    double* cal_D = new double[D_colunm * D_row];

    int n_increase = 0;
    int n_above_min = 0;
    errs.clear();

    _C  = (C != nullptr  ? copyMatrix(C,targetMatrix_colunm,targetMatrix_row,new double[targetMatrix_colunm * targetMatrix_row]) : new double[D_colunm *  D_row]);
    _ST = (ST != nullptr ? copyMatrix(ST,targetMatrix_colunm,targetMatrix_row,new double[targetMatrix_colunm * targetMatrix_row]) : new double[D_colunm *  D_row]);

    for(int i = 0; i < _maxIter; i++)
    {
        if(_ST != nullptr)
        {
            //std::cout << "ST" << std::endl;
            //showMatrix(_ST,targetMatrix_colunm,targetMatrix_row);
            transpose(D,D_colunm,D_row,D_T);
            transpose(_ST,targetMatrix_colunm,targetMatrix_row,Matrix_T);
            
            //showMatrix(D_T,D_row,D_colunm);
            //showMatrix(Matrix_T,targetMatrix_colunm,targetMatrix_row);
            regr->fit(Matrix_T,targetMatrix_row,targetMatrix_colunm,D_T,D_row,D_colunm);

            //showMatrix(Matrix_T,targetMatrix_row,targetMatrix_colunm);
            //showMatrix(D_T,D_row,D_colunm);

            double* C_tmp = regr->getCoef();
            //showMatrix(C_tmp,D_row,D_colunm);

            constraints(C_tmp,D_row,D_colunm);
            //showMatrix(C_tmp,D_row,D_colunm);


            cal_D = dot(C_tmp,D_row,D_colunm,_ST,targetMatrix_colunm,targetMatrix_row,cal_D);
            double tmp_err = meanSquareError(cal_D,D,D_colunm,D_row);

            if(ismin_err(tmp_err))n_above_min = 0;
            else n_above_min += 1;

            if(n_above_min > tol_n_above_min)break;

            if(errs.size() == 0)
            {
                errs.push_back(tmp_err);
                _C = copyMatrix(C_tmp,D_row,D_colunm,_C);
            }
            else if (tmp_err <= errs.end()[-1] * (1 + tol_increase))
            {
                errs.push_back(tmp_err);
                _C = copyMatrix(C_tmp,D_row,D_colunm,_C);;
            }
            else break;

            if(errs.size() > 1)
            {
                if(errs.end()[-1] > errs.end()[-2])n_increase += 1;
                else n_increase *= 0;
            }

            if(n_increase > tol_n_increase)break;
        }

        if (_C != nullptr)
        {
            //std::cout << "C" << std::endl;
            //showMatrix(_C,D_row,D_colunm);
            //showMatrix(_C,D_row,D_colunm);
            //showMatrix(D,D_colunm,D_row);
            regr->fit(_C,D_row,D_colunm,D,D_colunm,D_row);
            //showMatrix(_C,4,4);
            double* ST_tmp = regr->getCoef(false);

            constraints(ST_tmp,D_row,D_row);
            //showMatrix(ST_tmp,4,4);

            cal_D = dot(_C,D_row,D_colunm,ST_tmp,targetMatrix_colunm,targetMatrix_row,cal_D);
            double tmp_err = meanSquareError(cal_D,D,D_colunm,D_row);

            if(ismin_err(tmp_err))n_above_min = 0;
            else n_above_min += 1;

            if(n_above_min > tol_n_above_min)break;

            //std::cout << tmp_err << std::endl << errs.end()[-1] << std::endl << (1 + tol_increase) << std::endl;
            if(errs.size() == 0)
            {
                errs.push_back(tmp_err);
                _ST = copyMatrix(ST_tmp,targetMatrix_colunm,targetMatrix_row,_ST);
            }
            else if (tmp_err <= errs.end()[-1] * (1 + tol_increase))
            {
                errs.push_back(tmp_err);
                _ST = copyMatrix(ST_tmp,targetMatrix_colunm,targetMatrix_row,_ST);
            }
            else break;

            if(errs.size() > 1)
            {
                if(errs.end()[-1] > errs.end()[-2])n_increase += 1;
                else n_increase *= 0;
            }

            if(n_increase > tol_n_increase)break;
        }
    }
    delete[] D_T;
    delete[] Matrix_T;
    delete[] cal_D;

    return *this;
}

bool MCR_ALS::ismin_err(double val)
{
    if(errs.size() == 0)return true;
    else
    {
        bool tmp = false;
        for(double err : errs)tmp |= (val > err);
        return !tmp;
    }
}

double* MCR_ALS::get_C()
{
    return _C;
}


double* MCR_ALS::get_St()
{
    return _ST;
}