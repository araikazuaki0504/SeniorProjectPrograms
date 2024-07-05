#include "MCR-ALS.hpp"

MCR_ALS::MCR_ALS(int maxIter):_maxIter{maxIter}
{
    _C_regr = new OLS();
    _St_regr = new OLS();
}

MCR_ALS::~MCR_ALS()
{
    delete[] _C;
    delete[] _ST;
    delete _C_regr;
    delete _St_regr;
}

MCR_ALS& MCR_ALS::changeRegressorType_forALL(std::string C_RegressorType,std::string St_RegressorType)
{
    if(C_RegressorType == "OLS")
    {
        delete _C_regr;
        _C_regr = new OLS();
    }
    else if (C_RegressorType == "NNLS")
    {
        delete _C_regr;
        _C_regr = new NNLS();
    }

    if(St_RegressorType == "OLS")
    {
        delete _St_regr;
        _St_regr = new OLS();
    }
    else if (St_RegressorType == "NNLS")
    {
        delete _St_regr;
        _St_regr = new NNLS();
    }

    std::cout << "Change To C_Regressor -> " << C_RegressorType << " & St_Regressor -> " << St_RegressorType <<  std::endl;

    return *this;
}

MCR_ALS& MCR_ALS::changeRegressorType_forC(std::string RegressorType)
{
    if(RegressorType == "OLS")
    {
        delete _C_regr;
        _C_regr = new OLS();
    }
    else if (RegressorType == "NNLS")
    {
        delete _C_regr;
        _C_regr = new NNLS();
    }

    std::cout << "changeTo C_Regressor -> " << RegressorType <<  std::endl;

    return *this;
}

MCR_ALS& MCR_ALS::changeRegressorType_forSt(std::string RegressorType)
{
    if(RegressorType == "OLS")
    {
        delete _St_regr;
        _St_regr = new OLS();
    }
    else if (RegressorType == "NNLS")
    {
        delete _St_regr;
        _St_regr = new NNLS();
    }

    std::cout << "change To St_Regressor -> " << RegressorType <<  std::endl;

    return *this;
}

MCR_ALS& MCR_ALS::fit(double* D, int D_colunm, int D_row, double* C, double* ST, int targetMatrix_colunm, int targetMatrix_row)
{
    if(C == nullptr && ST == nullptr)return *this;
    if(_C != nullptr)delete[] _C;
    if(_ST != nullptr)delete[] _ST;

    _C_Colunm = targetMatrix_colunm;
    _C_Row = (C != nullptr ? targetMatrix_row : targetMatrix_colunm);;

    _ST_Colunm = (ST != nullptr ? targetMatrix_colunm : targetMatrix_row);
    _ST_Row = targetMatrix_row;

    double* D_T = new double[D_colunm * D_row];
    double* Matrix_T = new double[_ST_Colunm * _ST_Row];
    double* cal_D = new double[D_colunm * D_row];

    int n_increase = 0;
    int n_above_min = 0;
    errs.clear();

    _C  = (C != nullptr  ? copyMatrix(C,_C_Colunm,_C_Row,new double[_C_Colunm * _C_Row]) : new double[_C_Colunm *  _C_Row]);
    _ST = (ST != nullptr ? copyMatrix(ST,_ST_Colunm,_ST_Row,new double[_ST_Colunm * _ST_Row]) : new double[_ST_Colunm *  _ST_Row]);

    for(int i = 0; i < _maxIter; i++)
    {
        if(_ST != nullptr)
        {
            transpose(D,D_colunm,D_row,D_T);
            transpose(_ST,_ST_Colunm,_ST_Row,Matrix_T);

            _C_regr->fit(Matrix_T,_ST_Row,_ST_Colunm,D_T,D_row,D_colunm);

            double* C_tmp = _C_regr->getCoef();

            constraints(C_tmp,_C_Colunm,_C_Row);

            cal_D = product(C_tmp,_C_Colunm,_C_Row,_ST,_ST_Colunm,_ST_Row,cal_D);
            double tmp_err = meanSquareError(cal_D,D,D_colunm,D_row);

            if(ismin_err(tmp_err))n_above_min = 0;
            else n_above_min += 1;

            if(n_above_min > tol_n_above_min)break;

            if(errs.size() == 0)
            {
                errs.push_back(tmp_err);
                _C = copyMatrix(C_tmp,_C_Colunm,_C_Row,_C);
            }
            else if (tmp_err <= errs.end()[-1] * (1 + tol_increase))
            {
                errs.push_back(tmp_err);
                _C = copyMatrix(C_tmp,_C_Colunm,_C_Row,_C);;
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
            _St_regr->fit(_C,_C_Colunm,_C_Row,D,D_colunm,D_row);
            double* ST_tmp = _St_regr->getCoef(false);

            constraints(ST_tmp,_ST_Colunm,_ST_Row);

            cal_D = product(_C,_C_Colunm,_C_Row,ST_tmp,_ST_Colunm,_ST_Row,cal_D);
            double tmp_err = meanSquareError(cal_D,D,D_colunm,D_row);

            if(ismin_err(tmp_err))n_above_min = 0;
            else n_above_min += 1;

            if(n_above_min > tol_n_above_min)break;

            if(errs.size() == 0)
            {
                errs.push_back(tmp_err);
                _ST = copyMatrix(ST_tmp,_ST_Colunm,_ST_Row,_ST);
            }
            else if (tmp_err <= errs.end()[-1] * (1 + tol_increase))
            {
                errs.push_back(tmp_err);
                _ST = copyMatrix(ST_tmp,_ST_Colunm,_ST_Row,_ST);
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