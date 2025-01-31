#include "MCR-ALS.hpp"

MCR_ALS::MCR_ALS(int maxIter):_maxIter{maxIter}
{

}

MCR_ALS::~MCR_ALS()
{
    delete[] _C;
    delete[] _ST;
    delete _C_regr;
    delete _St_regr;
}

MCR_ALS& MCR_ALS::changeRegressorType_forALL(const char* C_RegressorType, const char* St_RegressorType)
{

    if(std::strcmp(C_RegressorType, "OLS") == 0)
    {
        if(_C_regr != nullptr)delete _C_regr;
        _C_regr = new OLS();
    }
    else if (std::strcmp(C_RegressorType, "NNLS") == 0)
    {
        if (_C_regr != nullptr)delete _C_regr;
        _C_regr = new NNLS();
    }

    if(std::strcmp(St_RegressorType, "OLS") == 0)
    {
        if (_St_regr != nullptr)delete _St_regr;
        _St_regr = new OLS();
    }
    else if (std::strcmp(St_RegressorType, "NNLS") == 0)
    {
        if (_St_regr != nullptr)delete _St_regr;
        _St_regr = new NNLS();
    }

    return *this;
}

MCR_ALS& MCR_ALS::changeRegressorType_forC(const char* RegressorType)
{
    if (std::strcmp(RegressorType, "OLS") == 0)
    {
        delete _C_regr;
        _C_regr = new OLS();
    }
    else if (std::strcmp(RegressorType, "NNLS") == 0)
    {
        delete _C_regr;
        _C_regr = new NNLS();
    }

    return *this;
}

MCR_ALS& MCR_ALS::changeRegressorType_forSt(const char* RegressorType)
{
    if (std::strcmp(RegressorType, "OLS") == 0)
    {
        delete _St_regr;
        _St_regr = new OLS();
    }
    else if (std::strcmp(RegressorType, "NNLS") == 0)
    {
        delete _St_regr;
        _St_regr = new NNLS();
    }

    return *this;
}

MCR_ALS& MCR_ALS::fit(double* D, int D_colunm, int D_row, int n_component, double* C, double* ST)
{
    if(C == nullptr && ST == nullptr)return *this;
    if(_C != nullptr)delete[] _C;
    if(_ST != nullptr)delete[] _ST;

    _C_Colunm = D_colunm;
    _C_Row = n_component;

    _ST_Colunm = n_component;
    _ST_Row = D_row;

    double* D_T = new double[D_colunm * D_row];
    double* Matrix_T = new double[_ST_Colunm * _ST_Row];
    double* cal_D = new double[D_colunm * D_row];

    int n_increase = 0;
    int n_above_min = 0;
    int n_iter = 0;
    double tol_err_change = 1e-10;

    bool loop_flag = true;

    errs.clear();

    _C  = (C != nullptr  ? copyMatrix(C,_C_Colunm,_C_Row,new double[_C_Colunm * _C_Row]) : new double[_C_Colunm *  _C_Row]);
    _ST = (ST != nullptr ? copyMatrix(ST,_ST_Colunm,_ST_Row,new double[_ST_Colunm * _ST_Row]) : new double[_ST_Colunm *  _ST_Row]);

    for(int i = 0; i < _maxIter; i++)
    {
        n_iter = i + 1;
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

            if(n_above_min > tol_n_above_min)
            {
                std::cout << "Half-iterated" << std::endl;
                break;
            }

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
            else
            {
                std::cout << "Error increased above fractional tol_increase (C iter). Exiting 2" << std::endl;
                break;
            }

            if(errs.size() > 1)
            {
                if(errs.end()[-1] > errs.end()[-2])n_increase += 1;
                else n_increase *= 0;
            }

            if(n_increase > tol_n_increase)
            {
                std::cout << "Maximum error increases reached (C iter)" << std::endl;
                break;
            }
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

            if(n_above_min > tol_n_above_min)
            {
                std::cout << "Half-iterated" << std::endl;
                break;
            }

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
            else if (loop_flag)
            {
                loop_flag = true;
            }
            else
            {
                std::cout << "Error increased above fractional tol_increase (St iter). Exiting" << std::endl;
                break;
            }

            if(errs.size() > 1)
            {
                if(errs.end()[-1] > errs.end()[-2])n_increase += 1;
                else n_increase *= 0;
            }

            if(n_increase > tol_n_increase)
            {
                std::cout << "Maximum error increases reached (St iter)" << std::endl;
                break;
            }
        }

        if (n_iter >= _maxIter)
        {
            std::cout << "Max iterations reached" << std::endl;
            break;
        }

        n_iter = i + 1;

        if (errs.size() > 2)
        {
            double error_diff = abs(errs.end()[-1] - errs.end()[-3]);
            if (error_diff < abs(tol_err_change))
            {
                std::cout << "Change in err below tol_err_change(" << std::scientific << std::setprecision(5) << error_diff << ")" << std::endl;
                std::cout << "Residue = " << std::scientific << std::setprecision(5) << errs.end()[-1] << std::endl;
                break;
            }
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