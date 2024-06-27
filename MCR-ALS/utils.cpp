#include "utils.hpp"

lstsq_result::lstsq_result(int solution_colunm_size, int solution_row_size):_solution_colunm_size{solution_colunm_size},_solution_row_size{solution_row_size}
{
    _lstsq_solution = new double[solution_colunm_size * solution_row_size];
    _residues = new double[solution_row_size];
    _single_value   = new double[solution_colunm_size];
};

lstsq_result::~lstsq_result()
{
    delete[] _lstsq_solution;
    delete[] _residues;
    delete[] _single_value;
};

void lstsq_result::showData()
{
    std::cout << "最小二乗解：" << std::endl << "{" << std::endl;
    for(int i = 0; i < _solution_colunm_size; i++)
    {
        for(int j = 0; j < _solution_row_size; j++)
        {
            cout << std::fixed << std::setprecision(5) << _lstsq_solution[i * _solution_row_size + j] << " ";
        }
        std::cout << std::endl;
    }
    cout << "}" << std::endl;
    std::cout << "残差：{";
    for(int i = 0; i < _solution_colunm_size; i++)std::cout << std::fixed << std::setprecision(16) <<_residues[i] << ",";
    std::cout << "}" << std::endl;
    if(_rank != -1)
    {
        std::cout << "rank：" << _rank << std::endl << "特異値：{";
        for(int i = 0; i < (_solution_colunm_size > _solution_row_size ? _solution_row_size : _solution_colunm_size);i++)std::cout << _single_value[i] << ",";
        std::cout << "}" << endl;
    }
}

bool check_P(bool* P, int P_size)
{
    bool all_P_bool = true;

    for(int i = 0; i < P_size; i++)all_P_bool &= all_P_bool;

    return all_P_bool;
}

bool check_w(double* w, bool* P, int P_size, double tol)
{
    bool isConv = false;
    for(int i = 0; i < P_size; i++)
    {
        if(P[i])continue;
        isConv |= (w[i] > tol);
    }
    return isConv;
}

bool check_s(double* s, bool* P, int P_size)
{
    bool isConv = false;
    double minValue = numeric_limits<double>::max();
    for(int i = 0; i < P_size; i++)
    {
        if(!P[i])continue;
        minValue = min(minValue,s[i]);
    }

    return minValue < 0;
}

int maxValueIndex(double* w, bool* P, int P_size)
{
    double maxValue = 0;
    int maxValueIndex = 0;
    for(int i = 0; i < P_size; i++)
    {
        if(P[i])continue;
        if(w[i] > maxValue)
        {
            maxValueIndex = i;
            maxValue = w[i];
        }
    }

    return maxValueIndex;
}

double get_miniumAlphaWithInds(double* x, double* s, bool* inds, int inds_size)
{
    double minValue = numeric_limits<double>::max();
    for(int i = 0; i < inds_size; i++)
    {
        if(!inds[i])continue;

        double alpha = x[i] / (x[i] - s[i]);
        minValue = min(minValue,alpha);
    }

    return minValue;
}

int lstsq(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_ ,lstsq_result* cal_result)
{
    std::cout << "start Least square..." << std::endl;

    double* Pinv_Matrix = new double[Matrix_A_colunm_ * Matrix_A_row_];
    double* vector_B_from_Matrix_B = new double[Matrix_B_colunm_];
    double* solution_buffer = new double[Matrix_A_row_];
    double* residues_buffer = new double[Matrix_B_colunm_];

    //疑似逆行列を求める&特異値を求める
    Pseudo_inverse(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Pinv_Matrix);
    //showMatrix(Vt_Sstar_Ut,Matrix_A_row_,Matrix_A_colunm_);

    for(int i = 0; i < Matrix_B_row_; i++)
    {
        for(int j = 0; j< Matrix_B_colunm_; j++)
        {
            vector_B_from_Matrix_B[j] = Matrix_B[j * Matrix_B_row_ + i];
        }

        //最小二乗解を求める
        
        //SETDEBUGFLAG(true);
        dot(Pinv_Matrix,Matrix_A_row_,Matrix_A_colunm_,vector_B_from_Matrix_B,Matrix_B_colunm_,1,solution_buffer);
        //SETDEBUGFLAG(false);
        for(int j = 0; j < Matrix_A_colunm_; j++)
        {
            cal_result->_lstsq_solution[j * cal_result->_solution_row_size + i] = solution_buffer[j];
        }

        //残差(二乗誤差)
        double residues = 0;
        //SETDEBUGFLAG(true);
        dot(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,solution_buffer,Matrix_A_row_,1,residues_buffer);
        //SETDEBUGFLAG(false);
        sub(vector_B_from_Matrix_B,residues_buffer,Matrix_B_colunm_,1,residues_buffer);
        residues = norm(residues_buffer,Matrix_B_colunm_);
        cal_result->_residues[i] = sqrt(residues);
    }

    //rankを求める
    //cal_result->_rank = eig_rank(cal_result->_single_value,(Matrix_A_colunm_ > Matrix_A_row_ ? Matrix_A_row_ : Matrix_A_colunm_));

    delete[] Pinv_Matrix;
    delete[] vector_B_from_Matrix_B;
    delete[] solution_buffer;
    delete[] residues_buffer;


    //double* Matrix_A_t = new double[Matrix_A_colunm_ * Matrix_A_row_];
    //double* MtM = new double[Matrix_A_row_ * Matrix_A_row_];
    //double* buffer = new double[Matrix_A_colunm_ * Matrix_A_row_];
    //double* vector = new double[Matrix_A_colunm_];
    //double* vector_B_from_Matrix_B = new double[Matrix_B_colunm_];
//
//
    //
//
    //for(int i = 0; i < cal_result->_solution_row_size; i++)
    //{
    //    for(int j = 0; j< Matrix_B_colunm_; j++)
    //    {
    //        vector_B_from_Matrix_B[j] = Matrix_B[j * Matrix_B_row_ + i];
    //    }
//
    //    double* current_resultPointer = cal_result->_lstsq_solution + i * cal_result->_solution_colunm_size;
//
    //     //最小二乗解を求める
    //    transpose(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_A_t);
    //    showMatrix(Matrix_A_t,Matrix_A_row_,Matrix_A_colunm_);
    //    dot(Matrix_A_t,Matrix_A_row_,Matrix_A_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,MtM);
    //    showMatrix(MtM,Matrix_A_row_,Matrix_A_row_);
    //    inverse(MtM,Matrix_A_row_,Matrix_A_row_,MtM);
    //    showMatrix(MtM,Matrix_A_row_,Matrix_A_row_);
    //    dot(MtM,Matrix_A_row_,Matrix_A_row_,Matrix_A_t,Matrix_A_row_,Matrix_A_colunm_,buffer);
    //    showMatrix(buffer,Matrix_A_colunm_,Matrix_A_row_);
    //    showMatrix(vector_B_from_Matrix_B,1,Matrix_B_colunm_);
    //    dot(buffer,Matrix_A_row_,Matrix_A_colunm_,vector_B_from_Matrix_B,Matrix_B_colunm_,1,current_resultPointer);
//
    //    showMatrix(current_resultPointer,1,cal_result->_solution_colunm_size);
//
    //    //残差(二乗誤差)
    //    int residues = 0;
//
    //    dot(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,current_resultPointer,cal_result->_solution_colunm_size,1,vector);
    //    sub(vector,Matrix_B,Matrix_B_colunm_,1,vector);
//
    //    residues = norm(vector,Matrix_A_row_);
//
    //    cal_result->_residues[i] = sqrt(residues);
    //}
//
    //transpose(cal_result->_lstsq_solution,cal_result->_solution_row_size,cal_result->_solution_colunm_size,cal_result->_lstsq_solution);
//
    //delete[] Matrix_A_t;
    //delete[] MtM;
    //delete[] buffer;
    //delete[] vector;
    //delete[] vector_B_from_Matrix_B;

    //特異値分解
    //int Matrix_buffer_size = (Matrix_A_colunm_ > Matrix_A_row_ ? Matrix_A_row_ : Matrix_A_colunm_);
    //double* Matrix_t = new double[Matrix_A_colunm_ * Matrix_A_row_];
    //double* Matrix_buffer = new double[Matrix_A_colunm_ * Matrix_A_colunm_];
    //double* eigVal = new double[Matrix_buffer_size];
//
    //Matrix_t = transpose(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_t);
//
    //if(Matrix_A_colunm_ > Matrix_A_row_)//U
    //{
    //    Matrix_buffer = dot(Matrix_t,Matrix_A_row_,Matrix_A_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_buffer);
    //}
    //else//V
    //{
    //    Matrix_buffer = dot(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_t,Matrix_A_row_,Matrix_A_colunm_,Matrix_buffer);
    //}
//
    //eigVal = eig(Matrix_buffer,Matrix_A_colunm_,Matrix_A_colunm_,eigVal);

    //Matrix_Aのrank計算
    //cal_result->_rank = Matrix_rank(eigVal,Matrix_A_colunm_,Matrix_buffer_size);

    //singleValueの計算
    //int j = 0;
    //for(int i = 0; i < cal_result->_solution_colunm_size;i++)
    //{
    //    if(eigVal[i * Matrix_A_row_ + i] == 0)continue;
    //    cal_result->_single_value[j] = std::sqrt(eigVal[i * Matrix_A_row_ + i]);
    //    j++;
    //}

    //delete[] Matrix_t;
    //delete[] Matrix_buffer;
    //delete[] eigVal;
    std::cout << "end Least square..." << std::endl;
    return 1;
}

int nnls(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_,lstsq_result* cal_result)
{
    std::cout << "start None Negative Least square..." << std::endl;
    double* Matrix_A_T = new double[Matrix_A_colunm_ * Matrix_A_row_];
    double* AtA = new double[Matrix_A_row_ * Matrix_A_row_];
    double* x = new double[Matrix_A_row_];
    double* s = new double[Matrix_A_row_];
    bool* P = new bool[Matrix_A_row_];
    double* w = new double[Matrix_A_row_];
    double* vector_B_from_MatrixB = new double[Matrix_B_colunm_];

    for(int i = 0; i < Matrix_B_row_; i++)
    {
        for(int j = 0; j < Matrix_B_colunm_; j++)
        {
            vector_B_from_MatrixB[j] = Matrix_B[j * Matrix_B_row_ + i];
        }

        double* AtB = cal_result->_lstsq_solution + i * cal_result->_solution_colunm_size;
        transpose(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_A_T);//A^T
        dot(Matrix_A_T,Matrix_A_row_,Matrix_A_row_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,AtA);//A^T * A
        dot(vector_B_from_MatrixB,1,Matrix_B_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,AtB);//Atb

        int maxiter = 3 * Matrix_A_row_;

        double tol = 10 * max(Matrix_A_colunm_,Matrix_A_row_) * numeric_limits<double>::epsilon();

        for(int i = 0; i < Matrix_A_row_; i++)
        {
            x[i] = 0;
            s[i] = 0;
            P[i] = false;
            w[i] = AtB[i];
        }

        int iter = 0;
        int number_of_true = 0;

        while(check_P(P,Matrix_A_row_) && check_w(w,P,Matrix_A_row_,tol))
        {
            int k = maxValueIndex(w,P,Matrix_A_row_);
            P[k] = true;
            number_of_true += 1;

            for(int i = 0; i < Matrix_A_row_; i++)s[i] = 0;

            double* new_Matrix = new double[number_of_true * number_of_true];
            double* new_vector = new double[number_of_true];
            double* tmp        = new double[number_of_true];

            vector<int> indexList;
            for(int i = 0; i < Matrix_A_row_; i++)
            {
                if(!P[i])continue;
                indexList.push_back(i);
            }

            int n = 0;
            for(int Colunm_index : indexList)
            {
                int m = 0;
                for(int Row_index : indexList)
                {
                    new_Matrix[n * number_of_true + m] = AtA[Colunm_index * Matrix_A_row_ + Row_index];
                    m++;
                }
                new_vector[n] = AtB[Colunm_index];
                n++;
            }

            indexList.clear();

            equal_solve(new_Matrix,number_of_true,number_of_true,new_vector,number_of_true,tmp);

            for(int i = 0, j = 0; i < Matrix_A_row_; i++)
            {
                if(P[i])
                {
                    s[i] = tmp[j];
                    j++;
                }
            }

            while(iter < maxiter && check_s(s,P,Matrix_A_row_))
            {
                //showMatrix(P,1,Matrix_A_row_);
                iter += 1;

                bool* inds = new bool[Matrix_A_row_];

                for(int i = 0; i < Matrix_A_row_; i++)inds[i] = P[i] && (s[i] < 0);

                double alpha = get_miniumAlphaWithInds(x,s,inds,Matrix_A_row_);

                //showMatrix(x,1,Matrix_A_row_);
                //showMatrix(s,1,Matrix_A_row_);
                //showMatrix(inds,1,Matrix_A_row_);
                //std::cout << alpha << std::endl;
                //showMatrix(P,1,Matrix_A_row_);
                //showMatrix(s,1,Matrix_A_row_);
                //std::cout << tol << std::endl;

                for(int i = 0; i < Matrix_A_row_; i++)
                {
                    x[i] *= 1 - alpha;
                    x[i] += alpha * s[i];
                }

                //showMatrix(x,1,Matrix_A_row_);
                //showMatrix(P,1,Matrix_A_row_);

                for(int i = 0; i < Matrix_A_row_; i++)
                {
                    if(x[i] <= tol && P[i] == true)
                    {
                        P[i] = false;
                        number_of_true -= 1;
                    }
                }
                //std::cout << tol << std::endl;
                //showMatrix(x,1,Matrix_A_row_);
                //showMatrix(P,1,Matrix_A_row_);
                //std::cout << number_of_true << std::endl;

                double* new_sub_Matrix = new double[number_of_true * number_of_true];
                double* new_sub_vector = new double[number_of_true];
                double* tmp_sub        = new double[number_of_true];


                vector<int> sub_indexList;

                for(int i = 0; i < Matrix_A_row_; i++)
                {
                    if(!P[i])continue;
                    sub_indexList.push_back(i);
                }

                int n = 0;
                for(int Colunm_index : sub_indexList)
                {
                    int m = 0;
                    for(int Row_index : sub_indexList)
                    {
                        new_sub_Matrix[n * number_of_true + m] = AtA[Colunm_index * Matrix_A_row_ + Row_index];
                        m++;
                    }
                    new_sub_vector[n] = AtB[Colunm_index];
                    n++;
                }

                sub_indexList.clear();

                //showMatrix(new_sub_Matrix,number_of_true,number_of_true);
                //showMatrix(new_sub_vector,1,number_of_true);

                equal_solve(new_sub_Matrix,number_of_true,number_of_true,new_sub_vector,number_of_true,tmp_sub);

                for(int i = 0, j = 0; i < Matrix_A_row_; i++)
                {
                    if(P[i])
                    {
                        s[i] = tmp_sub[j];
                        j++;
                    }
                    else
                    {
                        s[i] = 0;
                    }
                    //x[i] = s[i];       
                }

                delete[] new_sub_Matrix;
                delete[] new_sub_vector;
                delete[] tmp_sub;
                delete[] inds;
            }

            for(int i = 0; i < Matrix_A_row_; i++)x[i] = s[i];

            double* tmp_vector = new double[Matrix_B_colunm_];
            dot(AtA,Matrix_A_colunm_,Matrix_A_row_,x,Matrix_A_row_,1,tmp_vector);
            sub(AtB,tmp_vector,Matrix_B_colunm_,1,w);

            //showMatrix(x,1,Matrix_A_row_);
            //showMatrix(w,1,Matrix_B_colunm_);

            delete[] tmp_vector;

            delete[] new_Matrix;
            delete[] new_vector;
            delete[] tmp;

            if(iter == maxiter)
            {
                for(int i = 0; i < cal_result->_solution_colunm_size; i++)AtB[i] = x[i];
                cal_result->_residues = nullptr;
                delete[] Matrix_A_T;
                delete[] AtA;
                delete[] x;
                delete[] s;
                delete[] P;
                delete[] w;
                delete[] vector_B_from_MatrixB;
                return -1;
            }
        }

        //最小二乗解
        for(int i = 0; i < cal_result->_solution_colunm_size; i++)AtB[i] = x[i];
        double* tmp_vector = new double[Matrix_B_colunm_];
        dot(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,x,Matrix_A_row_,1,tmp_vector);
        sub(tmp_vector,vector_B_from_MatrixB,1,Matrix_B_colunm_,tmp_vector);
        //残差計算
        cal_result->_residues[i] = norm(tmp_vector,Matrix_B_colunm_);

        delete[] tmp_vector;

        //showMatrix(AtB,1,cal_result->_solution_colunm_size);

    }
    
    transpose(cal_result->_lstsq_solution,cal_result->_solution_row_size,cal_result->_solution_colunm_size,cal_result->_lstsq_solution);

    delete[] Matrix_A_T;
    delete[] AtA;
    delete[] x;
    delete[] s;
    delete[] P;
    delete[] w;
    delete[] vector_B_from_MatrixB;

    std::cout << "end None Negative Least square..." << std::endl;
    return 1;
}

double meanSquareError(double* acutuallyMatrix, double* calMatrix, int Matrix_colunm, int Matrix_row)
{
    double sum = 0;
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)sum += (acutuallyMatrix[i * Matrix_row + j] - calMatrix[i * Matrix_row + j]) * (acutuallyMatrix[i * Matrix_row + j] - calMatrix[i * Matrix_row + j]);
    return sum / (Matrix_colunm * Matrix_row);
}

void constraints(double* Matrix, int Matrix_colunm, int Matrix_row)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)if(Matrix[i * Matrix_row + j] < 0)Matrix[i * Matrix_row + j] = 0;
}

double* copyMatrix(double* Matrix, int Matrix_colunm, int Matrix_row, double* copiedMatrix)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)copiedMatrix[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
    return copiedMatrix;
}