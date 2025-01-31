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
    std::cout << "最小二乗解：" << std::endl;
    showMatrix(_lstsq_solution, _solution_colunm_size, _solution_row_size);
    std::cout << std::endl << "残差：{";
    for(int i = 0; i < _solution_row_size; i++)std::cout << std::fixed << std::setprecision(16) <<_residues[i] << ",";
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
    double maxValue = std::numeric_limits<double>::min();
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
    double* Pinv_Matrix = new double[Matrix_A_colunm_ * Matrix_A_row_];
    double* vector_B_from_Matrix_B = new double[Matrix_B_colunm_];
    double* residues_buffer = new double[Matrix_A_colunm_ * Matrix_B_row_];
    double* Buffer = new double[Matrix_B_colunm_];

    //疑似逆行列を求める
    Pseudo_inverse(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Pinv_Matrix);

    product(Pinv_Matrix,Matrix_A_row_,Matrix_A_colunm_,Matrix_B,Matrix_B_colunm_,Matrix_B_row_,cal_result->_lstsq_solution);
    
    product(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,cal_result->_lstsq_solution,Matrix_A_row_,Matrix_B_row_,residues_buffer);
    sub(residues_buffer,Matrix_B,Matrix_B_colunm_,Matrix_B_row_,residues_buffer);

    for(int i = 0; i < Matrix_B_row_; i++)
    {
        for(int j = 0; j < Matrix_B_colunm_; j++)Buffer[j] = residues_buffer[j * Matrix_B_row_ + i];

        double residues = 0;
        residues = norm(Buffer,Matrix_B_colunm_);
        cal_result->_residues[i] = sqrt(residues);
    }

    delete[] Pinv_Matrix;
    delete[] vector_B_from_Matrix_B;
    delete[] residues_buffer;
    delete[] Buffer;

    return 1;
}

int nnls(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_,lstsq_result* cal_result)
{
    double* Matrix_A_T = new double[Matrix_A_row_ * Matrix_A_colunm_];
    double* AtA = new double[Matrix_A_row_ * Matrix_A_row_];
    double* AtB = new double[Matrix_A_row_];
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

        transpose(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_A_T);//A^T
        product(Matrix_A_T,Matrix_A_row_,Matrix_A_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,AtA);//A^T * A
        product(vector_B_from_MatrixB,1,Matrix_B_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,AtB);//Atb

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
                iter += 1;

                bool* inds = new bool[Matrix_A_row_];

                for(int i = 0; i < Matrix_A_row_; i++)inds[i] = P[i] && (s[i] < 0);

                double alpha = get_miniumAlphaWithInds(x,s,inds,Matrix_A_row_);

                for(int i = 0; i < Matrix_A_row_; i++)
                {
                    x[i] *= 1 - alpha;
                    x[i] += alpha * s[i];
                }

                for(int i = 0; i < Matrix_A_row_; i++)
                {
                    if(x[i] <= tol && P[i] == true)
                    {
                        P[i] = false;
                        number_of_true -= 1;
                    }
                }

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
                }

                if(new_sub_Matrix != nullptr)delete[] new_sub_Matrix;
                if(new_sub_vector != nullptr)delete[] new_sub_vector;
                if(tmp_sub != nullptr)delete[] tmp_sub;
                if(inds != nullptr)delete[] inds;
            }


            for(int i = 0; i < Matrix_A_row_; i++)x[i] = s[i];

            double* tmp_vector = new double[Matrix_A_row_];
            product(AtA, Matrix_A_row_,Matrix_A_row_,x,Matrix_A_row_,1,tmp_vector);
            sub(AtB,tmp_vector, Matrix_A_row_,1,w);

            if(tmp_vector != nullptr)delete[] tmp_vector;

            if(new_Matrix != nullptr)delete[] new_Matrix;
            if(new_vector != nullptr)delete[] new_vector;
            if(tmp != nullptr)delete[] tmp;

            if(iter == maxiter)
            {
                for(int j = 0; j < cal_result->_solution_colunm_size; j++)cal_result->_lstsq_solution[j * cal_result->_solution_row_size + i] = x[j];
                cal_result->_residues = nullptr;
                if (Matrix_A_T != nullptr)delete[] Matrix_A_T;
                if (AtA != nullptr)delete[] AtA;
                if (AtB != nullptr)delete[] AtB;
                if (x != nullptr)delete[] x;
                if (s != nullptr)delete[] s;
                if (P != nullptr)delete[] P;
                if (w != nullptr)delete[] w;
                if (vector_B_from_MatrixB != nullptr)delete[] vector_B_from_MatrixB;
                return -1;
            }
        }

        //最小二乗解
        for(int j = 0; j < cal_result->_solution_colunm_size; j++)cal_result->_lstsq_solution[j * cal_result->_solution_row_size + i] = x[j];
        double* tmp_vector = new double[Matrix_A_row_];
        product(AtA,Matrix_A_row_,Matrix_A_row_,x,Matrix_A_row_,1,tmp_vector);
        sub(tmp_vector,vector_B_from_MatrixB,1, Matrix_A_row_,tmp_vector);
        //残差計算
        cal_result->_residues[i] = norm(tmp_vector, Matrix_A_row_);

        if(tmp_vector != nullptr)delete[] tmp_vector;
    }

    if (Matrix_A_T != nullptr)delete[] Matrix_A_T;
    if (AtA != nullptr)delete[] AtA;
    if (AtB != nullptr)delete[] AtB;
    if (x != nullptr)delete[] x;
    if (s != nullptr)delete[] s;
    if (P != nullptr)delete[] P;
    if (w != nullptr)delete[] w;
    if (vector_B_from_MatrixB != nullptr)delete[] vector_B_from_MatrixB;

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