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
    //std::cout << "start Least square..." << std::endl;
    //showMatrix(Matrix_A,Matrix_A_colunm_,Matrix_A_row_);
    //showMatrix(Matrix_B,Matrix_B_colunm_,Matrix_B_row_);

    double* Pinv_Matrix = new double[Matrix_A_colunm_ * Matrix_A_row_];
    double* vector_B_from_Matrix_B = new double[Matrix_B_colunm_];
    double* residues_buffer = new double[Matrix_A_colunm_ * Matrix_B_row_];
    double* Buffer = new double[Matrix_B_colunm_];

    //疑似逆行列を求める
    Pseudo_inverse(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Pinv_Matrix);
    //showMatrix(Pinv_Matrix,Matrix_A_row_,Matrix_A_colunm_);
    //showMatrix(Matrix_B, Matrix_B_colunm_, Matrix_B_row_);

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

    //for(int i = 0; i < Matrix_B_row_; i++)
    //{
    //    for(int j = 0; j< Matrix_B_colunm_; j++)
    //    {
    //        vector_B_from_Matrix_B[j] = Matrix_B[j * Matrix_B_row_ + i];
    //    }
    //    //showMatrix(Pinv_Matrix,Matrix_A_row_,Matrix_A_colunm_);
    //    //std::cout << "x" << std::endl;
    //    //showMatrix(vector_B_from_Matrix_B,1,Matrix_B_colunm_);
//
    //    //最小二乗解を求める
    //    //SETDEBUGFLAG(true);
    //    //product(Pinv_Matrix,Matrix_A_row_,Matrix_A_colunm_,vector_B_from_Matrix_B,Matrix_B_colunm_,1,solution_buffer);
    //    //std::cout << "=" << std::endl;
    //    //showMatrix(solution_buffer,1,Matrix_A_row_);
    //    //SETDEBUGFLAG(false);
    //    //for(int j = 0; j < Matrix_A_colunm_; j++)
    //    //{
    //    //    cal_result->_lstsq_solution[j * cal_result->_solution_row_size + i] = solution_buffer[j];
    //    //}
//
    //    //残差(二乗誤差)
    //    double residues = 0;
    //    //SETDEBUGFLAG(true);
    //    product(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,solution_buffer,Matrix_A_row_,1,residues_buffer);
    //    //SETDEBUGFLAG(false);
    //    sub(vector_B_from_Matrix_B,residues_buffer,Matrix_B_colunm_,1,residues_buffer);
    //    residues = norm(residues_buffer,Matrix_B_colunm_);
    //    cal_result->_residues[i] = sqrt(residues);
    //}

    //rankを求める
    //cal_result->_rank = eig_rank(cal_result->_single_value,(Matrix_A_colunm_ > Matrix_A_row_ ? Matrix_A_row_ : Matrix_A_colunm_));

    delete[] Pinv_Matrix;
    delete[] vector_B_from_Matrix_B;
    //delete[] solution_buffer;
    delete[] residues_buffer;
    delete[] Buffer;


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
    //    product(Matrix_A_t,Matrix_A_row_,Matrix_A_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,MtM);
    //    showMatrix(MtM,Matrix_A_row_,Matrix_A_row_);
    //    inverse(MtM,Matrix_A_row_,Matrix_A_row_,MtM);
    //    showMatrix(MtM,Matrix_A_row_,Matrix_A_row_);
    //    product(MtM,Matrix_A_row_,Matrix_A_row_,Matrix_A_t,Matrix_A_row_,Matrix_A_colunm_,buffer);
    //    showMatrix(buffer,Matrix_A_colunm_,Matrix_A_row_);
    //    showMatrix(vector_B_from_Matrix_B,1,Matrix_B_colunm_);
    //    product(buffer,Matrix_A_row_,Matrix_A_colunm_,vector_B_from_Matrix_B,Matrix_B_colunm_,1,current_resultPointer);
//
    //    showMatrix(current_resultPointer,1,cal_result->_solution_colunm_size);
//
    //    //残差(二乗誤差)
    //    int residues = 0;
//
    //    product(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,current_resultPointer,cal_result->_solution_colunm_size,1,vector);
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
    //    Matrix_buffer = product(Matrix_t,Matrix_A_row_,Matrix_A_colunm_,Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_buffer);
    //}
    //else//V
    //{
    //    Matrix_buffer = product(Matrix_A,Matrix_A_colunm_,Matrix_A_row_,Matrix_t,Matrix_A_row_,Matrix_A_colunm_,Matrix_buffer);
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
    //std::cout << "end Least square..." << std::endl;
    return 1;
}

//void _nnls(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* b_vector, int b_vector_size, double* Solution_vector)
//{
//    double* At = new double[Matrix_row_ * Matrix_colunm_];
//    double* AtA = new double[Matrix_row_, Matrix_row_];
//    double* Atb = new double[Matrix_row_];
//
//    transpose(Matrix, Matrix_colunm_, Matrix_row_, At);
//    product(At, Matrix_row_, Matrix_colunm_, Matrix, Matrix_colunm_, Matrix_row_, AtA);
//    product(At, Matrix_row_, Matrix_colunm_, b_vector, b_vector_size, 1, Atb);
//
//    int max_iter = 3 * Matrix_row_;
//    int tol = 10 * max(Matrix_colunm_, Matrix_row_) * numeric_limits<double>::epsilon();
//
//    double* x = new double[Matrix_row_];
//    double* s = new double[Matrix_row_];
//    double* w = new double[Matrix_row_];
//    bool* P = new bool[Matrix_row_];
//
//    for (int i = 0; i < Matrix_row_; i++)
//    {
//        x[i] = 0;
//        s[i] = 0;
//        P[i] = false;
//        w[i] = Atb[i];
//    }
//
//    int iter = 0;
//
//    while (!check_P(P,Matrix_row_) && check_w(w, P, Matrix_row_, tol))
//    {
//        int k = maxValueIndex(w, P, Matrix_row_);
//        P[k] = true;
//
//        for (int i = 0; i < Matrix_row_; i++)s[i] = 0;
//
//        vector<int> index_list;
//
//        for (int i = 0; i < Matrix_row_; i++)
//        {
//            if (!P[i])continue;
//            index_list.push_back(i);
//        }
//
//        int new_matrix_size = index_list.size();
//
//        double* new_matrix = new double[new_matrix_size * new_matrix_size];
//        double* new_vector = new double[new_matrix_size];
//        double* tmp = new double[new_matrix_size];
//
//        int n = 0;
//        for(int Colunm_index : index_list)
//        {
//            int m = 0;
//            for(int Row_index : index_list)
//            {
//                new_matrix[n * new_matrix_size + m] = AtA[Colunm_index * Matrix_row_ + Row_index];
//                m++;
//            }
//            new_vector[n] = Atb[Colunm_index];
//            n++;
//        }
//
//        equal_solve(new_matrix, new_matrix_size, new_matrix_size, new_vector, new_matrix_size, tmp);
//
//        n = 0;
//
//        for (int index : index_list)
//        {
//            s[index] = tmp[n];
//            n++;
//        }
//
//        while ((max_iter > iter) && check_s(s,P,Matrix_row_))
//        {
//            iter++;
//            bool* inds = new bool[Matrix_row_];
//            for (int i = 0; i < Matrix_row_; i++)inds[i] = P[i] && (s[i] < 0);
//            double alpha = get_miniumAlphaWithInds(x, s, inds, Matrix_row_);
//
//            for (int i = 0; i < Matrix_row_; i++)
//            {
//                x[i] *= 1 - alpha;
//                x[i] += alpha * s[i];
//            }
//        }
//    }
//
//}
//
//void nnls(double* Matrix_A, int Matrix_A_colunm_, int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_, lstsq_result* cal_result)
//{
//    
//}

int nnls(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_,lstsq_result* cal_result)
{
    //std::cout << "start None Negative Least square..." << std::endl;
    double* Matrix_A_T = new double[Matrix_A_row_ * Matrix_A_colunm_];
    double* AtA = new double[Matrix_A_row_ * Matrix_A_row_];
    double* AtB = new double[Matrix_A_row_];
    double* x = new double[Matrix_A_row_];
    double* s = new double[Matrix_A_row_];
    bool* P = new bool[Matrix_A_row_];
    double* w = new double[Matrix_A_row_];
    double* vector_B_from_MatrixB = new double[Matrix_B_colunm_];

    //showMatrix(Matrix_A,Matrix_A_colunm_,Matrix_A_row_);
    //showMatrix(vector_B_from_MatrixB,1,Matrix_B_colunm_);

    //showMatrix(Matrix_B,Matrix_B_colunm_,Matrix_B_row_);

    for(int i = 0; i < Matrix_B_row_; i++)
    {
        //std::cout << "i：" << i << std::endl;
        for(int j = 0; j < Matrix_B_colunm_; j++)
        {
            vector_B_from_MatrixB[j] = Matrix_B[j * Matrix_B_row_ + i];
        }

        //showMatrix(Matrix_A,Matrix_A_colunm_,Matrix_A_row_);
        //showMatrix(vector_B_from_MatrixB,1,Matrix_B_colunm_);

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
            //std::cout << "number of true：" << number_of_true << std::endl;
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
            //showMatrix(new_Matrix,number_of_true,number_of_true);
            //showMatrix(new_vector,1,number_of_true);

            indexList.clear();

            equal_solve(new_Matrix,number_of_true,number_of_true,new_vector,number_of_true,tmp);
            //showMatrix(tmp,1,number_of_true);
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
                //std::cout << "iter：" << iter << std::endl;
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

                if(new_sub_Matrix != nullptr)delete[] new_sub_Matrix;
                if(new_sub_vector != nullptr)delete[] new_sub_vector;
                if(tmp_sub != nullptr)delete[] tmp_sub;
                if(inds != nullptr)delete[] inds;
            }


            for(int i = 0; i < Matrix_A_row_; i++)x[i] = s[i];

            double* tmp_vector = new double[Matrix_A_row_];
            product(AtA, Matrix_A_row_,Matrix_A_row_,x,Matrix_A_row_,1,tmp_vector);
            sub(AtB,tmp_vector, Matrix_A_row_,1,w);

            //showMatrix(x,1,Matrix_A_row_,true);
            //showMatrix(w,1,Matrix_B_colunm_);

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

        //showMatrix(x,1,Matrix_A_row_,true);

        //最小二乗解
        for(int j = 0; j < cal_result->_solution_colunm_size; j++)cal_result->_lstsq_solution[j * cal_result->_solution_row_size + i] = x[j];
        double* tmp_vector = new double[Matrix_A_row_];
        product(AtA,Matrix_A_row_,Matrix_A_row_,x,Matrix_A_row_,1,tmp_vector);
        sub(tmp_vector,vector_B_from_MatrixB,1, Matrix_A_row_,tmp_vector);
        //残差計算
        cal_result->_residues[i] = norm(tmp_vector, Matrix_A_row_);

        if(tmp_vector != nullptr)delete[] tmp_vector;

        //showMatrix(AtB,1,cal_result->_solution_colunm_size);

    }

    if (Matrix_A_T != nullptr)delete[] Matrix_A_T;
    if (AtA != nullptr)delete[] AtA;
    if (AtB != nullptr)delete[] AtB;
    if (x != nullptr)delete[] x;
    if (s != nullptr)delete[] s;
    if (P != nullptr)delete[] P;
    if (w != nullptr)delete[] w;
    if (vector_B_from_MatrixB != nullptr)delete[] vector_B_from_MatrixB;

    //std::cout << "end None Negative Least square..." << std::endl;
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