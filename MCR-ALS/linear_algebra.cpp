#include "linear_algebra.hpp"

void showMatrix(double* Matrix, int colunm, int row)
{
    std::cout << "{" << std::endl;
    for(int i = 0; i < colunm; i++)
    {
        for(int j = 0; j < row; j++)
        {
            if(std::abs(Matrix[i * row + j]) < 1e-5)std::cout << 0 << ",";
            else std::cout << Matrix[i * row + j] << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}

void showMatrix(bool* Matrix, int colunm, int row)
{
    std::cout << "{" << std::endl;
    for(int i = 0; i < colunm; i++)
    {
        for(int j = 0; j < row; j++)
        {
            std::cout << std::boolalpha << Matrix[i * row + j] << ",";
        }
        std::cout << std::endl;
    }
    std::cout << "}" << std::endl;
}

double norm(double* vector, int vector_size)
{
    double norm = 0;

    for(int i = 0; i < vector_size; i++)
    {
        norm += vector[i] * vector[i];
    }

    return std::sqrt(norm);
}

double determinant(double* Matrix, int Matrix_colunm_, int Matrix_row_)
{
    double determinant_solution = 0.0;
    int sign = 1;
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    if(Matrix_colunm == 2)return Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2];
    
    double* child_Matrix = new double[Matrix_colunm * Matrix_row];

    for(int i = 0; i < Matrix_row; i++)
    {
        int n = 0;
        for(int j = 0; j < Matrix_colunm; j++)
        {
            if(j == 0)continue;
            int m = 0;
            for(int k = 0; k < Matrix_row; k++)
            {
                if(k == i)continue;
                child_Matrix[n * (Matrix_row - 1) + m] = Matrix[j * Matrix_row + k];
                m++;
            }
            n++;
        }

        double a = determinant( child_Matrix, Matrix_colunm - 1, Matrix_row - 1);

        determinant_solution += sign * Matrix[i] * a;
        sign *= -1;
    }

    delete child_Matrix;

    return determinant_solution;
}

double* transpose(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* Matrix_t)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)Matrix_t[j * Matrix_colunm + i] = Matrix[i * Matrix_row + j];
    
    return Matrix_t;
}

double* sub(double* Matrix_A, double* Matrix_B, int Matrix_colunm, int Matrix_row, double* Matrix)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)Matrix[i * Matrix_row + j] = Matrix_A[i * Matrix_row + j] - Matrix_B[i * Matrix_row + j];
    return Matrix;
}

double* add(double* Matrix_A, double* Matrix_B, int Matrix_colunm, int Matrix_row, double* Matrix)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)Matrix[i * Matrix_row + j] = Matrix_A[i * Matrix_row + j] + Matrix_B[i * Matrix_row + j];
    return Matrix;
}

//double Strassen(double* Matrix_A, double* Matrix_B, double* Matrix_C, int A_col, int A_row, int B_col, int B_row, int C_col, int C_row)
//{
//    if(N == 1)return {Matrix_A[0] * Matrix_B[0]};
//
//    for()
//}

double* dot(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_, double* Matrix)
{
    const unsigned int Matrix_A_colunm = Matrix_A_colunm_;
    const unsigned int Matrix_A_row    = Matrix_A_row_;
    const unsigned int Matrix_B_colunm = Matrix_B_colunm_;
    const unsigned int Matrix_B_row    = Matrix_B_row_;

    if(Matrix_A_row != Matrix_B_colunm)return nullptr;

    double* tmpMatrix = new double[Matrix_A_colunm * Matrix_B_row];
    for(int i = 0; i < Matrix_A_colunm; i++)
    {
        for(int j = 0; j < Matrix_B_row; j++)
        {
            tmpMatrix[i * Matrix_B_row + j] = 0;
        }
    }

    for(int i = 0; i < Matrix_A_colunm; i++)
    {
        for(int j = 0; j < Matrix_B_row; j++)
        {
            for(int k = 0; k < Matrix_A_row; k++)
            {
                tmpMatrix[i * Matrix_B_row + j] += Matrix_A[i * Matrix_A_row + k] * Matrix_B[k * Matrix_B_row + j]; 
            }
        }
    }
    
    for(int i = 0; i < Matrix_A_colunm; i++)
    {
        for(int j = 0; j < Matrix_B_row; j++)
        {
            Matrix[i * Matrix_B_row + j] = tmpMatrix[i * Matrix_B_row + j];
        }
    }

    delete tmpMatrix;
    return Matrix;
}

//QR分解を行う(Q : 規格直交行列、R : 上三角行列)
bool HouseHolderTransform(double* Matrix, int Matrix_colunm, int Matrix_row, double* R_Matrix,double* Q_Matrix)
{
    //std::cout << "HouseHolder Transform..." << std::endl;
    double* u_vector = new double[Matrix_colunm];
    double* H_Matrix = new double[Matrix_colunm * Matrix_colunm];
    int m = Matrix_colunm;

    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < Matrix_colunm; j++)
        {
            Q_Matrix[i + j * Matrix_colunm] = (i == j ? 1 : 0);
            R_Matrix[i + j * Matrix_colunm] = Matrix[i + j * Matrix_colunm];
        }
    }
    for(int i = 0; i < Matrix_row - 1; i++)
    {
        double absx = 0;
        for(int j = i; j < Matrix_colunm; j++)absx += R_Matrix[i + j * Matrix_row] * R_Matrix[i + j * Matrix_row];
        absx = sqrt(absx);
        if(absx == 0)continue;

        u_vector[i] = R_Matrix[i + i * Matrix_row] + (R_Matrix[i + i * Matrix_row] >= 0 ?  1 : -1 ) * absx;
        double absu = u_vector[i] * u_vector[i];
        for(int j = i + 1; j < Matrix_row; j++)
        {
            u_vector[j] = R_Matrix[i + j * Matrix_row];
            absu += u_vector[j] * u_vector[j];
        }

        for(int k = 0; k < Matrix_colunm; k++)
        {
            for(int j = 0; j < Matrix_colunm; j++)
            {
                H_Matrix[k + j * Matrix_row] = (k == j ? 1 : 0);
            }   
        }
        for(int n = i; n < Matrix_colunm; n++)
        {
            for(int m = i; m < Matrix_colunm; m++)
            {
                H_Matrix[n * Matrix_colunm + m] -= 2 * u_vector[n] * u_vector[m] / absu;
            }
        }

        R_Matrix = dot(H_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix);
        Q_Matrix = dot(Q_Matrix,Matrix_colunm,Matrix_colunm,H_Matrix,Matrix_colunm,Matrix_colunm,Q_Matrix);
    }
    //for(int i = 0; i < Matrix_colunm; i++)
    //{
    //    for(int j = Matrix_row - 1; j > 0; j--)
    //    {
    //        if(Q_Matrix[i * Matrix_row + j] == 0)continue;
    //        else
    //        {
    //            Q_Matrix[i * Matrix_row + j] *= -1;
    //            break;
    //        }
    //    }
    //}
    //std::cout << "HouseHolder Transform end..." << std::endl;
    delete u_vector;
    delete H_Matrix;
    return true;
}

bool checkDiagonal(double* Matrix, int Matrix_colunm, int Matrix_row, double epsilon)
{
    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < Matrix_row; j++)
        {
            if(i == j)continue;
            if(std::abs(Matrix[i * Matrix_row + j]) > epsilon)return false;
        }
    }
    return true;
}

double* eig(double* Matrix, int Matrix_colunm, int Matrix_row,double* eigValue)
{
    std::cout << "finding eigenvalues..." << std::endl;
    const int maxIter = 1000;
    int iter = 0;
    double S = 0;

    double* Matrix_buff = new double[Matrix_colunm * Matrix_row];
    double* Q = new double[Matrix_colunm * Matrix_row];
    double* R = new double[Matrix_colunm * Matrix_row];
    
    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j< Matrix_row; j++)
        {
            Matrix_buff[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
        }
    }

    while(true)
    {
        S = Matrix_buff[Matrix_colunm * Matrix_row - 1];
        for(int i = 0; i < Matrix_colunm; i++)Matrix_buff[i * Matrix_row + i] -= S;
        HouseHolderTransform(Matrix_buff,Matrix_colunm,Matrix_row,R,Q);
        //if(0 <= iter && iter < 10)showMatrix(R,Matrix_colunm,Matrix_row);
        //if(0 <= iter && iter < 10)showMatrix(Q,Matrix_colunm,Matrix_row);
        dot(R,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_row,Matrix_buff);
        for(int i = 0; i < Matrix_colunm; i++)Matrix_buff[i * Matrix_row + i] += S;

        //if(0 <= iter && iter < 10)showMatrix(Matrix_buff,Matrix_colunm,Matrix_row);

        if(checkDiagonal(Matrix_buff,Matrix_colunm,Matrix_row))break;
        else iter+= 1;
        
        if(iter > maxIter)break;
    }

    for(int i = 0; i < Matrix_colunm; i++)eigValue[i] = Matrix_buff[i * Matrix_row + i];

    delete Matrix_buff;
    delete Q;
    delete R;
    std::cout << "found eigenvalues..." << std::endl;
    return eigValue;
}

//void jacobi(double* Matrix, int Matrix_colunm, int Matrix_row, double* eigVector, double* eigValue)
//{
//    const double epsilon = 1e-10;
//    int p = 0;
//    int q = 0;
//    double theta = 0;
//    double cos_value = 0;
//    double sin_value = 0;
//    double t = 0;
//
//    double* Matrix_buffer = new double[Matrix_colunm * Matrix_row];
//
//    for(int i = 0; i < Matrix_row; i++)
//    {
//        for(int j = 0; j < Matrix_row; j++)
//        {
//            eigVector[i * Matrix_row + j] = (i == j ? 1 : 0);
//            Matrix_buffer[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
//        }
//    }
//
//    while(true)
//    {
//        double max = 0;
//        for(int i = 0; i < Matrix_colunm - 1; i++)
//        {
//            for(int j = i + 1; j < Matrix_row; j++)
//            {
//                if(fabs(Matrix_buffer[i * Matrix_row + j] > max))
//                {
//                    max = fabs(Matrix_buffer[i * Matrix_row + j]);
//                    p = i;
//                    q = i;
//                }
//            }
//        }
//
//        if(max < epsilon)break;
//
//        if(Matrix_buffer[p * Matrix_row + p] == Matrix_buffer[q * Matrix_row + q])
//        {
//            theta = M_PI / 4;
//        }
//        else
//        {
//            theta = 0.5 * atan(2 * Matrix_buffer[p * Matrix_row + q] / (Matrix_buffer[p * Matrix_row + p] - Matrix_buffer[q * Matrix_row + q]));
//        }
//
//        cos_value = cos(theta);
//        sin_value = sin(theta);
//
//        for(int i = 0; i < Matrix_colunm; i++)
//        {
//            t = eigVector[i * Matrix_row + p];
//            eigVector[i * Matrix_row + p] = cos_value * t - sin_value * eigVector[i * Matrix_row + q];
//            eigVector[i * Matrix_row + q] = sin_value * t + cos_value * eigVector[i * Matrix_row + q];
//        }
//
//        for(int i = 0; i < Matrix_colunm; i++)
//        {
//            if(i != p && i != q)
//            {
//                t = Matrix_buffer[i * Matrix_row + p];
//                Matrix_buffer[i * Matrix_row + p] = cos_value * t - sin_value * Matrix_buffer[i * Matrix_row + q];
//                Matrix_buffer[i * Matrix_row + q] = sin_value * t + cos_value * Matrix_buffer[i * Matrix_row + p];
//                Matrix_buffer[p * Matrix_row + i] = Matrix_buffer[i * Matrix_row + p];
//                Matrix_buffer[q * Matrix_row + i] = Matrix_buffer[i * Matrix_row + q];
//             }
//        }
//        
//        t = Matrix[p * Matrix_row + p];
//        Matrix_buffer[p * Matrix_row + p] = cos_value * cos_value * t - 2 * sin_value * cos_value * Matrix_buffer[p * Matrix_row + q] + sin_value * sin_value * Matrix_buffer[q * Matrix_row + q];
//        Matrix_buffer[q * Matrix_row + q] = sin_value * sin_value * t + 2 * sin_value * cos_value * Matrix_buffer[p * Matrix_row + q] + cos_value * cos_value * Matrix_buffer[q * Matrix_row + q];
//        Matrix_buffer[p * Matrix_row + q] = 0;
//        Matrix_buffer[q * Matrix_row + p] = 0;
//    }
//        for(int i = 0; i < Matrix_colunm; i++)
//        {
//            eigValue[i] = Matrix_buffer[i * Matrix_row + i];
//        }
//
//    delete Matrix_buffer;
//}



double* inverse_lower_triangular_matrix(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    for(int i = 0; i < Matrix_colunm; i++)
    {
        double* x_vector = new double[Matrix_colunm];
        
        for(int j = 0; j < Matrix_colunm; j++)x_vector[j] = 0;
        x_vector[i] = 1.0;

        for(int j = i + 1; j < Matrix_row; j++)
        {
            double tmp = 0;
            for(int k = 0; k < j; k++)
            {
                tmp += Matrix[j * Matrix_row + k] * x_vector[k];
            }
            x_vector[j] -= tmp;
        }
        for(int j = 0; j < Matrix_row; j++)
        {
            inv_Matrix[j * Matrix_row + i] = x_vector[j];
        }
        delete x_vector;
    }
    return inv_Matrix;
}

double* inverse_upper_triangular_matrix(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    for(int i = Matrix_row - 1; i >= 0; i--)
    {
        double* x_vector = new double[Matrix_colunm];
        
        for(int j = 0; j < Matrix_colunm; j++)x_vector[j] = 0;
        x_vector[i] = 1.0 / Matrix[i * Matrix_row + i];

        for(int j = i - 1; j >= 0; j--)
        {
            double tmp = 0;
            for(int k = Matrix_row - 1; k > j; k--)
            {
                tmp += Matrix[j * Matrix_row + k] * x_vector[k];
            }
            x_vector[j] -= tmp / Matrix[j * Matrix_row + j];
        }
        for(int j = 0; j < Matrix_row; j++)
        {
            inv_Matrix[j * Matrix_row + i] = x_vector[j];
        }

        delete x_vector;
    }
    return inv_Matrix;
}

//LU分解を行う(mode = 1は逆行列用、mode = 0は連立方程式用)
bool LU_Decompose(int mode, double* Matrix, int Matrix_colunm_, int Matrix_row_, double* L, double* U, double* Q)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    for(int i = 0; i < Matrix_colunm; i++)    
    {
        for(int j = 0; j < Matrix_row; j++)
        {
            U[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
            if(mode == 1)Q[i * Matrix_row + j] = (i == j ? 1.0 : 0.0);;
        }
    }

    for(int i = 0; i < Matrix_colunm - 1; i++)
    {
        int holdValueIndex = 0;
        double holdValue = 0.0;
        for(int j = i; j < Matrix_row; j++)
        {
            if(holdValue < abs(U[i + Matrix_row * j]))
            {
                holdValue = abs(U[i + Matrix_row * j]);
                holdValueIndex = j;
            }
        }

        if(holdValueIndex != i)
        {
            for(int j = 0; j < Matrix_colunm; j++)std::swap(U[j + Matrix_row * holdValueIndex],U[j + Matrix_row * i]);
            for(int j = 0; j < Matrix_colunm; j++)std::swap(L[j + Matrix_row * holdValueIndex],L[j + Matrix_row * i]);
            if(mode == 1)for(int j = 0; j < Matrix_colunm; j++)std::swap(Q[j + Matrix_row * holdValueIndex],Q[j + Matrix_row * i]);
            else if(mode == 0)std::swap(Q[holdValueIndex],Q[i]);
        }

        for(int j = 0; j < i; j++)L[j * Matrix_row + i] = 0.0;
        L[i * Matrix_row + i] = 1.0;

        for(int j = i + 1; j < Matrix_row; j++)
        {
            double cat = U[j * Matrix_row + i] / U[i * Matrix_row + i];
            L[j * Matrix_row + i] = cat;
            for(int k = 0; k < Matrix_row; k++)U[j * Matrix_row + k] -= U[i * Matrix_row + k] * cat;
        }

        if(mode == 1)L[Matrix_colunm * Matrix_row - 1] = 1.0;
    }
    return true;
}

double* inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    double* L = new double[Matrix_colunm * Matrix_row];
    double* U = new double[Matrix_colunm * Matrix_row];
    double* Q = new double[Matrix_colunm * Matrix_row];
    double* inverse_L = new double[Matrix_colunm * Matrix_row];
    double* inverse_U = new double[Matrix_colunm * Matrix_row];
    double* P = new double[Matrix_colunm * Matrix_row];

    for(int i = 0; i < Matrix_colunm * Matrix_row; i++)
    {
        L[i] = 0;
        U[i] = 0;
        Q[i] = 0;
        inverse_L[i] = 0;
        inverse_U[i] = 0;
        P[i] = 0;
    }

    LU_Decompose(1,Matrix,Matrix_colunm,Matrix_row,L,U,Q);

    inverse_L = inverse_lower_triangular_matrix(L,Matrix_colunm,Matrix_row,inverse_L);
    inverse_U = inverse_upper_triangular_matrix(U,Matrix_colunm,Matrix_row,inverse_U);

    //showMatrix(inverse_L,Matrix_colunm,Matrix_row);
    //showMatrix(inverse_U,Matrix_colunm,Matrix_row);
    //showMatrix(Q,Matrix_colunm,Matrix_row);

    P = dot(inverse_U,Matrix_colunm,Matrix_row,inverse_L,Matrix_colunm,Matrix_row,P);
    inv_Matrix = dot(P,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_row,inv_Matrix);
    //showMatrix(inv_Matrix,Matrix_colunm,Matrix_row);

    delete L;
    delete U;
    delete Q;
    delete inverse_L;
    delete inverse_U;
    delete P;

    return inv_Matrix;
}

double* equal_solve(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* vector, int vector_size_, double* solution_vector)
{
    std::cout << "equal solve.." << std::endl;
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    const unsigned int vector_size = vector_size_;

    double* L = new double[Matrix_colunm * Matrix_row];
    double* U = new double[Matrix_colunm * Matrix_row];
    double* y = new double[Matrix_colunm];

    for(int i = 0; i < Matrix_colunm * Matrix_row; i++)
    {
        L[i] = 0;
        U[i] = 0;
    }

    LU_Decompose(0,Matrix,Matrix_colunm,Matrix_row,L,U,vector);
    y[0] = vector[0];

    for(int i = 1; i < vector_size; i++)
    {
        y[i] = vector[i];
        for(int j = 0; j < i; j++)y[i] -= L[i * Matrix_row + j] * y[j];
    }

    solution_vector[vector_size - 1] = y[vector_size - 1] / U[Matrix_colunm * Matrix_row - 1];

    for(int i = vector_size - 2; i >= 0; i--)
    {
        solution_vector[i] = y[i];
        for(int j = i + 1; j < vector_size; j++)solution_vector[i] -= U[i * Matrix_row + j] * solution_vector[j]; 
        solution_vector[i] /= U[i * Matrix_row + i];
    }

    delete L;
    delete U;
    delete y;

    std::cout << "equal solve end..." << std::endl;
    return solution_vector;
}

void SVD(double* Matrix, int Matrix_colunm, int Matrix_row, double* U, double* S, double* V,int maxIter)
{
    std::cout << "start SVD..." << std::endl;

    int smaller_Size = std::min(Matrix_colunm,Matrix_row);
    int larger_Size = std::max(Matrix_colunm,Matrix_row);
    int diffrence_of_col_and_row = std::abs(Matrix_colunm - Matrix_row);
    bool eigValueSolution_flag = (Matrix_colunm > Matrix_row ? true : false);

    double* Matrix_t = new double[Matrix_colunm * Matrix_row];
    double* MMt = new double[Matrix_colunm * Matrix_colunm];
    double* MtM = new double[Matrix_row * Matrix_row];
    double* u_eigValues = new double[Matrix_colunm];
    double* v_eigValues = new double[Matrix_row];
    double* tmp_for_U = new double[Matrix_colunm * Matrix_colunm];
    double* tmp_for_V = new double[Matrix_row * Matrix_row];
    double* u_vector = new double[Matrix_colunm];
    double* v_vector = new double[Matrix_row];
    bool* isExist = new bool[larger_Size];

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_colunm; j++)U[i * Matrix_colunm + j] = 0;
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)S[i * Matrix_row + j] = 0;
    for(int i = 0; i < Matrix_row; i++)for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = 0;
    for(int i = 0; i < larger_Size; i++)isExist[i] = false;

    Matrix_t = transpose(Matrix,Matrix_colunm,Matrix_row,Matrix_t);
    MMt = dot(Matrix,Matrix_colunm,Matrix_row,Matrix_t,Matrix_row,Matrix_colunm,MMt);
    MtM = dot(Matrix_t,Matrix_row,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,MtM);

    u_eigValues = eig(MMt,Matrix_colunm,Matrix_colunm,u_eigValues);
    v_eigValues = eig(MtM,Matrix_row,Matrix_row,v_eigValues);

    double* eigValues = (Matrix_colunm > Matrix_row ? v_eigValues : u_eigValues);
    double* another_eigValues = (Matrix_colunm > Matrix_row ? u_eigValues : v_eigValues);

    for(int i = 0; i < smaller_Size; i++)
    {
        double eigValue = eigValues[i];

        for(int j = 0; j < larger_Size; j++)if(std::abs(another_eigValues[j] - eigValue) < 1e-10)isExist[j] = true;

        if(eigValue < 0)eigValue = 0;

        S[i * Matrix_row + i] = std::sqrt(eigValue);
        int iter = 0;

        for(int n = 0; n < Matrix_colunm; n++)for(int m = 0; m < Matrix_colunm; m++)tmp_for_U[n * Matrix_colunm + m] = (n == m ? eigValue : 0);
        for(int n = 0; n < Matrix_row; n++)for(int m = 0; m < Matrix_row; m++)tmp_for_V[n * Matrix_row + m] = (n == m ? eigValue : 0);
        double* sub_tmp_for_U = sub(MMt,tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
        double* inv_tmp_for_U = inverse(sub_tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
        double* sub_tmp_for_V = sub(MtM,tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
        double* inv_tmp_for_V = inverse(sub_tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
        for(int j = 0; j < Matrix_colunm; j++)u_vector[j] = 1;
        for(int j = 0; j < Matrix_row; j++)v_vector[j] = 1;
        while (iter < maxIter)
        {
            u_vector = dot(inv_tmp_for_U,Matrix_colunm,Matrix_colunm,u_vector,Matrix_colunm,1,u_vector);
            double norm_u = norm(u_vector,Matrix_colunm);
            for(int n = 0; n < Matrix_colunm; n++)u_vector[n] /= norm_u;

            v_vector = dot(inv_tmp_for_V,Matrix_row,Matrix_row,v_vector,Matrix_row,1,v_vector);
            double norm_v = norm(v_vector,Matrix_row);
            for(int n = 0; n < Matrix_row; n++)v_vector[n] /= norm_v;

            iter++;
        }
        for(int j = 0; j < Matrix_colunm; j++)U[j * Matrix_colunm + i] = u_vector[j];
        for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = v_vector[j];
    }

    for(int i = 0, index = 0; i < larger_Size; i++)
    {
        if(isExist[i])continue;
        double eigValue = another_eigValues[i];
        if(Matrix_colunm == larger_Size)
        {
            int iter_for_U = 0;
            for(int n = 0; n < Matrix_colunm; n++)for(int m = 0; m < Matrix_colunm; m++)tmp_for_U[n * Matrix_colunm + m] = (n == m ? eigValue : 0);
            double* sub_tmp_for_U = sub(MMt,tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
            double* inv_tmp_for_U = inverse(sub_tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
            for(int j = 0; j < Matrix_colunm; j++)u_vector[j] = 1;
            while (iter_for_U < maxIter)
            {
                u_vector = dot(inv_tmp_for_U,Matrix_colunm,Matrix_colunm,u_vector,Matrix_colunm,1,u_vector);
                double norm_u = norm(u_vector,Matrix_colunm);
                for(int n = 0; n < Matrix_colunm; n++)u_vector[n] /= norm_u;
                iter_for_U++;
            }
            for(int j = 0; j < Matrix_colunm; j++)
            {
                U[j * Matrix_colunm + (index + smaller_Size)] = u_vector[j];
            }
        }
        if(Matrix_row == larger_Size)
        {
            int iter_for_V = 0;
            for(int n = 0; n < Matrix_row; n++)for(int m = 0; m < Matrix_row; m++)tmp_for_V[n * Matrix_row + m] = (n == m ? eigValue : 0);
            double* sub_tmp_for_V = sub(MtM,tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
            double* inv_tmp_for_V = inverse(sub_tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
            for(int j = 0; j < Matrix_row; j++)v_vector[j] = 1;
            while (iter_for_V < maxIter)
            {
                v_vector = dot(inv_tmp_for_V,Matrix_row,Matrix_row,v_vector,Matrix_row,1,v_vector);
                double norm_v = norm(v_vector,Matrix_row);
                for(int n = 0; n < Matrix_row; n++)v_vector[n] /= norm_v;
                iter_for_V++;
            }
            for(int j = 0; j < Matrix_row; j++)
            {
                V[(index + smaller_Size) * Matrix_row + j] = v_vector[j];
            }
        }
        index++;
    }

    //固有値ベクトルの向きを合わせる
    for(int k = 0; k < Matrix_row; k++)
    {
        double A_v_multi = 0;
        double sigma_u_multi = U[k] * S[k * Matrix_row + k];
        for(int i = 0; i < Matrix_row; i++)
        {
            A_v_multi += Matrix[i] * V[k * Matrix_row + i];
        }

        if(std::abs(A_v_multi - sigma_u_multi) < 1e-10)continue;

        for(int i = 0; i < Matrix_row; i++)V[k * Matrix_row + i] *= -1;
    }

    delete Matrix_t;
    delete MMt;
    delete MtM;
    delete u_eigValues;
    delete v_eigValues;
    delete tmp_for_U;
    delete tmp_for_V;
    delete u_vector;
    delete v_vector;
    delete isExist;

    std::cout << "SVD end..." << std::endl;
}

void NMF(double* Matrix, int Matrix_colunm, int Matrix_row, int n_components, double* W, double* H, double epsilon)
{
    std::mt19937 mt{ std::random_device{}() };
    double* H_T = new double[Matrix_row * n_components];
    double* W_T = new double[Matrix_colunm * n_components];
    double* YH = new double[n_components * Matrix_colunm];
    double* WH_TH = new double[n_components * Matrix_colunm];
    double* W_TY = new double[n_components * Matrix_row];
    double* W_TWH_T_T = new double[n_components * Matrix_row];
    double* WH_T = new double[Matrix_colunm * Matrix_row];
    double maxValue = std::numeric_limits<double>::min();
    double minValue = std::numeric_limits<double>::max();

    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < Matrix_row; j++)
        {
            maxValue = std::max(Matrix[i * Matrix_row + j], maxValue);
            minValue = std::min(Matrix[i * Matrix_row + j], minValue);
        }
    }

    std::uniform_real_distribution<double> random(minValue, maxValue);

    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < Matrix_row; j++)
        {
            for(int k = 0; k < n_components; k++)
            {
                W[i * n_components + k] = random(mt);
                H_T[k * Matrix_row + j] = random(mt);
            }
        }
    }


    while(true)
    {
        YH = dot(Matrix,Matrix_colunm,Matrix_row,H_T,Matrix_row,n_components,YH);
        H = transpose(H_T,Matrix_row,n_components,H);
        WH_T = dot(W,Matrix_colunm,n_components,H,n_components,Matrix_row,WH_T);
        WH_TH = dot(WH_T,Matrix_colunm,Matrix_row,H_T,Matrix_row,n_components,WH_TH);
        W_T = transpose(W,Matrix_colunm,n_components,W_T);
        W_TY = dot(W_T,n_components,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,W_TY);
        W_TWH_T_T = dot(W_T,n_components,Matrix_colunm,WH_T,Matrix_colunm,Matrix_row,W_TWH_T_T);
        W_TWH_T_T = transpose(W_TWH_T_T,n_components,Matrix_row,W_TWH_T_T);

        for(int i = 0; i < Matrix_colunm; i++)
        {
            for(int j = 0; j < Matrix_row; j++)
            {
                for(int k = 0; k < n_components; k++)
                {
                    W[i * n_components + k] = W[i * n_components + k] * YH[i * n_components + k] / WH_TH[i * n_components + k];
                    H_T[k * Matrix_row + j] = H_T[k * Matrix_row + j] * W_TY[k * Matrix_row + j] / W_TWH_T_T[k * Matrix_row + j];
                }
            }
        }
        
        double residues = 0;
        WH_T = dot(W,Matrix_colunm,n_components,H,n_components,Matrix_row,WH_T);
        for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)residues += (WH_T[i * Matrix_row + j] - Matrix[i * Matrix_row + j]) * (WH_T[i * Matrix_row + j] - Matrix[i * Matrix_row + j]);

        if(residues < epsilon)break;
    }

    H = transpose(H_T,Matrix_row,n_components,H);

    delete H_T;
    delete YH;
    delete WH_TH;
    delete W_TY;
    delete W_TWH_T_T;
    delete WH_T;

}

int Matrix_rank(double* Matrix, int Matrix_colunm, int Matrix_row, bool DoSVD, double* U, double* S, double* V)
{
    int rank = 0;
    if(DoSVD)
    {
        if(U == nullptr || S == nullptr || V == nullptr)return -1;
        SVD(Matrix,Matrix_colunm,Matrix_row,U,S,V);
    }

    for(int i = 0; i < (Matrix_colunm > Matrix_row ? Matrix_row : Matrix_colunm); i++)if((DoSVD ? S[i * Matrix_row + i] : Matrix[i * Matrix_row + i]) != 0)rank += 1;

    return rank;
}