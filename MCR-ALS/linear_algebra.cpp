#include "linear_algebra.hpp"

bool DEBUGFLAG = false;

void showMatrix(double* Matrix, int colunm, int row, bool allShow)
{
    double epsilon = std::numeric_limits<double>::epsilon();

    if((colunm <= 10 && row <= 10)|| allShow)
    {
        std::cout << "{" << std::endl;
        for(int i = 0; i < colunm; i++)
        {
            std::cout << "{";
            for(int j = 0; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }
        std::cout << "}" << std::endl;
    }
    else if ((colunm <= 10 && row > 10) || allShow)
    {
        std::cout << "{" << std::endl;
        for(int i = 0; i < colunm; i++)
        {
            std::cout << "{";
            for(int j = 0; j < 3; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "...";
            for(int j = row - 3; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }
        std::cout << "}" << std::endl;
    }
    else if ((colunm > 10 && row <= 10)|| allShow)
    {
        std::cout << "{" << std::endl;
        for(int i = 0; i < 3; i++)
        {
            std::cout << "{";
            for(int j = 0; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }

        std::cout << "..." << std::endl;

        for(int i = colunm - 3; i < colunm; i++)
        {
            std::cout << "{";
            for(int j = 0; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }
        std::cout << "}" << std::endl;
    }
    else
    {
        std::cout << "{" << std::endl;
        for(int i = 0; i < 3; i++)
        {
            std::cout << "{";
            for(int j = 0; j < 3; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "...";
            for(int j = row - 3; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }

        std::cout << "..." << std::endl;

        for(int i = colunm - 3; i < colunm; i++)
        {
            std::cout << "{";
            for(int j = 0; j < 3; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "...";
            for(int j = row - 3; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < epsilon)std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }
        std::cout << "}" << std::endl;
    }

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

double frobenius_norm(double* Matrix, int Matrix_colunm, int Matrix_row)
{
    double fro_norm_sum = 0;

    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < Matrix_row; j++)
        {
            fro_norm_sum += Matrix[i * Matrix_row + j] * Matrix[i * Matrix_row + j];
        }
    }

    return std::sqrt(fro_norm_sum);
}

double determinant(double* Matrix, int Matrix_colunm_, int Matrix_row_)
{
    double determinant_solution = 0.0;
    int sign = 1;
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    if(Matrix_colunm == 2)return Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2];
    
    double* child_Matrix = new double[Matrix_colunm * Matrix_row];

    for(int i = 0; i < Matrix_row_; i++)
    {
        int n = 0;
        for(int j = 0; j < Matrix_colunm_; j++)
        {
            if(j == 0)continue;
            int m = 0;
            for(int k = 0; k < Matrix_row_; k++)
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

    delete[] child_Matrix;

    return determinant_solution;
}

double* transpose(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* Matrix_t)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    if(Matrix == Matrix_t)
    {
        double* buffer = new double[Matrix_colunm * Matrix_row];

        for(int i = 0; i < Matrix_colunm_; i++)for(int j = 0; j < Matrix_row_; j++)buffer[j * Matrix_colunm + i] = Matrix[j * Matrix_colunm + i];

        for(int i = 0; i < Matrix_colunm_; i++)for(int j = 0; j < Matrix_row_; j++)Matrix_t[j * Matrix_colunm + i] = buffer[i * Matrix_row + j];

        delete[] buffer;
    }
    else
    {
        for(int i = 0; i < Matrix_colunm_; i++)for(int j = 0; j < Matrix_row_; j++)Matrix_t[j * Matrix_colunm + i] = Matrix[i * Matrix_row + j];
    }

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

double dot(double* vector_A, double* vector_B, int vectorSize)
{
    double sum = 0;
    for(int i = 0; i < vectorSize; i++)sum += vector_A[i] * vector_B[i];
    return sum;
}

//double Strassen(double* Matrix_A, double* Matrix_B, double* Matrix_C, int A_col, int A_row, int B_col, int B_row, int C_col, int C_row)
//{
//    if(N == 1)return {Matrix_A[0] * Matrix_B[0]};
//
//    for()
//}

double* squareMatrix_product(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_, double* Matrix)
{
    const std::size_t double_byte = sizeof(double);
    const std::size_t L1_cache_size = 704 * 1024;
    const std::size_t L2_cache_size = 7 * 1024 * 1024;
    const std::size_t L3_cache_size = 12 * 1024 * 1024;

    const int L1_cache_Matrix_size = (int)(L1_cache_size / double_byte);
    const int L2_cache_Matrix_size = (int)(L2_cache_size / double_byte);
    const int L3_cache_Matrix_size = (int)(L3_cache_size / double_byte);



    double* tmpMatrix = new double[Matrix_A_colunm_ * Matrix_B_row_];

    for(int i = 0; i < Matrix_A_colunm_; i++)
    {
        for(int j = 0; j < Matrix_B_row_; j++)
        {
            tmpMatrix[i * Matrix_B_row_ + j] = 0;
        }
    }

    auto L1_cache_prod = [Matrix_A_row_,Matrix_B_row_,L1_cache_size](double* MatrixA,double* MatrixB,double* MatrixC)
    {
        #pragma omp parallel for
        for(int i = 0; i < L1_cache_size; i++)
        {
            #pragma omp parallel for
            for(int k = 0; k < L1_cache_size; k++)
            {
                #pragma omp parallel for
                for(int j = 0; j < L1_cache_size; j++)
                {
                    MatrixC[i * Matrix_B_row_ + j] += MatrixA[i * Matrix_A_row_ + k] * MatrixB[k * Matrix_B_row_ + j];
                }
            }
        }
    };

    auto L2_cache_prod = [Matrix_A_row_,Matrix_B_row_,L1_cache_size,L2_cache_size,L1_cache_prod](double* MatrixA,double* MatrixB,double* MatrixC)
    {
        for(int i = 0; i < L2_cache_size; i+=L1_cache_size)
        {
            for(int j = 0; j < L2_cache_size; j+=L1_cache_size)
            {
                for(int k = 0; k < L2_cache_size; k+=L1_cache_size)
                {
                    double* A = &MatrixA[i * Matrix_A_row_ + k];
                    double* B = &MatrixB[k * Matrix_B_row_ + j];
                    double* C = &MatrixC[i * Matrix_B_row_ + j];
                    L1_cache_prod(A,B,C);
                }
            }
        }
    };

    for(int i = 0; i < L3_cache_size; i+=L2_cache_size)
    {
        for(int j = 0; j < L3_cache_size; j+=L2_cache_size)
        {
            for(int k = 0; k < L3_cache_size; k+=L2_cache_size)
            {
                double* A = &Matrix_A[i * Matrix_A_row_ + k];
                double* B = &Matrix_B[k * Matrix_B_row_ + j];
                double* C = &tmpMatrix[i * Matrix_B_row_ + j];
                L2_cache_prod(A,B,C);
            }
        }
    }
        
    for(int i = 0; i < Matrix_A_colunm_; i++)
    {
        for(int j = 0; j < Matrix_B_row_; j++)
        {
            Matrix[i * Matrix_B_row_ + j] = tmpMatrix[i * Matrix_B_row_ + j];
        }
    }
    
    delete[] tmpMatrix;
    return Matrix;
}

double* product(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_, double* Matrix)
{
    //if(DEBUGFLAG)showMatrix(Matrix_A,Matrix_A_colunm_,Matrix_A_row_);
    //if(DEBUGFLAG)showMatrix(Matrix_B,Matrix_B_colunm_,Matrix_B_row_);
    if(Matrix_A_row_ != Matrix_B_colunm_)return nullptr;

    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans,  Matrix_A_colunm_, Matrix_B_row_, Matrix_A_row_, 1, Matrix_A, Matrix_A_row_, Matrix_B, Matrix_B_row_, 0, Matrix, Matrix_B_row_);

    //double* tmpMatrix = new double[Matrix_A_colunm_ * Matrix_B_row_];
    //for(int i = 0; i < Matrix_A_colunm_; i++)
    //{
    //    for(int j = 0; j < Matrix_B_row_; j++)
    //    {
    //        tmpMatrix[i * Matrix_B_row_ + j] = 0;
    //    }
    //}
    //
    //#pragma omp parallel for
    //for(int i = 0; i < Matrix_A_colunm_; i++)
    //{
    //    #pragma omp parallel for
    //    for(int k = 0; k < Matrix_A_row_; k++)
    //    {
    //        #pragma omp parallel for
    //        for(int j = 0; j < Matrix_B_row_; j++)
    //        {
    //            tmpMatrix[i * Matrix_B_row_ + j] += Matrix_A[i * Matrix_A_row_ + k] * Matrix_B[k * Matrix_B_row_ + j]; 
    //        }
    //    }
    //}
    //    
    //for(int i = 0; i < Matrix_A_colunm_; i++)
    //{
    //    for(int j = 0; j < Matrix_B_row_; j++)
    //    {
    //        Matrix[i * Matrix_B_row_ + j] = tmpMatrix[i * Matrix_B_row_ + j];
    //    }
    //}
    //if(DEBUGFLAG)std::cout << "(" << Matrix_A_colunm_ << "," << Matrix_A_row_ << ")" << " (" << Matrix_B_colunm_ << "," << Matrix_B_row_ << ")" << std::endl;
    //
    //delete[] tmpMatrix;
    return Matrix;
}

//ギブンス回転行列を用いたQR分解(Q : 規格直交行列、R : 上三角行列)
bool QR_Decompose(double* Matrix, int Matrix_colunm, int Matrix_row, double* Q_Matrix,double* R_Matrix)
{
    std::vector<std::pair<int,int>> lower_tril;

    double* G_Matrix = new double[Matrix_colunm * Matrix_colunm];

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_colunm; j++)Q_Matrix[i * Matrix_colunm + j] = (i == j ? 1 : 0);
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)R_Matrix[i * Matrix_row + j] = Matrix[i * Matrix_row + j];

    for(int i = 1; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < (Matrix_row > i ? i : Matrix_row); j++)
        {
            lower_tril.push_back({j,i});
        }
    }

    for(std::pair<int,int> tril : lower_tril)
    {
        int j = tril.first;
        int i = tril.second;

        //std::cout << "(" << i  << " " << j << ")" << std::endl;

        if(R_Matrix[i * Matrix_row + j] != 0)
        {
            for(int n = 0; n < Matrix_colunm; n++)
            {
                for(int m = 0; m < Matrix_colunm; m++)
                {
                    G_Matrix[n * Matrix_colunm + m] = (n == m ? 1 : 0);
                }
            }

            double r = std::sqrt(R_Matrix[j * Matrix_row + j] * R_Matrix[j * Matrix_row + j] + R_Matrix[i * Matrix_row + j] * R_Matrix[i * Matrix_row + j]);
            double c = R_Matrix[j * Matrix_row + j] / r;
            double s = - R_Matrix[i * Matrix_row + j] / r;
            //std::cout << "(" << r << "," << c << "," << s << ")" << std::endl;
            
            G_Matrix[i * Matrix_colunm + i] = c;
            G_Matrix[j * Matrix_colunm + j] = c;
            G_Matrix[i * Matrix_colunm + j] = s;
            G_Matrix[j * Matrix_colunm + i] = -s;

            //showMatrix(G_Matrix,Matrix_colunm,Matrix_colunm);

            R_Matrix = product(G_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix,Matrix_colunm,Matrix_row,R_Matrix);
            double* G_Matrix_t = transpose(G_Matrix,Matrix_colunm,Matrix_colunm,G_Matrix);
            Q_Matrix= product(Q_Matrix,Matrix_colunm,Matrix_colunm,G_Matrix_t,Matrix_colunm,Matrix_colunm,Q_Matrix);

            //showMatrix(Q_Matrix,Matrix_colunm,Matrix_colunm);
            //showMatrix(R_Matrix,Matrix_colunm,Matrix_row);
        }
    }
    
    delete[] G_Matrix;

    return true;
}

//正方行列にのみ対応
void QRDecompotion_with_SR(double* Matrix, int Matrix_colunm, int Matrix_row, double* Q_Matrix,double* R_Matrix)
{
    double* Q_i = new double[Matrix_colunm];
    double* Q_k = new double[Matrix_colunm];

    Q_Matrix = copyMatrix(Matrix,Matrix_colunm,Matrix_row,Q_Matrix);
    for(int i = 0; i < Matrix_row; i++)for(int j = 0; j < Matrix_row; j++)R_Matrix[i * Matrix_row + j] = 0;

    #pragma omp parallel for
    for(int k = 0; k < Matrix_row; k++)
    {   
        #pragma omp parallel for
        for(int j = 0; j < Matrix_colunm; j++)Q_k[j] = Q_Matrix[j * Matrix_row + k];

        #pragma omp parallel for
        for(int i = 0; i < k; i++)
        {
            #pragma omp parallel for
            for(int j = 0; j < Matrix_colunm; j++)Q_i[j] = Q_Matrix[j * Matrix_row + i];
                
            R_Matrix[i * Matrix_row + k] = dot(Q_i,Q_k,Matrix_colunm);

            #pragma omp parallel for
            for(int j = 0; j < Matrix_colunm; j++)Q_k[j] -= R_Matrix[i * Matrix_row + k] * Q_Matrix[j * Matrix_row + i];
        }

        R_Matrix[k * Matrix_row + k] = norm(Q_k,Matrix_colunm);

        #pragma omp parallel for
        for(int j = 0; j < Matrix_colunm; j++)Q_Matrix[j * Matrix_row + k] = Q_k[j] / R_Matrix[k * Matrix_row + k];
    }

    #pragma omp parallel for
    for(int i = 0; i < Matrix_colunm; i++)
    {
        #pragma omp parallel for
        for(int j = 0; j < Matrix_row; j++)
        {
            R_Matrix[i * Matrix_row + j] *= -1;
            Q_Matrix[i * Matrix_row + j] *= -1;
        }
    }

    delete[] Q_i;
    delete[] Q_k;
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
        absx = std::sqrt(absx);
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

        R_Matrix = product(H_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix);
        Q_Matrix = product(Q_Matrix,Matrix_colunm,Matrix_colunm,H_Matrix,Matrix_colunm,Matrix_colunm,Q_Matrix);
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
    delete[] u_vector;
    delete[] H_Matrix;
    return true;
}

double* tridiagnalization_With_HouseHolder(double* Matrix, int Matrix_colunm, int Matrix_row, double* tridiagnalMatrix,double* HouseHolderMatrix)
{
    bool needHouseHolderMatrix_flag = (HouseHolderMatrix != nullptr ? true : false);

    double* tmp = new double[Matrix_colunm];
    double* x = new double[Matrix_colunm];
    double* y = new double[Matrix_colunm];
    double* T = tridiagnalMatrix;
    double* H = new double[Matrix_colunm * Matrix_row];
    double* Ht = new double[Matrix_row * Matrix_colunm];

    if(needHouseHolderMatrix_flag)for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)HouseHolderMatrix[i * Matrix_row + j] = (i == j ? 1 : 0);

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)T[i * Matrix_row + j] = Matrix[i * Matrix_row + j];

    if(needHouseHolderMatrix_flag)
    {
        for(int i = 0; i < Matrix_row - 2; i++)
        {
            for(int j = 0; j < Matrix_colunm; j++)
            {
                tmp[j] = 0;
                x[j] = 0;
                y[j] = 0;
            }

            for(int j = i; j < Matrix_colunm; j++)tmp[j] = T[j * Matrix_row + i];
            tmp[i] = 0;
            double s = norm(tmp,Matrix_colunm);

            for(int j = i; j < Matrix_colunm; j++)x[j] = T[j * Matrix_row + i];
            y[i] = T[i * Matrix_row + i];
            y[(i != Matrix_row ? i + 1 : i - 1)] = -s;

            double* v = sub(x,y,Matrix_colunm,1,x);
            double v_norm = norm(v,Matrix_colunm);

            product(v,Matrix_colunm,1,v,1,Matrix_colunm,H);

            for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)H[i * Matrix_row + j] = (i == j ? 1 : 0) - 2 * H[i * Matrix_row + j] / (v_norm * v_norm);

            product(HouseHolderMatrix,Matrix_colunm,Matrix_row,H,Matrix_colunm,Matrix_row,HouseHolderMatrix);

            transpose(H,Matrix_colunm,Matrix_row,Ht);

            product(Ht,Matrix_colunm,Matrix_row,T,Matrix_colunm,Matrix_row,T);
            product(T,Matrix_colunm,Matrix_row,H,Matrix_colunm,Matrix_row,T);
        }
    }
    else
    {
        for(int i = 0; i < Matrix_row - 2; i++)
        {
            for(int j = 0; j < Matrix_colunm; j++)
            {
                tmp[j] = 0;
                x[j] = 0;
                y[j] = 0;
            }

            for(int j = i; j < Matrix_colunm; j++)tmp[j] = T[j * Matrix_row + i];
            tmp[i] = 0;
            double s = norm(tmp,Matrix_colunm);

            for(int j = i; j < Matrix_colunm; j++)x[j] = T[j * Matrix_row + i];
            y[i] = T[i * Matrix_row + i];
            y[(i != Matrix_row ? i + 1 : i - 1)] = -s;

            double* v = sub(x,y,Matrix_colunm,1,x);
            double v_norm = norm(v,Matrix_colunm);

            product(v,Matrix_colunm,1,v,1,Matrix_colunm,H);

            for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)H[i * Matrix_row + j] = (i == j ? 1 : 0) - 2 * H[i * Matrix_row + j] / (v_norm * v_norm);

            transpose(H,Matrix_colunm,Matrix_row,Ht);

            product(Ht,Matrix_colunm,Matrix_row,T,Matrix_colunm,Matrix_row,T);
            product(T,Matrix_colunm,Matrix_row,H,Matrix_colunm,Matrix_row,T);
        }
    }

    delete[] tmp;
    delete[] x;
    delete[] y;
    delete[] H;
    delete[] Ht;

    return tridiagnalMatrix;
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

double* eig(double* Matrix, int Matrix_colunm, int Matrix_row,double* eigenValue,double* eigenVector)
{
    std::cout << "finding eigenvalues..." << std::endl;

        int iter = 0;
    int maxIter = 1000;

    double epsilon = std::numeric_limits<double>::epsilon();

    double* T = new double[Matrix_colunm * Matrix_row];
    double* H = new double[Matrix_colunm * Matrix_row];
    double* Q = new double[Matrix_colunm * Matrix_row];
    double* R = new double[Matrix_colunm * Matrix_row];

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)eigenVector[i * Matrix_row + j] = (i == j ? 1 : 0);

    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j< Matrix_row; j++)
        {
            T[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
        }
    }

    //tridiagnalization_With_HouseHolder(Matrix,Matrix_colunm,Matrix_row,T,H);

    //for(int i = 0; i < Matrix_colunm; i++)
    //{
    //    if(std::abs(T[i * Matrix_row + i]) > epsilon && std::abs(T[i * Matrix_row + i + 1]) > epsilon && i > Matrix_colunm - 1)break;
    //    if(std::abs(T[i * Matrix_row + i]) > epsilon && std::abs(T[i * Matrix_row + i - 1]) > epsilon && i == Matrix_colunm - 1)break;
//
    //    if(i == Matrix_colunm - 1)return false;
    //}

    while (true)
    {
        HouseHolderTransform(T,Matrix_colunm,Matrix_row,R,Q);

        product(R,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_colunm,T);

        product(eigenVector,Matrix_colunm,Matrix_colunm,Q,Matrix_colunm,Matrix_colunm,eigenVector);
        
        if(checkDiagonal(T,Matrix_colunm,Matrix_row,epsilon))break;
        else iter++;
        if(iter > maxIter)break;
    }

    for(int i = 0; i < Matrix_row; i++)eigenValue[i] = T[i * Matrix_row + i];

    //product(H,Matrix_colunm,Matrix_row,eigenVector,Matrix_colunm,Matrix_row,eigenVector);

    delete[] T;
    delete[] H;
    delete[] Q;
    delete[] R;
    //showMatrix(Matrix,Matrix_colunm,Matrix_row);
    //const int maxIter = 1000;
    //int iter = 0;
    //double S = 0;
//
    //double* Matrix_buff = new double[Matrix_colunm * Matrix_row];
    //double* Q = new double[Matrix_colunm * Matrix_row];
    //double* R = new double[Matrix_colunm * Matrix_row];
    //
    //for(int i = 0; i < Matrix_colunm; i++)
    //{
    //    for(int j = 0; j< Matrix_row; j++)
    //    {
    //        Matrix_buff[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
    //    }
    //}
//
    ////tridiagnalization_With_HouseHolder(Matrix,Matrix_colunm,Matrix_row,Matrix_buff);
//
    //while(true)
    //{
    //    S = Matrix_buff[Matrix_colunm * Matrix_row - 1];
    //    for(int i = 0; i < Matrix_colunm; i++)Matrix_buff[i * Matrix_row + i] -= S;
    //    HouseHolderTransform(Matrix_buff,Matrix_colunm,Matrix_row,R,Q);
    //    //if(0 <= iter && iter < 10)showMatrix(R,Matrix_colunm,Matrix_row);
    //    //if(0 <= iter && iter < 10)showMatrix(Q,Matrix_colunm,Matrix_row);
    //    product(R,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_row,Matrix_buff);
    //    for(int i = 0; i < Matrix_colunm; i++)Matrix_buff[i * Matrix_row + i] += S;
//
    //    //if(0 <= iter && iter < 10)showMatrix(Matrix_buff,Matrix_colunm,Matrix_row);
//
    //    if(checkDiagonal(Matrix_buff,Matrix_colunm,Matrix_row))break;
    //    else iter+= 1;
    //    
    //    if(iter > maxIter)break;
    //}
//
    //for(int i = 0; i < Matrix_colunm; i++)eigValue[i] = Matrix_buff[i * Matrix_row + i];
    //delete[] Matrix_buff;
    //delete[] Q;
    //delete[] R;
    std::cout << "found eigenvalues..." << std::endl;
    return eigenValue;
}

void eig_with_divide(double* Matrix, int Matrix_colunm, int Matrix_row,double* eigenValue,double* eigenVector)
{
    //showMatrix(Matrix,Matrix_colunm,Matrix_row);
    const int divide_number = 2;
    const int divide_timesNumber = Matrix_colunm / 2;
    const int subMatrix_size = Matrix_colunm / divide_number;
    const int surplusMatrix_size = subMatrix_size + Matrix_colunm % divide_number;
    const int maxIter = 1000; 
    const double epsilon = std::numeric_limits<double>::epsilon();

    int iter = 0;
    double s = 0;

    double* T = new double[Matrix_colunm * Matrix_row];
    double* H = new double[Matrix_colunm * Matrix_row];
    double* Q_stack = new double[Matrix_colunm * Matrix_row];
    double* Qt_stack = new double[Matrix_row * Matrix_colunm];
    double* sub_Q = new double[subMatrix_size * subMatrix_size];
    double* sub_R = new double[subMatrix_size * subMatrix_size];
    double* sub_Q_stack = new double[subMatrix_size * subMatrix_size];
    double* eigenValueMatrix = new double[Matrix_colunm * Matrix_row];
    double* adjustmentMatrix = new double[Matrix_colunm * Matrix_row];
    double* subMatrix = new double[subMatrix_size * subMatrix_size];
    double* surplusMatrix = new double[surplusMatrix_size * surplusMatrix_size];
    double* surplus_Q = new double[surplusMatrix_size * surplusMatrix_size];
    double* surplus_R = new double[surplusMatrix_size * surplusMatrix_size];
    double* surplus_Q_stack = new double[Matrix_colunm * Matrix_row];

    for(int i = 0; i < Matrix_colunm; i++)
    {
        for(int j = 0; j < Matrix_row; j++)
        {
            eigenValueMatrix[i * Matrix_row + j] = 0;
            adjustmentMatrix[i * Matrix_row + j] = 0;
        }
    }

    tridiagnalization_With_HouseHolder(Matrix,Matrix_colunm,Matrix_row,T,H);
    //showMatrix(T,Matrix_colunm,Matrix_row,true);

    for(int k = 0; k < divide_number - 1; k++)
    {
        for(int i = 0; i < subMatrix_size; i++)
        {
            for(int j = 0; j < subMatrix_size; j++)
            {
                subMatrix[i * subMatrix_size + j] = T[(subMatrix_size * k + i) * Matrix_row + (subMatrix_size * k + j)];
                sub_Q_stack[i * subMatrix_size + j] = (i == j ? 1 : 0);
            }
        }

        double b_k = T[(subMatrix_size * (k + 1) - 1) * Matrix_row + (subMatrix_size * (k + 1))];

        if(k != 0)subMatrix[0] -= std::abs(T[(subMatrix_size * k) * Matrix_row + (subMatrix_size * k - 1)]);
        subMatrix[subMatrix_size * subMatrix_size - 1] -= std::abs(b_k);


        adjustmentMatrix[(subMatrix_size * (k + 1) - 1) * Matrix_row + subMatrix_size * (k + 1) - 1] = std::abs(b_k);
        adjustmentMatrix[(subMatrix_size * (k + 1)) * Matrix_row + subMatrix_size * (k + 1) - 1] = b_k;
        adjustmentMatrix[(subMatrix_size * (k + 1) - 1) * Matrix_row + subMatrix_size * (k + 1)] = b_k;
        adjustmentMatrix[(subMatrix_size * (k + 1)) * Matrix_row + subMatrix_size * (k + 1)] = std::abs(b_k);

        //showMatrix(subMatrix,subMatrix_size,subMatrix_size,true);
        //showMatrix(adjustmentMatrix,Matrix_colunm,Matrix_row,true);

        iter = 0;
        while (true)
        {
            QRDecompotion_with_SR(subMatrix,subMatrix_size,subMatrix_size,sub_Q,sub_R);
            product(sub_R,subMatrix_size,subMatrix_size,sub_Q,subMatrix_size,subMatrix_size,subMatrix);

            product(sub_Q_stack,subMatrix_size,subMatrix_size,sub_Q,subMatrix_size,subMatrix_size,sub_Q_stack);
            
            if(checkDiagonal(subMatrix,subMatrix_size,subMatrix_size,epsilon))break;
            else iter++;

            if(iter > maxIter)break;

        }

        //showMatrix(sub_Q_stack,subMatrix_size,subMatrix_size,true);

        for(int i = 0; i < subMatrix_size; i++)eigenValueMatrix[(subMatrix_size * k + i) * Matrix_row + (subMatrix_size * k + i)] = subMatrix[i * subMatrix_size + i];
        for(int i = 0; i < subMatrix_size; i++)for(int j = 0; j < subMatrix_size; j++)Q_stack[(subMatrix_size * k + i) * Matrix_row + (subMatrix_size * k + j)] = sub_Q_stack[i * subMatrix_size + j];

        //showMatrix(Q_stack,Matrix_colunm,Matrix_row);
    }

    int last_divide_top_index = Matrix_colunm - surplusMatrix_size;

    for(int i = 0; i < surplusMatrix_size; i++)
    {
        for(int j = 0; j < surplusMatrix_size; j++)
        {
            surplusMatrix[i * surplusMatrix_size + j] = T[(last_divide_top_index + i) * Matrix_row + (last_divide_top_index + j)];
            surplus_Q_stack[i * surplusMatrix_size + j] =  (i == j ? 1 : 0);
        }
    }

    surplusMatrix[0] -= std::abs(T[last_divide_top_index * Matrix_row + last_divide_top_index - 1]);

    iter = 0;
    while (true)
    {
        QRDecompotion_with_SR(surplusMatrix,surplusMatrix_size,surplusMatrix_size,surplus_Q,surplus_R);
        product(surplus_R,surplusMatrix_size,surplusMatrix_size,surplus_Q,surplusMatrix_size,surplusMatrix_size,surplusMatrix);

        product(surplus_Q_stack,surplusMatrix_size,surplusMatrix_size,surplus_Q,surplusMatrix_size,surplusMatrix_size,surplus_Q_stack);
        
        if(checkDiagonal(surplusMatrix,surplusMatrix_size,surplusMatrix_size,epsilon))break;
        else iter++;
        if(iter > maxIter)break;
    }

    for(int i = 0; i < surplusMatrix_size; i++)eigenValueMatrix[(last_divide_top_index + i) * Matrix_row + (last_divide_top_index + i)] = surplusMatrix[i * surplusMatrix_size + i];    
    for(int i = 0; i < surplusMatrix_size; i++)for(int j = 0; j < surplusMatrix_size; j++)Q_stack[(last_divide_top_index + i) * Matrix_row + (last_divide_top_index + j)] = sub_Q_stack[i * surplusMatrix_size + j];

    transpose(Q_stack,Matrix_colunm,Matrix_row,Qt_stack);

    product(Qt_stack,Matrix_row,Matrix_colunm,adjustmentMatrix,Matrix_colunm,Matrix_row,Qt_stack);
    
    //showMatrix(Qt_stack,Matrix_colunm,Matrix_row);
    product(Qt_stack,Matrix_row,Matrix_colunm,Q_stack,Matrix_colunm,Matrix_row,Qt_stack);

    //showMatrix(Qt_stack,Matrix_colunm,Matrix_row);

    add(eigenValueMatrix,Qt_stack,Matrix_colunm,Matrix_row,eigenValueMatrix);

    //showMatrix(Q_stack,Matrix_colunm,Matrix_row,true);
    //showMatrix(H,Matrix_colunm,Matrix_row,true);
    product(H,Matrix_colunm,Matrix_row,Q_stack,Matrix_colunm,Matrix_row,eigenVector);
    for(int i = 0; i < Matrix_colunm; i++)eigenValue[i] = eigenValueMatrix[i * Matrix_row + i];

    delete[] T;
    delete[] H;
    delete[] Q_stack;
    delete[] Qt_stack;
    delete[] sub_Q;
    delete[] sub_R;
    delete[] sub_Q_stack;
    delete[] adjustmentMatrix;
    delete[] subMatrix;
    delete[] surplusMatrix;
    delete[] surplus_Q;
    delete[] surplus_R;
    delete[] surplus_Q_stack;
}

bool eig_with_tridiag(double* Matrix, int Matrix_colunm, int Matrix_row, double* eigenValue, double* eigenVector)
{
    //showMatrix(Matrix,Matrix_colunm,Matrix_row,true);
    int iter = 0;
    int maxIter = 1000;

    double epsilon = std::numeric_limits<double>::epsilon();

    double* T = new double[Matrix_colunm * Matrix_row];
    double* H = new double[Matrix_colunm * Matrix_row];
    double* Q = new double[Matrix_colunm * Matrix_row];
    double* R = new double[Matrix_colunm * Matrix_row];

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)eigenVector[i * Matrix_row + j] = (i == j ? 1 : 0);

    //for(int i = 0; i < Matrix_colunm; i++)
    //{
    //    for(int j = 0; j< Matrix_row; j++)
    //    {
    //        T[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
    //    }
    //}

    tridiagnalization_With_HouseHolder(Matrix,Matrix_colunm,Matrix_row,T,H);

    //for(int i = 0; i < Matrix_colunm; i++)
    //{
    //    if(std::abs(T[i * Matrix_row + i]) > epsilon && std::abs(T[i * Matrix_row + i + 1]) > epsilon && i > Matrix_colunm - 1)break;
    //    if(std::abs(T[i * Matrix_row + i]) > epsilon && std::abs(T[i * Matrix_row + i - 1]) > epsilon && i == Matrix_colunm - 1)break;
//
    //    if(i == Matrix_colunm - 1)return false;
    //}
    while (true)
    {
        QRDecompotion_with_SR(T,Matrix_colunm,Matrix_row,Q,R);

        product(R,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_colunm,T);

        product(eigenVector,Matrix_colunm,Matrix_colunm,Q,Matrix_colunm,Matrix_colunm,eigenVector);
        
        if(checkDiagonal(T,Matrix_colunm,Matrix_row,epsilon))break;
        else iter++;
        if(iter > maxIter)break;
    }

    for(int i = 0; i < Matrix_row; i++)eigenValue[i] = T[i * Matrix_row + i];

    product(H,Matrix_colunm,Matrix_row,eigenVector,Matrix_colunm,Matrix_row,eigenVector);

    delete[] T;
    delete[] H;
    delete[] Q;
    delete[] R;

    return true;
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
//    delete[] Matrix_buffer;
//}

double* inverse_lower_triangular_matrix(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    for(int i = 0; i < Matrix_colunm_; i++)
    {
        double* x_vector = new double[Matrix_colunm];
        
        for(int j = 0; j < Matrix_colunm_; j++)x_vector[j] = 0;
        x_vector[i] = 1.0;

        for(int j = i + 1; j < Matrix_row_; j++)
        {
            double tmp = 0;
            for(int k = 0; k < j; k++)
            {
                tmp += Matrix[j * Matrix_row + k] * x_vector[k];
            }
            x_vector[j] -= tmp;
        }
        for(int j = 0; j < Matrix_row_; j++)
        {
            inv_Matrix[j * Matrix_row + i] = x_vector[j];
        }
        delete[] x_vector;
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
        
        for(int j = 0; j < Matrix_colunm_; j++)x_vector[j] = 0;
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
        for(int j = 0; j < Matrix_row_; j++)
        {
            inv_Matrix[j * Matrix_row + i] = x_vector[j];
        }

        delete[] x_vector;
    }
    return inv_Matrix;
}

//LU分解を行う(mode = 1は逆行列用、mode = 0は連立方程式用)
bool LU_Decompose(int mode, double* Matrix, int Matrix_colunm_, int Matrix_row_, double* L, double* U, double* Q)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;

    for(int i = 0; i < Matrix_colunm_; i++)    
    {
        for(int j = 0; j < Matrix_row_; j++)
        {
            U[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
            if(mode == 1)Q[i * Matrix_row + j] = (i == j ? 1.0 : 0.0);;
        }
    }

    for(int i = 0; i < Matrix_colunm_ - 1; i++)
    {
        int holdValueIndex = 0;
        double holdValue = 0.0;
        for(int j = i; j < Matrix_row_; j++)
        {
            if(holdValue < abs(U[i + Matrix_row * j]))
            {
                holdValue = abs(U[i + Matrix_row * j]);
                holdValueIndex = j;
            }
        }

        if(holdValueIndex != i)
        {
            for(int j = 0; j < Matrix_colunm_; j++)std::swap(U[j + Matrix_row * holdValueIndex],U[j + Matrix_row * i]);
            for(int j = 0; j < Matrix_colunm_; j++)std::swap(L[j + Matrix_row * holdValueIndex],L[j + Matrix_row * i]);
            if(mode == 1)for(int j = 0; j < Matrix_colunm_; j++)std::swap(Q[j + Matrix_row * holdValueIndex],Q[j + Matrix_row * i]);
            else if(mode == 0)std::swap(Q[holdValueIndex],Q[i]);
        }

        for(int j = 0; j < i; j++)L[j * Matrix_row + i] = 0.0;
        L[i * Matrix_row + i] = 1.0;

        for(int j = i + 1; j < Matrix_row_; j++)
        {
            double cat = U[j * Matrix_row + i] / U[i * Matrix_row + i];
            L[j * Matrix_row + i] = cat;
            for(int k = 0; k < Matrix_row_; k++)U[j * Matrix_row + k] -= U[i * Matrix_row + k] * cat;
        }

        if(mode == 1)L[Matrix_colunm * Matrix_row - 1] = 1.0;
    }
    return true;
}

double* inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix)
{
    //showMatrix(Matrix,Matrix_colunm_,Matrix_row_);
    //double* Q = new double[Matrix_colunm_ * Matrix_colunm_];
    //double* R = new double[Matrix_colunm_ * Matrix_row_];
    //double* R_inverse = new double[Matrix_row_ * Matrix_colunm_];
//
    //for(int i = 0; i < Matrix_colunm_; i++)for(int j = 0; j < Matrix_row_; j++)R_inverse[i * Matrix_row_ + j] = (i == j ? 1 : 0);
//
//
    //QR_Decompose(Matrix,Matrix_colunm_,Matrix_row_,Q,R);
    //showMatrix(Q,Matrix_colunm_,Matrix_colunm_);
    //showMatrix(R,Matrix_colunm_,Matrix_row_);
//
    //inverse_upper_triangular_matrix(R,Matrix_colunm_,Matrix_row_,R_inverse);
    //
    //showMatrix(R_inverse,Matrix_row_,Matrix_colunm_);
//
    //transpose(Q,Matrix_colunm_,Matrix_colunm_,Q);
//
    //product(R_inverse,Matrix_row_,Matrix_colunm_,Q,Matrix_colunm_,Matrix_colunm_,inv_Matrix);
    //showMatrix(inv_Matrix,Matrix_row_,Matrix_colunm_);
//
    //delete[] Q;
    //delete[] R; 

    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    const double epsilon = 1e-15;

    //2次以下の行列では、LU分解がうまくいかないため
    //if(Matrix_colunm == 1 && Matrix_row == 1)
    //{
    //    inv_Matrix[0] = 1 / Matrix[0];
    //    return inv_Matrix;
    //}
    //else if(Matrix_colunm == 2 && Matrix_row == 2)
    //{
    //    double* tmp = new double[2 * 2];
//
    //    double det = Matrix[0] * Matrix[3] - Matrix[1] * Matrix[2];
    //    std::cout << Matrix[0] << "*" << Matrix[3] << "-" << Matrix[1] << "*" << Matrix[2] << std::endl;
    //    tmp[0] = Matrix[3] / det;
    //    tmp[1] = - Matrix[1] / det;
    //    tmp[2] = - Matrix[2] / det;
    //    tmp[3] = Matrix[0] / det;
//
    //    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)inv_Matrix[i * Matrix_row + j] = tmp[i * Matrix_row + j];
//
    //    delete[] tmp;
    //    return inv_Matrix;
    //}

    double* L = new double[Matrix_colunm * Matrix_row];
    double* U = new double[Matrix_colunm * Matrix_row];
    double* Q = new double[Matrix_colunm * Matrix_row];
    double* inverse_L = new double[Matrix_colunm * Matrix_row];
    double* inverse_U = new double[Matrix_colunm * Matrix_row];
    double* P = new double[Matrix_colunm * Matrix_row];
    
    for(int i = 0; i < Matrix_colunm_; i++)
    {
        for(int j = 0; j < Matrix_row_; j++)
        {
            L[i * Matrix_row + j] = 0;
            U[i * Matrix_row + j] = 0;
            inverse_L[i * Matrix_row + j] = (i == j ? 1 : 0);
            inverse_U[i * Matrix_row + j] = (i == j ? 1 : 0);
        }
    }
//
//    ////pivoting
//
//    //for(int i = 0; i < Matrix_colunm; i++)
//    //{
//    //    int maxindex = 0;
//    //    for(int j = 0; j < Matrix_row; j++)
//    //    {
//    //        if(Matrix[i * Matrix_colunm + j] < Matrix[i * Matrix_colunm + maxindex])continue;
//    //        maxindex = j;
    //    }
    //    
    //    for(int j = 0; j < Matrix_row; j++)Matrix[i * Matrix_row + j] /= Matrix[i * Matrix_row + maxindex];
    //}

    

    //for(int i = 0 ;i < Matrix_colunm; i++)
    //{
    //    L[i * Matrix_row] = Matrix[i * Matrix_row];
    //    U[i] = Matrix[i] / Matrix[0];
    //    U[i * Matrix_row + i] = 1;
    //    //std::cout << "L：" << std::endl; 
    //    //showMatrix(L,3,3);
    //    //std::cout << "U：" << std::endl;
    //    //showMatrix(U,3,3);
    //    //std::cout << std::endl;
    //}
//
    ////LU分解
    //for(int i = 1; i < Matrix_colunm; i++)
    //{
    //    for(int j = i; j < Matrix_colunm; j++)
    //    {
    //        double L_pivot = Matrix[j * Matrix_row + i];
    //        //std::cout << L_pivot << std::endl;
//
    //        for(int k = i - 1; k >= 0; k--)
    //        {
    //            L_pivot -= L[j * Matrix_row + k] * U[k * Matrix_row + i];
    //            //std::cout << L[j * Matrix_row + k] << "x" << U[k * Matrix_row + i] << std::endl;
    //        }
    //        L[j* Matrix_row + i] = L_pivot;
    //        //std::cout << "L[" << j << "," << i << "] = " << L_pivot << std::endl; 
    //        //showMatrix(L,3,3);
    //    }
    //    for(int j = i + 1; j < Matrix_row; j++)
    //    {
    //        double U_pivot = Matrix[i * Matrix_row + j];
    //        for(int k = i - 1; k >= 0; k--)
    //        {
    //            U_pivot -= L[i * Matrix_row + k] * U[k * Matrix_row + j];
    //            //std::cout << L[i * Matrix_row + k] << "x" << U[k * Matrix_row + j] << std::endl;
    //        }
    //        U[i * Matrix_row + j] = U_pivot / L[i * Matrix_row + i];
    //        //std::cout << "U：" << std::endl; 
    //        //showMatrix(U,3,3);
    //        //std::cout << std::endl;
    //    }
    //}
//
    //int maxIndex = Matrix_colunm - 1;
//
    //for(int n = 0; n < Matrix_colunm; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(L[n * Matrix_colunm + m]) > 1e-10 && std::abs(L[m * Matrix_colunm + n]) > 1e-10)break;
    //        if(m == maxIndex)
    //        {
    //            L[n * Matrix_colunm + n] = 1;
    //            break;
    //        }
    //    }
    //}
//
    //for(int n = 0; n < Matrix_colunm; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(U[n * Matrix_colunm + m]) > 1e-10 && std::abs(U[m * Matrix_colunm + n]) > 1e-10)break;
    //        if(m == maxIndex)
    //        {
    //            U[n * Matrix_colunm + n] = 1;
    //            break;
    //        }
    //    }
    //}
//
    ////Lの逆行列計算
    //for(int i = 0; i < Matrix_row; i++)
    //{
    //    for(int j = i; j < Matrix_colunm; j++)
    //    {
    //        if(i == j)
    //        {
    //            for(int k = 0; k <= i; k++)
    //            {
    //                inverse_L[j * Matrix_row + k] /= L[i * Matrix_row + i];
    //            }
    //        }
    //        else
    //        {
    //            for(int k = 0; k <= i; k++)
    //            {
    //                inverse_L[j * Matrix_row + k] -= inverse_L[i * Matrix_row + k] * L[j * Matrix_row + i];
    //            }
    //        }
    //    }
    //}
//
    ////Uの逆行列計算
    //for(int i = Matrix_row - 1; i >= 0; i--)
    //{
    //    for(int j = i - 1; j >= 0; j--)
    //    {
    //        for(int k = Matrix_row - 1; k >= i; k--)
    //        {
    //            inverse_U[j * Matrix_colunm + k] -= inverse_U[i * Matrix_colunm + k] * U[j * Matrix_colunm + i];
    //        }
    //    }
    //}
    //double* bufferMatrix = new double[Matrix_colunm * Matrix_row];
//
    //std::cout << "L_inv" << std::endl;
    //showMatrix(inverse_L,Matrix_colunm,Matrix_row);
//
    //product(L,Matrix_colunm,Matrix_row,inverse_L,Matrix_colunm,Matrix_row,bufferMatrix);
    //std::cout << "L x invL" << std::endl;
    //showMatrix(bufferMatrix,Matrix_colunm,Matrix_row);
//
    //product(U,Matrix_colunm,Matrix_row,inverse_U,Matrix_colunm,Matrix_row,bufferMatrix);
    //std::cout << "U x invU" << std::endl;
    //showMatrix(bufferMatrix,Matrix_colunm,Matrix_row);
//
    //std::cout << "inputMatrix" << std::endl;
    //showMatrix(Matrix,Matrix_colunm,Matrix_row);
    //product(L,Matrix_colunm,Matrix_row,U,Matrix_colunm,Matrix_row,bufferMatrix);
    //std::cout << "outputMatrix" << std::endl;
    //showMatrix(bufferMatrix,Matrix_colunm,Matrix_row);
//
    //std::cout << std::endl;
    //delete[] bufferMatrix;

    //inv_Matrix = product(inverse_U,Matrix_colunm,Matrix_row,inverse_L,Matrix_colunm,Matrix_row,inv_Matrix);


    LU_Decompose(1,Matrix,Matrix_colunm,Matrix_row,L,U,Q);
//
    //for(int i = 0; i < Matrix_colunm; i++)
    //{
    //    if(L[i * Matrix_row + i] == 0)L[i * Matrix_row + i] = epsilon;
    //    if(U[i * Matrix_row + i] == 0)U[i * Matrix_row + i] = epsilon;
    //}


    //showMatrix(L,Matrix_colunm,Matrix_row);
    //showMatrix(U,Matrix_colunm,Matrix_row);
    //int maxIndex = Matrix_colunm - 1;
//
    //for(int n = 0; n < Matrix_colunm; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(L[n * Matrix_colunm + m]) > epsilon && std::abs(L[m * Matrix_colunm + n]) > epsilon)break;
    //        if(m == maxIndex)
    //        {
    //            L[n * Matrix_colunm + n] = 1;
    //            break;
    //        }
    //    }
    //}
//
    //for(int n = 0; n < Matrix_colunm; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(U[n * Matrix_colunm + m]) > epsilon && std::abs(U[m * Matrix_colunm + n]) > epsilon)break;
    //        if(m == maxIndex)
    //        {
    //            U[n * Matrix_colunm + n] = 1;
    //            break;
    //        }
    //    }
    //}

    inverse_L = inverse_lower_triangular_matrix(L,Matrix_colunm,Matrix_row,inverse_L);
    inverse_U = inverse_upper_triangular_matrix(U,Matrix_colunm,Matrix_row,inverse_U);
    //showMatrix(inverse_L,Matrix_colunm,Matrix_row);
    //showMatrix(inverse_U,Matrix_colunm,Matrix_row);
    //showMatrix(Q,Matrix_colunm,Matrix_row);

    //std::cout << std::endl;


    P = product(inverse_U,Matrix_colunm,Matrix_row,inverse_L,Matrix_colunm,Matrix_row,P);
    inv_Matrix = product(P,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_row,inv_Matrix);
    
    //double* test = new double[Matrix_colunm * Matrix_colunm];
//
    //product(inv_Matrix,Matrix_colunm,Matrix_row,Matrix,Matrix_colunm,Matrix_row,test);
    //showMatrix(test,Matrix_colunm,Matrix_row);
//
    //delete[] test;

    delete[] L;
    delete[] U;
    delete[] Q;
    delete[] inverse_L;
    delete[] inverse_U;
    delete[] P;

    return inv_Matrix;
}

//まだRank落ちの正方行列にしか対応していない
double* Pseudo_inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* pinv_Matrix)
{
    //showMatrix(Matrix,Matrix_colunm_,Matrix_row_);

    int eig_size = (Matrix_colunm_ > Matrix_row_ ? Matrix_row_ : Matrix_colunm_);

    double* U = new double[Matrix_colunm_ * Matrix_colunm_];
    double* S = new double[Matrix_colunm_ * Matrix_row_];
    double* V = new double[Matrix_row_ * Matrix_row_];
    double* tmpMatrix = new double[Matrix_row_ * Matrix_colunm_];

    double* S_inverse = new double[Matrix_colunm_ * Matrix_row_];

    for (int i = 0; i < Matrix_colunm_ * Matrix_row_; i++)S_inverse[i] = 0;

    SVD_double(Matrix, Matrix_colunm_, Matrix_row_, U, S, V);

    for (int i = 0; i < eig_size; i++)
    {
        S_inverse[i * Matrix_colunm_ + i] = 1 / S[i * Matrix_row_ + i];
    }

    //showMatrix(U, Matrix_colunm_, Matrix_colunm_);
    //showMatrix(S_inverse, Matrix_row_, Matrix_colunm_);
    //showMatrix(V, Matrix_row_, Matrix_row_);

    cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, Matrix_row_, Matrix_colunm_, Matrix_row_, 1, V, Matrix_row_, S_inverse, Matrix_colunm_, 0, tmpMatrix, Matrix_colunm_);
    cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Matrix_row_, Matrix_colunm_, Matrix_colunm_, 1, tmpMatrix, Matrix_colunm_, U, Matrix_colunm_, 0, pinv_Matrix, Matrix_colunm_);

    delete[] U;
    delete[] S;
    delete[] V;
    delete[] tmpMatrix;
    delete[] S_inverse;

    //const int smaller_size = (Matrix_colunm_ > Matrix_row_ ? Matrix_row_ : Matrix_colunm_);
    //const double epsilon = std::numeric_limits<double>::epsilon();
    //
    //int R_rank = 0;
    //int new_R_colunm = 0;
    //int new_R_row = 0;
    //std::vector<int> noneZero_colunm_Index;
    //std::vector<int> noneZero_row_Index;
    //
    //double* Q = new double[Matrix_colunm_ * Matrix_colunm_];
    //double* R = new double[Matrix_colunm_ * Matrix_row_];

    //QR_Decompose(Matrix,Matrix_colunm_,Matrix_row_,Q,R);

    //showMatrix(Q,Matrix_colunm_,Matrix_colunm_);
    //showMatrix(R,Matrix_colunm_,Matrix_row_);
    //double* test = new double[Matrix_colunm_ * Matrix_row_];
//
    //product(Q,Matrix_colunm_,Matrix_colunm_,R,Matrix_colunm_,Matrix_row_,test);
//
    //showMatrix(test,Matrix_colunm_,Matrix_row_);
//
    //sub(Matrix,test,Matrix_colunm_,Matrix_row_,test);
//
    //std::cout << frobenius_norm(test,Matrix_colunm_,Matrix_row_) << std::endl;

    //delete[] test;
    //for(int i = 0; i < Matrix_row_; i++)
    //{
    //    int zero_Counter = 0; 
    //    for(int j = 0; j <= (i > Matrix_colunm_ - 1 ? Matrix_colunm_ - 1 : i); j++)
    //    {
    //        //std::cout << R[j * Matrix_row_ + i] << std::endl;
    //        if(std::abs(R[j * Matrix_row_ + i]) < 1e-12)zero_Counter++;
    //    }
    //    //std::cout << (i > Matrix_colunm_ - 1 ? Matrix_colunm_ : i + 1) << ":" << zero_Counter << std::endl;
    //    if(zero_Counter != (i > Matrix_colunm_ - 1 ? Matrix_colunm_ : i + 1))noneZero_row_Index.push_back(i);
    //}
    //
    //for(int i = 0; i < Matrix_colunm_; i++)
    //{
    //    int zero_Counter = 0;
    //    for(int j = i; j < Matrix_row_; j++)
    //    {
    //        //std::cout << R[i * Matrix_row_ + j] << std::endl;
    //        if(std::abs(R[i * Matrix_row_ + j]) < 1e-12)zero_Counter++;
    //    }
    //    if(zero_Counter != (Matrix_row_ - i))noneZero_colunm_Index.push_back(i);
    //}
    //
    //new_R_colunm = noneZero_colunm_Index.size();
    //new_R_row = noneZero_row_Index.size();
    //R_rank = std::min(noneZero_colunm_Index.size(), noneZero_row_Index.size());
    //
    //if(R_rank == Matrix_colunm_ && R_rank == Matrix_row_)
    //{
    //    delete[] Q;
    //    delete[] R;
    //
    //    return inverse(Matrix,Matrix_colunm_,Matrix_row_,pinv_Matrix);
    //}
    //else if (R_rank == Matrix_colunm_ || R_rank == Matrix_row_)
    //{   
    //    int new_Q_row = new_R_colunm;
    //
    //    double* new_Q = new double[Matrix_colunm_ * new_Q_row];
    //    double* new_R = new double[new_R_colunm * Matrix_row_];
    //
    //    int n = 0;
    //    for(int i : noneZero_colunm_Index)
    //    {
    //        for(int j = 0; j < Matrix_colunm_; j++)
    //        {
    //            new_Q[j * new_Q_row + n] = Q[j * Matrix_colunm_ + i];
    //        }
    //        n++;
    //    }
    //
    //    n = 0;
    //    for(auto i : noneZero_colunm_Index)
    //    {
    //        for(int j = 0; j < Matrix_row_; j++)
    //        {
    //            new_R[n * Matrix_row_ + j] = R[i * Matrix_row_ + j];
    //        }
    //        n++;
    //    }
    //
    //    double* new_Rt = new double[Matrix_row_ * new_R_colunm];
    //    double* new_RRt = new double[new_R_colunm * new_R_colunm];
    //
    //    double* new_Qt = transpose(new_Q,Matrix_colunm_,new_Q_row,new_Q); 
    //    new_Rt = transpose(new_R,new_R_colunm,Matrix_row_,new_Rt);
    //    //showMatrix(new_R,new_R_colunm,Matrix_row_);
    //
    //    new_RRt = product(new_R,new_R_colunm,Matrix_row_,new_Rt,Matrix_row_,new_R_colunm,new_RRt);
    //    //showMatrix(new_RRt,new_R_colunm,new_R_colunm);
    //    double* inv_RRt = inverse(new_RRt,new_R_colunm,new_R_colunm,new_RRt);
    //    //showMatrix(new_RRt,new_R_colunm,new_R_colunm);
    //
    //    double* Rt_inv_RRt = product(new_Rt,Matrix_row_,new_R_colunm,inv_RRt,new_R_colunm,new_R_colunm,new_Rt);
    //    //showMatrix(Rt_inv_RRt,Matrix_row_,new_R_colunm);
    //    product(Rt_inv_RRt,Matrix_row_,new_R_colunm,new_Qt,new_Q_row,Matrix_colunm_,pinv_Matrix);
    //
    //    delete[] new_Q;
    //    delete[] new_R;
    //    delete[] new_Rt;
    //    delete[] new_RRt;
    //}
    //else
    //{
    //    double* U = new double[Matrix_colunm_ * Matrix_colunm_];
    //    double* S = new double[Matrix_colunm_ * Matrix_row_];
    //    double* V = new double[Matrix_row_ * Matrix_row_];
    //
    //    //showMatrix(Matrix,Matrix_colunm_,Matrix_row_,true);
    //    SVD_double(Matrix,Matrix_colunm_,Matrix_row_,U,S,V);
    //
    //    //showMatrix(U,Matrix_colunm_,Matrix_colunm_);
    //    //showMatrix(S,Matrix_colunm_,Matrix_row_);
    //    //showMatrix(V,Matrix_row_,Matrix_row_);
    //
    //    transpose(U,Matrix_colunm_,Matrix_colunm_,U);
    //    transpose(S,Matrix_colunm_,Matrix_row_,S);
    //    transpose(V,Matrix_row_,Matrix_row_,V);
    //
    //    for(int i = 0; i < smaller_size; i++)S[i * Matrix_row_ + i] = (std::abs(S[i * Matrix_row_ + i]) > epsilon ? 1 / S[i * Matrix_row_ + i] : 0);
    //
    //    product(V,Matrix_row_,Matrix_row_,S,Matrix_row_,Matrix_colunm_,S);
    //    product(S,Matrix_row_,Matrix_colunm_,U,Matrix_colunm_,Matrix_colunm_,pinv_Matrix); 
    //
    //    delete[] U;
    //    delete[] S;
    //    delete[] V;
    //}
    //
    //delete[] Q;
    //delete[] R;
    //
    //double* fullRank_R = new double[new_R_colunm * new_R_row];
    ////for(int i : noneZero_row_Index)std::cout << i << " ";
//  //
    //n = 0;
    //for(int i : noneZero_row_Index)
    //{
    //    for(int j = 0; j < new_R_colunm; j++)
    //    {
    //        //std::cout << "(" << j << "," << i << ")" << std::endl;
    //        //std::cout << "(" << j << "," << n << ")" << std::endl;
    //        fullRank_R[j * new_R_row + n] = new_R[j * Matrix_row_ + i];
    //        //showMatrix(fullRank_R,new_R_colunm,new_R_row);
    //    }
    //    n++;
    //}
//
    //showMatrix(new_Q,Matrix_colunm_,new_Q_row);
    ////showMatrix(new_R,new_R_colunm,Matrix_row_);
//
    ////showMatrix(fullRank_R,new_R_colunm,new_R_row);
//
    //if(R_rank == Matrix_colunm_ || R_rank == Matrix_row_)
    //{
//
    //    double* fullRank_Rt = new double[new_R_row * new_R_colunm];
    //    fullRank_Rt = transpose(fullRank_R,new_R_colunm,new_R_row,fullRank_Rt);
//
    //    double* RRt = new double[new_R_colunm * new_R_colunm];
    //    RRt = product(fullRank_R,new_R_colunm,new_R_row,fullRank_Rt,new_R_row,new_R_colunm,RRt);
//
    //    double* RRt_inv = inverse(RRt,new_R_colunm,new_R_colunm,RRt);
    //    //showMatrix(RRt_inv,new_R_colunm,new_R_colunm);
//
    //    double* Rt = transpose(new_R,new_R_colunm,Matrix_row_,new_R);
    //    //showMatrix(Rt,Matrix_row_,new_R_colunm);
    //    double* Qt = transpose(new_Q,Matrix_colunm_,new_Q_row,new_Q);
    //    //showMatrix(Qt,new_Q_row,Matrix_colunm_);
//
    //    double* calBuffer = new double[Matrix_row_ * new_R_colunm];
    //    double* Rt_RRt_inv = product(Rt,Matrix_row_,new_R_colunm,RRt_inv,new_R_colunm,new_R_colunm,calBuffer);
//
    //    pinv_Matrix = product(Rt_RRt_inv,Matrix_row_,new_R_colunm,Qt,new_Q_row,Matrix_colunm_,pinv_Matrix);
//
    //    delete[] fullRank_Rt;
    //    delete[] RRt;
    //    delete[] calBuffer;
    //}
    //else
    //{
    //    //DEBUGFLAG = true;
    //    double* fullRank_R_t = new double[new_R_row * new_R_colunm];
    //    fullRank_R_t = transpose(fullRank_R,new_R_colunm,new_R_row,fullRank_R_t);
//
    //    double* fullRank_RtR = new double[new_R_row * new_R_row];
    //    fullRank_RtR = product(fullRank_R_t,new_R_row,new_R_colunm,fullRank_R,new_R_colunm,new_R_row,fullRank_RtR);
    //    //showMatrix(fullRank_RtR,new_R_row,new_R_row);
//
    //    //DEBUGFLAG = false;
    //    double* RtR_inv = inverse(fullRank_RtR,new_R_row,new_R_row,fullRank_RtR);
    //    showMatrix(RtR_inv,new_R_row,new_R_row);
//
    //    double* calBuffer = new double[Matrix_row_ * new_R_colunm];
    //    double* RMRt = new double[new_R_colunm * new_R_colunm];
//
    //    //DEBUGFLAG = true;
    //    double* M = product(RtR_inv,new_R_row,new_R_row,RtR_inv,new_R_row,new_R_row,fullRank_RtR);
    //    //showMatrix(M,new_R_row,new_R_row);
//
    //    double* MRt = product(M,new_R_row,new_R_row,fullRank_R_t,new_R_row,new_R_colunm,fullRank_R_t);
    //    RMRt = product(fullRank_R,new_R_colunm,new_R_row,MRt,new_R_row,new_R_colunm,RMRt);
    //    //showMatrix(RMRt,new_R_colunm,new_R_colunm);
//
    //    double* Rt = transpose(new_R,new_R_colunm,Matrix_row_,new_R);
    //    //showMatrix(Rt,Matrix_row_,new_R_colunm);
    //    double* Qt = transpose(new_Q,Matrix_colunm_,new_Q_row,new_Q);
    //    //showMatrix(Qt,new_Q_row,Matrix_colunm_);
//
    //    double* Rt_A = product(Rt,Matrix_row_,new_R_colunm,RMRt,new_R_colunm,new_R_colunm,calBuffer);
    //    //showMatrix(calBuffer,Matrix_row_,new_R_colunm);
    //    double* Rt_A_Qt = product(Rt_A,Matrix_row_,new_R_colunm,Qt,new_Q_row,Matrix_colunm_,pinv_Matrix);
    //    //showMatrix(Rt_A_Qt,Matrix_row_,Matrix_colunm_);
    //    //DEBUGFLAG = false;
//
    //    delete[] fullRank_R_t;
    //    delete[] fullRank_RtR;
    //    delete[] calBuffer;
    //    delete[] RMRt;
    //}
    

    //if(Matrix_colunm_ > Matrix_row_)
    //{
    //    double* fullRank_R = new double[Matrix_colunm_ * R_rank];
//
    //    int fullRankR_index = 0;
    //    for(int i : noneZeroIndex)
    //    {
    //        for(int j = 0; j < Matrix_row_; j++)
    //        {
    //            fullRank_R[j * R_rank + fullRankR_index] = R[j * Matrix_row_ + i];
    //        }
    //        fullRankR_index++;
    //    }
//
    //    //showMatrix(fullRank_R,Matrix_row_,R_rank);
//
    //    double* fullRank_R_t = new double[R_rank * Matrix_row_];
    //    double* fullRank_RtR = new double[R_rank * R_rank];
//
    //    fullRank_R_t = transpose(fullRank_R,Matrix_row_,R_rank,fullRank_R_t);
    //    //showMatrix(fullRank_R_t,R_rank,Matrix_row_);
    //    fullRank_RtR = product(fullRank_R_t,R_rank,Matrix_row_,fullRank_R,Matrix_row_,R_rank,fullRank_RtR);
    //    //showMatrix(fullRank_RtR,R_rank,R_rank);
//
    //    double* RtR_inv = inverse(fullRank_RtR,R_rank,R_rank,fullRank_RtR);
    //    //showMatrix(RtR_inv,R_rank,R_rank);
//
    //    double* M = product(RtR_inv,R_rank,R_rank,RtR_inv,R_rank,R_rank,fullRank_RtR);
    //    double* MRt = product(M,R_rank,R_rank,fullRank_R_t,R_rank,Matrix_row_,fullRank_R_t);
    //    double* RMRt = product(fullRank_R,Matrix_row_,R_rank,MRt,R_rank,Matrix_row_,pinv_Matrix);
//
    //    //showMatrix(RMRt,Matrix_row_,Matrix_row_);
//
    //    double* Rt = transpose(R,Matrix_colunm_,Matrix_row_,R);
//
    //    double* Rt_A = product(Rt,Matrix_row_,Matrix_colunm_,RMRt,Matrix_colunm_,Matrix_row_,pinv_Matrix);
    //    double* Rt_A_Qt = product(Rt_A,Matrix_row_,Matrix_row_,Qt,Matrix_colunm_,Matrix_colunm_,pinv_Matrix);
//
    //    delete[] fullRank_R;
    //    delete[] fullRank_R_t;
    //    delete[] fullRank_RtR;
    //}
    //else
    //{
    //    double* fullRank_R = new double[(Matrix_colunm_ - (Matrix_row_ - R_rank)) * ];
//
    //    int fullRankR_index = 0;
    //    for(int i : noneZeroIndex)
    //    {
    //        for(int j = 0; j < Matrix_colunm_; j++)
    //        {
    //            fullRank_R[j * Matrix_row_ + fullRankR_index] = R[j * Matrix_row_ + i];
    //        }
    //        fullRankR_index++;
    //    }
//
    //    showMatrix(fullRank_R,Matrix_colunm_,R_rank);
//
    //    double* fullRank_R_t = new double[R_rank * Matrix_colunm_];
    //    double* fullRank_RtR = new double[R_rank * R_rank];
//
    //    fullRank_R_t = transpose(fullRank_R,Matrix_colunm_,R_rank,fullRank_R_t);
    //    showMatrix(fullRank_R_t,R_rank,Matrix_colunm_);
    //    fullRank_RtR = product(fullRank_R_t,R_rank,Matrix_colunm_,fullRank_R,Matrix_colunm_,R_rank,fullRank_RtR);
    //    showMatrix(fullRank_RtR,R_rank,R_rank);
//
    //    double* RtR_inv = inverse(fullRank_RtR,R_rank,R_rank,fullRank_RtR);
    //    double* M = product(RtR_inv,R_rank,R_rank,RtR_inv,R_rank,R_rank,fullRank_RtR);
    //    double* MRt = product(M,R_rank,R_rank,fullRank_R_t,R_rank,Matrix_colunm_,fullRank_R_t);
    //    double* RMRt = product(fullRank_R,Matrix_colunm_,R_rank,MRt,R_rank,Matrix_colunm_,pinv_Matrix);
//
    //    double* Rt = transpose(R,Matrix_colunm_,Matrix_row_,R);
//
    //    double* Rt_A = product(Rt,Matrix_row_,Matrix_colunm_,RMRt,Matrix_colunm_,Matrix_row_,pinv_Matrix);
    //    double* Rt_A_Qt = product(Rt_A,Matrix_row_,Matrix_row_,Qt,Matrix_colunm_,Matrix_colunm_,pinv_Matrix);
    //    std::cout << "a" << std::endl;
//
    //    delete[] fullRank_R;
    //    delete[] fullRank_R_t;
    //    delete[] fullRank_RtR;
    //}

    //std::cout << "Q : " << Q << std::endl;
    //delete[] Q;
    //std::cout << "R : " << R << std::endl;
    //delete[] R;
    //std::cout << "new_Q : " << new_Q << std::endl;
    //delete[] new_Q;
    //std::cout << "new_R : " << new_R << std::endl;
    //delete[] new_R;
    //std::cout << "fullRank_R : " << fullRank_R << std::endl;
    //delete[] fullRank_R;
    ////std::cout << "fullRank_R_t : " << fullRank_R_t << std::endl;
    //
    ////std::cout << "fullRank_RtR : " << fullRank_RtR << std::endl;
    //
    ////std::cout << "calBuffer : " << calBuffer << std::endl;
    //
    ////std::cout << "RMRt : " << RMRt << std::endl;
    //

    return pinv_Matrix;
}

double* equal_solve(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* vector, int vector_size_, double* solution_vector)
{
    double* piv_Matrix = new double[Matrix_row_ * Matrix_colunm_];

    Pseudo_inverse(Matrix, Matrix_colunm_, Matrix_row_, piv_Matrix);

    product(piv_Matrix, Matrix_row_, Matrix_colunm_, vector, vector_size_, 1, solution_vector);

    delete[] piv_Matrix;

    //std::cout << "equal solve.." << std::endl;
    //const unsigned int Matrix_colunm = Matrix_colunm_;
    //const unsigned int Matrix_row    = Matrix_row_;
    //const unsigned int vector_size = vector_size_;
    //
    //double* pinv_Matrix = new double[Matrix_row * Matrix_row];
//
    //Pseudo_inverse(Matrix,Matrix_colunm,Matrix_row,pinv_Matrix);
    //product(pinv_Matrix,Matrix_row,Matrix_colunm,vector,vector_size,1,solution_vector);

  // // double* L = new double[Matrix_colunm * Matrix_row];
   // double* U = new double[Matrix_colunm * Matrix_row];
   // double* y = new double[Matrix_colunm];

  // // for(int i = 0; i < Matrix_colunm_ * Matrix_row_; i++)
   // {
   //     L[i] = 0;
   //     U[i] = 0;
   // }

  // // LU_Decompose(0,Matrix,Matrix_colunm,Matrix_row,L,U,vector);
   // y[0] = vector[0];

  // // for(int i = 1; i < vector_size; i++)
   // {
   //     y[i] = vector[i];
   //     for(int j = 0; j < i; j++)y[i] -= L[i * Matrix_row + j] * y[j];
   // }

  // // solution_vector[vector_size - 1] = y[vector_size - 1] / U[Matrix_colunm * Matrix_row - 1];

  // // for(int i = vector_size - 2; i >= 0; i--)
   // {
   //     solution_vector[i] = y[i];
   //     for(int j = i + 1; j < vector_size; j++)solution_vector[i] -= U[i * Matrix_row + j] * solution_vector[j]; 
   //     solution_vector[i] /= U[i * Matrix_row + i];
   // }

  // // delete[] L;
   // delete[] U;
   // delete[] y;
    //delete[] pinv_Matrix;

    //std::cout << "equal solve end..." << std::endl;
    return solution_vector;
}

void SVD_float(float* Matrix, int Matrix_colunm, int Matrix_row, float* U, float* S, float* V)
{
    float* tmp_Matrix = new float[Matrix_colunm * Matrix_row];
    float* tmp_S = new float[Matrix_colunm * Matrix_row];

    int target_size = (Matrix_colunm > Matrix_row ? Matrix_row : Matrix_colunm);

    for (int i = 0; i < Matrix_colunm; i++)
    {
        for (int j = 0; j < Matrix_row; j++)
        {
            tmp_Matrix[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
            S[i * Matrix_row + j] = 0;
        }
    }

    LAPACKE_sgesdd(LAPACK_COL_MAJOR, 'A', Matrix_row, Matrix_colunm, tmp_Matrix, Matrix_row, tmp_S, V, Matrix_row, U, Matrix_colunm);

    for (int i = 0; i < target_size; i++)
    {
        S[i * Matrix_row + i] = tmp_S[i];
    }

    delete[] tmp_Matrix;
    delete[] tmp_S;
}

void SVD_double(double* Matrix, int Matrix_colunm, int Matrix_row, double* U, double* S, double* V)
{
    //std::cout << "start SVD..." << std::endl;

    double* tmp_Matrix = new double[Matrix_colunm * Matrix_row];
    double* tmp_S = new double[Matrix_colunm * Matrix_row];

    int target_size = (Matrix_colunm > Matrix_row ? Matrix_row : Matrix_colunm);

    for (int i = 0; i < Matrix_colunm; i++)
    {
        for (int j = 0; j < Matrix_row; j++)
        {
            tmp_Matrix[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
            S[i * Matrix_row + j] = 0;
        }
    }

    LAPACKE_dgesdd(LAPACK_COL_MAJOR, 'A', Matrix_row, Matrix_colunm, tmp_Matrix, Matrix_row, tmp_S, V, Matrix_row, U, Matrix_colunm);

    for (int i = 0; i < target_size; i++)
    {
        S[i * Matrix_row + i] = tmp_S[i];
    }

    delete[] tmp_Matrix;
    delete[] tmp_S;

    //const double epsilon = std::numeric_limits<double>::epsilon();

    //std::vector<int> MMt_index;
    //std::vector<int> MtM_index;
    //std::vector<int> MMt_erase_index;
    //std::vector<int> MtM_erase_index;

    //double* Mt = new double[Matrix_row * Matrix_colunm];
    //double* MMt = new double[Matrix_colunm * Matrix_colunm];
    //double* MtM = new double[Matrix_row * Matrix_row];
    //double* eigenValue_U = new double[Matrix_colunm];
    //double* eigenValue_V = new double[Matrix_row];
    //
    //for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_colunm; j++)U[i * Matrix_colunm + j] = 0;
    //for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)S[i * Matrix_row + j] = (i == j ? 1 : 0);
    //for(int i = 0; i < Matrix_row; i++)for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = 0;
    //
    //for(int i = 0; i < Matrix_colunm; i++)eigenValue_U[i] = 0;
    //for(int i = 0; i < Matrix_row; i++)eigenValue_V[i] = 0;
    //
    //transpose(Matrix,Matrix_colunm,Matrix_row,Mt);
    //
    //product(Matrix,Matrix_colunm,Matrix_row,Mt,Matrix_row,Matrix_colunm,MMt);
    //product(Mt,Matrix_row,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,MtM);
    //
    //eig(MMt,Matrix_colunm,Matrix_colunm,eigenValue_U,U);
    //eig(MtM,Matrix_row,Matrix_row,eigenValue_V,V);

    //showMatrix(eigenValue_U,1,Matrix_colunm);
    //showMatrix(eigenValue_V,1,Matrix_row);

    //for(int i = 0; i < Matrix_colunm; i++)MMt_index.push_back(i);
    //for(int i = 0; i < Matrix_row; i++)MtM_index.push_back(i);
//
    //int maxIndex = Matrix_colunm - 1;
    //for(int n = 0; n < Matrix_colunm; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(MMt[n * Matrix_colunm + m]) > epsilon && std::abs(MMt[m * Matrix_colunm + n]) > epsilon)break;
    //        if(m == maxIndex)
    //        {
    //            MMt_index.erase(std::remove(MMt_index.begin(), MMt_index.end(), n));
    //            MMt_erase_index.push_back(n);
    //            break;
    //        }
    //    }
    //}
//
    //maxIndex = Matrix_row - 1;
    //for(int n = 0; n < Matrix_row; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(MtM[n * Matrix_row + m]) > epsilon && std::abs(MtM[m * Matrix_row + n]) > epsilon)break;
    //        if(m == maxIndex)
    //        {
    //            MtM_index.erase(std::remove(MtM_index.begin(), MtM_index.end(), n));
    //            MtM_erase_index.push_back(n);
    //            break;
    //        }
    //    }
    //}
//
    //if(MMt_index.size() == Matrix_colunm)
    //{
    //    eig_with_tridiag(MMt,Matrix_colunm,Matrix_colunm,eigenValue_U,U);
    //}
    //else
    //{
    //    double* tmpMatrix = new double[MMt_index.size() * MMt_index.size()];
    //    double* tmpEigenValue = new double[MMt_index.size()];
    //    double* tmpEigenVector = new double[MMt_index.size() * MMt_index.size()];
    //    
    //    int n = 0;
//
    //    n = 0;
    //    for(int col : MMt_index)
    //    {
    //        int m = 0;
    //        for(int row : MMt_index)
    //        {
    //            tmpMatrix[n * MMt_index.size() + m] = MMt[col * Matrix_colunm + row];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    eig_with_tridiag(tmpMatrix,MMt_index.size(),MMt_index.size(),tmpEigenValue,tmpEigenVector);
    // 
    //    n = 0;
    //    for(int col : MMt_index)
    //    {   std::cout << col << std::endl;
    //        eigenValue_U[col] = tmpEigenValue[n];
    //        int m = 0;
    //        for(int row : MMt_index)
    //        {
    //            U[col * Matrix_colunm + row] = tmpEigenVector[n * MMt_index.size() + m];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    delete[] tmpMatrix;
    //    delete[] tmpEigenValue;
    //    delete[] tmpEigenVector;
    //}
//
    //if(MtM_index.size() == Matrix_row)
    //{
    //    eig_with_tridiag(MtM,Matrix_row,Matrix_row,eigenValue_V,V);
    //}
    //else
    //{
    //    double* tmpMatrix = new double[MtM_index.size() * MtM_index.size()];
    //    double* tmpEigenValue = new double[MtM_index.size()];
    //    double* tmpEigenVector = new double[MtM_index.size() * MtM_index.size()];
    //    
    //    int n = 0;
//
    //    n = 0;
    //    for(int col : MtM_index)
    //    {
    //        int m = 0;
    //        for(int row : MtM_index)
    //        {
    //            tmpMatrix[n * MtM_index.size() + m] = MtM[col * Matrix_row + row];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    eig_with_tridiag(tmpMatrix,MtM_index.size(),MtM_index.size(),tmpEigenValue,tmpEigenVector);
//
    //    n = 0;
    //    for(int col : MtM_index)
    //    {
    //        std::cout << col << std::endl;
    //        eigenValue_V[col] = tmpEigenValue[n];
    //        int m = 0;
    //        for(int row : MtM_index)
    //        {
    //            V[col * Matrix_row + row] = tmpEigenVector[n * MtM_index.size() + m];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    delete[] tmpMatrix;
    //    delete[] tmpEigenValue;
    //    delete[] tmpEigenVector;
    //}

    //if(Matrix_colunm > Matrix_row)
    //{
    //    int smaller_Size = Matrix_row;
    //    int larger_Size = Matrix_colunm;
    //
    //    double* eigValues = eigenValue_V;
    //    double* another_eigValues = eigenValue_U;
    //
    //    double* eigVector = V;
    //    double* another_eigVector = U;
    //
    //    //固有値の大きい順に並べ、それに従ってU,Vの固有値ベクトルも並び替える。
    //    std::vector<int> smaller_size_index;
    //    std::vector<int> larger_size_index;
    //
    //    for(int i = 0; i < smaller_Size; i++)
    //    {
    //        smaller_size_index.push_back(i);
    //    }
    //
    //    for(int i = 0; i < larger_Size; i++)
    //    {
    //        larger_size_index.push_back(i);
    //    }
    //
    //    std::sort(smaller_size_index.begin(), smaller_size_index.end(), [eigValues](int l_idx, int r_idx){return eigValues[l_idx] > eigValues[r_idx];});
    //    std::vector<std::pair<int,int>> smaller_size_transposition = get_transposition(smaller_size_index);
    //
    //    for(std::pair<int,int> index_pair : smaller_size_transposition)
    //    {
    //        std::swap(eigValues[index_pair.first], eigValues[index_pair.second]);
    //        for(int i = 0; i < smaller_Size; i++)std::swap(eigVector[index_pair.first * smaller_Size + i],eigVector[index_pair.second * smaller_Size + i]);
    //    }
    //
    //    for(int i = 0; i < smaller_Size; i++)
    //    {
    //        double eigvalue = eigValues[i];
    //        S[i * Matrix_row + i] = (eigvalue > epsilon ? std::sqrt(eigvalue) : 0);
    //    }
    //
    //    for(int i : smaller_size_index)
    //    {
    //        double tmpepsilon = 1e-10;
    //        int index = larger_Size;
    //        for(int j = 0; j < larger_Size; j++)
    //        {
    //            double diff = std::abs(eigValues[i] - another_eigValues[j]);
    //            if(diff < tmpepsilon)
    //            {
    //                index = j;
    //                tmpepsilon = diff;
    //            }
    //        }
    //
    //        if(index != larger_Size)
    //        {
    //            larger_size_index.erase(std::remove(larger_size_index.begin(), larger_size_index.end(), index));
    //            larger_size_index.insert(larger_size_index.begin(), index);
    //        }
    //
    //    }
    //
    //    std::sort(larger_size_index.begin(), larger_size_index.begin() + smaller_Size, [another_eigValues](int l_idx, int r_idx){return another_eigValues[l_idx] > another_eigValues[r_idx];});
    //    std::sort(larger_size_index.begin() + smaller_Size, larger_size_index.end(), [another_eigValues](int l_idx, int r_idx){return std::abs(another_eigValues[l_idx]) > std::abs(another_eigValues[r_idx]);});
    //
    //    std::vector<std::pair<int,int>> larger_size_transposition = get_transposition(larger_size_index);
    //
    //    for(std::pair<int,int> index_pair : larger_size_transposition)
    //    {
    //        std::swap(another_eigValues[index_pair.first], another_eigValues[index_pair.second]);
    //        for(int i = 0; i < larger_Size; i++)std::swap(another_eigVector[i * larger_Size + index_pair.first],another_eigVector[i * larger_Size + index_pair.second]);
    //    }
    //}
    //else
    //{
    //    int smaller_Size = Matrix_colunm;
    //    int larger_Size = Matrix_row;
    //
    //    double* eigValues = eigenValue_U;
    //    double* another_eigValues = eigenValue_V;
    //
    //    double* eigVector = U;
    //    double* another_eigVector = V;
    //
    //    //固有値の大きい順に並べ、それに従ってU,Vの固有値ベクトルも並び替える。
    //    std::vector<int> smaller_size_index;
    //    std::vector<int> larger_size_index;
    //
    //    for(int i = 0; i < smaller_Size; i++)
    //    {
    //        smaller_size_index.push_back(i);
    //    }
    //
    //    for(int i = 0; i < larger_Size; i++)
    //    {
    //        larger_size_index.push_back(i);
    //    }
    //
    //    std::sort(smaller_size_index.begin(), smaller_size_index.end(), [eigValues](int l_idx, int r_idx){return eigValues[l_idx] > eigValues[r_idx];});
    //    std::vector<std::pair<int,int>> smaller_size_transposition = get_transposition(smaller_size_index);
    //
    //    for(std::pair<int,int> index_pair : smaller_size_transposition)
    //    {
    //        std::swap(eigValues[index_pair.first], eigValues[index_pair.second]);
    //        for(int i = 0; i < smaller_Size; i++)std::swap(eigVector[i * smaller_Size + index_pair.first],eigVector[i * smaller_Size + index_pair.second]);
    //    }
    //
    //    for(int i = 0; i < smaller_Size; i++)
    //    {
    //        double eigvalue = eigValues[i];
    //        S[i * Matrix_row + i] = (eigvalue > 1e-10 ? std::sqrt(eigvalue) : 0);
    //    }
    //
    //    int loopCounter = 0;
    //    for(int i : smaller_size_index)
    //    {
    //        double tmpepsilon = 1e-10;
    //        int index = larger_Size;
    //        for(int j = 0; j < larger_Size; j++)
    //        {
    //            double diff = std::abs(eigValues[i] - another_eigValues[j]);
    //            if(diff < tmpepsilon)
    //            {
    //                index = j;
    //                tmpepsilon = diff;
    //            }
    //        }
    //
    //        auto result = std::find(larger_size_index.begin(),larger_size_index.end(),index);
    //        std::swap(larger_size_index[loopCounter],larger_size_index[*result]);
    //
    //        loopCounter++;
    //    }
    //
    //    std::sort(larger_size_index.begin(), larger_size_index.begin() + smaller_Size, [another_eigValues](int l_idx, int r_idx){return another_eigValues[l_idx] > another_eigValues[r_idx];});
    //    //std::cout << *(larger_size_index.begin() + smaller_Size) << std::endl;
    //    //std::sort(larger_size_index.begin() + smaller_Size, larger_size_index.end(), [another_eigValues,larger_size_index](int l_idx, int r_idx){std::cout << "{";for(int i : larger_size_index)std::cout << i << ",";std::cout << "}" << std::endl;return std::abs(another_eigValues[l_idx]) > std::abs(another_eigValues[r_idx]);});
    //
    //
    //    std::vector<std::pair<int,int>> larger_size_transposition = get_transposition(larger_size_index);
    //
    //    for(std::pair<int,int> index_pair : larger_size_transposition)
    //    {
    //        std::swap(another_eigValues[index_pair.first], another_eigValues[index_pair.second]);
    //        for(int i = 0; i < larger_Size; i++)std::swap(another_eigVector[index_pair.first * larger_Size + i],another_eigVector[index_pair.second * larger_Size + i]);
    //    }
    //}
    //
    
    
    //if(MMt_erase_index.size() != 0)
    //{
    //    for(int i : MMt_erase_index)U[i * Matrix_colunm + i] = 1;
    //}
//  
    //if(MtM_erase_index.size() != 0)
    //{
    //    for(int i : MtM_erase_index)V[i * Matrix_row + i] = 1;
    //}
    //showMatrix(V,Matrix_row,Matrix_row,true);
    //
    //transpose(V,Matrix_row,Matrix_row,V);
    //
    //for(int k = 0; k < Matrix_row; k++)
    //{
    //    double A_v_multi = 0;
    //    double sigma_u_multi = U[k] * S[k * Matrix_row + k];
    //    for(int i = 0; i < Matrix_row; i++)
    //    {
    //        A_v_multi += Matrix[i] * V[k * Matrix_row + i];
    //    }
    //
    //    if(std::abs(A_v_multi - sigma_u_multi) < 1e-10)continue;
    //
    //    for(int i = 0; i < Matrix_row; i++)V[k * Matrix_row + i] *= -1;
    //}

    //double* tmp = new double[Matrix_colunm * Matrix_row];
//
    //product(U,Matrix_colunm,Matrix_colunm,S,Matrix_colunm,Matrix_row,tmp);
    //product(tmp,Matrix_colunm,Matrix_row,V,Matrix_row,Matrix_row,tmp);
//
    //showMatrix(tmp,Matrix_colunm,Matrix_row,true);

    //delete[] tmp;

    //delete[] Mt;
    //delete[] MMt;
    //delete[] MtM;
    //delete[] eigenValue_U;
    //delete[] eigenValue_V;
    //
    //showMatrix(Matrix,Matrix_colunm,Matrix_row);

    //const double epsilon = std::numeric_limits<double>::epsilon();

    //int smaller_Size = std::min(Matrix_colunm,Matrix_row);
    //int larger_Size = std::max(Matrix_colunm,Matrix_row);
    //int diffrence_of_col_and_row = std::abs(Matrix_colunm - Matrix_row);
    //bool eigValueSolution_flag = (Matrix_colunm > Matrix_row ? true : false);
//
    //double* Matrix_t = new double[Matrix_colunm * Matrix_row];
    //double* MMt = new double[Matrix_colunm * Matrix_colunm];
    //double* MtM = new double[Matrix_row * Matrix_row];
    //double* u_eigValues = new double[Matrix_colunm];
    //double* v_eigValues = new double[Matrix_row];
    //double* tmp_for_U = new double[Matrix_colunm * Matrix_colunm];
    //double* tmp_for_V = new double[Matrix_row * Matrix_row];
    //double* u_vector = new double[Matrix_colunm];
    //double* v_vector = new double[Matrix_row];
    //bool* isExist = new bool[larger_Size];
//
    //for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_colunm; j++)U[i * Matrix_colunm + j] = 0;
    //for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)S[i * Matrix_row + j] = 0;
    //for(int i = 0; i < Matrix_row; i++)for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = 0;
    //for(int i = 0; i < larger_Size; i++)isExist[i] = false;
//
    //Matrix_t = transpose(Matrix,Matrix_colunm,Matrix_row,Matrix_t);
//
    ////showMatrix(Matrix_t,Matrix_row,Matrix_colunm);
//
    ////DEBUGFLAG = true;
    //MMt = product(Matrix,Matrix_colunm,Matrix_row,Matrix_t,Matrix_row,Matrix_colunm,MMt);
    //MtM = product(Matrix_t,Matrix_row,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,MtM);
    ////DEBUGFLAG = false;
//
    ////showMatrix(MMt,Matrix_colunm,Matrix_colunm);
    ////showMatrix(MtM,Matrix_row,Matrix_row);
    //
    //u_eigValues = eig(MMt,Matrix_colunm,Matrix_colunm,u_eigValues);
    //v_eigValues = eig(MtM,Matrix_row,Matrix_row,v_eigValues);
//
    ////showMatrix(u_eigValues,1,Matrix_colunm);
    ////showMatrix(v_eigValues,1,Matrix_row);
//
    //double* eigValues = (Matrix_colunm > Matrix_row ? v_eigValues : u_eigValues);
    //double* another_eigValues = (Matrix_colunm > Matrix_row ? u_eigValues : v_eigValues);
//
    ////固有値の大きい順に並べ、それに従ってU,Vの固有値ベクトルも並び替える。
    //std::vector<int> smaller_size_index;
    //std::vector<int> larger_size_index;
//
    //for(int i = 0; i < smaller_Size; i++)
    //{
    //    smaller_size_index.push_back(i);
    //}
//
    //for(int i = 0; i < larger_Size; i++)
    //{
    //    larger_size_index.push_back(i);
    //}
//
    //std::sort(smaller_size_index.begin(), smaller_size_index.end(), [eigValues](int l_idx, int r_idx){return eigValues[l_idx] > eigValues[r_idx];});
    //std::vector<std::pair<int,int>> smaller_size_transposition = get_transposition(smaller_size_index);
//
    //for(std::pair<int,int> index_pair : smaller_size_transposition)
    //{
    //    std::swap(eigValues[index_pair.first], eigValues[index_pair.second]);
    //}
//
    //for(int i : smaller_size_index)
    //{
    //    double epsilon = 1e-10;
    //    int index = larger_Size;
    //    for(int j = 0; j < larger_Size; j++)
    //    {
    //        double diff = std::abs(eigValues[i] - another_eigValues[j]);
    //        if(diff < epsilon)
    //        {
    //            index = j;
    //            epsilon = diff;
    //        }
    //    }
//
    //    if(index != larger_Size)
    //    {
    //        larger_size_index.erase(std::remove(larger_size_index.begin(), larger_size_index.end(), index));
    //        larger_size_index.insert(larger_size_index.begin(), index);
    //        //std::cout << index << std::endl;
    //    }
    //}
//
    //std::sort(larger_size_index.begin(), larger_size_index.begin() + smaller_Size, [another_eigValues](int l_idx, int r_idx){return another_eigValues[l_idx] > another_eigValues[r_idx];});
    //std::sort(larger_size_index.begin() + smaller_Size, larger_size_index.end(), [another_eigValues](int l_idx, int r_idx){return std::abs(another_eigValues[l_idx]) > std::abs(another_eigValues[r_idx]);});
    ////std::cout << another_eigValues[*larger_size_index.begin()] << another_eigValues[*(larger_size_index.begin() + smaller_Size)] << std::endl;
    ////std::cout << another_eigValues[*(larger_size_index.begin() + smaller_Size)] << another_eigValues[*larger_size_index.end()] << std::endl;
    //std::vector<std::pair<int,int>> larger_size_transposition = get_transposition(larger_size_index);
//
    //for(std::pair<int,int> index_pair : larger_size_transposition)
    //{
    //    std::swap(another_eigValues[index_pair.first], another_eigValues[index_pair.second]);
    //}
//
//
    //for(int i = 0; i < smaller_Size; i++)
    //{
    //    double eigValue = (std::abs(eigValues[i]) > 1e-10 ? eigValues[i] : 0);
//
    //    if(eigValue < 0)eigValue = 0;
    //    S[i * Matrix_row + i] = std::sqrt(eigValue);
//
    //    int iter = 0;
//
    //    for(int n = 0; n < Matrix_colunm; n++)for(int m = 0; m < Matrix_colunm; m++)tmp_for_U[n * Matrix_colunm + m] = (n == m ? eigValue : 0);
    //    for(int n = 0; n < Matrix_row; n++)for(int m = 0; m < Matrix_row; m++)tmp_for_V[n * Matrix_row + m] = (n == m ? eigValue : 0);
//
//
    //    double* sub_tmp_for_U = sub(MMt,tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
    //    double* sub_tmp_for_V = sub(MtM,tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
//
    //    int maxIndex = Matrix_colunm - 1;
//
    //    for(int n = 0; n < Matrix_colunm; n++)
    //    {
    //        for(int m = 0;;m++)
    //        {
    //            if(std::abs(sub_tmp_for_U[n * Matrix_colunm + m]) > 1e-10 && std::abs(sub_tmp_for_U[m * Matrix_colunm + n]) > 1e-10)break;
    //            if(m == maxIndex)
    //            {
    //                sub_tmp_for_U[n * Matrix_colunm + n] = 1;
    //                break;
    //            }
    //        }
    //    }
//
    //    maxIndex = Matrix_row - 1;
//
    //    for(int n = 0; n < Matrix_row; n++)
    //    {
    //        for(int m = 0;;m++)
    //        {
    //            if(std::abs(sub_tmp_for_V[n * Matrix_row + m]) > 1e-10 && std::abs(sub_tmp_for_V[m * Matrix_row + n]) > 1e-10)break;
    //            if(m == maxIndex)
    //            {
    //                sub_tmp_for_V[n * Matrix_row + n] = 1;
    //                break;
    //            }
    //        }
    //    }
//
    //    showMatrix(sub_tmp_for_U,Matrix_colunm,Matrix_colunm);
    //    showMatrix(sub_tmp_for_V,Matrix_row,Matrix_row);
//
    //    double* inv_tmp_for_U = inverse(sub_tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
    //    double* inv_tmp_for_V = inverse(sub_tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
//
    //    showMatrix(sub_tmp_for_U,Matrix_colunm,Matrix_colunm);
    //    showMatrix(sub_tmp_for_V,Matrix_row,Matrix_row);
    //    for(int j = 0; j < Matrix_colunm; j++)u_vector[j] = 1;
    //    for(int j = 0; j < Matrix_row; j++)v_vector[j] = 1;
//
    //    while (iter < maxIter)
    //    {
    //        u_vector = product(inv_tmp_for_U,Matrix_colunm,Matrix_colunm,u_vector,Matrix_colunm,1,u_vector);
    //        double norm_u = norm(u_vector,Matrix_colunm);
    //        for(int n = 0; n < Matrix_colunm; n++)
    //        {
    //            if(norm_u != 0)u_vector[n] /= norm_u;
    //            else u_vector[n] = 1;
    //        }
    //        //showMatrix(inv_tmp_for_U,Matrix_colunm,Matrix_colunm);
    //        //showMatrix(u_vector,Matrix_colunm,1);
//
    //        v_vector = product(inv_tmp_for_V,Matrix_row,Matrix_row,v_vector,Matrix_row,1,v_vector);
    //        double norm_v = norm(v_vector,Matrix_row);
    //        for(int n = 0; n < Matrix_row; n++)
    //        {
    //            if(norm_v != 0)v_vector[n] /= norm_v;
    //            else v_vector[n] = 1;
    //        }
    //        //showMatrix(inv_tmp_for_V,Matrix_row,Matrix_row);
    //        //showMatrix(v_vector,Matrix_row,1);
//
    //        iter++;
    //    }
    //    for(int j = 0; j < Matrix_colunm; j++)U[j * Matrix_colunm + i] = u_vector[j];
    //    for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = v_vector[j];
    //    //showMatrix(U,Matrix_colunm,Matrix_colunm);
    //    //showMatrix(S,Matrix_colunm,Matrix_row);
    //    //showMatrix(V,Matrix_row,Matrix_row);
    //    //std::cout << std::endl;
    //}
//
    ////showMatrix(U,Matrix_colunm,Matrix_colunm);
    ////showMatrix(S,Matrix_colunm,Matrix_row);
    ////showMatrix(V,Matrix_row,Matrix_row);
//
    //for(int i = smaller_Size, index = 0; i < larger_Size; i++)
    //{
    //    double eigValue = another_eigValues[i];
    //    if(Matrix_colunm == larger_Size)
    //    {
    //        int iter_for_U = 0;
    //        for(int n = 0; n < Matrix_colunm; n++)for(int m = 0; m < Matrix_colunm; m++)tmp_for_U[n * Matrix_colunm + m] = (n == m ? eigValue : 0);
    //        double* sub_tmp_for_U = sub(MMt,tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
    //        int maxIndex = Matrix_colunm - 1;
    //        for(int n = 0; n < Matrix_colunm; n++)
    //        {
    //            for(int m = 0;;m++)
    //            {
    //                if(std::abs(sub_tmp_for_U[n * Matrix_colunm + m]) > 1e-10 && std::abs(sub_tmp_for_U[m * Matrix_colunm + n]) > 1e-10)break;
    //                if(m == maxIndex)
    //                {
    //                    sub_tmp_for_U[n * Matrix_colunm + n] = 1;
    //                    break;
    //                }
    //            }
    //        }
    //        double* inv_tmp_for_U = inverse(sub_tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
    //        for(int j = 0; j < Matrix_colunm; j++)u_vector[j] = 1;
    //        while (iter_for_U < maxIter)
    //        {
    //            u_vector = product(inv_tmp_for_U,Matrix_colunm,Matrix_colunm,u_vector,Matrix_colunm,1,u_vector);
    //            double norm_u = norm(u_vector,Matrix_colunm);
    //            for(int n = 0; n < Matrix_colunm; n++)u_vector[n] /= norm_u;
    //            iter_for_U++;
    //        }
    //        for(int j = 0; j < Matrix_colunm; j++)
    //        {
    //            U[j * Matrix_colunm + (index + smaller_Size)] = u_vector[j];
    //        }
    //    }
//
    //    if(Matrix_row == larger_Size)
    //    {
    //        int iter_for_V = 0;
    //        for(int n = 0; n < Matrix_row; n++)for(int m = 0; m < Matrix_row; m++)tmp_for_V[n * Matrix_row + m] = (n == m ? eigValue : 0);
    //        double* sub_tmp_for_V = sub(MtM,tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
    //        int maxIndex = Matrix_row - 1;
//
    //        for(int n = 0; n < Matrix_row; n++)
    //        {
    //            for(int m = 0;;m++)
    //            {
    //                if(std::abs(sub_tmp_for_V[n * Matrix_row + m]) > 1e-10 && std::abs(sub_tmp_for_V[m * Matrix_row + n]) > 1e-10)break;
    //                if(m == maxIndex)
    //                {
    //                    sub_tmp_for_V[n * Matrix_row + n] = 1;
    //                    break;
    //                }
    //            }
    //        }
    //        double* inv_tmp_for_V = inverse(sub_tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
    //        for(int j = 0; j < Matrix_row; j++)v_vector[j] = 1;
    //        while (iter_for_V < maxIter)
    //        {
    //            v_vector = product(inv_tmp_for_V,Matrix_row,Matrix_row,v_vector,Matrix_row,1,v_vector);
    //            double norm_v = norm(v_vector,Matrix_row);
    //            for(int n = 0; n < Matrix_row; n++)v_vector[n] /= norm_v;
    //            iter_for_V++;
    //        }
    //        for(int j = 0; j < Matrix_row; j++)
    //        {
    //            V[(index + smaller_Size) * Matrix_row + j] = v_vector[j];
    //        }
    //    }
    //    index++;
    //}
//
    ////showMatrix(U,Matrix_colunm,Matrix_colunm);
    ////showMatrix(V,Matrix_row,Matrix_row);
//
    ////固有値ベクトルの向きを合わせる(もう少し高速化できる)
    ////double* Ut = new double[Matrix_colunm * Matrix_colunm];
    ////double* Vt = new double[Matrix_row * Matrix_row];
////
    ////transpose(U,Matrix_colunm,Matrix_colunm,Ut);
    ////product(Ut,Matrix_colunm,Matrix_colunm,MMt,Matrix_colunm,Matrix_colunm,Ut);
    ////product(Ut,Matrix_colunm,Matrix_colunm,U,Matrix_colunm,Matrix_colunm,Ut);
////
    ////showMatrix(Ut,Matrix_colunm,Matrix_colunm);
////
    ////if(std::abs(Ut[0]) > 1e-10)for(int j = 0; j < Matrix_colunm; j++)U[j * Matrix_colunm] *= -1;
////
    ////for(int i = 1; i < Matrix_colunm; i++)
    ////{
    ////    if(std::abs(Ut[i]) < 1e-10)continue;
    ////    for(int j = 0; j < Matrix_colunm; j++)U[j * Matrix_colunm + i] *= -1;
    ////}
////
    ////transpose(V,Matrix_row,Matrix_row,Vt);
    ////product(Vt,Matrix_row,Matrix_row,MtM,Matrix_row,Matrix_row,Vt);
    ////product(Vt,Matrix_row,Matrix_row,V,Matrix_row,Matrix_row,Vt);
////
    ////showMatrix(Vt,Matrix_row,Matrix_row);
////
    ////if(std::abs(Vt[1]) > 1e-10)for(int j = 0; j < Matrix_row; j++)V[j] *= -1;
////
    ////for(int i = 1; i < Matrix_row; i++)
    ////{
    ////    if(std::abs(Vt[i]) < 1e-10)continue;
    ////    for(int j = 0; j < Matrix_row; j++)U[i * Matrix_colunm + j] *= -1;
    ////}
////
    ////delete[] Ut;
    ////delete[] Vt;
//
    ////固有値ベクトルの向きを合わせる
    //for(int k = 0; k < Matrix_row; k++)
    //{
    //    double A_v_multi = 0;
    //    double sigma_u_multi = U[k] * S[k * Matrix_row + k];
    //    for(int i = 0; i < Matrix_row; i++)
    //    {
    //        A_v_multi += Matrix[i] * V[k * Matrix_row + i];
    //    }
//
    //    if(std::abs(A_v_multi - sigma_u_multi) < 1e-10)continue;
//
    //    for(int i = 0; i < Matrix_row; i++)V[k * Matrix_row + i] *= -1;
    //}
//
    //delete[] Matrix_t;
    //delete[] MMt;
    //delete[] MtM;
    //delete[] u_eigValues;
    //delete[] v_eigValues;
    //delete[] tmp_for_U;
    //delete[] tmp_for_V;
    //delete[] u_vector;
    //delete[] v_vector;
    //delete[] isExist;

    //std::cout << "SVD end..." << std::endl;
}

void svd_with_tridiag(double* Matrix, int Matrix_colunm, int Matrix_row, double* U, double* S, double* V)
{
    const double epsilon = std::numeric_limits<double>::epsilon();

    //std::vector<int> MMt_index;
    //std::vector<int> MtM_index;
    //std::vector<int> MMt_erase_index;
    //std::vector<int> MtM_erase_index;

    double* Mt = new double[Matrix_row * Matrix_colunm];
    double* MMt = new double[Matrix_colunm * Matrix_colunm];
    double* MtM = new double[Matrix_row * Matrix_row];
    double* eigenValue_U = new double[Matrix_colunm];
    double* eigenValue_V = new double[Matrix_row];

    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_colunm; j++)U[i * Matrix_colunm + j] = 0;
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)S[i * Matrix_row + j] = (i == j ? 1 : 0);
    for(int i = 0; i < Matrix_row; i++)for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = 0;

    for(int i = 0; i < Matrix_colunm; i++)eigenValue_U[i] = 0;
    for(int i = 0; i < Matrix_row; i++)eigenValue_V[i] = 0;

    transpose(Matrix,Matrix_colunm,Matrix_row,Mt);

    product(Matrix,Matrix_colunm,Matrix_row,Mt,Matrix_row,Matrix_colunm,MMt);
    product(Mt,Matrix_row,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,MtM);

    eig_with_divide(MMt,Matrix_colunm,Matrix_colunm,eigenValue_U,U);
    eig_with_divide(MtM,Matrix_row,Matrix_row,eigenValue_V,V);

    //showMatrix(eigenValue_U,1,Matrix_colunm);
    //showMatrix(eigenValue_V,1,Matrix_row);

    //for(int i = 0; i < Matrix_colunm; i++)MMt_index.push_back(i);
    //for(int i = 0; i < Matrix_row; i++)MtM_index.push_back(i);
//
    //int maxIndex = Matrix_colunm - 1;
    //for(int n = 0; n < Matrix_colunm; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(MMt[n * Matrix_colunm + m]) > epsilon && std::abs(MMt[m * Matrix_colunm + n]) > epsilon)break;
    //        if(m == maxIndex)
    //        {
    //            MMt_index.erase(std::remove(MMt_index.begin(), MMt_index.end(), n));
    //            MMt_erase_index.push_back(n);
    //            break;
    //        }
    //    }
    //}
//
    //maxIndex = Matrix_row - 1;
    //for(int n = 0; n < Matrix_row; n++)
    //{
    //    for(int m = 0;;m++)
    //    {
    //        if(std::abs(MtM[n * Matrix_row + m]) > epsilon && std::abs(MtM[m * Matrix_row + n]) > epsilon)break;
    //        if(m == maxIndex)
    //        {
    //            MtM_index.erase(std::remove(MtM_index.begin(), MtM_index.end(), n));
    //            MtM_erase_index.push_back(n);
    //            break;
    //        }
    //    }
    //}
//
    //if(MMt_index.size() == Matrix_colunm)
    //{
    //    eig_with_tridiag(MMt,Matrix_colunm,Matrix_colunm,eigenValue_U,U);
    //}
    //else
    //{
    //    double* tmpMatrix = new double[MMt_index.size() * MMt_index.size()];
    //    double* tmpEigenValue = new double[MMt_index.size()];
    //    double* tmpEigenVector = new double[MMt_index.size() * MMt_index.size()];
    //    
    //    int n = 0;
//
    //    n = 0;
    //    for(int col : MMt_index)
    //    {
    //        int m = 0;
    //        for(int row : MMt_index)
    //        {
    //            tmpMatrix[n * MMt_index.size() + m] = MMt[col * Matrix_colunm + row];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    eig_with_tridiag(tmpMatrix,MMt_index.size(),MMt_index.size(),tmpEigenValue,tmpEigenVector);
    // 
    //    n = 0;
    //    for(int col : MMt_index)
    //    {   std::cout << col << std::endl;
    //        eigenValue_U[col] = tmpEigenValue[n];
    //        int m = 0;
    //        for(int row : MMt_index)
    //        {
    //            U[col * Matrix_colunm + row] = tmpEigenVector[n * MMt_index.size() + m];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    delete[] tmpMatrix;
    //    delete[] tmpEigenValue;
    //    delete[] tmpEigenVector;
    //}
//
    //if(MtM_index.size() == Matrix_row)
    //{
    //    eig_with_tridiag(MtM,Matrix_row,Matrix_row,eigenValue_V,V);
    //}
    //else
    //{
    //    double* tmpMatrix = new double[MtM_index.size() * MtM_index.size()];
    //    double* tmpEigenValue = new double[MtM_index.size()];
    //    double* tmpEigenVector = new double[MtM_index.size() * MtM_index.size()];
    //    
    //    int n = 0;
//
    //    n = 0;
    //    for(int col : MtM_index)
    //    {
    //        int m = 0;
    //        for(int row : MtM_index)
    //        {
    //            tmpMatrix[n * MtM_index.size() + m] = MtM[col * Matrix_row + row];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    eig_with_tridiag(tmpMatrix,MtM_index.size(),MtM_index.size(),tmpEigenValue,tmpEigenVector);
//
    //    n = 0;
    //    for(int col : MtM_index)
    //    {
    //        std::cout << col << std::endl;
    //        eigenValue_V[col] = tmpEigenValue[n];
    //        int m = 0;
    //        for(int row : MtM_index)
    //        {
    //            V[col * Matrix_row + row] = tmpEigenVector[n * MtM_index.size() + m];
    //            m++;
    //        }
    //        n++;
    //    }
//
    //    delete[] tmpMatrix;
    //    delete[] tmpEigenValue;
    //    delete[] tmpEigenVector;
    //}

    if(Matrix_colunm > Matrix_row)
    {
        int smaller_Size = Matrix_row;
        int larger_Size = Matrix_colunm;

        double* eigValues = eigenValue_V;
        double* another_eigValues = eigenValue_U;

        double* eigVector = V;
        double* another_eigVector = U;

        //固有値の大きい順に並べ、それに従ってU,Vの固有値ベクトルも並び替える。
        std::vector<int> smaller_size_index;
        std::vector<int> larger_size_index;

        for(int i = 0; i < smaller_Size; i++)
        {
            smaller_size_index.push_back(i);
        }

        for(int i = 0; i < larger_Size; i++)
        {
            larger_size_index.push_back(i);
        }

        std::sort(smaller_size_index.begin(), smaller_size_index.end(), [eigValues](int l_idx, int r_idx){return eigValues[l_idx] > eigValues[r_idx];});
        std::vector<std::pair<int,int>> smaller_size_transposition = get_transposition(smaller_size_index);

        for(std::pair<int,int> index_pair : smaller_size_transposition)
        {
            std::swap(eigValues[index_pair.first], eigValues[index_pair.second]);
            for(int i = 0; i < smaller_Size; i++)std::swap(eigVector[index_pair.first * smaller_Size + i],eigVector[index_pair.second * smaller_Size + i]);
        }

        for(int i = 0; i < smaller_Size; i++)
        {
            double eigvalue = eigValues[i];
            S[i * Matrix_row + i] = (eigvalue > epsilon ? std::sqrt(eigvalue) : 0);
        }

        for(int i : smaller_size_index)
        {
            double tmpepsilon = 1e-10;
            int index = larger_Size;
            for(int j = 0; j < larger_Size; j++)
            {
                double diff = std::abs(eigValues[i] - another_eigValues[j]);
                if(diff < tmpepsilon)
                {
                    index = j;
                    tmpepsilon = diff;
                }
            }

            if(index != larger_Size)
            {
                larger_size_index.erase(std::remove(larger_size_index.begin(), larger_size_index.end(), index));
                larger_size_index.insert(larger_size_index.begin(), index);
            }

        }

        std::sort(larger_size_index.begin(), larger_size_index.begin() + smaller_Size, [another_eigValues](int l_idx, int r_idx){return another_eigValues[l_idx] > another_eigValues[r_idx];});
        std::sort(larger_size_index.begin() + smaller_Size, larger_size_index.end(), [another_eigValues](int l_idx, int r_idx){return std::abs(another_eigValues[l_idx]) > std::abs(another_eigValues[r_idx]);});

        std::vector<std::pair<int,int>> larger_size_transposition = get_transposition(larger_size_index);

        for(std::pair<int,int> index_pair : larger_size_transposition)
        {
            std::swap(another_eigValues[index_pair.first], another_eigValues[index_pair.second]);
            for(int i = 0; i < larger_Size; i++)std::swap(another_eigVector[i * larger_Size + index_pair.first],another_eigVector[i * larger_Size + index_pair.second]);
        }
    }
    else
    {
        int smaller_Size = Matrix_colunm;
        int larger_Size = Matrix_row;

        double* eigValues = eigenValue_U;
        double* another_eigValues = eigenValue_V;

        double* eigVector = U;
        double* another_eigVector = V;

        //固有値の大きい順に並べ、それに従ってU,Vの固有値ベクトルも並び替える。
        std::vector<int> smaller_size_index;
        std::vector<int> larger_size_index;

        for(int i = 0; i < smaller_Size; i++)
        {
            smaller_size_index.push_back(i);
        }

        for(int i = 0; i < larger_Size; i++)
        {
            larger_size_index.push_back(i);
        }

        std::sort(smaller_size_index.begin(), smaller_size_index.end(), [eigValues](int l_idx, int r_idx){return eigValues[l_idx] > eigValues[r_idx];});
        std::vector<std::pair<int,int>> smaller_size_transposition = get_transposition(smaller_size_index);

        for(std::pair<int,int> index_pair : smaller_size_transposition)
        {
            std::swap(eigValues[index_pair.first], eigValues[index_pair.second]);
            for(int i = 0; i < smaller_Size; i++)std::swap(eigVector[i * smaller_Size + index_pair.first],eigVector[i * smaller_Size + index_pair.second]);
        }

        for(int i = 0; i < smaller_Size; i++)
        {
            double eigvalue = eigValues[i];
            S[i * Matrix_row + i] = (eigvalue > 1e-10 ? std::sqrt(eigvalue) : 0);
        }

        //int loopCounter = 0;
        for(int i : smaller_size_index)
        {
            double tmpepsilon = 1e-10;
            int index = larger_Size;
            for(int j = 0; j < larger_Size; j++)
            {
                double diff = std::abs(eigValues[i] - another_eigValues[j]);
                if(diff < tmpepsilon)
                {
                    index = j;
                    tmpepsilon = diff;
                }
            }

            if(index != larger_Size)
            {
                larger_size_index.erase(std::remove(larger_size_index.begin(), larger_size_index.end(), index));
                larger_size_index.insert(larger_size_index.begin(), index);
            }

            //auto result = std::find(larger_size_index.begin(),larger_size_index.end(),index);
            //std::swap(larger_size_index[loopCounter],larger_size_index[*result]);

            //loopCounter++;
        }

        std::sort(larger_size_index.begin(), larger_size_index.begin() + smaller_Size, [another_eigValues](int l_idx, int r_idx){return another_eigValues[l_idx] > another_eigValues[r_idx];});
        std::sort(larger_size_index.begin() + smaller_Size, larger_size_index.end(), [another_eigValues,larger_size_index](int l_idx, int r_idx){return std::abs(another_eigValues[l_idx]) > std::abs(another_eigValues[r_idx]);});


        std::vector<std::pair<int,int>> larger_size_transposition = get_transposition(larger_size_index);

        for(std::pair<int,int> index_pair : larger_size_transposition)
        {
            std::swap(another_eigValues[index_pair.first], another_eigValues[index_pair.second]);
            for(int i = 0; i < larger_Size; i++)std::swap(another_eigVector[index_pair.first * larger_Size + i],another_eigVector[index_pair.second * larger_Size + i]);
        }
    }



    //if(MMt_erase_index.size() != 0)
    //{
    //    for(int i : MMt_erase_index)U[i * Matrix_colunm + i] = 1;
    //}
//
    //if(MtM_erase_index.size() != 0)
    //{
    //    for(int i : MtM_erase_index)V[i * Matrix_row + i] = 1;
    //}
    //showMatrix(V,Matrix_row,Matrix_row,true);
    
    transpose(V,Matrix_row,Matrix_row,V);

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

    //double* tmp = new double[Matrix_colunm * Matrix_row];
//
    //product(U,Matrix_colunm,Matrix_colunm,S,Matrix_colunm,Matrix_row,tmp);
    //product(tmp,Matrix_colunm,Matrix_row,V,Matrix_row,Matrix_row,tmp);
//
    //showMatrix(tmp,Matrix_colunm,Matrix_row,true);

    //delete[] tmp;

    delete[] Mt;
    delete[] MMt;
    delete[] MtM;
    delete[] eigenValue_U;
    delete[] eigenValue_V;
}

void NMF(double* Matrix, int Matrix_colunm, int Matrix_row, int n_components, double* W, double* H, double epsilon)
{
    std::mt19937 mt{ std::random_device{}() };

    double* MH_T = new double[Matrix_colunm * n_components];
    double* WH = new double[Matrix_colunm * Matrix_row];
    double* WHH_T = new double[Matrix_colunm * n_components];
    double* W_TM = new double[Matrix_colunm * n_components];
    double* W_TWH = new double[Matrix_colunm * n_components];

    double maxValue = std::numeric_limits<double>::min();
    double minValue = std::numeric_limits<double>::max();

    for (int i = 0; i < Matrix_colunm; i++)
    {
        for (int j = 0; j < Matrix_row; j++)
        {
            maxValue = std::max(Matrix[i * Matrix_row + j], maxValue);
            minValue = std::min(Matrix[i * Matrix_row + j], minValue);
        }
    }

    std::uniform_real_distribution<double> random(minValue, maxValue);

    for (int i = 0; i < Matrix_colunm; i++)
    {
        for (int j = 0; j < n_components; j++)
        {
            W[i * n_components + j] = random(mt);
        }
    }

    for (int i = 0; i < n_components; i++)
    {
        for (int j = 0; j < Matrix_row; j++)
        {
            H[i * Matrix_row + j] = random(mt);
        }
    }

    showMatrix(Matrix, Matrix_colunm, Matrix_row);
    showMatrix(W, Matrix_colunm, n_components);
    showMatrix(H, n_components, Matrix_row);


    while (true)
    {
        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n_components, Matrix_row, Matrix_colunm, 1, W, n_components, Matrix, Matrix_row, 0, W_TM, Matrix_row);

        //showMatrix(W_TM, n_components, Matrix_row);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Matrix_colunm, Matrix_row, n_components, 1, W, n_components, H, Matrix_row, 0, WH, Matrix_row);

        //showMatrix(WH, Matrix_colunm, Matrix_row);

        cblas_dgemm(CblasRowMajor, CblasTrans, CblasNoTrans, n_components, Matrix_row, Matrix_colunm, 1, W, n_components, WH, Matrix_row, 0, W_TWH, Matrix_row);

        //showMatrix(W_TWH, Matrix_row, n_components);

        for (int i = 0; i < n_components; i++)
        {
            for (int j = 0; j < Matrix_row; j++)
            {
                H[i * Matrix_row + j] = H[i * Matrix_row + j] * W_TM[i * Matrix_row + j] / W_TWH[i * Matrix_row + j];
            }
        }

        //showMatrix(H, Matrix_row, n_components);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Matrix_colunm, n_components, Matrix_row, 1, Matrix, Matrix_row, H, Matrix_row, 0, MH_T, n_components);

        //showMatrix(MH_T, Matrix_colunm, n_components);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Matrix_colunm, Matrix_row, n_components, 1, W, n_components, H, Matrix_row, 0, WH, Matrix_row);

        //showMatrix(WH, Matrix_colunm, Matrix_row);

        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasTrans, Matrix_colunm, n_components, Matrix_row, 1, WH, Matrix_row, H, Matrix_row, 0, WHH_T, n_components);

        //showMatrix(WHH_T, Matrix_colunm, n_components);

        for (int i = 0; i < Matrix_colunm; i++)
        {
            for (int j = 0; j < n_components; j++)
            {
                W[i * n_components + j] = W[i * n_components + j] * MH_T[i * n_components + j] / WHH_T[i * n_components + j];
            }
        }

        //showMatrix(W, Matrix_colunm, n_components);

        double residues = 0;
        cblas_dgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, Matrix_colunm, Matrix_row, n_components, 1, W, n_components, H, Matrix_row, 0, WH, Matrix_row);

        for (int i = 0; i < Matrix_colunm; i++)for (int j = 0; j < Matrix_row; j++)residues += (WH[i * Matrix_row + j] - Matrix[i * Matrix_row + j]) * (WH[i * Matrix_row + j] - Matrix[i * Matrix_row + j]);

        std::cout << residues << std::endl;

        if (residues < epsilon)break;
    }


    delete[] MH_T;
    delete[] WH;
    delete[] WHH_T;
    delete[] W_TM;
    delete[] W_TWH;
}

int Matrix_rank(double* Matrix, int Matrix_colunm, int Matrix_row, bool DoSVD, double* U, double* S, double* V)
{
    int rank = 0;
    if(DoSVD)
    {
        if(U == nullptr || S == nullptr || V == nullptr)return -1;
        SVD_double(Matrix,Matrix_colunm,Matrix_row,U,S,V);
    }

    for(int i = 0; i < (Matrix_colunm > Matrix_row ? Matrix_row : Matrix_colunm); i++)if((DoSVD ? S[i * Matrix_row + i] : Matrix[i * Matrix_row + i]) != 0)rank += 1;

    return rank;
}

int eig_rank(double* eig_value, int eig_value_size)
{
    int rank = 0;
    for(int i = 0; i < eig_value_size; i++)if(std::abs(eig_value[i]) > 1e-10)rank++;
    return rank;
}

std::vector<std::pair<int,int>> get_transposition(std::vector<int> index)
{
    std::vector<bool> Waspassed(index.size(),false);
    std::vector<std::pair<int,int>> transposition;

    for(int i = 0; i < index.size(); i++)
    {
        if(Waspassed[i])continue;
        int start_index = i;
        int next_index = i;

        while(start_index != index[next_index]){
            Waspassed[next_index] = true;
            transposition.push_back({next_index,index[next_index]});
            next_index = index[next_index];
        }
        Waspassed[next_index] = true;
    }

    return transposition;
}

void SETDEBUGFLAG(bool flag)
{
    DEBUGFLAG = flag;
}

double* copyMatrix(double* Matrix, int Matrix_colunm, int Matrix_row, double* copiedMatrix)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)copiedMatrix[i * Matrix_row + j] = Matrix[i * Matrix_row + j];
    return copiedMatrix;
}
