#include "linear_algebra.hpp"

//Debug用
void showMatrix(double* Matrix, int colunm, int row, bool allShow)
{
    if((colunm <= 10 && row <= 10)|| allShow)
    {
        std::cout << "{" << std::endl;
        for(int i = 0; i < colunm; i++)
        {
            std::cout << "{";
            for(int j = 0; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
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
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "...";
            for(int j = row - 3; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
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
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
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
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
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
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "...";
            for(int j = row - 3; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
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
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "...";
            for(int j = row - 3; j < row; j++)
            {
                if(std::abs(Matrix[i * row + j]) < std::numeric_limits<double>::epsilon())std::cout << "0.0000000" << ",";
                else std::cout << std::scientific << std::setprecision(3) << Matrix[i * row + j] << ",";
            }
            std::cout << "}," << std::endl;
        }
        std::cout << "}" << std::endl;
    }

}

//Debug用
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

//ノルム計算
double norm(double* vector, int vector_size)
{
    double norm = 0;

    for(int i = 0; i < vector_size; i++)
    {
        norm += vector[i] * vector[i];
    }

    return std::sqrt(norm);
}

//フベルニウスノルム計算
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

//行列式を求める(QR分解を用いてないのでとても遅いし、メモリを食う、改善の余地あり)
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

//転置行列を求める
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

//行列同士の和
double* add(double* Matrix_A, double* Matrix_B, int Matrix_colunm, int Matrix_row, double* Matrix)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)Matrix[i * Matrix_row + j] = Matrix_A[i * Matrix_row + j] + Matrix_B[i * Matrix_row + j];
    return Matrix;
}

//行列同士の差
double* sub(double* Matrix_A, double* Matrix_B, int Matrix_colunm, int Matrix_row, double* Matrix)
{
    for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)Matrix[i * Matrix_row + j] = Matrix_A[i * Matrix_row + j] - Matrix_B[i * Matrix_row + j];
    return Matrix;
}

//内積計算
double dot(double* vector_A, double* vector_B, int vectorSize)
{
    double sum = 0;
    for(int i = 0; i < vectorSize; i++)sum += vector_A[i] * vector_B[i];
    return sum;
}

//行列積計算(改善の余地あり、手計算と同じ方法なのでとても遅い)
double* product(double* Matrix_A, int Matrix_A_colunm_,  int Matrix_A_row_, double* Matrix_B, int Matrix_B_colunm_, int Matrix_B_row_, double* Matrix)
{
    if(Matrix_A_row_ != Matrix_B_colunm_)return nullptr;

    double* tmpMatrix = new double[Matrix_A_colunm_ * Matrix_B_row_];
    for(int i = 0; i < Matrix_A_colunm_; i++)
    {
        for(int j = 0; j < Matrix_B_row_; j++)
        {
            tmpMatrix[i * Matrix_B_row_ + j] = 0;
        }
    }


    for(int i = 0; i < Matrix_A_colunm_; i++)
    {
        for(int j = 0; j < Matrix_B_row_; j++)
        {
            for(int k = 0; k < Matrix_A_row_; k++)
            {
                tmpMatrix[i * Matrix_B_row_ + j] += Matrix_A[i * Matrix_A_row_ + k] * Matrix_B[k * Matrix_B_row_ + j]; 
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
            
            G_Matrix[i * Matrix_colunm + i] = c;
            G_Matrix[j * Matrix_colunm + j] = c;
            G_Matrix[i * Matrix_colunm + j] = s;
            G_Matrix[j * Matrix_colunm + i] = -s;

            R_Matrix = product(G_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix,Matrix_colunm,Matrix_row,R_Matrix);
            double* G_Matrix_t = transpose(G_Matrix,Matrix_colunm,Matrix_colunm,G_Matrix);
            Q_Matrix= product(Q_Matrix,Matrix_colunm,Matrix_colunm,G_Matrix_t,Matrix_colunm,Matrix_colunm,Q_Matrix);

        }
    }
    
    delete[] G_Matrix;

    return true;
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

        R_Matrix = product(H_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix,Matrix_colunm,Matrix_colunm,R_Matrix);
        Q_Matrix = product(Q_Matrix,Matrix_colunm,Matrix_colunm,H_Matrix,Matrix_colunm,Matrix_colunm,Q_Matrix);
    }

    delete[] u_vector;
    delete[] H_Matrix;
    return true;
}

//固有値以外がepsilon以下かどうか
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

//QR分解を用いた固有値計算
double* eig(double* Matrix, int Matrix_colunm, int Matrix_row,double* eigValue)
{
    std::cout << "finding eigenvalues..." << std::endl;

    const int maxIter = 2000;
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
        product(R,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_row,Matrix_buff);

        for(int i = 0; i < Matrix_colunm; i++)Matrix_buff[i * Matrix_row + i] += S;

        if(checkDiagonal(Matrix_buff,Matrix_colunm,Matrix_row))break;
        else iter+= 1;
        
        if(iter > maxIter)break;
    }

    for(int i = 0; i < Matrix_colunm; i++)eigValue[i] = Matrix_buff[i * Matrix_row + i];
    delete[] Matrix_buff;
    delete[] Q;
    delete[] R;
    std::cout << "found eigenvalues..." << std::endl;
    return eigValue;
}

//下三角行列の逆行列計算
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

//上三角行列の逆行列計算
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

//LU分解による逆行列計算(正則のみ)
double* inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* inv_Matrix)
{
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    const double epsilon = 1e-15;

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

    LU_Decompose(1,Matrix,Matrix_colunm,Matrix_row,L,U,Q);

    inverse_L = inverse_lower_triangular_matrix(L,Matrix_colunm,Matrix_row,inverse_L);
    inverse_U = inverse_upper_triangular_matrix(U,Matrix_colunm,Matrix_row,inverse_U);

    P = product(inverse_U,Matrix_colunm,Matrix_row,inverse_L,Matrix_colunm,Matrix_row,P);
    inv_Matrix = product(P,Matrix_colunm,Matrix_row,Q,Matrix_colunm,Matrix_row,inv_Matrix);
    
    delete[] L;
    delete[] U;
    delete[] Q;
    delete[] inverse_L;
    delete[] inverse_U;
    delete[] P;

    return inv_Matrix;
}

//LU分解、QR分解またはSVDを用いた疑似逆行列計算
double* Pseudo_inverse(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* pinv_Matrix)
{
    const int smaller_size = (Matrix_colunm_ > Matrix_row_ ? Matrix_row_ : Matrix_colunm_);
    const double epsilon = std::numeric_limits<double>::epsilon();
    
    int R_rank = 0;
    int new_R_colunm = 0;
    int new_R_row = 0;
    std::vector<int> noneZero_colunm_Index;
    std::vector<int> noneZero_row_Index;

    double* Q = new double[Matrix_colunm_ * Matrix_colunm_];
    double* R = new double[Matrix_colunm_ * Matrix_row_];

    QR_Decompose(Matrix,Matrix_colunm_,Matrix_row_,Q,R);

    for(int i = 0; i < Matrix_row_; i++)
    {
        int zero_Counter = 0; 
        for(int j = 0; j <= (i > Matrix_colunm_ - 1 ? Matrix_colunm_ - 1 : i); j++)
        {
            if(std::abs(R[j * Matrix_row_ + i]) < 1e-12)zero_Counter++;
        }
        if(zero_Counter != (i > Matrix_colunm_ - 1 ? Matrix_colunm_ : i + 1))noneZero_row_Index.push_back(i);
    }

    for(int i = 0; i < Matrix_colunm_; i++)
    {
        int zero_Counter = 0;
        for(int j = i; j < Matrix_row_; j++)
        {
            if(std::abs(R[i * Matrix_row_ + j]) < 1e-12)zero_Counter++;
        }
        if(zero_Counter != (Matrix_row_ - i))noneZero_colunm_Index.push_back(i);
    }

    new_R_colunm = noneZero_colunm_Index.size();
    new_R_row = noneZero_row_Index.size();
    R_rank = std::min(noneZero_colunm_Index.size(), noneZero_row_Index.size());

    if(R_rank == Matrix_colunm_ && R_rank == Matrix_row_)
    {
        delete[] Q;
        delete[] R;

        return inverse(Matrix,Matrix_colunm_,Matrix_row_,pinv_Matrix);
    }
    else if (R_rank == Matrix_colunm_ || R_rank == Matrix_row_)
    {
        int new_Q_row = new_R_colunm;

        double* new_Q = new double[Matrix_colunm_ * new_Q_row];
        double* new_R = new double[new_R_colunm * Matrix_row_];

        int n = 0;
        for(int i : noneZero_colunm_Index)
        {
            for(int j = 0; j < Matrix_colunm_; j++)
            {
                new_Q[j * new_Q_row + n] = Q[j * Matrix_colunm_ + i];
            }
            n++;
        }

        n = 0;
        for(auto i : noneZero_colunm_Index)
        {
            for(int j = 0; j < Matrix_row_; j++)
            {
                new_R[n * Matrix_row_ + j] = R[i * Matrix_row_ + j];
            }
            n++;
        }

        double* new_Rt = new double[Matrix_row_ * new_R_colunm];
        double* new_RRt = new double[new_R_colunm * new_R_colunm];

        double* new_Qt = transpose(new_Q,Matrix_colunm_,new_Q_row,new_Q); 
        new_Rt = transpose(new_R,new_R_colunm,Matrix_row_,new_Rt);

        new_RRt = product(new_R,new_R_colunm,Matrix_row_,new_Rt,Matrix_row_,new_R_colunm,new_RRt);

        double* inv_RRt = inverse(new_RRt,new_R_colunm,new_R_colunm,new_RRt);

        double* Rt_inv_RRt = product(new_Rt,Matrix_row_,new_R_colunm,inv_RRt,new_R_colunm,new_R_colunm,new_Rt);

        product(Rt_inv_RRt,Matrix_row_,new_R_colunm,new_Qt,new_Q_row,Matrix_colunm_,pinv_Matrix);

        delete[] new_Q;
        delete[] new_R;
        delete[] new_Rt;
        delete[] new_RRt;
    }
    else
    {
        double* U = new double[Matrix_colunm_ * Matrix_colunm_];
        double* S = new double[Matrix_colunm_ * Matrix_row_];
        double* V = new double[Matrix_row_ * Matrix_row_];

        SVD(Matrix,Matrix_colunm_,Matrix_row_,U,S,V);

        transpose(U,Matrix_colunm_,Matrix_colunm_,U);
        transpose(S,Matrix_colunm_,Matrix_row_,S);
        transpose(V,Matrix_row_,Matrix_row_,V);

        for(int i = 0; i < smaller_size; i++)S[i * Matrix_row_ + i] = (std::abs(S[i * Matrix_row_ + i]) > epsilon ? 1 / S[i * Matrix_row_ + i] : 0);

        product(V,Matrix_row_,Matrix_row_,S,Matrix_row_,Matrix_colunm_,S);
        product(S,Matrix_row_,Matrix_colunm_,U,Matrix_colunm_,Matrix_colunm_,pinv_Matrix); 

        delete[] U;
        delete[] S;
        delete[] V;
    }

    delete[] Q;
    delete[] R;

    return pinv_Matrix;
}

//LU分解を用いた線形solver
double* equal_solve(double* Matrix, int Matrix_colunm_, int Matrix_row_, double* vector, int vector_size_, double* solution_vector)
{
    //std::cout << "equal solve.." << std::endl;
    const unsigned int Matrix_colunm = Matrix_colunm_;
    const unsigned int Matrix_row    = Matrix_row_;
    const unsigned int vector_size = vector_size_;

    double* L = new double[Matrix_colunm * Matrix_row];
    double* U = new double[Matrix_colunm * Matrix_row];
    double* y = new double[Matrix_colunm];

    for(int i = 0; i < Matrix_colunm_ * Matrix_row_; i++)
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

    delete[] L;
    delete[] U;
    delete[] y;

    //std::cout << "equal solve end..." << std::endl;
    return solution_vector;
}

//SVD計算(いろいろ遅い、改善の余地あり)
void SVD(double* Matrix, int Matrix_colunm, int Matrix_row, double* U, double* S, double* V,int maxIter)
{
    std::cout << "start SVD..." << std::endl;

    const double epsilon = std::numeric_limits<double>::epsilon();

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

    MMt = product(Matrix,Matrix_colunm,Matrix_row,Matrix_t,Matrix_row,Matrix_colunm,MMt);
    MtM = product(Matrix_t,Matrix_row,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,MtM);

    u_eigValues = eig(MMt,Matrix_colunm,Matrix_colunm,u_eigValues);
    v_eigValues = eig(MtM,Matrix_row,Matrix_row,v_eigValues);

    double* eigValues = (Matrix_colunm > Matrix_row ? v_eigValues : u_eigValues);
    double* another_eigValues = (Matrix_colunm > Matrix_row ? u_eigValues : v_eigValues);

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
    }

    for(int i : smaller_size_index)
    {
        double epsilon = 1e-10;
        int index = larger_Size;
        for(int j = 0; j < larger_Size; j++)
        {
            double diff = std::abs(eigValues[i] - another_eigValues[j]);
            if(diff < epsilon)
            {
                index = j;
                epsilon = diff;
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
    }

    for(int i = 0; i < smaller_Size; i++)
    {
        double eigValue = (std::abs(eigValues[i]) > epsilon ? eigValues[i] : 0);

        if(eigValue < 0)eigValue = 0;
        S[i * Matrix_row + i] = std::sqrt(eigValue);

        int iter = 0;

        for(int n = 0; n < Matrix_colunm; n++)for(int m = 0; m < Matrix_colunm; m++)tmp_for_U[n * Matrix_colunm + m] = (n == m ? eigValue : 0);
        for(int n = 0; n < Matrix_row; n++)for(int m = 0; m < Matrix_row; m++)tmp_for_V[n * Matrix_row + m] = (n == m ? eigValue : 0);


        double* sub_tmp_for_U = sub(MMt,tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
        double* sub_tmp_for_V = sub(MtM,tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);

        int maxIndex = Matrix_colunm - 1;

        for(int n = 0; n < Matrix_colunm; n++)
        {
            for(int m = 0;;m++)
            {
                if(std::abs(sub_tmp_for_U[n * Matrix_colunm + m]) > 1e-10 && std::abs(sub_tmp_for_U[m * Matrix_colunm + n]) > 1e-10)break;
                if(m == maxIndex)
                {
                    sub_tmp_for_U[n * Matrix_colunm + n] = 1;
                    break;
                }
            }
        }

        maxIndex = Matrix_row - 1;

        for(int n = 0; n < Matrix_row; n++)
        {
            for(int m = 0;;m++)
            {
                if(std::abs(sub_tmp_for_V[n * Matrix_row + m]) > 1e-10 && std::abs(sub_tmp_for_V[m * Matrix_row + n]) > 1e-10)break;
                if(m == maxIndex)
                {
                    sub_tmp_for_V[n * Matrix_row + n] = 1;
                    break;
                }
            }
        }

        double* inv_tmp_for_U = inverse(sub_tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
        double* inv_tmp_for_V = inverse(sub_tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);

        for(int j = 0; j < Matrix_colunm; j++)u_vector[j] = 1;
        for(int j = 0; j < Matrix_row; j++)v_vector[j] = 1;

        while (iter < maxIter)
        {
            u_vector = product(inv_tmp_for_U,Matrix_colunm,Matrix_colunm,u_vector,Matrix_colunm,1,u_vector);
            double norm_u = norm(u_vector,Matrix_colunm);
            for(int n = 0; n < Matrix_colunm; n++)
            {
                if(norm_u != 0)u_vector[n] /= norm_u;
                else u_vector[n] = 1;
            }

            v_vector = product(inv_tmp_for_V,Matrix_row,Matrix_row,v_vector,Matrix_row,1,v_vector);
            double norm_v = norm(v_vector,Matrix_row);
            for(int n = 0; n < Matrix_row; n++)
            {
                if(norm_v != 0)v_vector[n] /= norm_v;
                else v_vector[n] = 1;
            }

            iter++;
        }
        for(int j = 0; j < Matrix_colunm; j++)U[j * Matrix_colunm + i] = u_vector[j];
        for(int j = 0; j < Matrix_row; j++)V[i * Matrix_row + j] = v_vector[j];
    }

    for(int i = smaller_Size, index = 0; i < larger_Size; i++)
    {
        double eigValue = another_eigValues[i];
        if(Matrix_colunm == larger_Size)
        {
            int iter_for_U = 0;
            for(int n = 0; n < Matrix_colunm; n++)for(int m = 0; m < Matrix_colunm; m++)tmp_for_U[n * Matrix_colunm + m] = (n == m ? eigValue : 0);
            double* sub_tmp_for_U = sub(MMt,tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
            int maxIndex = Matrix_colunm - 1;
            for(int n = 0; n < Matrix_colunm; n++)
            {
                for(int m = 0;;m++)
                {
                    if(std::abs(sub_tmp_for_U[n * Matrix_colunm + m]) > 1e-10 && std::abs(sub_tmp_for_U[m * Matrix_colunm + n]) > 1e-10)break;
                    if(m == maxIndex)
                    {
                        sub_tmp_for_U[n * Matrix_colunm + n] = 1;
                        break;
                    }
                }
            }
            double* inv_tmp_for_U = inverse(sub_tmp_for_U,Matrix_colunm,Matrix_colunm,tmp_for_U);
            for(int j = 0; j < Matrix_colunm; j++)u_vector[j] = 1;
            while (iter_for_U < maxIter)
            {
                u_vector = product(inv_tmp_for_U,Matrix_colunm,Matrix_colunm,u_vector,Matrix_colunm,1,u_vector);
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
            int maxIndex = Matrix_row - 1;

            for(int n = 0; n < Matrix_row; n++)
            {
                for(int m = 0;;m++)
                {
                    if(std::abs(sub_tmp_for_V[n * Matrix_row + m]) > 1e-10 && std::abs(sub_tmp_for_V[m * Matrix_row + n]) > 1e-10)break;
                    if(m == maxIndex)
                    {
                        sub_tmp_for_V[n * Matrix_row + n] = 1;
                        break;
                    }
                }
            }
            double* inv_tmp_for_V = inverse(sub_tmp_for_V,Matrix_row,Matrix_row,tmp_for_V);
            for(int j = 0; j < Matrix_row; j++)v_vector[j] = 1;
            while (iter_for_V < maxIter)
            {
                v_vector = product(inv_tmp_for_V,Matrix_row,Matrix_row,v_vector,Matrix_row,1,v_vector);
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

    delete[] Matrix_t;
    delete[] MMt;
    delete[] MtM;
    delete[] u_eigValues;
    delete[] v_eigValues;
    delete[] tmp_for_U;
    delete[] tmp_for_V;
    delete[] u_vector;
    delete[] v_vector;
    delete[] isExist;

    std::cout << "SVD end..." << std::endl;
}

//NMF計算
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
        YH = product(Matrix,Matrix_colunm,Matrix_row,H_T,Matrix_row,n_components,YH);
        H = transpose(H_T,Matrix_row,n_components,H);
        WH_T = product(W,Matrix_colunm,n_components,H,n_components,Matrix_row,WH_T);
        WH_TH = product(WH_T,Matrix_colunm,Matrix_row,H_T,Matrix_row,n_components,WH_TH);
        W_T = transpose(W,Matrix_colunm,n_components,W_T);
        W_TY = product(W_T,n_components,Matrix_colunm,Matrix,Matrix_colunm,Matrix_row,W_TY);
        W_TWH_T_T = product(W_T,n_components,Matrix_colunm,WH_T,Matrix_colunm,Matrix_row,W_TWH_T_T);
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
        WH_T = product(W,Matrix_colunm,n_components,H,n_components,Matrix_row,WH_T);
        for(int i = 0; i < Matrix_colunm; i++)for(int j = 0; j < Matrix_row; j++)residues += (WH_T[i * Matrix_row + j] - Matrix[i * Matrix_row + j]) * (WH_T[i * Matrix_row + j] - Matrix[i * Matrix_row + j]);

        if(residues < epsilon)break;
    }

    H = transpose(H_T,Matrix_row,n_components,H);

    delete[] H_T;
    delete[] YH;
    delete[] WH_TH;
    delete[] W_TY;
    delete[] W_TWH_T_T;
    delete[] WH_T;

}

//rank計算
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

//使われていないが今後使われるかも?
int eig_rank(double* eig_value, int eig_value_size)
{
    int rank = 0;
    for(int i = 0; i < eig_value_size; i++)if(std::abs(eig_value[i]) > 1e-10)rank++;
    return rank;
}

//互換の積計算(sortした行列をかっこよく並び替えるためだけに作った)
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