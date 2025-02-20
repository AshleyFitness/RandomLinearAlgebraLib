#pragma once 
#ifndef AUGMENTED_MATRIX_H
#define AUGMENTED_MATRIX_H

#include "Matrix.h"
#include "Vector.h"
#include "MatrixSpace.h"
enum class SolutionType {
    UNIQUE,
    INFINITE,
    NONE,
    INVALID_DATA,
    UNDEFINED = -1,
};

template <typename T>
class AugmentedMatrix {
public:
    AugmentedMatrix();
    AugmentedMatrix(const Matrix<T>& matrix, const Vector<T>& vec);
    AugmentedMatrix(const Matrix<T>& matrix, const Matrix<T>& target);
    ~AugmentedMatrix();
    
    void gaussian_jordan_elimination();
    void gaussian_elimination();
    
    Vector<T> get_solution();
    Matrix<T> get_matrix_solution();

    int lines() const;

    SolutionType get_solution_type() ;

    static MatrixSpace<T> get_null_space(const Matrix<T>& matrix);
    static MatrixSpace<T> get_column_space(const Matrix<T>& matrix);
    static MatrixSpace<T> get_row_space(const Matrix<T>& matrix);
    static MatrixSpace<T> get_left_null_space(const Matrix<T>& matrix);

    Matrix<T> left_side;
    Matrix<T> right_side;

private:
    void row_addition(int at, int add, T n_time);
    bool normalise(int line);
    void jordan_elimination();

    bool rref;
    int N;

    SolutionType type;
};



template <typename T>
int AugmentedMatrix<T>::lines() const { return this->N; }

// Default Constructor
template <typename T>
AugmentedMatrix<T>::AugmentedMatrix() : N(0), rref(false), type(SolutionType::UNDEFINED) {}

// Constructor with Matrix and Vector
template <typename T>
AugmentedMatrix<T>::AugmentedMatrix(const Matrix<T>& matrix, const Vector<T>& vec) : type(SolutionType::UNDEFINED) {
    this->rref=false;
    if (matrix.lines() != vec.dimension()) {
        this->N = 0;
        return;
    }
    this->left_side = matrix;
    this->right_side = Matrix<T>(vec);
    this->N = matrix.lines();
}

// Constructor with Two Matrices
template <typename T>
AugmentedMatrix<T>::AugmentedMatrix(const Matrix<T>& matrix, const Matrix<T>& target) : type(SolutionType::UNDEFINED) {
    this->rref=false;
    if (matrix.lines() != target.lines()) {
        this->N = 0;
        return;
    }
    this->left_side = matrix;
    this->right_side = target;
    this->N = matrix.lines();
}


template <typename T> 
Vector<T> AugmentedMatrix<T>::get_solution() {
    if(!this->rref)
        this->gaussian_jordan_elimination();
    if(this->type == SolutionType::UNDEFINED) this->get_solution_type();
    if(this->type == SolutionType::UNIQUE) return Vector<T>(this->right_side);
    else return Vector<T>();
}

template <typename T>
Matrix<T> AugmentedMatrix<T>::get_matrix_solution() {
    if(!this->rref)
        this->gaussian_jordan_elimination();
    return this->right_side;
}

template <typename T>
AugmentedMatrix<T>::~AugmentedMatrix() {}


template <typename T>
void AugmentedMatrix<T>::gaussian_jordan_elimination() {

    this->gaussian_elimination();
    this->jordan_elimination();
    this->rref=true;
}

template <typename T> 
void AugmentedMatrix<T>::gaussian_elimination() {
    for(int i = 0; i < this->N;i++ ) {
        if(!this->normalise(i)) continue;
        for(int j = i+1; j < this->N;j++) {
            T to_eliminate = this->left_side[j][i];
            this->row_addition(j,i, to_eliminate * (-1) );
        }
    }
}

template <typename T>
void AugmentedMatrix<T>::jordan_elimination() {
    for(int i = this->N-1;  i >=0;i--) {
        for(int j = i-1;j >=0;j--) {
            T to_eliminate = this->left_side[j][i];
            this->row_addition(j,i, to_eliminate * (-1) );
        }
    }
}

template <typename T>
void AugmentedMatrix<T>::row_addition(int at, int add, T n_time) {
    for(int i = 0; i < this->left_side.cols(); i++) {
        T right =  this->left_side[add][i] * n_time;
        this->left_side[at][i] = this->left_side[at][i] + right;
    }
    for(int i = 0; i < this->right_side.cols();i++) {
        T right =  this->right_side[add][i] * n_time;
        this->right_side[at][i] = this->right_side[at][i] + right;
    }
}

template <typename T>
bool AugmentedMatrix<T>::normalise(int line) {
    bool found_factor = false;
    int factor = -1;
    for(int i = 0; i < this->left_side[line].size(); i++) {
        if(found_factor) 
            this->left_side[line][i]/=factor;
        else if(this->left_side[line][i] != 0) {
            found_factor=true;
            factor = this->left_side[line][i];
            this->left_side[line][i]/=factor;
        }
    }
    if(!found_factor) return false;
    for(int i = 0; i < this->right_side[line].size();i++)
        this->right_side[line][i]/=factor;
    return true;
}
template <typename T>
std::ostream& operator<<(std::ostream& stream, const AugmentedMatrix<T> matrix) {
    for(int i = 0; i < matrix.lines(); i++) {
        stream << "| ";
        for(int j = 0; j < matrix.left_side.cols(); j++) {
            stream << matrix.left_side[i][j];
            if(j < matrix.left_side.cols()-1)
                stream << " , ";
        }
        stream << " | ";

        for(int j = 0; j < matrix.right_side.cols(); j++) {
            stream << matrix.right_side[i][j];
            if(j < matrix.right_side.cols()-1)
                stream << " , ";
        }


        stream << " |\n";
    }
    return stream;
}

template <typename T>
SolutionType AugmentedMatrix<T>::get_solution_type()  {
    
    if(this->type != SolutionType::UNDEFINED) return this->type;
    if(this->right_side.cols() > 1) return SolutionType::INVALID_DATA;
    if(!this->rref)
        this->gaussian_jordan_elimination();
    SolutionType type = SolutionType::UNIQUE;
    for(int i = 0; i < this->N; i++) {
         bool found_factor = false;
         for(int j = 0; j < this->left_side.cols(); j++) { 
            if(!found_factor && this->left_side[i][j] == 1) {
                found_factor=true;
            }
            else if(found_factor && this->left_side[i][j] != 0 ) {
                type = SolutionType::INFINITE;
            } 
        }
        if(!found_factor && this->right_side[i][0] != 0) {
            this->type = SolutionType::NONE;
            return SolutionType::NONE;
        }
    }
    this->type = type;
    return type;
}

template <typename T>
MatrixSpace<T> AugmentedMatrix<T>::get_null_space(const Matrix<T>& matrix) 
{
    AugmentedMatrix<T> working_matrix = AugmentedMatrix<T>(matrix,Vector<T>(matrix.lines()));
    working_matrix.gaussian_jordan_elimination();
    MatrixSpace<T> null_space;

    int N = matrix.lines();
    int M = matrix.cols();
    for(int i = 0; i < M; i++) {
       bool vector_solution = true;
       DimensionalSolution<T> potential_solution;
       potential_solution.col=i;
       for(int j = 0; j < M; j++) {
            if( j >= N ) {
                if(i == j)
                    potential_solution.vector_solution.push_back(1);
                else
                    potential_solution.vector_solution.push_back(0);
            } else {
                potential_solution.vector_solution.push_back(working_matrix.left_side[j][i]);
                if( j != i && working_matrix.left_side[j][i] != 0)
                    vector_solution=false;
            }
        }
        if(!vector_solution)
            null_space.solutions.push_back(potential_solution);
    }
    null_space.rank=null_space.solutions.size();
    if(!null_space.rank) {
        DimensionalSolution<T> null_solution;
        null_solution.col=0;
        null_solution.vector_solution = std::vector<T>(M);
        null_space.solutions.push_back(null_solution);
    } 
    return null_space;
}

template <typename T>
MatrixSpace<T> AugmentedMatrix<T>::get_column_space(const Matrix<T>& matrix) { 
    AugmentedMatrix<T> working_matrix = AugmentedMatrix<T>(matrix,Vector<T>(matrix.lines()));
    working_matrix.gaussian_jordan_elimination();
    MatrixSpace<T> column_space;
    
    int N = matrix.lines();
    int M = matrix.cols();
    for(int i = 0; i < M; i++) {
        bool vector_solution = true; 
        DimensionalSolution<T> potential_solution;
        potential_solution.col=i;
        for(int j = 0; j < N; j++) {
            potential_solution.vector_solution.push_back(matrix[j][i]);
            if(j != i && working_matrix.left_side[j][i] != 0) {
                vector_solution=false;
                break;
            }
        }
        if(vector_solution)
            column_space.solutions.push_back(potential_solution);
    }
    column_space.rank = column_space.solutions.size();
    if(!column_space.rank) {
        DimensionalSolution<T> null_solution;
        null_solution.col = 0; 
        null_solution.vector_solution = std::vector<T>(N);
        column_space.solutions.push_back(null_solution);
    }
    return column_space;
}

template <typename T>
MatrixSpace<T> AugmentedMatrix<T>::get_row_space(const Matrix<T>& matrix) {
    Matrix<T> alt = matrix.get_transpose();
    return AugmentedMatrix<T>::get_null_space(alt);
}

template <typename T>
MatrixSpace<T> AugmentedMatrix<T>::get_left_null_space(const Matrix<T>& matrix) {
    Matrix<T> alt = matrix.get_transpose();
    return AugmentedMatrix<T>::get_column_space(alt);
}

#endif 