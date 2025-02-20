#pragma once
#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>
#include <cmath>
#include "Fraction.h"
#include "Vector.h"
template <typename T>
class Matrix {
public:
    Matrix();
    Matrix(int N, int M);
    Matrix(const Matrix<T>& matrix);
    Matrix(const Vector<T>& vec);
    Matrix(std::vector<Vector<T>>& span);
    ~Matrix();
    Matrix<T> operator+(const Matrix<T>& matrix);
    Matrix<T> operator-(const Matrix<T>& matrix);
    Matrix<T> operator*(const Matrix<T>& matrix);
    Matrix<T> operator*(const T& scalar);
    Matrix<T> operator+=(const Matrix<T>& matrix);
    Matrix<T> operator-=(const Matrix<T>& matrix);
    Matrix<T> operator*=(const T& scalar);
    Matrix<T> operator=(const Matrix<T>& matrix);
    
    static Matrix<T> identity(int N);

    Vector<T> operator*(const Vector<T>& vec);
    
    bool operator==(const Matrix<T>& matrix);

    std::vector<T>& operator[](int index);
    const std::vector<T>& operator[](int index) const;
    Fraction get_discriminant();
    Matrix<T> get_comatrix();

    Matrix<T> get_transpose() const;
    Matrix<T> get_inverse_matrix();



    int lines() const;
    int cols() const;
private:
    Matrix<T> get_minor(int row, int col);
    int N;
    int M;
    std::vector<std::vector<T>> matrix;
};
template<typename T> 
Matrix<T> Matrix<T>::identity(int N) {
    Matrix<T> identity = Matrix<T>(N,N);
    for(int i = 0; i < N; i++)
        identity[i][i] = 1;
    return identity;
}

template <typename T>
Matrix<T>::Matrix(const Vector<T>& vec) {
    this->N=vec.dimension();
    this->M=1;
    for (int i = 0; i < this->N; i++) {
        this->matrix.push_back(std::vector<T>());
        this->matrix[i].push_back(vec[i]);
    }
}
template<typename T>
Matrix<T>::Matrix(std::vector<Vector<T>>& span) {
    if (span.size() < 1) {
        this->N = 0;
        this->M = 0;
        return;
    }
    this->N = span.size();
    this->M = span[0].dimension();
    for (int i = 0; i < this->N; i++) {
        this->matrix.push_back(std::vector<T>());
        for (int j = 0; j < this->M; j++) {
            this->matrix[i].push_back(span[i][j]);
        }
    }
}   

template<typename T>
Matrix<T>::Matrix(int N, int M) {
    for (int i = 0; i < N; i++) {
        this->matrix.push_back(std::vector<T>());
        for (int j = 0; j < M; j++) {
            this->matrix[i].push_back(T());
        }
    }
    this->N = N;
    this->M = M;
}

/*template <typename T>
Matrix<T>::Matrix(int N, int M) {
    this->matrix.resize(N, std::vector<T>(M));
}
*/
template <typename T>
std::vector<T>& Matrix<T>::operator[](int index) {
    return this->matrix[index];
}
template <typename T>
const std::vector<T>& Matrix<T>::operator[](int index) const {
    return this->matrix[index];
}

template <typename T>
Fraction Matrix<T>::get_discriminant() {
    if (this->N < 2 || this->N != this->M)
        return 0;
    if (this->N == 2)
        return this->matrix[0][0] * this->matrix[1][1] - this->matrix[0][1] * this->matrix[1][0];
    Fraction res = 0;
    for (int j = 0; j < this->M; j++)
        res += this->matrix[0][j] * (int)std::pow(-1, 1 + (j + 1)) * (this->get_minor(0, j).get_discriminant());
    return res;
}

template<typename T>
int Matrix<T>::lines() const { return this->N; }
template<typename T>
int Matrix<T>::cols() const { return this->M; }

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Matrix<T>& matrix) {
    for (int i = 0; i < matrix.lines(); i++) {
        stream << "| ";
        for (int j = 0; j < matrix.cols(); j++) {
            stream << matrix[i][j];
            if (j < matrix.cols() - 1)
                stream << " , ";
        }
        stream << " |\n";
    }
    return stream;
}

template <typename T>
Matrix<T>::Matrix() {
    this->N=0;
    this->M=0;
}
template<typename T>
Matrix<T>::Matrix(const Matrix<T>& matrix) {
    for (int i = 0; i < matrix.N; i++) {
        this->matrix.push_back(std::vector<T>());
        for (int j = 0; j < matrix.M; j++)
            this->matrix[i].push_back(matrix[i][j]);
    }
    this->N = matrix.N;
    if (this->N > 1)
        this->M = matrix.M;
    else
        this->M = 0;
}

template <typename T>
Matrix<T>::~Matrix() {}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T>& matrix) {
    if (this->N != matrix.N || this->M != matrix.M)
        return NULL;
    Matrix<T> result(this->N, this->M);
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            result[i][j] = this->matrix[i][j] + matrix[i][j];
    return result;
}


template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T>& matrix) {
    if (this->N != matrix.N || this->M != matrix.M)
        return NULL;
    Matrix<T> result(this->N, this->M);
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            result[i][j] = this->matrix[i][j] - matrix[i][j];
    result.N = this->N;
    return result;
}
template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T>& matrix) {
    //dont help me
    if (this->M != matrix.N)
        throw "The number of lines must match the number of columns of the target ";
    Matrix<T> result = Matrix<T>(this->N,matrix.M); 
    for (int i = 0; i < this->N; i++) {
        for (int j = 0; j < matrix.M; j++) {
            T val = 0;
            for (int k = 0; k < this->M; k++) 
                val += this->matrix[i][k] * matrix[k][j];
            result[i][j] = val;
        }
    }
    return result;
}
template <typename T>
Matrix<T> Matrix<T>::operator*(const T& scalar) {
    Matrix<T> result = Matrix<T>(*this);
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            result[i][j] *= scalar;
    return result;
}

template <typename T>
Vector<T> Matrix<T>::operator*(const Vector<T>& vec) {
    if (this->M != vec.dimension())
        return NULL;
    Vector<T> res = Vector<T>(this->N);
    for (int i = 0; i < this->N; i++) {
        T sigma = T();
        for (int j = 0; j < this->M; j++) 
            sigma += this->matrix[i][j] * vec[j];
        res[i] = sigma;
    }
    return res;
}

template <typename T>
Matrix<T> Matrix<T>::operator+=(const Matrix<T>& matrix) {
    if (this->N != matrix.N || this->M != matrix.M)
        return *this;
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            this->matrix[i][j] += matrix[i][j];
    return *this;
}
template<typename T>
Matrix<T> Matrix<T>::operator-=(const Matrix<T>& matrix) {
    if (this->N != matrix.N || this->M != matrix.M)
        return *this;
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            this->matrix[i][j] -= matrix[i][j];
    return *this;
}
template<typename T>
Matrix<T> Matrix<T>::operator*=(const T& scalar) {
    for (int i = 0; i < this->N; i++) {
        for (int j = 0; i < this->M; j++)
            this->matrix[i][j] *= scalar;
    }
    return *this;
}

template<typename T>
Matrix<T> Matrix<T>::operator=(const Matrix<T>& matrix) {
    if (this != &matrix) {
        while (this->N < matrix.N) {
            this->matrix.push_back(std::vector<T>(matrix[this->N]));
            this->N++;
        }
        while (this->N > matrix.N) {
            this->matrix.pop_back();
            this->N--;
        }
        int diff = matrix.M - this->M;
        for (int i = 0; i < matrix.N; i++) {
            if (diff < 0) {
                for (int j = 0; j < abs(diff); j++)
                    this->matrix[i].pop_back();
            } else if (diff > 0) {
                for (int j = 0; j < diff; j++)
                    this->matrix[i].push_back(0);
            }
            for (int j = 0; j < matrix.M; j++)
                this->matrix[i][j] = matrix[i][j];
        }
        this->M = matrix.M;
    }
    return *this;
}

template <typename T>
bool Matrix<T>::operator==(const Matrix<T>& matrix) {
    if (this->N != matrix.N || this->M != matrix.M)
        return false;
    for (int i = 0; i < this->N; i++)
        for (int j = 0; i < this->M; j++)
            if (this->matrix[i][j] != matrix[i][j])
                return false;
    return true;
}


template <typename T>
Matrix<T> Matrix<T>::get_minor(int row, int col) {
    if (this->N < 2 || this->N != this->M)
        throw "Matrix is not square or is too small";
    Matrix<T> minor = Matrix<T>(this->N - 1, this->M - 1);

    int min_i = 0;
    int min_j = 0;
    for (int i = 0; i < this->N; i++) {
        if (i == row)
            continue;
        for (int j = 0; j < this->M; j++) {
            if (j == col)
                continue;
            minor[min_i][min_j] = this->matrix[i][j];
            min_j++;
        }
        min_i++;
        min_j = 0;
    }
    return minor;
}



template<typename T>
Matrix<T> Matrix<T>::get_transpose() const {
    Matrix<T> transpose = Matrix<T>(this->M, this->N);
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            transpose[j][i] = this->matrix[i][j];
    return transpose;
}

template <typename T>
Matrix<T> Matrix<T>::get_comatrix() {
    Matrix<T> comatrix = Matrix<T>(this->N, this->M);
    for (int i = 0; i < this->N; i++)
        for (int j = 0; j < this->M; j++)
            comatrix[i][j] = (int)std::pow(-1, (i + 1) + (j + 1)) * this->get_minor(i, j).get_discriminant();
    return comatrix;
}


#endif