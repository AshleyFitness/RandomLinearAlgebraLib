#pragma once
#ifndef DISCRETEMATH_H
#define DISCRETEMATH_H
#include "Vector.h"
#include "AugmentedMatrix.h"
#include <cmath>
#include <array>
template <typename T>
class Matrix;




template <typename T> 
inline Matrix<T> gram_schmidt(Matrix<T>& span) {

    //I'll assume that Vectors are set in a column way.
    //Just to respect Mathematical standards, even tho i agree it ain't the most convenient way to work with
    //vectors.

    if (span.lines() < 1)
        return Matrix<T>();
    Vector<T> res_vec = Vector<T>(span.cols());

    Matrix<T> working_span = span.get_transpose();
    //To work more easily we will work with the tranpose

    working_span[0] = Vector<T>(working_span[0]).normalized().to_array();

    for (int i = 1; i < working_span.lines(); i++) {
        Vector<T> res(working_span.cols());
        for (int j = i - 1; j >= 0; j--) {
            T dot_product = Vector<T>(working_span[j]) * Vector<T>(working_span[i]);
            Vector<T> normalized_vector = ( Vector<T>(working_span[j]) * dot_product);
            res += normalized_vector;
        }
        res = Vector<T>(working_span[i]) - res;
        working_span[i] = (res * (double)(1.0 / res.len())).to_array();
    }
    return working_span.get_transpose();
}

template <typename T>
inline Vector<T> least_square_solution(Matrix<T>& transformation, Vector<T>&  target) {
    Matrix<T> transpose = transformation.get_transpose();
    Matrix<T> adapted_transformation = transpose * transformation;
    Vector<T> adapted_vector = transpose * target;

    AugmentedMatrix<T> working_matrix = AugmentedMatrix<T>(adapted_transformation,adapted_vector);
    working_matrix.gaussian_jordan_elimination();
    return working_matrix.get_solution();
}

inline Matrix<double> rotate_2d(int degree) {
    Matrix<double> res = Matrix<double>(2,2);
    double PI = 	3.14159265358979323846;
    double radians = degree * PI / 180.0;
    res[0][0] = cos(radians);
    res[0][1]=-1*sin(radians);
    res[1][0]=sin(radians);
    res[1][1]=cos(radians);
    return res;
}
inline Matrix<double> rotate_3d(int degree,int point) {
    Matrix<double> res = Matrix<double>(3,3);
    double PI = 	3.14159265358979323846;
    double radians = degree * PI / 180.0;
    switch(point) {
        case 1:
            res[0][0]=1; res[0][1]=0; res[0][2]=0; res[1][0]=0; res[2][0]=0;
            res[1][1]=cos(radians); res[1][2]=-1*sin(radians); res[2][1]=sin(radians); res[2][2]=cos(radians);
            return res;
        case 2:
            res[1][1]=1; res[1][0]=0; res[1][2]=0; res[0][1]=0; res[2][1]=0;
            res[0][0]=cos(radians); res[0][2]=-1*sin(radians); res[2][0]=sin(radians); res[2][2]=cos(radians);
            return res;
        case 3:
            res[2][2]=1; res[2][0]=0; res[2][1]=0; res[0][2]=0; res[1][2]=0;
            res[0][0]=cos(radians); res[0][1]=-1*sin(radians); res[1][0]=sin(radians); res[1][1]=cos(radians);
            return res;
        default:
            return res;
    }
}
template <typename T>
inline std::array<Matrix<T>,2> lu_factorisation(const Matrix<T>& target) {
    if(target.lines() != target.cols())
        throw std::runtime_error("Matrix is not square");
    int n = target.lines();
    Matrix<T> L = Matrix<T>::identity(target.lines());
    Matrix<T> U = Matrix<T>(target);

    for(int i = 0; i < n; i++) {
        
        if(U[i][i] == 0)
            throw std::runtime_error("Matrix is not invertible");
        
        for(int j = i + 1; j < n; j++) {
            T factor = U[j][i] / U[i][i];
            L[j][i] = factor;
            for(int k = i; k < n; k++)
                U[j][k] -= factor * U[i][k];
        }

    }
    return {L,U};
}

#endif