#pragma once
#ifndef MATRIX_SPACE_H
#define MATRIX_SPACE_H
#include <vector>
template <typename T>
struct DimensionalSolution {
    int col;
    std::vector<T> vector_solution;
};
template <typename T>
struct MatrixSpace {
    std::vector<DimensionalSolution<T>> solutions;
    int rank;
};

template <typename T>
std::ostream& operator<<(std::ostream& stream, const MatrixSpace<T>& space) {
    if(space.rank==0) {
        stream << "[ ";
        for(int i = 0 ; i < space.solutions[0].vector_solution.size(); i++) {
            stream << space.solutions[0].vector_solution[i];
            if(i < space.solutions[0].vector_solution.size()-1)
                stream << ", ";
        }
        stream << " ]\n";
        return stream;
    }
    for(int i = 0;  i < space.solutions.size();i++) {
        stream << "x(" << space.solutions[i].col+1 << ") : [ "; 
        for(int j = 0; j < space.solutions[i].vector_solution.size(); j++) {
            stream << space.solutions[i].vector_solution[j];
            if( j < space.solutions[i].vector_solution.size()-1)
                stream << " , ";
        }
        stream << "]";
    }
    stream << "\n"; 
    return stream;
}
#endif