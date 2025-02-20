#pragma once
#ifndef VECTOR_H
#define VECTOR_H
#include <ostream>
#include <vector>
#include <cmath>

template <typename T>
class Matrix;

template <typename T>
class Vector
{
public:
	Vector();
	Vector(int n);
	Vector(const Vector<T>& vec);
	Vector(std::vector<T>& vec_arr);
	Vector(const Matrix<T>& matrix);
	~Vector();

	Vector<T> operator+(const Vector<T>& vec) const;
	Vector<T> operator+=(const Vector<T>& vec);
	
	Vector<T> operator-(const Vector<T>& vec) const;
	Vector<T> operator-=(const Vector<T>& vec);

	Vector<T> cross_product(const Vector<T>& vec) const;
	T operator*(const Vector<T>& vec) const;


    Vector<T> operator*(double x);
    Vector<T> operator*(double x) const;


	bool operator==(const Vector<T>& vec) const;
	Vector<T> operator=(const Vector<T>& vec);
	
	Vector<T> normalized();

	T& operator[](int index) ;
	const T& operator[](int index) const;

	std::vector<T> to_array();
	
	double len();
	int dimension() const;
private:
	int N;
	std::vector<T> container;
};


template <typename T>
Vector<T>::Vector(const Matrix<T>& matrix) {
	if(matrix.cols() > 1) {
		this->N=0;
		return;
	}
	this->N = matrix.lines();
	for(int i = 0; i < this->N; i++) 
		this->container.push_back(matrix[i][0]);
}


template <typename T>
std::vector<T> Vector<T>::to_array() {
	return this->container;
}
template <typename T>
Vector<T>::Vector(std::vector<T>& vec_arr) {
	this->N = vec_arr.size();
	this->container = vec_arr;
}

template <typename T>
Vector<T>::Vector() { this->N = 0; }

template <typename T>
Vector<T>::Vector(int n) {
	for (int i = 0; i < n; i++)
		this->container.push_back(T());
	this->N = n;
}

template <typename T>
Vector<T>::Vector(const Vector<T>& vec) {
	for (int i = 0; i < vec.N; i++)
		this->container.push_back(vec[i]);
	this->N = vec.N;
}

template <typename T>
Vector<T>::~Vector() {}

template <typename T>
Vector<T> Vector<T>::operator+(const Vector<T>& vec) const {
	if (this->N != vec.N)
		return NULL;
	Vector<T> res = Vector<T>(*this);
	for (int i = 0; i < this->N; i++)
		res.container[i] += vec[i];
	return res;
}
template <typename T>
Vector<T> Vector<T>::operator+=(const Vector<T>& vec) {
    if (this->N != vec.N)
		return NULL;
    
	for (int i = 0; i < this->N; i++)
		this->container[i] += vec[i];
	return *this;
}

template <typename T>
Vector<T> Vector<T>::operator-(const Vector<T>& vec) const {
    if (this->N != vec.N)
		return NULL;
	Vector<T> res = Vector<T>(*this);
	for (int i = 0; i < this->N; i++)
		res.container[i] -= vec[i];
	return res;
}
template <typename T>
Vector<T> Vector<T>::operator-=(const Vector<T>& vec) {
	if (this->N != vec.N)
		return NULL;
	for (int i = 0; i < this->N; i++)
		this->container[i] -= vec[i];
	return *this;
}


template<typename T>
Vector<T> Vector<T>::normalized() {
	double scalar = 1 / this->len();
	Vector<T> res = Vector<T>(*this);
	for (int i = 0; i < this->N; i++) 
		res[i] *= scalar;
	return res;
}

template <typename T>
Vector<T> Vector<T>::cross_product(const Vector<T>& vec) const {
	if (vec.N != 3 || this->N != 3)
		return NULL;
	Vector<T> res = Vector<T>(3);
	res[0] = ((*this)[1] * vec[2]) - ((*this)[2] * vec[1]);
	res[1] = -1 * (((*this)[0] * vec[2]) - ((*this)[2] * vec[0]));
	res[2] = ((*this)[0] * vec[1]) - ((*this)[1] * vec[0]);
	return res;
}
template<typename T>
T Vector<T>::operator*(const Vector<T>& vec) const {
	if (this->N != vec.N)
		return NULL;
	T res = T();
	for (int i = 0; i < this->N; i++)
		res += (*this)[i] * vec[i];
	return res;
}


template<typename T>
Vector<T> Vector<T>::operator*(double x) {
    Vector<T> res = Vector(*this);
    for(int i = 0; i < res.N;i++) 
        res[i]*=x;
    return res;
}
template<typename T>
Vector<T> Vector<T>::operator*(double x) const {
    Vector<T> res = Vector(*this);
    for(int i = 0; i < res.N;i++) 
        res[i]*=x;
    return res;
}
template <typename T>
bool Vector<T>::operator==(const Vector<T>& vec) const {
	if (this->N != vec.N)
		return false;
	for (int i = 0; i < this->N; i++)
		if ((*this)[i] != vec[i])
			return false;
	return true;
}

template <typename T>
Vector<T> Vector<T>::operator=(const Vector<T>& vec) {
	if (this != &vec) {
		this->N = vec.N;
		this->container = vec.container;
	}
	return *this;
}

template <typename T>
T& Vector<T>::operator[](int index) {
	return this->container[index];
}

template <typename T>

const T& Vector<T>::operator[](int index) const {	
	return this->container[index];
}

template <typename T>
double Vector<T>::len() {
	T sigma = T();
	for (int i = 0; i < this->N; i++)
		sigma += std::pow((*this)[i],2);
	return std::sqrt((double)sigma);
}

template <typename T>
int Vector<T>::dimension() const { return this->N; }

template <typename T>
std::ostream& operator<<(std::ostream& stream, const Vector<T>& vec) {
	stream << "| ";
	int n = vec.dimension();
	for (int i = 0; i < n; i++) {
		stream << vec[i];
		if (i < n - 1)
			stream << " , ";
	}
	stream << " |";
	return stream;
}


#endif