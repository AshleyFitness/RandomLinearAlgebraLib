#pragma once
#ifndef FRACTION_H
#define FRACTION_H
#include <memory>
#include <type_traits>
#include <iostream>
class Fraction {
public:
    Fraction(int x, int y);
    Fraction(int x);
    Fraction(const Fraction& fraction);
    Fraction();

    ~Fraction();

    Fraction operator+(const Fraction& fraction);
    Fraction operator+(const int x);
    
    Fraction operator+=(const int x);
    Fraction operator+=(const Fraction& fraction);

    Fraction operator-(const Fraction& fraction);
    Fraction operator-(const int x);

    Fraction operator-=(const int x);
    Fraction operator-=(const Fraction& fraction);

    Fraction operator*(const Fraction& fraction);
    Fraction operator*(const int x);

    Fraction operator*=(const Fraction& fraction);
    Fraction operator*=(const int x);


    Fraction operator/(const Fraction& fraction);
    Fraction operator/(const int x);

    Fraction operator/=(const Fraction& fraction);
    Fraction operator/=(const int x);

    Fraction operator=(const Fraction& fraction);
    
    bool operator==(int x);
    bool operator==(double x);
    bool operator==(float x);

    bool operator==(const Fraction& fraction);

    operator double() const {
        return static_cast<double>(numerator) / denumerator;
    }

    void simplify();
    int numerator;
    int denumerator;
    private:
    int gcd(int x , int y);
};
std::ostream& operator<<(std::ostream& stream, const Fraction& fraction);
#endif