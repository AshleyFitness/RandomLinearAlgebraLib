#include "Fraction.h"

Fraction::Fraction(int x, int y) {
    this->numerator = x;
    this->denumerator = y;
    this->simplify();
}

Fraction::Fraction(int x) {
    this->numerator = x;
    this->denumerator = 1;
}
Fraction::Fraction() {
    this->numerator = 0;
    this->denumerator = 1;
}

Fraction::Fraction(const Fraction& fraction) {
    this->numerator = fraction.numerator;
    this->denumerator = fraction.denumerator;
    this->simplify();
}

Fraction::~Fraction() {}

Fraction Fraction::operator+(const Fraction& fraction) {
    Fraction res = Fraction(*this);
    if (fraction.denumerator != this->denumerator) {
        int temp_denumerator = res.denumerator;
        res.denumerator *= fraction.denumerator;
        res.numerator *= fraction.denumerator;
        res.numerator += fraction.numerator * temp_denumerator;
    }
    else {
        res.numerator += fraction.numerator;
    }
    res.simplify();
    return res;
}

Fraction Fraction::operator+(const int x) {
    Fraction res = Fraction(*this);
    res.numerator += x * res.denumerator;
    res.simplify();
    return res;
}

Fraction Fraction::operator+=(const Fraction& fraction) {
    if (fraction.denumerator != this->denumerator) {
        int temp_denumerator = this->denumerator;
        this->denumerator *= fraction.denumerator;
        this->numerator *= fraction.denumerator;
        this->numerator += fraction.numerator * temp_denumerator;
    }
    else
        this->numerator += fraction.numerator;
    this->simplify();
    return *this;
}

Fraction Fraction::operator+=(const int x) {
    this->numerator += x * this->denumerator;
    this->simplify();
    return *this;
}

Fraction Fraction::operator-(const Fraction& fraction) {
    Fraction res = Fraction(*this);
    if (fraction.denumerator != this->denumerator) {
        int temp_denumerator = res.denumerator;
        res.denumerator *= fraction.denumerator;
        res.numerator *= fraction.denumerator;
        res.numerator -= fraction.numerator * temp_denumerator;
    }
    else
        res.numerator -= fraction.numerator;
    res.simplify();
    return res;
}

Fraction Fraction::operator-(const int x) {
    Fraction res = Fraction(*this);
    res.numerator -= x * res.denumerator;
    res.simplify();
    return res;
}


Fraction Fraction::operator-=(const Fraction& fraction) {
    if (fraction.denumerator != this->denumerator) {
        int temp_denumerator = this->denumerator;
        this->denumerator *= fraction.denumerator;
        this->numerator *= fraction.denumerator;
        this->numerator -= fraction.numerator * temp_denumerator;
    }
    else
        this->numerator -= fraction.numerator;
    this->simplify();
    return *this;
}

Fraction Fraction::operator-=(const int x) {
    this->numerator -= x * this->denumerator;
    this->simplify();
    return *this;
}

Fraction Fraction::operator*(const int x) {
    Fraction res = Fraction(*this);
    res.numerator *= x;
    res.simplify();
    return res;
}

Fraction Fraction::operator*(const Fraction& fraction) {
    Fraction res = Fraction(*this);
    res.numerator *= fraction.numerator;
    res.denumerator *= fraction.denumerator;
    res.simplify();
    return res;
}

Fraction Fraction::operator*=(const int x) {
    this->numerator *= x;
    this->simplify();
    return *this;
}

Fraction Fraction::operator*=(const Fraction& fraction) {
    this->numerator *= fraction.numerator;
    this->denumerator *= fraction.denumerator;
    this->simplify();
    return *this;
}

Fraction Fraction::operator/(const int x) {
    Fraction res = Fraction(*this);
    res.denumerator *= x;
    res.simplify();
    return res;
}

Fraction Fraction::operator/(const Fraction& fraction) {
    Fraction res = Fraction(*this);
    res.numerator = res.numerator * fraction.denumerator;
    res.denumerator = res.denumerator * fraction.numerator;
    res.simplify();
    return res;
}

Fraction Fraction::operator/=(const int x) {
    this->denumerator *= x;
    this->simplify();
    return *this;
}

Fraction Fraction::operator/=(const Fraction& fraction) {
    this->numerator = this->numerator * fraction.denumerator;
    this->denumerator = this->denumerator * fraction.numerator;
    this->simplify();
    return *this;
}
bool Fraction::operator==(int x) {
    return (double)(this->numerator/this->denumerator) == (double)x;   
}

bool Fraction::operator==(double x) {
    return (double)(this->numerator/this->denumerator) == x;   
}

bool Fraction::operator==(float x) {
    return (float)(this->numerator/this->denumerator) == x;   
}

Fraction Fraction::operator=(const Fraction& fraction) {
    this->denumerator = fraction.denumerator;
    this->numerator = fraction.numerator;
    this->simplify();
    return *this;
}

bool Fraction::operator==(const Fraction& fraction) {
    Fraction right_op = Fraction(fraction);
    this->simplify();
    right_op.simplify();
    return (this->numerator == right_op.numerator && this->denumerator == right_op.denumerator);
}

void Fraction::simplify() {
    if (this->numerator == 0)
        return;
    int div = this->gcd(std::max(this->numerator, this->denumerator),
        std::min(this->denumerator, this->numerator));

    this->numerator /= div;
    this->denumerator /= div;

    if(this->denumerator < 0 || this->numerator < 0) {
        this->numerator = -abs(this->numerator);
        this->denumerator = abs(this->denumerator);
    }

    if (this->numerator < 0 && this->denumerator < 0) {
        this->numerator = abs(this->numerator);
        this->denumerator = abs(this->denumerator);
    }
}

std::ostream& operator<<(std::ostream& stream, const Fraction& fraction) {
    if (fraction.numerator % fraction.denumerator == 0) 
        stream << fraction.numerator / fraction.denumerator;
    else {
        //since the fraction is always simplified, we won't need to check if both are negatives
        if (fraction.numerator < 0 || fraction.denumerator < 0)
            stream << " - ";
        stream << abs(fraction.numerator) << " / " << abs(fraction.denumerator);

    }
    return stream;
} 
int Fraction::gcd(int x, int y) {
    int num = x;
    int div = y;
    int r = -1;
    while (r != 0) {
        r = num % div;
        if (r != 0) {
            num = div;
            div = r;
        }
    }
    return div;
}