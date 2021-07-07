#ifndef FRACTION_HEADER

#define FRACTION_HEADER
#include <iostream>


class Fraction{
friend std::ostream& operator<<(std::ostream& os, Fraction fr);
private:
    int numerator, denominator;
    void no_zero_division();

public:
    Fraction();
    Fraction(int n);
    Fraction(int n, int d);
    //why doesn't it work??
    // Fraction(const Fraction& fr);


    Fraction& operator=(Fraction fr2);

    Fraction& operator+=(Fraction& fr2);
    Fraction& operator-=(Fraction& fr2);
    Fraction& operator*=(Fraction& fr2);
    Fraction& operator/=(Fraction& fr2);

    Fraction operator+(Fraction fr2);
    Fraction operator-(Fraction fr2);
    Fraction operator*(Fraction fr2);
    Fraction operator/(Fraction fr2);


    Fraction& simplify();


};


#endif
