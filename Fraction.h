#ifndef FRACTION_HEADER

#define FRACTION_HEADER
#include <iostream>


class Fraction{
friend std::ostream& operator<<(std::ostream& os, Fraction fr);
friend Fraction& abs(Fraction&);
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

// put abs overload here instead of Fraction.cpp cause it needs to run in Matrix.h
// and i don't want to spend time on thinkink about why i can't
Fraction& abs(Fraction& fr){
    fr.numerator = abs(fr.numerator);
    fr.denominator = abs(fr.denominator);
    return fr;
}


#endif
