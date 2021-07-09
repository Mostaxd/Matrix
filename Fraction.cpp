#include "Fraction.h"
using std::cout;
using std::cerr;
using std::endl;

// function zur ermittlung groessten gemeinsamen teilers - aus 1. uebung
long ggT(long a, long b){

    a = abs(a);
    b = abs(b);

    while (a != b) {
        if (a > b)
            a -= b;
        else
            b -= a;
    }
    // both equal -> return either one
    return a;
}




////////////////////////////////////////////////////////////////////////////////
// private member functions
////////////////////////////////////////////////////////////////////////////////

void Fraction::no_zero_division(){ cerr << "--- ERROR: DIVISION BY 0 - REFUSED! ---" << endl; exit(0); }




////////////////////////////////////////////////////////////////////////////////
// public member functions
////////////////////////////////////////////////////////////////////////////////


// constructor functions:
////////////////////////////////////////////////////////////////////////////////

// constructor without parameter    --> initialized as 1/1
Fraction::Fraction(){
    numerator = 1;
    denominator = 1;
}

// constructor with one parameter   --> initialized as value/1
Fraction::Fraction(int value){
    numerator = value;
    denominator = 1;
}

//constructor with two parameters   --> initialized as p1/p2
Fraction::Fraction(int value1, int value2){
    if(value2 == 0) no_zero_division();
    numerator = value1;
    denominator = value2;
}

//constructor with existing fraction
// why doesn't it work??
// Fraction::Fraction(const Fraction& fr) : numerator(fr.numerator), denominator(fr.denominator){}



// overloading operators:
////////////////////////////////////////////////////////////////////////////////

std::ostream& operator<<(std::ostream& os, Fraction fr){
    os << "(" << fr.numerator << " / " << fr.denominator << ")";
    return os;
}


Fraction& Fraction::operator=(Fraction fr2){
    numerator = fr2.numerator;
    denominator = fr2.denominator;
    return *this;
}


Fraction& Fraction::operator+=(Fraction& fr2){
    numerator = numerator * fr2.denominator + fr2.numerator * denominator;
    denominator *= fr2.denominator;
    return *this;
}

Fraction& Fraction::operator-=(Fraction& fr2){
    numerator = numerator * fr2.denominator - fr2.numerator * denominator;
    denominator *= fr2.denominator;
    return *this;
}

Fraction& Fraction::operator*=(Fraction& fr2){
    numerator *= fr2.numerator;
    denominator *= fr2.denominator;
    return *this;
}

Fraction& Fraction::operator/=(Fraction& fr2){
    if(fr2.numerator == 0) Fraction::no_zero_division();
    numerator *= fr2.denominator;
    denominator *= fr2.numerator;
    return *this;
}


Fraction Fraction::operator+(Fraction fr2){
    Fraction fr3;
    fr3 = *this;
    fr3 += fr2;
    return fr3;
}

Fraction Fraction::operator-(Fraction fr2){
    Fraction fr3;
    fr3 = *this;
    fr3 -= fr2;
    return fr3;
}

Fraction Fraction::operator*(Fraction fr2){
    Fraction fr3;
    fr3 = *this;
    fr3 *= fr2;
    return fr3;
}

Fraction Fraction::operator/(Fraction fr2){
    Fraction fr3;
    fr3 = *this;
    fr3 /= fr2;
    return fr3;
}


Fraction& Fraction::operator+=(int value){
    numerator += (denominator * value);
    return *this;
}

Fraction& Fraction::operator-=(int value){
    numerator -= (denominator * value);
    return *this;
}

Fraction& Fraction::operator*=(int value){
    numerator *= value;
    return *this;
}
Fraction& Fraction::operator/=(int value){
    if (value == 0) no_zero_division();
    denominator/= value;
    return *this;
}


Fraction Fraction::operator+(int value){
    // this stays untouched!
    Fraction fr = *this;
    fr += value;
    return fr;
}

Fraction Fraction::operator-(int value){
    // this stays untouched!
    Fraction fr = *this;
    fr -= value;
    return fr;
}

Fraction Fraction::operator*(int value){
    // this stays untouched!
    Fraction fr = *this;
    fr *= value;
    return fr;
}
Fraction Fraction::operator/(int value){
    // this stays untouched!
    Fraction fr = *this;
    fr /= value;
    return fr;
}
// other public members
////////////////////////////////////////////////////////////////////////////////

Fraction& Fraction::simplify(){
    if (numerator == 0) return *this;
    int divisor = ggT(numerator, denominator);
    numerator /= divisor;
    denominator /= divisor;
    return *this;
}



// abs overload to return absolute value of fraction
Fraction& abs(Fraction& fr){
    fr.numerator = abs(fr.numerator);
    fr.denominator = abs(fr.denominator);
    return fr;
}
