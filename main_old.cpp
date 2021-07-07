#include <iostream>
#include <random>
#include "Matrix.h"
#include "Fraction.h"
#include <complex>
// using namespace std;
using std::cout;
using std::endl;
using std::complex;


int main(int argc, char* argv[])
{



    std::vector<std::vector<int>> ArrayA =
    {
        {1,2,3},
        {4,5,6},
        {7,8,9}
    };

    std::vector<std::vector<int>> ArrayB =
    {
        {1,2,3},
        {4,5,6},
        {7,8,9}
    };


    std::vector<std::vector<double>> ArrayC =
    {
        {1,2},
        {4,5},
        {7,8}
    };

    std::vector<std::vector<double>> ArrayD =
    {
        {1,2,3,10},
        {4,5,6,10},


    };

    std::vector<std::vector<int>> ArrayE =
    {
        {1,1,1},
        {1,1,1},
        {1,1,1}
    };

    std::vector<std::vector<int>> ArrayF =
    {
        {6,6,6},
        {6,6,6},
        {6,6,6}
    };

    complex<float> c1(2.5, 5.2);
    complex<float> c2(22.2, -2);
    complex<float> c3(-1.125,3.8);
    complex<float> c4(-1.25, -1.25);

    std::vector<std::vector<complex<float>>> ArrayG =
    {
        {c1,c2},
        {c3,c4}
    };






//    Creating Matrices with two-dimensional constructor (columns , rows, Array) ---> Matrix

    Matrix<int> m_a(3, 3, ArrayA);
    Matrix<int> m_b(3, 3, ArrayB);

    Matrix<double> m_c(3, 2, ArrayC);
    Matrix<double> m_d(2, 4, ArrayD);

    Matrix<int> m_e(3, 3, ArrayE);
    Matrix<int> m_f(3, 3, ArrayF);

    Matrix<complex<float>> m_g(2,2, ArrayG);

    cout << 2 * m_e * 3 * m_f << endl; //Scalar and vector product

    cout << "Matrix A: \n" << m_a <<endl;
    cout << "Matrix B: \n" << m_b <<endl;
    cout << "MatrixA + MatrixB: \n" << m_a + m_b  << endl;
    cout << "MatrixA after changes \n" << m_b << endl;
    cout << "MatrixA - MatrixB: \n" << m_a - m_b  << endl;

    cout << "MatrixA * 3 \n" << m_a * 3;// Skalar-Multiplication - matrix<int> * int 3 type
    cout << "\n\n";
    cout << "Matrix C: \n" << m_c <<endl;
    cout << "Matrix D: \n" << m_d <<endl;

    cout << "Matrixa * Matrixb: \n" << m_a * m_b << endl; //Working matrix multiplication vector multiplication
    cout << "MatrixD * MatrixC: \n" << m_c * m_d << endl; //Working matrix multiplication

    cout << "complex matrixG:" << endl << m_g << endl;
    // this doesn't work yet:
   // cout << "gaussian_elemination of MatrixD:" << endl << m_d.gaussian_elemination() << endl; 

    
   //example of the new () operator:
    cout << "Matrix e: \n";
    cout << m_e;
    cout << "accessing elemnts one by one: \n";
    cout << m_e(0,0) << " " << m_e(0,1) << " " << m_e(0,2) << endl;
    cout << "changing an element manually: \n";
    m_e(0,0) = 10;
    cout << m_e;
    
    return 0;
}
