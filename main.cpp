#include <iostream>
#include <random>
#include "Matrix.h"
#include "Fraction.h"
#include <complex>
// using namespace std;
using std::cout;
using std::endl;
using std::complex;

// we don't need command line arguments in this program
int main()
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
        {1,-2,1,0},
        {2,1,-3,5},
        {4,-7,1,-1}

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
    complex<float> c5(12, 2);
    complex<float> c6(-2, 4);
    complex<float> c7(1, 0);
    complex<float> c8(6, -1);;
    complex<float> c9(1,1);
    complex<float> c10(-1,0);
    complex<float> c11(0, 1);
    complex<float> c12(0, -1);
    std::vector<std::vector<complex<float>>> ArrayG =
    {
        {c1,c2,c3},
        {c5,c6,c7},
        {c9,c10,c11}
    };
    std::vector<std::vector<complex<float>>> ArrayG1 =
    {
        {c1,c2,c3},
        {c5,c6,c7},
        {c9,c10,c11}
    };


    // testing fractions
   Fraction fr1(6,3);
   fr1.simplify();
   Fraction fr2(2);
   Fraction fr3(3);
   Fraction fr4(4,2);
   Fraction fr5(1,3);
   Fraction fr6(3,2);
   Fraction fr7(9,2);
   Fraction fr8(1,4);
   Fraction fr9(3,4);
   Fraction fr10(1,8);
   Fraction fr11(-2);
   Fraction fr12(2);
   std::vector<std::vector<Fraction>> ArrayH = {
       {fr1, fr2, fr3},
       {fr5, fr6, fr7},
       {fr9, fr10, fr11}
   };







//    Creating Matrices with two-dimensional constructor (columns , rows, Array) ---> Matrix

    Matrix<int> m_a(3, 3, ArrayA);
    Matrix<int> m_b(3, 3, ArrayB);

    Matrix<double> m_c(3, 2, ArrayC);
    Matrix<double> m_d(3, 4, ArrayD);

    Matrix<int> m_e(3, 3, ArrayE);
    Matrix<int> m_f(3, 3, ArrayF);


    Matrix<complex<float>> m_g1(3,3, ArrayG);
    Matrix<complex<float>> m_g2(3,3, ArrayG1);

    Matrix<Fraction> m_h(3,3,ArrayH);
    Matrix<Fraction> m_h2(m_h);
  //  m_g1 = m_g1 + m_g;


  //  cout << 2 * m_e * 3 * m_f << endl; //Scalar and vector product

    cout << "Matrix A: \n" << m_a <<endl;
    cout << "Matrix B: \n" << m_b <<endl;
    cout << "MatrixA + MatrixB: \n" << m_a + m_b  << endl;
    cout << "MatrixA - MatrixB: \n" << m_a - m_b  << endl;
    cout << "Matrixa * Matrixb: \n" << m_a * m_b << endl; //Matrix multiplication
    cout << "MatrixA * 3 \n" << m_a * 3;// Skalar-Multiplication - (matrix<int> m_a * 3)
    //cout << "MatrixA after changes \n" << m_a << endl;
    cout << "\n\n";
    cout << "Matrix C: \n" << m_c <<endl;
    cout << "Matrix D: \n" << m_d <<endl;


    cout << "complex matrixG1:" << endl << m_g1 << endl;
    cout << "complex matrixG2:" << endl << m_g2 << endl;
    cout << "complex matrix G1 * G2:" << endl << m_g1 * m_g1 << endl;

    cout << "fraction matrix1:" << endl << m_h << endl;
    cout << "fraction matrix1 + matrix1:" << endl << m_h + m_h2 << endl;

    m_d.gauss();
/*
    Matrix<double> mat(3,4);
    std::cin>>mat;
    mat.gauss();
*/
}
