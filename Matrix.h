#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <complex>
#include<iomanip>
//using namespace std;
using std::complex;
using std::vector;

template  <typename T>
class Matrix
{
private:
	int cols; // Number of columns
	int rows; // number of rows
	int size; // array size
	T* data; // Element array




public:
	Matrix(); // Default constructor
	Matrix(int rows, int cols, T val); // constructor that fills matrix with val
	Matrix(int rows, int cols);	//empty constructor initialized with 0
	Matrix(int rows, int cols, const vector<vector<T>>& Array); // A Two-dimensional array to construct a matrix
	Matrix(const Matrix& matrix);                      // Use an existing matrix object to construct a matrix

	~Matrix() { delete[]data; }; //deconstructor (deletes arrays after usage, manually emptying heap memory)


	int  getCols() const { return cols; };                 // Get the number of columns
	int  getRows() const { return rows; };                 // Get the number of rows
	int  getSize() const { return rows * cols; };          // Get the size of the array of Matrix



    Matrix<T>& gauss();

	T& operator()(int row, int col);
	template  <typename ElemType>
	friend Matrix<ElemType>          operator +(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2); // Matrix and matrix addition (+ operator overload)
	template  <typename ElemType>
	friend Matrix<ElemType>          operator -(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2); // Matrix and matrix subtraction(-operator overloading)
	template  <typename ElemType>
	friend Matrix<ElemType>          operator *(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2); // Matrix and matrix multiplication (* operator overload)


	// stream operators
	template <typename COMPLEXTYPE>
	friend std::ostream&            operator <<(std::ostream& os, const Matrix<complex<COMPLEXTYPE>>& matrix); // Output matrix of type Complex to output stream
	template  <typename ElemType>
<<<<<<< HEAD
	friend std::ostream& operator <<(std::ostream& os, const Matrix<ElemType>& matrix); // Output matrix to output stream
    template  <typename ElemType>
	friend std::istream& operator >>(std::istream& os, const Matrix<ElemType>& matrix); // Output matrix to output stream
=======
	friend std::ostream&            operator <<(std::ostream& os, const Matrix<ElemType>& matrix); // Output matrix to output stream
    template  <typename ElemType>
	friend std::istream&            operator >>(std::istream& os, const Matrix<ElemType>& matrix); // Output matrix to output stream
>>>>>>> defcb6fc27c0d0c76d6531a47697fbf637fdca9f

	//Skalar Funktionen:
	template  <typename ElemType>
	friend Matrix<ElemType>         operator *(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number multiplication * operator overload (1)
	template  <typename ElemType>
<<<<<<< HEAD
	friend Matrix<ElemType>   operator *(const ElemType value, const Matrix<ElemType>& matrix);       // Matrix and number multiplication * operator overload (2)


	Matrix<T>&   operator/(const T value);       // Matrix and number division/operator overload (1)
=======
	friend Matrix<ElemType>         operator *(const ElemType value, const Matrix<ElemType>& matrix);       // Matrix and number multiplication * operator overload (2)
    template  <typename ElemType>
	friend Matrix<ElemType>         operator /(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number division/operator overload (1)
>>>>>>> defcb6fc27c0d0c76d6531a47697fbf637fdca9f

};










////////////////////////////////////////////////////////////////////////////////
// class member initialization

// Default constructor
template  <typename T>
Matrix<T>::Matrix()
{
	std::cout << "Matrix constructor without paramerters called" << std::endl;
	cols = 0;
	rows = 0;
	size = 0;
	data = nullptr;
}

// This is one of the main constructors that 1- creates a matrix 2-fills it with values (USED in * operations!!)
template  <typename T>
Matrix<T>::Matrix(int rows, int cols, T val)
{
	this -> cols = cols;
	this -> rows = rows;
	size = cols * rows;
	data = new T[size];

	for (int i = 0; i < size; i++)
		data[i] = val;
}






// constructor to generate empty matrix with initializeation.
template <typename T>
Matrix<T>::Matrix(int rows, int cols){
	this -> cols = cols;
	this -> rows = rows;
	size = cols * rows;
	data = new T[size];
}


// Constructor: consists of a two-dimensional array  -----   Main Constructor!
template  <typename T>
Matrix<T>::Matrix(int rows, int cols, const vector<vector<T>>& Array)
{
	this -> cols = cols;
	this -> rows = rows;
	size = cols * rows;
	data = new T[size];
	for (int i = 0; i < cols; i++)
		for (int j = 0; j < rows; j++)
			data[i * rows + j] = Array[j][i];
}

// Constructor: constructed by the class
template  <typename T>
Matrix<T>::Matrix(const Matrix& matrix)
{
	cols = matrix.cols;
	rows = matrix.rows;
	size = cols * rows;
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = matrix.data[i];
}

//-----------------------------------------------------------------------------------

// Matrix output
template  <typename ElemType>
std::ostream& operator<<(std::ostream& os, const Matrix<ElemType>& matrix)
{
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{
			os << matrix.data[j * matrix.rows + i];
			if (j != matrix.cols - 1)
				os << " , ";
		}
		os << ";" << std::endl;
	}
	return os;
}
// Matrix output
template  <typename ElemType>
std::istream& operator>>(std::istream& is, const Matrix<ElemType>& matrix)
{
    std::cout << "Enter the elements of the " << matrix.rows << "X" <<matrix.cols << " matrix: \n";
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{   std::cout << "[" <<i<<"]"<<"["<<j<<"] : ";
			is >> matrix.data[j * matrix.rows + i] ;
			//if (j != matrix.cols - 1 )
				//std::cout << " , " ;
		}
		std::cout << "\t;\n";
	}
	std::cout << "The matrix " << matrix.rows << "X" <<matrix.cols << " is: \n";
	std::cout << matrix;
	return is;
}


// Matrix output
template  <typename ElemType>
std::istream& operator>>(std::istream& is, const Matrix<ElemType>& matrix)
{
    std::cout << "Enter the elements of the " << matrix.rows << "X" <<matrix.cols << " matrix: \n";
	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{   std::cout << "[" <<i<<"]"<<"["<<j<<"] : ";
			is >> matrix.data[j * matrix.rows + i] ;
			//if (j != matrix.cols - 1 )
				//std::cout << " , " ;
		}
		std::cout << "\t;\n";
	}
	std::cout << "The matrix " << matrix.rows << "X" <<matrix.cols << " is: \n";
	std::cout << matrix;
	return is;
}


template  <typename COMPLEXTYPE>
std::ostream& operator<<(std::ostream& os, const Matrix<complex<COMPLEXTYPE>>& matrix)
{
	COMPLEXTYPE a, b;

	for (int i = 0; i < matrix.rows; i++)
	{
		for (int j = 0; j < matrix.cols; j++)
		{
			a = real(matrix.data[j * matrix.rows + i]);
			b = imag(matrix.data[j * matrix.rows + i]);

		 	os << "(" << a;
			(b >= 0 )? os << "+" << b  << "i" : os << b << "i";
			os << ")";

			if (j != matrix.cols - 1)
				os << " , ";
		}
		os << ";" << std::endl;
	}
	return os;
}

//overloading the () operator, makes it easier to access data inside the matrix,  example: m_x(0,1) = 1; this code can change the data inside a matrix or print it.
template  <class T>
T& Matrix<T>::operator()(int row, int col)
{
	return data[col * rows + row];
}


// Operator + overload (matrix addition)
template  <typename ElemType>
Matrix<ElemType>  operator+(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2)
{
	Matrix<ElemType> res(matrix1);
	if (matrix1.cols != matrix2.cols || matrix1.rows != matrix2.rows)
	{
		std::cerr << "Fehler:Die Anzahl der Zeilen oder Spalten ist nicht gleich(matrix::operator+)" << std::endl;
	}
	else if (typeid(matrix1.data).name() != typeid(matrix2.data).name())
	{
		std::cerr << "Fehler:Verschiedene Arten von Matrixdaten (matrix::operator+)" << std::endl;
	}
	else
	{
		for (int i = 0; i < matrix1.size; i++)
			res.data[i] = matrix1.data[i] + matrix2.data[i];
	}
	return res;
}

// Operator-Overload (Matrix Subtraction)
template  <typename ElemType>
Matrix<ElemType>  operator-(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2)
{
	Matrix<ElemType> res(matrix1);
	if (matrix1.cols != matrix2.cols || matrix1.rows != matrix2.rows)
	{
		std::cerr << "Fehler:Die Anzahl der Zeilen oder Spalten ist nicht gleich(matrix::operator-)" << std::endl;
	}
	else if (typeid(matrix1.data).name() != typeid(matrix2.data).name())
	{
		std::cerr << "Fehler:Verschiedene Arten von Matrixdaten(matrix::operator-)" << std::endl;
	}
	else
	{
		for (int i = 0; i < matrix1.size; i++)
			res.data[i] = matrix1.data[i] - matrix2.data[i];
	}
	return res;
}

// Operator* overload (matrix multiplication)
template  <typename ElemType>
Matrix<ElemType>  operator*(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2)
{

	//Re-FIXED
	Matrix<ElemType> res(matrix1.rows, matrix2.cols);

	// i used them for debugging
	// std::cout << "\nresult matrix dimensions:: " << res.rows << " X " << res.cols << std::endl <<std::endl;
	//	std::cout << "data multiplication: \n" << matrix1 << "\n * \n" << matrix2<< " = " <<std::endl;


	 if (matrix1.cols != matrix2.rows) {
	 	std::cerr << "Fehler:Die Spalten von Matrix1 sind nicht gleich den Zeilen von Matrix2, \
und können nicht multipliziert werden.(matrix::operator*)" << std::endl;
	 	return res;
	 }


	if (typeid(matrix1.data).name() != typeid(matrix2.data).name()) {
		std::cerr << "Fehler:Verschiedene Typen von Matrixdaten(matrix::operator*)" << std::endl;
		return res;
	}
	else
	{
        std::cout << "Vector multiplication result: \n";
		for (int i = 0; i < matrix1.rows; i++){
			for (int j = 0; j < matrix2.cols; j++)
			{
				ElemType temp = 0.0;
				for (int k = 0; k < matrix1.cols; k++)
				{
					temp += matrix1.data[k * matrix1.rows + i] * matrix2.data[j * matrix2.rows + k];
					res.data[j * matrix1.rows + i] = temp;
				}
			}

		}

	}
	return res;
}







// Operator* overload (Skalar-Multiplizieren) (1)
template  <typename ElemType>
Matrix<ElemType>  operator*(const Matrix<ElemType>& matrix, const ElemType value)
{   float res = static_cast<double>(value);
	for (int i = 0; i < matrix.size; i++)
		matrix.data[i] *= res;
	return matrix;
}
// Operator* overload (Skalar-Multiplizieren) (2)
template  <typename ElemType>
Matrix<ElemType>  operator*(const ElemType value, const Matrix<ElemType>& matrix)
{
      return matrix * value;
}


// Operator/ overload (Matrix division by number)
template  <typename ElemType>
Matrix<ElemType>  operator/(const Matrix<ElemType>& matrix, const ElemType value)
{   float res = static_cast<double>(value);
	for (int i = 0; i < matrix.size; i++)
		matrix.data[i] /= res;
	return matrix;
}

template <typename T>
Matrix<T>& Matrix<T>::gauss()
{
    std::cout<< "Matrix: \n";
    std::cout << *this;

    int n,i,j,k;
    Matrix<T> mat(rows, cols);
    n = rows;

    for (int i = 0; i < cols; i++)
		for (int j = 0; j < rows; j++)
			{mat(j,i) = data[i * rows + j];}   //put data of matrix in an array;

    for (i=0;i<n;i++)                    //Pivotisation
        for (k=i+1;k<n;k++)
            if (abs(mat(i,i))<abs(mat(k,i)))
                for (j=0;j<=n;j++)
                {
                    double temp=mat(i,j);
                    mat(i,j)=mat (k,j);
                    mat(k,j)=temp;
                }
    std::cout<<"\nDie Matrix nach der Pivotisierung ist:\n";
    std::cout << mat;
    for (i=0;i<n-1;i++)            //do gauss elimination
        for (k=i+1;k<n;k++)
            {
                double t=mat(k,i)/mat(i,i);
                for (j=0;j<=n;j++)
                    mat(k,j)=mat(k,j)-t*mat(i,j);    //elimnate variables that are under pivot
            }
             std::cout<<"\n\nThe matrix after gauss-elimination is:\n";
    std::cout << mat;
    float lsg[n];
    for (i=n-1;i>=0;i--)                //back-substitution
    {                        //lsg is an array whose values correspond to the values of x[1],x[2],x[3]..
        lsg[i] = mat(i,n);                //make the variable to be calculated equal to the rhs of the last equation
        for (j=i+1;j<n;j++)
            if (j!=i)            //then subtract all the lhs values except the coefficient of the variable whose value                                   is being calculated
                lsg[i]=lsg[i]-mat(i,j)*lsg[j];
        lsg[i]=lsg[i]/mat(i,i);            //now finally divide the rhs by the coefficient of the variable to be calculated
    }
    std::cout<<"\nDie Werte der Variablen sind:\n";
    for (i=0;i<n;i++)
    {
        std::cout<<"x["<<i<<"] = " << lsg[i]<<std::endl;
    }


        return *this;



}







#endif
