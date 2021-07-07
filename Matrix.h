#ifndef MATRIX
#define MATRIX
#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <array>
#include <typeinfo>
#include <complex>

using std::complex;

template  <class T>
class Matrix
{
private:
	int cols; // Number of columns
	int rows; // number of rows
	int size; // array size
	T* data; // Element array


public:
	Matrix(); // Default constructor
	Matrix(int rows_, int cols_, T val); // constructor that fills matrix with val
	Matrix(int rows_, int cols_);	//empty constructor initialized with 0
	Matrix(int rows_, int cols_, const std::vector<std::vector<T>>& Array); // A Two-dimensional array to construct a matrix
	Matrix(const Matrix& matrix);                      // Use an existing matrix object to construct a matrix

	~Matrix() { delete[]data; }; //deconstructor (deletes arrays after usage, emptying memory)


	int  getCols() const { return cols; };                 // Get the number of columns
	int  getRows() const { return rows; };                 // Get the number of rows
	int  getSize() const { return rows * cols; };          // Get the size of the array of Matrix


	Matrix<T> gaussian_elemination();
	T* substitute();




    T& operator()(int row, int col);
	template  <typename ElemType>
	friend Matrix<ElemType>   operator +(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2); // Matrix and matrix addition (+ operator overload)
	template  <typename ElemType>
	friend Matrix<ElemType>   operator -(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2); // Matrix and matrix subtraction(-operator overloading)
	template  <typename ElemType>
	friend Matrix<ElemType>   operator *(const Matrix<ElemType>& matrix1, const Matrix<ElemType>& matrix2); // Matrix and matrix multiplication (* operator overload)


	// stream operators
	template <typename COMPLEXTYPE>
	friend std::ostream& operator <<(std::ostream& os, const Matrix<complex<COMPLEXTYPE>>& matrix);
	template  <typename ElemType>
	friend std::ostream& operator <<(std::ostream& os, const Matrix<ElemType>& matrix); // Output matrix to output stream


	//Skalar Funktionen:
	template  <typename ElemType>
	friend Matrix<ElemType>   operator *(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number multiplication * operator overload (1)
	template  <typename ElemType>
	friend Matrix<ElemType>   operator *(const ElemType value, const Matrix<ElemType>& matrix);       // Matrix and number multiplication * operator overload (2)
   	template  <typename ElemType>
	friend Matrix<ElemType>   operator /(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number division/operator overload (1)
	template  <typename ElemType>
	friend Matrix<ElemType>   operator +(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number addition + operator overloading (1)
	template  <typename ElemType>
	friend Matrix<ElemType>   operator -(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number addition + operator overloading (1)
};

// Default constructor
template  <class T>
Matrix<T>::Matrix()
{
	cols = 0;
	rows = 0;
	size = 0;
	data = nullptr;
}

// This is one of the main constructors that 1- creates a matrix 2-fills it with values (USED in * operations!!)
template  <class T>
Matrix<T>::Matrix(int rows, int cols, T val)
{
	this -> cols = cols;
	this -> rows = rows;
	size = cols * rows;
	data = new T[size];
}



// constructor to generate empty matrix with initializeation.
template <class T>
Matrix<T>::Matrix(int rows, int cols){
	this -> cols = cols;
	this -> rows = rows;
	size = cols * rows;
	data = new T[size];
	for (int i = 0; i < size; i++)
		data[i] = 0;
}


// Constructor: consists of a two-dimensional array  -----   Main Constructor!
template  <class T>
Matrix<T>::Matrix(int rows_, int cols_, const std::vector<std::vector<T>>& Array)
{
	cols = cols_;
	rows = rows_;
	size = cols * rows;
	data = new T[size];
	for (int i = 0; i < cols; i++)
		for (int j = 0; j < rows; j++)
			data[i * rows + j] = Array[j][i];
}

// Constructor: constructed by the class
template  <class T>
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
				os << ",";
		}
		os << ";" << std::endl;
	}
	return os;
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






// Operator+ overload (Skalar-Addieren) (1)
template  <typename ElemType>
Matrix<ElemType>  operator+(const Matrix<ElemType>& matrix, const ElemType value)
{   float res = static_cast<double>(value);
	for (int i = 0; i < matrix.size; i++)
		matrix.data[i] += res;
	return matrix;
}


// Operator+ overload (Skalar-Addieren) (2)
 template  <typename ElemType>
 Matrix<ElemType>  operator+(const ElemType value, const Matrix<ElemType>& matrix)
 {
	 return matrix+value;
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
{
	for (int i = 0; i < matrix.size; i++)
		matrix.data[i] /= value;
	return matrix;
}






template <class T>
Matrix<T> Matrix<T>::gaussian_elemination(){
	// need at least cols-1 rows to perform successful elimination
	if (rows < cols -1){ return *this;}

	// very likely not an integer --- conversion to double needed!
	// obviously not working with complex and fractions
	// copy int / float  matrix to double matrix would be much cleaner
	// cleaner way needs different constructor but since we use templates it would interfere with existing one
	std::vector<std::vector<double>> lgs;
	// set expected dimensions so we dont have to use vector.push_back() to add new element
	lgs.reserve(rows);
	for ( int i = 0; i < rows; i++) lgs[i].reserve(cols);

	for (int i = 0; i < rows; i++){
		for (int j = 0; j < cols; j++){
			lgs[i][j] = data[i * rows + j];		// implicit typecast
		}
	}

	double divisor;
	// why is i = 1? --> first row stays same -- start with second :)
	for (int i = 1; i < rows; i++){
		divisor = lgs[i][i-1] / lgs[i-1][i-1];
		// why is j = i-1?  building a diagonal matrix --> value at [i][i] = 0
		for (int j = i-1; j < cols; j++){
			lgs[i][j] /= divisor;
			lgs[i][j] -= lgs[i-1][j];

		}
	}


	Matrix<double> lgs_matrix(rows, cols, lgs);

	return lgs_matrix;
}

template <class T>
T* Matrix<T>::substitute(){

}

template <class T>
T* solve(const Matrix<T>& matrix){
	// should work but still need to test this
	return matrix.gaussian_elemination().substitute();
}

#endif
