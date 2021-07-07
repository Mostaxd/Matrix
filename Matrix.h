
#ifndef MATRIX
#define MATRIX

#include <iostream>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <typeinfo>
#include <complex>

using std::complex;

template  <typename T>
class Matrix
{
private:
	int cols; // Number of columns
	int rows; // number of rows
	int size; // array size
	T* data; // Element array
	// member of child class lgs
	std::vector<T> x;


public:
	Matrix(); // Default constructor
	Matrix(int rows_, int cols_, T val); // constructor that fills matrix with val
	Matrix(int rows_, int cols_);	//empty constructor initialized with 0
	Matrix(int rows_, int cols_, const std::vector<std::vector<T>>& Array); // A Two-dimensional array to construct a matrix
	Matrix(const Matrix& matrix);                      // Use an existing matrix object to construct a matrix

	~Matrix() { delete[]data; }; //deconstructor (deletes arrays after usage, manually emptying heap memory)


	int  getCols() const { return cols; };                 // Get the number of columns
	int  getRows() const { return rows; };                 // Get the number of rows
	int  getSize() const { return rows * cols; };          // Get the size of the array of Matrix


	Matrix<T> gaussian_elemination();
	Matrix<T>& substitute();
	Matrix<T> solve();


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

	//much much cleaner and easyer to read -- this bugged me a long time
	// we need to simplify every operator overload like this!
	Matrix<T>&   operator/(const T value);       // Matrix and number division/operator overload (1)
	// template  <typename ElemType>
	// friend Matrix<ElemType>   operator +(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number addition + operator overloading (1)
	// template  <typename ElemType>
	// friend Matrix<ElemType>   operator -(const Matrix<ElemType>& matrix, const ElemType value);       // Matrix and number addition + operator overloading (1)
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





// this is adding the same scalar to every element in a matrix
// which is not adding an adding a scalar to a matrix operation
// to do that you would need to multiply a identity matrix with the scalar and then do matrix + matrix addition
// which would only permit adding a scalar to square matrices
// nevertheless i strongly think this wasn't requested and we shouldn't implement it
// // Operator+ overload (Skalar-Addieren)
// template  <typename ElemType>
// Matrix<ElemType>  operator+(const Matrix<ElemType>& matrix, const ElemType value)
// {   float res = static_cast<double>(value);
// 	for (int i = 0; i < matrix.size; i++)
// 		matrix.data[i] += res;
// 	return matrix;
// }

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
template  <typename T>
Matrix<T>& Matrix<T>::operator/(const T value)
{
	for (int i = 0; i < size; i++)
		data[i] /= value;
	return this;
}






template <typename T>
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
			lgs[i][j] = data[j * rows + i];		// implicit typecast
		}
	}

	double divisor;
	for (int latest_row = 0; latest_row < rows; latest_row++){		//latest row == latest column <-- diagonal
		for (int i = latest_row + 1; i < rows; i++){
			divisor = lgs[i][latest_row] / lgs[latest_row][latest_row];

			for (int j = latest_row; j < cols; j++){
				lgs[i][j] /= divisor;
				lgs[i][j] -= lgs[latest_row][j];
			}
		}

	}

	Matrix<double> lgs_matrix(rows, cols, lgs);
	return lgs_matrix;
}

// solve for x1 x2 x3 ... xn
template <typename T>
Matrix<T>& Matrix<T>::substitute(){
	int i, j;
	i = 0;
 	// std::vector<T> x;
	x.reserve(cols-1);
	while ( cols-1 >= i++) x.push_back(1);
	rows = cols - 1;
	T sum = 0;

	// why i = rows - 1?
	// counting starts at 0 -- rows is number of rows  --> if there is one row start at rows-1
	for (i = rows -1 ; i >= 0; i--){
		for (j = cols - 2; j >= i; j--){
			sum += (data[j * rows + i] * x[j]);

			// std::cout << "factor to mutliply cell value with (x[j]): " << x[j] << " at j: " << j << std::endl;
			// std::cout << "cell value: " << data[j * rows + i] << " sum: " << sum << std::endl;

		}
		// std::cout << "rows: " << cols -1 << std::endl;
		// cols - 1 = number of interesting rows  -- using cols-1 just in case rows is not set properly
		// std::cout << "cell value: " << data[(cols-1) * (cols-1) - (i*(cols-1))] << std::endl;
		// std::cout << "cell value: " << data[(rows)*i + i] << std::endl;

		x[i] = data[ (cols - 1) * rows + i] / sum;
		// std::cout << "value at last column: " << data[ (cols - 1) * rows + i] << std::endl;
		// std::cout << "value of i: " << i << std::endl;
		// std::cout << "solution of x[i]: " << x[i] << std::endl;
		// std::cout << "reset sum. sum = 0 -------------" << std::endl;
		sum = 0;
	}
	// this -> x = x;
	// just testing...
	std::cout << "lgs solutions: x1 " << x[0] << " , " << x[1] << " , " << x[2] << std::endl;
	return *this;
}



template <typename T>
Matrix<T> Matrix<T>::solve(){
	// nice chaining :)
	return this->gaussian_elemination().substitute();
}



////////////////////////////////////////////////////////////////////////////////
// new class: Lgs
////////////////////////////////////////////////////////////////////////////////
template <typename T>
class Lgs : public Matrix<T>{
private:
	T* x;
public:
	Lgs();
	Lgs<T>& test();
};


////////////////////////////////////////////////////////////////////////////////
/// lgs constructors
template <typename T>
Lgs<T>::Lgs(){
	x = nullptr;

}


#endif
