#pragma once
#include <vector>
#include <iostream>

using namespace std;

typedef std::vector<double> VecD;
typedef std::vector<std::vector<double>> VecVecD;

#define EPS ((double)(pow(10, -14)))

class Matrix
{
public:
	Matrix(size_t Rows = 0, size_t Cols = 0);
	Matrix(VecVecD& values);
	Matrix(VecD values, bool gorizontal = true);
	void Resize(size_t Rows, size_t Cols);
	size_t GetRows();
	size_t GetCols();
	//double GetDeterminant();
	Matrix GetTransposed();
	Matrix GetReverse();
	void SwitchRows(size_t Which, size_t Where);
	size_t GetIndexOfBiggestValue(size_t Row, size_t Col);
	double GetNorma();
	void PrintMatrix(ostream& out);

	//---------------------------//

	Matrix operator-(const Matrix& other);
	Matrix operator*(const Matrix& other);
	Matrix operator=(double& value);
	friend ostream& operator <<(ostream& out, const Matrix& mtr);
	VecD& operator[](size_t Index);

private:

	size_t MakeItTriangle();
	//---------------------------//
	VecVecD _Values;
	size_t _Rows;
	size_t _Cols;
};