#include "Matrix.h"
#include <iostream>

using namespace std;

Matrix::Matrix(size_t Rows, size_t Cols) : _Rows(Rows), _Cols(Cols)
{
	for (size_t i = 0; i < Rows; ++i) {
		this->_Values.push_back(VecD(Cols, 0.0));
	}
}

Matrix::Matrix(VecVecD& values) : _Values(values)
{
	this->_Rows = values.size();
	if (this->_Rows != 0) {
		this->_Cols = values[0].size();
	}
	else {
		this->_Cols = 0;
	}
}

Matrix::Matrix(VecD values, bool gorizontal)
{
	if (gorizontal)
	{
		this->_Rows = 1;
		this->_Cols = values.size();
		this->_Values.push_back(values);
	}
	else
	{
		this->_Cols = 1;
		this->_Rows = values.size();
		for (size_t i = 0; i < this->_Rows; i++) {
			this->_Values.push_back(VecD(1, values[i]));
		}
	}
}

void Matrix::Resize(size_t Rows, size_t Cols)
{
	this->_Rows = Rows;
	this->_Cols = Cols;
	for (size_t i = 0; i < Rows; ++i) {
		this->_Values.push_back(VecD(Cols, 0.0));
	}
}

size_t Matrix::GetRows()
{
	return this->_Rows;
}

size_t Matrix::GetCols()
{
	return this->_Cols;
}

//double Matrix::GetDeterminant()
//{
//	size_t Rows = this->_Rows, Cols = this->_Cols, NumberOfChanges = 0;
//	double ValueToMulty;
//	Matrix temp(*this);
//	if (Rows != Cols) {
//		throw std::exception("Imposible to calc determinant");
//	}
//	try
//	{
//		NumberOfChanges = temp.MakeItTriangle();
//		if (abs(temp._Values[Rows - 1][Cols - 1]) <= EPS)
//			return 0;
//		double Result = 1;
//		for (int i = 0; i < Rows; ++i)
//			Result *= temp._Values[i][i];
//		return (Result * pow(-1, NumberOfChanges));
//	}
//	catch (std::exception& exept)
//	{
//		return 0;
//	}
//}

Matrix Matrix::GetTransposed()
{
	Matrix Result(this->_Cols, this->_Rows);
	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[j][i] = this->_Values[i][j];
		}
	}
	return Result;
}

Matrix Matrix::GetReverse()
{
	int size = this->_Values.size();

	VecVecD mat(size, VecD(2 * size));
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			mat[i][j] = this->_Values[i][j];
		}
		mat[i][i + size] = 1.0;
	}

	for (size_t i = 0; i < size; i++) {
		double Divisor = mat[i][i];
		for (size_t j = 0; j < 2 * size; j++) {
			mat[i][j] /= Divisor;
		}
		for (size_t k = 0; k < size; k++) {
			if (k != i) {
				double Multiplier = mat[k][i];
				for (size_t j = 0; j < 2 * size; j++) {
					mat[k][j] -= Multiplier * mat[i][j];
				}
			}
		}
	}
	VecVecD Reveresed(size, std::vector<double>(size));
	for (size_t i = 0; i < size; i++) {
		for (size_t j = 0; j < size; j++) {
			Reveresed[i][j] = mat[i][j + size];
		}
	}
	return Matrix(Reveresed);
}

void Matrix::SwitchRows(size_t Which, size_t Where)
{
	if (Where >= this->_Values.size() || Where >= this->_Values.size()) {
		throw "Invalid Index";
	}
	VecD temp = this->_Values[Where];
	this->_Values[Where] = this->_Values[Which];
	this->_Values[Which] = temp;
}

size_t Matrix::GetIndexOfBiggestValue(size_t Row, size_t Col)
{
	double BiggestValue = 0;
	size_t IndexOfRow = Row;
	for (size_t i = Row + 1; i < this->_Rows; ++i)
	{
		if (abs(this->_Values[i][Col]) > BiggestValue)
		{
			BiggestValue = abs(this->_Values[i][Col]);
			IndexOfRow = i;
		}
	}
	return IndexOfRow;
}

double Matrix::GetNorma()
{
	double result = 0;
	for (size_t i = 0; i < this->_Rows; i++) {
		for (size_t j = 0; j < this->_Cols; j++) {
			result += (this->_Values[i][j] * this->_Values[i][j]);
		}
	}
	return sqrt(result);
}

Matrix Matrix::operator-(const Matrix& other)
{
	if (this->_Rows != other._Rows || this->_Cols != other._Cols) {
		throw std::exception("Imposible to subtract these matrices");
	}

	Matrix Result(this->_Rows, this->_Cols);
	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[i][j] = this->_Values[i][j] - other._Values[i][j];
		}
	}
	return Result;
}

Matrix Matrix::operator*(const Matrix& other)
{
	if (this->_Cols != other._Rows) {
		throw std::exception("Imposible to multiply these matrices");
	}

	Matrix Result(this->_Rows, other._Cols);
	for (size_t i = 0; i < Result._Rows; ++i) {
		for (size_t j = 0; j < Result._Cols; ++j) {
			for (size_t k = 0; k < this->_Cols; ++k) {
				Result._Values[i][j] += this->_Values[i][k] * other._Values[k][j];
			}
		}
	}
	return Result;
}

Matrix Matrix::operator=(double& value) // прописать остальные подобные операторы
{
	Matrix Result(this->_Rows, this->_Cols);
	for (size_t i = 0; i < this->_Rows; ++i) {
		for (size_t j = 0; j < this->_Cols; ++j) {
			Result._Values[i][j] = value;
		}
	}
	return Result;
}


VecD& Matrix::operator[](size_t Index)
{
	if (Index >= this->_Rows) {
		std::cout << Index << "  " << this->_Rows << endl;
		throw std::exception("Invalid index");
	}
	return this->_Values[Index];
}

/*void Matrix::PrintMatrix(ostream& out)
{

}*/

ostream& operator <<(ostream& out, const Matrix& mtr)
{
	//int Row = mtr.GetRows();
	//int Col = mtr.GetCols();
	for (int i = 0; i < mtr._Rows; i++)
	{
		for (int j = 0; j < mtr._Cols; j++)
		{
			//mtr._Values;
			out << mtr._Values[i][j] << "\t";
		}
		out << "\n";
	}
	return out;
}//*/

size_t Matrix::MakeItTriangle()
{
	size_t Rows = this->_Rows, Cols = this->_Cols, NumberOfChanges = 0;
	double ValueToMulty;
	for (size_t i = 0; i < Rows; ++i)
	{
		if (abs(this->_Values[i][i]) <= EPS)
		{
			this->SwitchRows(i, this->GetIndexOfBiggestValue(i, i));
			if (abs(this->_Values[i][i]) <= EPS) {
				throw std::exception("Det = 0");
			}
			NumberOfChanges++;
		}
		for (size_t j = i + 1; j < Rows; ++j)
		{
			ValueToMulty = this->_Values[j][i] / this->_Values[i][i];
			for (int k = i; k < Cols; k++) {
				this->_Values[j][k] -= (ValueToMulty * this->_Values[i][k]);
			}
		}
	}
	return NumberOfChanges;
}