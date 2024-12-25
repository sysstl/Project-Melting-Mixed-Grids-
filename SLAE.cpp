#include "SLAE.h"

SLAE::SLAE(size_t Rows, size_t Cols) : _Matrix(Rows, Cols), _FreeCol(VecD(Cols, 0.0))
{}

SLAE::SLAE(VecVecD& Values, VecD& FreeCol) : _Matrix(Values), _FreeCol(FreeCol)
{}

SLAE::SLAE(Matrix& Matrix, VecD& FreeCol) : _Matrix(Matrix), _FreeCol(FreeCol)
{}

SLAE::SLAE(VecD& Diagonal, VecD& FreeCol) : _Diagonal(Diagonal), _FreeCol(FreeCol)
{}

VecD SLAE::GetSolution()
{
	SLAE temp = *this;
	size_t size = temp._Matrix.GetRows(), size2 = temp._Matrix.GetCols();
	if ((size != size2) || size == 0 || size2 == 0) {
		throw std::exception("Imposible to solve this slae");
	}

	VecD result(size, 0.0);
	temp.MakeItTriangle();
	double Solution;
	for (int64_t i = size - 1; i >= 0; --i)
	{
		Solution = temp._FreeCol[i];
		for (size_t j = size - 1; j > i; --j)
		{
			Solution -= result[j] * temp._Matrix[i][j];
		}
		Solution /= temp._Matrix[i][i];
		if (abs(Solution) < EPS) {
			result[i] = 0;
		}
		else {
			result[i] = Solution;
		}
	}
	return result;
}

VecD SLAE::SolveProgon()
{
	size_t Size = this->_Matrix.GetRows();
	VecD alpha(Size, 0), beta(Size, 0);

	alpha[1] = (-this->_Matrix[0][1]) / this->_Matrix[0][0];
	beta[1] = this->_FreeCol[0] / this->_Matrix[0][0];

	for (size_t i = 2; i < Size; ++i)
	{
		alpha[i] = (-this->_Matrix[i - 1][i]) / (this->_Matrix[i - 1][i - 2] * alpha[i - 1] + this->_Matrix[i - 1][i - 1]);
		beta[i] = (this->_FreeCol[i - 1] - this->_Matrix[i - 1][i - 2] * beta[i - 1]) / (this->_Matrix[i - 1][i - 2] * alpha[i - 1] + this->_Matrix[i - 1][i - 1]);
	}

	VecD Result(Size, 0);

	Result[Size - 1] = (this->_FreeCol[Size - 1] - this->_Matrix[Size - 1][Size - 2] * beta[Size - 1]) / (this->_Matrix[Size - 1][Size - 1] + this->_Matrix[Size - 1][Size - 2] * alpha[Size - 1]);

	for (size_t i = Size - 1; i > 0; --i) {
		Result[i - 1] = alpha[i] * Result[i] + beta[i];
	}

	return Result;
}

VecD SLAE::SolveProgonLinearInterpolate()
{
	for (size_t i = 1; i < this->_Diagonal.size() - 2; i += 2) {
		this->_FreeCol[i] = (this->_FreeCol[i] - this->_FreeCol[i - 1]) / (this->_Diagonal[i]);
	}
	this->_FreeCol[this->_Diagonal.size() - 1] = (this->_FreeCol[this->_Diagonal.size() - 1] - this->_FreeCol[this->_Diagonal.size() - 2]) / (this->_Diagonal[this->_Diagonal.size() - 1]);
	return this->_FreeCol;
}

Matrix SLAE::GetMatrix()
{
	return this->_Matrix;
}

VecD SLAE::GetFreeCol()
{
	return this->_FreeCol;
}

void SLAE::SwitchRows(size_t Which, size_t Where)
{
	this->_Matrix.SwitchRows(Which, Where);
	double temp = this->_FreeCol[Where];
	this->_FreeCol[Where] = this->_FreeCol[Which];
	this->_FreeCol[Which] = temp;
}

void SLAE::MakeItTriangle()
{
	size_t Rows = this->_FreeCol.size();
	double ValueToMulty;
	for (size_t i = 0; i < Rows; ++i)
	{
		if (abs(this->_Matrix[i][i]) <= EPS) {
			this->SwitchRows(i, this->_Matrix.GetIndexOfBiggestValue(i, i));
			if (abs(this->_Matrix[i][i]) <= EPS) {
				throw std::exception("Determinant = 0");
			}
		}
		for (size_t j = i + 1; j < Rows; ++j)
		{
			ValueToMulty = this->_Matrix[j][i] / this->_Matrix[i][i];
			for (size_t k = i; k < Rows; ++k) {
				this->_Matrix[j][k] -= (ValueToMulty * this->_Matrix[i][k]);
			}
			this->_FreeCol[j] -= (ValueToMulty * this->_FreeCol[i]);
		}
	}
}
