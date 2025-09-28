#include "Splayn.h"
#include <iostream>

using namespace std;

Polinomial pow(Polinomial pol, size_t power)
{
	Polinomial Result = pol;
	for (size_t i = 1; i < power; i++) {
		Result *= pol;
	}
	return Result;
}

size_t GetCoeffDeriv(size_t value, size_t power)
{
	size_t Result = value;
	for (size_t i = 1; i < power; i++) {
		Result *= (--value);
	}
	return Result;
}

Polinomial::Polinomial()
{
	this->_Coefficients = VecD(1, 0.0);
}

Polinomial::Polinomial(VecD Coefficients) : _Coefficients(Coefficients)
{}

double Polinomial::GetY(double X)
{
	double Result = 0;
	for (size_t i = 0; i < this->_Coefficients.size(); i++) {
		Result += this->_Coefficients[i] * pow(X, i);
	}
	return Result;
}

size_t Polinomial::GetMaxPower()
{
	return (this->_Coefficients.size() - 1);
}

Polinomial Polinomial::operator*(Polinomial & other)
{
	int size = ((this->_Coefficients.size() - 1) + (other._Coefficients.size() - 1) + 1);
	VecD Result(size, 0.0);
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		for (size_t j = 0; j < other._Coefficients.size(); ++j) {
			Result[i + j] += this->_Coefficients[i] * other._Coefficients[j];
		}
	}
	return Polinomial(Result);
}

Polinomial Polinomial::operator*=(Polinomial other)
{
	int size = ((this->_Coefficients.size() - 1) + (other._Coefficients.size() - 1) + 1);
	VecD Result(size, 0.0);
	for (size_t i = 0; i < this->_Coefficients.size(); ++i) {
		for (size_t j = 0; j < other._Coefficients.size(); ++j) {
			Result[i + j] += this->_Coefficients[i] * other._Coefficients[j];
		}
	}
	this->_Coefficients = Result;
	return *this;
}

Polinomial Polinomial::operator+(Polinomial & other)
{
	VecD Result(this->_Coefficients), ToSum(other._Coefficients);
	(this->_Coefficients.size() > other._Coefficients.size()) ? ToSum.resize(Result.size(), 0.0) : Result.resize(ToSum.size(), 0.0);
	for (size_t i = 0; i < ToSum.size(); i++) {
		Result[i] += ToSum[i];
	}
	return Polinomial(Result);
}

Polinomial Polinomial::operator+=(Polinomial other)
{
	*this = *this + other;
	return *this;
}

Polinomial Polinomial::operator-(Polinomial & other)
{
	VecD Result(this->_Coefficients), ToSub(other._Coefficients);
	(this->_Coefficients.size() > other._Coefficients.size()) ? ToSub.resize(Result.size(), 0.0) : Result.resize(ToSub.size(), 0.0);
	for (size_t i = 0; i < ToSub.size(); i++) {
		Result[i] -= ToSub[i];
	}
	return Polinomial(Result);
}

Polinomial Polinomial::operator*=(double value)
{
	for (size_t i = 0; i < this->_Coefficients.size(); i++) {
		this->_Coefficients[i] *= value;
	}
	return *this;
}

Polinomial Polinomial::operator*(double value)
{
	Polinomial Result = *this;
	for (size_t i = 0; i < this->_Coefficients.size(); i++) {
		Result._Coefficients[i] = this->_Coefficients[i] * value;
	}
	return Result;
}

Polinomial Polinomial::operator/=(double value)
{
	for (size_t i = 0; i < this->_Coefficients.size(); i++) {
		this->_Coefficients[i] /= value;
	}
	return *this;
}

VecD& Polinomial::operator[](int Index)
{
	return this->_Coefficients;
}

void Polinomial::OutPut()
{
	size_t size = this->_Coefficients.size();
	for (size_t i = 0; i < size; i++) {
		if (i != size - 1) {
			std::cout << std::setprecision(15) << this->_Coefficients[i] << "*x^" << i << '+';
		}
		else {
			std::cout << std::setprecision(15) << this->_Coefficients[i] << "*x^" << i;
		}
	}
}

Splayn::Splayn()
{}

void Splayn::SetInitialData(VecD X, VecD Y, VecD Z, size_t power)
{
	this->_X = X;
	this->_Y = Y;
	this->_Z = Z;
	//this->_Polinoms.clear();
	//size_t NumberOfPoints = X.size();
	size_t size = (power + 1) * (Z.size() - 1);
	VecD Diagonal(size, 1.0);
	this->_Diagonal = Diagonal;
	double dx = this->_X[1] - this->_X[0];// равномерная сетка (если будет неравномерная сетка, то нужно одпраить в методе способ заполнения диагонали)
	for (size_t i = 1; i < this->_Diagonal.size(); i += 2) {
		this->_Diagonal[i] = dx;
	}
}

void Splayn::InterpolateFast1D(VecD X, VecD Y)
{
	this->_Polinoms.clear();
	this->_Polinoms.erase(this->_Polinoms.begin(), this->_Polinoms.end());

	this->SetInitialData(X, Y, Y, 1);

	//if (dr == x)
	{
		for (size_t i = 0; i < this->_X.size() - 1; i++)
		{
			this->_Polinoms.push_back(Polinomial(VecD{ Y[i] + (-this->_X[i]) * ((Y[i + 1] - Y[i]) / (this->_X[i + 1] - this->_X[i])), ((Y[i + 1] - Y[i]) / (this->_X[i + 1] - this->_X[i])) }));
			//this->_Polinoms.push_back(Polinomial(VecD{ Y[i] + (-this->_X[i]) * ((Y[i + 1] - Y[i]) / (this->_X[1] - this->_X[0])), ((Y[i + 1] - Y[i]) / (this->_X[1] - this->_X[0])) }));
		}
	}
}

void Splayn::InterpolateFast1D(string fiename)
{
	this->_Polinoms.clear();
	this->_Polinoms.erase(this->_Polinoms.begin(), this->_Polinoms.end());

	ifstream in(fiename); // окрываем файл для чтения
	VecD X;
	VecD Y;

	if (in.is_open())
	{
		double x, y;
		while (in >> x >> y)
		{
			/*	new_points.push_back(Point{ x, y });*/
			X.push_back(x);
			Y.push_back(y);
			//cout << x << "   " << y << endl;
		}
	}
	in.close();
	this->SetInitialData(X, Y, Y, 1);

	//if (dr == x)
	{
		for (size_t i = 0; i < this->_X.size() - 1; i++)
		{
			/*this->_Polinoms.push_back(Polinomial(VecD{ Y[i] + (-this->_X[i]) * ((Y[i + 1] - Y[i]) / (this->_X[1] - this->_X[0])), ((Y[i + 1] - Y[i]) / (this->_X[1] - this->_X[0])) } ));
			//cout << Y[i] + (-this->_X[i]) * ((Y[i + 1] - Y[i]) / (this->_X[1] - this->_X[0])) << "     " << ((Y[i + 1] - Y[i]) / (this->_X[1] - this->_X[0])) << endl;*/

			this->_Polinoms.push_back(Polinomial(VecD{ Y[i] + (-this->_X[i]) * ((Y[i + 1] - Y[i]) / (this->_X[i + 1] - this->_X[i])), ((Y[i + 1] - Y[i]) / (this->_X[i + 1] - this->_X[i])) }));
			//cout << Y[i] + (-this->_X[i]) * ((Y[i + 1] - Y[i]) / (this->_X[i + 1] - this->_X[i])) << "     " << ((Y[i + 1] - Y[i]) / (this->_X[i + 1] - this->_X[i])) << endl;
		}
	}
}

void Splayn::Interpolate(VecD X, VecD Y, size_t power)
{
	this->_X = X;
	this->_Polinoms.clear();
	this->_Polinoms.erase(this->_Polinoms.begin(), this->_Polinoms.end());
	size_t NumberOfPoints = X.size();
	size_t size = (power + 1) * (NumberOfPoints - 1);
	Matrix mat(size, size);
	VecD Free(size, 0.0);
	int num = 0;
	for (size_t i = 0; i < (NumberOfPoints - 1); i++)
	{
		for (size_t j = 0; j < 2; j++)
		{
			for (size_t k = 0; k < (power + 1); k++)
			{
				mat[(2 * i) + j][(i * (power + 1)) + k] = pow((X[i + j] - X[i]), k);
			}
			Free[(2 * i) + j] = Y[i + j];
			num++;
		}
	}

	for (size_t t = 1; t < power; t++)
	{
		for (size_t i = 0; i < (NumberOfPoints - 2); i++)
		{
			for (size_t j = 0; j < 2; j++)
			{
				for (size_t k = t; k < (power + 1); k++)
				{
					mat[2 * (NumberOfPoints - 1) - (NumberOfPoints - 2) + (t * (NumberOfPoints - 2)) + i][k + (j * (power + 1)) + (i * (power + 1))] = GetCoeffDeriv(k, t) * pow(X[i + 1] - X[i + j], (k - t)) * pow(-1, j);
				}
			}
			num++;
		}
	}

	for (size_t t = 2; t <= power / 2 + 1; t++)
	{
		for (size_t k = t; k < (power + 1); k++)
		{
			mat[num][k] = GetCoeffDeriv(k, t) * pow(X[0] - X[0], (k - t));
		}
		num++;
	}

	for (size_t t = 2; t <= (power)-(power / 2); t++)
	{
		for (size_t k = power; k >= t; k--)
		{
			mat[num][size - 1 - power + k] = GetCoeffDeriv(k, t) * pow(X[NumberOfPoints - 1] - X[NumberOfPoints - 2], (k - t));
		}
		num++;
	}

	/*for (size_t i = power - 1; i > 1; i--)
	{
		mat[size - 1][size - power + i] = GetCoeffDeriv(i, 2) * (X[NumberOfPoints - 1] - X[NumberOfPoints - 2]);
	}*/

	if (power == 2) {
		mat[size - 1][size - 1] = 2;
	}

	//////

	/*for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			cout << mat[i][j] << "   ";
		}
		cout << Free[i] << endl;
	}
	cout << endl;*/

	SLAE slae(mat, Free);
	VecD coeffs = slae.GetSolution();
	//VecD coeffs = slae.SolveProgon();

	/////////
	//cout << endl;
	/*for (size_t i = 0; i < size; i++)
	{
		for (size_t j = 0; j < size; j++)
		{
			cout << mat[i][j] << "   ";
		}
		cout << Free[i] << endl;
	}
	cout << endl;*/

	//for (size_t i = 0; i < size; i++)
	//{
	//	cout<< coeffs[i] << "   ";
	//}
	//cout << endl;

	//////////////////

	Polinomial temp;
	for (size_t i = 0; i < NumberOfPoints - 1; i++)
	{
		temp = Polinomial(VecD{ coeffs[(power + 1) * i] });
		for (size_t j = 1; j < (power + 1); j++) {
			temp += pow(Polinomial(VecD{ -X[i], 1 }), j) * coeffs[(power + 1) * i + j];
		}
		this->_Polinoms.push_back(temp);
	}
}

double Splayn::GetY(double X)
{
	for (size_t i = 0; i < this->_X.size(); i++)
	{
		if (X <= this->_X[i])
		{
			if (i != 0) {
				return (this->_Polinoms[i - 1].GetY(X));
			}
			else {
				return (this->_Polinoms[i].GetY(X));
			}
		}
	}
	return (this->_Polinoms[_Polinoms.size() - 1].GetY(X));
}

void Splayn::OutPut()
{
	for (size_t i = 0; i < this->_Polinoms.size(); i++)
	{
		//std::cout << i << '\n';
		std::cout << '{';
		_Polinoms[i].OutPut();
		std::cout << ',' << this->_X[i] << "<=x<=" << this->_X[i + 1] << "}," << '\n';
	}
}
