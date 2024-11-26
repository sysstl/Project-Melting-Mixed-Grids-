#pragma once
#include <vector>
#include <iostream>
#include <iomanip>
#include <fstream>
#include "SLAE.h"

using namespace std;

enum Direction { x, y, z };

struct Point
{
	Point(double x, double y) : x{ x }, y{ y } {}
	double x;
	double y;
};

class Polinomial;

const double PI = acos(-1);

Polinomial pow(Polinomial pol, size_t power);

size_t GetCoeffDeriv(size_t value, size_t power);

class Polinomial
{
public:
	Polinomial();
	Polinomial(VecD Coefficients);
	double GetY(double X);
	double GetYFast(double X);
	size_t GetMaxPower();

	Polinomial operator*(Polinomial& other);
	Polinomial operator*=(Polinomial other);
	Polinomial operator+(Polinomial& other);
	Polinomial operator+=(Polinomial other);
	Polinomial operator-(Polinomial& other);
	Polinomial operator*=(double value);
	Polinomial operator*(double value);
	Polinomial operator/=(double value);
	VecD& operator[](int Index);
	void OutPut();

private:

	VecD _Coefficients;
};

class Splayn
{
public:
	Splayn();
	void SetInitialData(VecD X, VecD Y, VecD Z, size_t power);
	void InterpolateFast(size_t power, double*** Tei, int fix_x, int fix_y, int fix_z, Direction dr);
	void Interpolate(VecD X, VecD Y, size_t power);
	double GetY(double X);
	//double GetYFast(double X);
	void OutPut();
	/*void OutPutFast(double X);*/
	void Calculation_Interpolation(string fiename);
	void InterpolateFast1D(VecD X, VecD Y);
	void InterpolateFast1D(string fiename);
	void Calculation_Interpolation(double*** Tei, int Nz_heat, double dz_heat, int fix_x, int fix_y);
	void Calculation_InterpolationFast(double*** Tei, int Nz_heat, double dz_heat, int fix_x, int fix_y);
	void Calculation_InterpolationFast(string fiename);

private:
	VecD _X, _Y, _Z;
	std::vector<Polinomial> _Polinoms;
	Matrix mat;
	VecD _Diagonal;
};