#pragma once

class Integral
{
private:
	double(*f)(double);
	double(*f2)(double, double);
	double min_x;
	double max_x;
	double answer;
	double gamma;
	double h;
	bool solved = false;

public:
	Integral(double(*f_)(double), double min_x_, double max_x_);
	Integral(double(*f2_)(double, double), double min_x_, double max_x_);

	double CentralRect(int n);
	double Trapezoid(int n);
	double Parabola(int n);

	double CentralRect(double eps);
	double Trapezoid(double eps);
	double Parabola(double eps);

	double GaussLegendre_IntervalVariety(int n);
	double GaussLegendre(int n);
	double GaussChebishev(int n);
	double GaussLiager(int n);
	double GaussErmit(int n);

	double GaussLegendre(int n, double i);
	double GaussLiager(int n, double i);
	double GaussChebishev(int n, double i);

	double GetH() { return h; };
	bool IfSolved() { return solved; };
};
