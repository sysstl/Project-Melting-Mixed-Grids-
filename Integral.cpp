#include "Integral.h"
#include <cmath>
#include <iostream>

Integral::Integral(double(*f_)(double), double min_x_, double max_x_)
{
	f = f_;
	min_x = min_x_;
	max_x = max_x_;
	answer = 0;
}

Integral::Integral(double(*f2_)(double, double), double min_x_, double max_x_)
{
	f2 = f2_;
	min_x = min_x_;
	max_x = max_x_;
	answer = 0;
}

double Integral::CentralRect(int n)
{
	answer = 0;
	h = (max_x - min_x) / n;
	for (int i = 1; i < n + 1; ++i)
		answer += h * f(min_x + h * (i - 0.5));
	gamma = 2;
	return answer;
}

double Integral::Trapezoid(int n)
{
	answer = 0;
	h = (max_x - min_x) / n;
	for (int i = 1; i < n + 1; ++i)
		answer += h * (f(min_x + h * i) + f(min_x + h * (i - 1))) / 2;
	gamma = 2;
	return answer;
}

double Integral::Parabola(int n)
{
	double pi = acos(-1);
	answer = 0;
	h = (max_x - min_x) / n;
	for (int i = 1; i < n + 1; ++i)
	{
		//std::cout << "{" << min_x + h * i << ",";
		answer += h * (f(min_x + h * i) + f(min_x + h * (i - 1)) + 4 * f(min_x + h * (i - 0.5))) / 6;
		//std::cout << 2/sqrt(pi)*answer << "},";
	}
	gamma = 4;
	return answer;
}

double Integral::CentralRect(double eps)
{
	double R0 = 1;
	double R1 = 0;
	double r = 2;
	double prev_h;
	double prev_ans;
	int n = 3;
	CentralRect(n);

	/*std::cout << "LiastPlot[{\n";
	int i = 25;*/
	do
	{
		R0 = R1;
		prev_h = (max_x - min_x) / n;
		prev_ans = answer;
		n *= r;
		CentralRect(n);
		R1 = (answer - prev_ans) / (pow(r, gamma) - 1);
		if (R0 / R1 < 0)
		{
			solved = false;
			return 0;
		}
		/*std::cout << "{" << log(n) << "," << log(abs(R1)) << "},";
		i--;*/
	} while (abs(R1) >= eps);
	/*while (i > 0);
	std::cout << "}]";*/
	solved = true;
	return answer;
}

double Integral::Trapezoid(double eps)
{
	double R0 = -1;
	double R1 = 0;
	double r = 2;
	double prev_h;
	double prev_ans;
	int n = 3;
	Trapezoid(n);

	//std::cout << "LiastPlot[{\n";
	do
	{
		R0 = R1;
		prev_h = (max_x - min_x) / n;
		prev_ans = answer;
		n *= r;
		Trapezoid(n);
		R1 = (answer - prev_ans) / (pow(r, gamma) - 1);
		if (R0 / R1 < 0)
		{
			solved = false;
			return 0;
		}
		//std::cout << "{" << log(n) << "," << log(abs(R1)) << "},";
	} while (abs(R1) >= eps);
	//std::cout << "}]";

	solved = true;
	return answer;
}

double Integral::Parabola(double eps)
{
	double R0 = -1;
	double R1 = 0;
	double r = 3;
	double prev_h;
	double prev_ans;
	int n = 2;
	Parabola(n);

	//std::cout << "LiastPlot[{\n";
	do
	{
		R0 = R1;
		prev_h = (max_x - min_x) / n;
		prev_ans = answer;
		n *= r;
		Parabola(n);
		R1 = (answer - prev_ans) / (pow(r, gamma) - 1);
		if (R0 / R1 < 0)
		{
			solved = false;
			return 0;
		}
		//std::cout << "{" << log(n) << "," << log(abs(R1)) << "},";
	} while (abs(R1) >= eps);

	//std::cout << "}]";
	answer += R1;
	solved = true;
	return answer;
}

double Integral::GaussLegendre_IntervalVariety(int n)
{
	switch (n)
	{
	case 2:
		solved = true;
		return ((this->max_x - this->min_x) / 2) * (f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *(1. / sqrt(3))) + f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) *(1. / sqrt(3))));
	case 3:
		solved = true;
		return ((this->max_x - this->min_x) / 2) * ((5./9.)*f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *(sqrt(3.) / sqrt(5))) + (5./9.)*f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) *(sqrt(3.) / sqrt(5)))  +  (8.)/(9.) * f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) * 0.));
	case 4:
		solved = true;
		return ((this->max_x - this->min_x) / 2) * ((1.0 / 2 - sqrt(5.0 / 6) / 6.0) *  f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *(sqrt((2.0 * sqrt(30) + 15.0) / 35.0))) + (1.0 / 2 - sqrt(5.0 / 6) / 6.0) * f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) *(sqrt((2.0 * sqrt(30) + 15.0) / 35.0))) +
			((sqrt(30) + 18) / 36)* f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *(sqrt((15 - 2 * sqrt(30)) / 35))) + ((sqrt(30) + 18) / 36) * f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) *(sqrt((15 - 2 * sqrt(30)) / 35))));
	case 5:
		solved = true;
		return ((this->max_x - this->min_x) / 2) * (((322.0 - 13.0 * sqrt(70)) / 900.0) * f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *(1.0 / 3 * sqrt(1.0 / 7 * (35.0 + 2.0 * sqrt(70))))) + ((322.0 - 13.0 * sqrt(70)) / 900.0) * f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) *(1.0 / 3 * sqrt(1.0 / 7 * (35.0 + 2.0 * sqrt(70))))) +
			((13.0 * sqrt(70) + 322.0) / 900.0) * f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *((1.0 / 3) * sqrt((1.0 / 7) * (35.0 - 2.0 * sqrt(70))))) + ((13.0 * sqrt(70) + 322.0) / 900.0) * f(((this->max_x + this->min_x) / 2) + ((this->max_x - this->min_x) / 2) *((1.0 / 3) * sqrt((1.0 / 7) * (35.0 - 2.0 * sqrt(70))))) +
			(128.0 / 225.) * f(((this->max_x + this->min_x) / 2) - ((this->max_x - this->min_x) / 2) *(0.)));
	default:
		solved = false;
		return 0;
	}
}

double Integral::GaussLegendre(int n)
{
	switch (n)
	{
	case 2:
		solved = true;
		return f(-sqrt(1.0 / 3)) + f(sqrt(1.0 / 3));
	case 3:
		solved = true;
		return 5.0 / 9.0 * (f(-sqrt(3.0 / 5)) + f(sqrt(3.0 / 5))) + 8.0 / 9 * f(0);
	case 4:
		solved = true;
		return (1.0 / 2 - sqrt(5.0 / 6) / 6.0) * (f(-sqrt((2.0 * sqrt(30) + 15.0) / 35.0)) + f(sqrt((2.0 * sqrt(30) + 15.0) / 35.0))) +
			(sqrt(30) + 18) / 36 * (f(-sqrt((15 - 2 * sqrt(30)) / 35)) + f(sqrt((15 - 2 * sqrt(30)) / 35)));
	case 5:
		solved = true;
		return (322.0 - 13.0 * sqrt(70)) / 900.0 * (f(1.0 / 3 * sqrt(1.0 / 7 * (35.0 + 2.0 * sqrt(70)))) + f(-1.0 / 3 * sqrt(1.0 / 7 * (35.0 + 2.0 * sqrt(70))))) +
			(13.0 * sqrt(70) + 322.0) / 900.0 * (f((1.0 / 3) * sqrt((1.0 / 7) * (35.0 - 2.0 * sqrt(70)))) + f(-(1.0 / 3) * sqrt((1.0 / 7) * (35.0 - 2.0 * sqrt(70))))) +
			128.0 / 225 * f(0);
	default:
		solved = false;
		return 0;
	}
}

double Integral::GaussChebishev(int n)
{
	double ans = 0;
	double PI = acos(-1);
	for (int i = 0; i < n; ++i)
		ans += PI / n * f(cos((2.0*(i + 1.0) - 1.0)*PI / (2.0*n)));
	solved = true;
	return ans;
}

double Integral::GaussLiager(int n)
{
	switch (n)
	{
	case 2:
		solved = true;
		return 0.853553 * f(0.585786) + 0.146447 * f(3.414214);
	case 3:
		solved = true;
		return f(0.415775) * 0.711093 + f(2.294280) * 0.278518 + f(6.289945) * 0.0103893;
	case 4:
		solved = true;
		return f(0.322548) * 0.603154 + f(1.745761) * 0.357419 + f(4.536620) * 0.0388879 + f(9.395071) * 0.000539295;
	case 5:
		solved = true;
		return f(0.263560) * 0.521756 + f(1.413403) * 0.398667 + f(3.596426) * 0.0759424 + f(7.085810) * 0.00361176 + f(12.640801) * 0.0000233700;
	default:
		solved = false;
		return 0;
	}
}

double Integral::GaussErmit(int n)
{
	switch (n)
	{
	case 2:
		solved = true;
		return 0.886227 * (f(0.707107) + f(-0.707107));
	case 3:
		solved = true;
		return f(0) * 1.181636 + (f(1.224745) + f(-1.224745)) * 0.295409;
	case 4:
		solved = true;
		return (f(0.524648) + f(-0.524648)) * 0.804914 + (f(1.650680) + f(-1.650680)) * 0.0813128;
	case 5:
		solved = true;
		return f(0) * 0.945309 + (f(0.958572) + f(-0.958572)) * 0.393619 + (f(2.020183) + f(-2.020183)) * 0.0199532;
	default:
		solved = false;
		return 0;
	}
}

double Integral::GaussLegendre(int n, double i)
{
	switch (n)
	{
	case 2:
		solved = true;
		return f2(-sqrt(1.0 / 3), i) + f2(sqrt(1.0 / 3), i);
	case 3:
		solved = true;
		return 5.0 / 9.0 * (f2(-sqrt(3.0 / 5), i) + f2(sqrt(3.0 / 5), i)) + 8.0 / 9 * f2(0, i);
	case 4:
		solved = true;
		return (1.0 / 2 - sqrt(5.0 / 6) / 6.0) * (f2(-sqrt((2.0 * sqrt(30) + 15.0) / 35.0), i) + f2(sqrt((2.0 * sqrt(30) + 15.0) / 35.0), i)) +
			(sqrt(30) + 18) / 36 * (f2(-sqrt((15 - 2 * sqrt(30)) / 35), i) + f2(sqrt((15 - 2 * sqrt(30)) / 35), i));
	case 5:
		solved = true;
		return (322.0 - 13.0 * sqrt(70)) / 900.0 * (f2(1.0 / 3.0 * sqrt(1.0 / 7.0 * (35.0 + 2.0 * sqrt(70))), i) + f2(-1.0 / 3.0 * sqrt(1.0 / 7.0 * (35.0 + 2.0 * sqrt(70))), i)) +
			(13.0 * sqrt(70) + 322.0) / 900.0 * (f2((1.0 / 3.0) * sqrt((1.0 / 7.0) * (35.0 - 2.0 * sqrt(70))), i) + f2(-(1.0 / 3.0) * sqrt((1.0 / 7.0) * (35.0 - 2.0 * sqrt(70))), i)) +
			128.0 / 225.0 * f2(0, i);
	default:
		solved = false;
		return 0;
	}
}

double Integral::GaussLiager(int n, double i)
{
	switch (n)
	{
	case 2:
		solved = true;
		return 0.853553 * f2(0.585786, i) + 0.146447 * f2(3.414214, i);
	case 3:
		solved = true;
		return f2(0.415775, i) * 0.711093 + f2(2.294280, i) * 0.278518 + f2(6.289945, i) * 0.0103893;
	case 4:
		solved = true;
		return f2(0.322548, i) * 0.603154 + f2(1.745761, i) * 0.357419 + f2(4.536620, i) * 0.0388879 + f2(9.395071, i) * 0.000539295;
	case 5:
		solved = true;
		return f2(0.263560, i) * 0.521756 + f2(1.413403, i) * 0.398667 + f2(3.596426, i) * 0.0759424 + f2(7.085810, i) * 0.00361176 + f2(12.640801, i) * 0.0000233700;
	default:
		solved = false;
		return 0;
	}
}

double Integral::GaussChebishev(int n, double i)
{
	double ans = 0;
	double PI = acos(-1);
	for (int j = 0; j < n; ++j)
		ans += PI / n * f2(cos((2.0 * (j + 1.0) - 1.0) * PI / (2.0 * n)), i);
	solved = true;
	return ans;
}

