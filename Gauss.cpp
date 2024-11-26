#include "Gauss.h"
#include <random>
#include <time.h>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;
Gauss::Gauss(int n_)
{
	n = n_;
	srand(time(NULL));
	A = vector<vector<double>>(n_, vector<double>(n_));
	B = vector<double>(n_);
	for (int i = 0; i < n_; ++i)
		for (int j = 0; j < n_; ++j)
			A[i][j] = (rand() % 2001 - 1000);
	for (int i = 0; i < n_; ++i)
		B[i] = (rand() % 2001 - 1000);
}

Gauss::Gauss(vector< vector<double> > & A_, vector<double> & B_)
{
	if (A_.size() != A_[0].size())
	{
		errFlag = true;
		err = "A is not square\n";
	}
	else if (A_.size() != B_.size())
	{
		errFlag = true;
		err = "A and B have different scale\n";
	}
	else
	{
		errFlag = false;
		n = A_.size();
		A = A_;
		B = B_;
	}
}

ostream& operator <<(ostream & out, const Gauss & sistem)
{
	if (sistem.errFlag)
		out << sistem.err;
	else if (sistem.solved)
	{
		for (int i = 0; i < sistem.n; ++i)
		{
			out << 'x' << i + 1 << " = " << setprecision(20) << sistem.B[i] << "\n";
		}
	}
	else
	{
		out << setw(7) << 'A';
		for (int i = 0; i < sistem.n - 1; ++i)
			out << setw(7) << ' ';
		out << setw(7) << 'B' << '\n';
		for (int i = 0; i < sistem.n; ++i)
		{
			for (int j = 0; j < sistem.n; ++j)
				out << setw(7) << sistem.A[i][j];
			out << setw(7) << sistem.B[i];
			out << '\n';
		}
	}
	return out;
}

void Gauss::PrintA(ostream & out)
{
	if (errFlag)
		out << err;
	else
	{
		for (int i = 0; i < n; ++i)
		{
			for (int j = 0; j < n; ++j)
				out << setw(7) << setprecision(3) << A[i][j];
			out << '\n';
		}
	}
}

bool Gauss::SolveSLAEGeneral()
{
	if (errFlag)
	{
		return false;
	}
	for (int i = 0; i < n; ++i)
	{
		for (int j = i + 1; j < n; ++j)
			A[i][j] /= A[i][i];
		B[i] /= A[i][i];

		for (int j = i + 1; j < n; ++j)
		{
			for (int k = i + 1; k < n; ++k)
				A[j][k] -= A[i][k] * A[j][i];
			B[j] -= B[i] * A[j][i];
		}
	}

	for (int i = n - 2; i >= 0; --i)
		for (int j = i + 1; j < n; ++j)
			B[i] -= B[j] * A[i][j];
	solved = true;
	return true;
}
bool Gauss::SolveSLAEMainElement()
{
	if (errFlag)
	{
		return false;
	}
	for (int i = 0; i < n; ++i)
	{
		int m = i;
		for (int j = i; j < n; ++j)
		{
			if (abs(A[m][i]) < abs(A[j][i]))
				m = j;
		}
		A[i].swap(A[m]);
		swap(B[i], B[m]);
		if (abs(A[i][i]) <= 0.00000001)
		{
			errFlag = true;
			err = "Determinant equals 0";
			return false;
		}

		for (int j = i + 1; j < n; ++j)
			A[i][j] /= A[i][i];
		B[i] /= A[i][i];

		for (int j = i + 1; j < n; ++j)
		{
			for (int k = i + 1; k < n; ++k)
				A[j][k] -= A[i][k] * A[j][i];
			B[j] -= B[i] * A[j][i];
		}
	}

	for (int i = n - 2; i >= 0; --i)
		for (int j = i + 1; j < n; ++j)
			B[i] -= B[j] * A[i][j];
	solved = true;
	return true;
}