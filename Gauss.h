#pragma once
#include <vector>
#include <string>
#include <iostream>

using namespace std;
class Gauss
{
private:
	int n;
	vector< vector<double> > A;
	vector<double> B;
	string err = "";
	bool solved = false;

public:
	bool errFlag = false;
	Gauss(int n_);
	Gauss(vector< vector<double> >& A_, vector<double>& B_);
	bool SolveSLAEGeneral();
	bool SolveSLAEMainElement();
	vector<double> ReturnAnswer() { return B; };
	friend ostream& operator <<(ostream&, const Gauss&);
	void PrintA(ostream& out);
	string GetError() { return err; };
};
