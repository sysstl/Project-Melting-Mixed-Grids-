// Master's work.cpp : Этот файл содержит функцию "main". Здесь начинается и заканчивается выполнение программы.
//

#include <iostream>
#include <thread>
#include <fstream>
#include <vector>
#include <windows.h>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <random>
#include <chrono>
#include <algorithm>
#include <map>
#include <exception>
#include <string>
#include <math.h>
#include "GnuPlot.h"
#include "Integral.h"
#include "Splayn.h"

using namespace std;

double pi = 3.141592654;

enum Metal { Au, Al, Cu, Ni };
enum TypeBeam { Gauss, Vortex };

struct Parametrs {
	double t0;
	double r0;
	double kabs;
	double P0;
	double T00;
	double beta;
};

struct Melting
{
	Metal mt;
	double Q_fusion;
	double T_melting;
	double Density;
	double DensityLiquid;
	double gi;
	double ge;
	double u0;
	double u0_Liquid;
};

struct Interval
{
	Interval(int begin, int end) : begin{ begin }, end{ end } {}
	int begin;
	int end;
};

struct P0V
{
	double p0v_s;
	double p0v_sl;
	double p0v_l;
};

struct Point3D // структура необходимая для определение точек, где произошел разрыв
{
	Point3D(int index_x, int index_y, int index_z) : index_x{ index_x }, index_y{ index_y }, index_z{ index_z } {}
	int index_x;
	int index_y;
	int index_z;
};

string ConvertNumToStringdouble(double i)
{
	string tmp_num_line;

	int length = snprintf(NULL, 0, "%lf", i);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%lf", i);
	for (int j = 0; j < length; j++)
	{
		tmp_num_line.push_back(str[j]); // номер столбца
	}

	return tmp_num_line;
}

string ConvertNumToString(int i)
{
	string tmp_num_line;

	int length = snprintf(NULL, 0, "%i", i);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%i", i);
	for (int j = 0; j < length; j++)
	{
		tmp_num_line.push_back(str[j]); // номер столбца
	}

	return tmp_num_line;
}

void MyCreateFile(int number_of_files_, vector<string> filename)
{
	vector<ofstream> file;
	file.reserve(number_of_files_);

	for (int i = 0; i < number_of_files_; i++)
	{
		file.emplace_back(ofstream{ filename[i] });
	}
}

void MyCreateFile_Animate(int number_of_plots, int count_frame_, int start_number_plot, int start_frame)
{
	vector<ofstream> file;
	file.reserve(number_of_plots * count_frame_);
	for (int i = start_number_plot; i < number_of_plots + start_number_plot; ++i)// цикл для создания системных файлов и файлов для данных для Анимации
	{
		string type_plot_Gif;
		string file_name;
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			type_plot_Gif.push_back(str[j]);
		}

		for (int j = start_frame; j < count_frame_ + start_frame; ++j)
		{
			string number_frame_Gif /*= "GIF"*/;

			int length = snprintf(NULL, 0, "%i", j);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", j);
			for (int k = 0; k < length; k++)
			{
				number_frame_Gif.push_back(str[k]);
			}
			// z_str		// mom_time
//file_name = "e(xy) z(index) = " + type_plot_Gif + " + " + number_frame_Gif + ".txt";
			file_name = "Data_Gif" + type_plot_Gif + " + " + number_frame_Gif + ".txt";
			file.emplace_back(ofstream{ file_name });
		}
	}
}

void Null_Array3D(double*** Array, int Nx, int Ny, int Nz)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
				Array[i][j][k] = 1e-16;
}

void Null_Array(double** Array, int Nx, int Ny)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			Array[i][j] = 1e-16;
}

void SetDataArray(double** Array, int Nx, int Ny, double a)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			Array[i][j] = a;
}

void SetDataArray3D(double*** Array, int Nx, int Ny, int Nz, double a)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
				Array[i][j][k] = a;
}

void CopyDataArray3D(double*** ArrayTo, double*** ArrayFrom, int Nx, int Ny, int Nz)
{
	for (int i = 0; i < Nx; i++)
		for (int j = 0; j < Ny; j++)
			for (int k = 0; k < Nz; k++)
				ArrayTo[i][j][k] = ArrayFrom[i][j][k];
}

void SelectDataForPlot2D(vector<string> & Initial_Data, vector<string> & Final_Data, double x_left_boundary, double x_right_boundary)
{
	// Initial_Data - вектор названий файлов с данными
	vector<string>::iterator PP_it, PP_it_res;
	int i = 1;
	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		string text;
		text = ConvertNumToString(i);
		string file_name = "AData" + text + ".txt";
		Final_Data.push_back(file_name);
		text.clear();
		i++;
	}

	MyCreateFile(Initial_Data.size(), Final_Data);
	PP_it_res = Final_Data.begin();

	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		ifstream in((*PP_it)); // окрываем файл для чтения

		VecD X;
		VecD Y1;
		VecD::iterator it_x, it_y1;

		if (in.is_open())
		{
			double x, y1;
			while (in >> x >> y1)
			{
				/*	new_points.push_back(Point{ x, y });*/
				X.push_back(x);
				Y1.push_back(y1);
				//cout << x << "   " << y << endl;
			}
		}
		in.close();

		ofstream fout((*PP_it_res));
		it_y1 = Y1.begin();
		for (it_x = X.begin(); it_x != X.end(); it_x++)
		{
			if ((*it_x) >= x_left_boundary && (*it_x) <= x_right_boundary)
			{
				fout << (*it_x) << "   " << (*it_y1) << endl;
			}
			it_y1++;
		}

		PP_it_res++;
		X.clear();
		Y1.clear();
		fout.close();
	}
}

void SelectDataForPlotColor(vector<string> & Initial_Data, vector<string> & Final_Data, double x_left_boundary, double x_right_boundary, double y_left_boundary, double y_right_boundary)
{
	// Initial_Data - вектор названий файлов с данными
	vector<string>::iterator PP_it, PP_it_res;
	int i = 1;
	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		string text;
		text = ConvertNumToString(i);
		string file_name = "ADataColor" + text + ".txt";
		Final_Data.push_back(file_name);
		text.clear();
		i++;
	}

	MyCreateFile(Initial_Data.size(), Final_Data);
	PP_it_res = Final_Data.begin();

	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		ifstream in((*PP_it)); // окрываем файл для чтения

		VecD X;
		VecD Y;
		VecD Z;
		VecD::iterator it_x, it_y, it_z;

		if (in.is_open())
		{
			double x, y, z;
			while (in >> x >> y >> z)
			{
				/*	new_points.push_back(Point{ x, y });*/
				X.push_back(x);
				Y.push_back(y);
				Z.push_back(z);
				//cout << x << "   " << y << endl;
			}
		}
		in.close();

		ofstream fout((*PP_it_res));
		it_y = Y.begin();
		it_z = Z.begin();
		for (it_x = X.begin(); it_x != X.end(); it_x++)
		{
			if ((*it_x) >= x_left_boundary && (*it_x) <= x_right_boundary && (*it_y) >= y_left_boundary && (*it_y) <= y_right_boundary)
			{
				fout << (*it_x) << "   " << (*it_y) << "   " << (*it_z) << endl;
			}
			it_y++;
			it_z++;
		}

		PP_it_res++;
		X.clear();
		Y.clear();
		Z.clear();
		fout.close();
	}
}

void SelectDataFromPlotColorForPlot2D(vector<string> & Initial_Data, vector<string> & Final_Data, double x_left_boundary, double x_right_boundary, double y_fixed)
{
	// Initial_Data - вектор названий файлов с данными
	vector<string>::iterator PP_it, PP_it_res;
	int i = 1;
	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		string text;
		text = ConvertNumToString(i);
		string file_name = "ADataColorFor2D" + text + ".txt";
		Final_Data.push_back(file_name);
		text.clear();
		i++;
	}

	MyCreateFile(Initial_Data.size(), Final_Data);
	PP_it_res = Final_Data.begin();

	for (PP_it = Initial_Data.begin(); PP_it != Initial_Data.end(); PP_it++)
	{
		ifstream in((*PP_it)); // окрываем файл для чтения

		VecD X;
		VecD Y;
		VecD Z;
		VecD::iterator it_x, it_y, it_z;

		if (in.is_open())
		{
			double x, y, z;
			while (in >> x >> y >> z)
			{
				/*	new_points.push_back(Point{ x, y });*/
				X.push_back(x);
				Y.push_back(y);
				Z.push_back(z);
				//cout << x << "   " << y << endl;
			}
		}
		in.close();

		ofstream fout((*PP_it_res));
		it_y = Y.begin();
		it_z = Z.begin();
		for (it_x = X.begin(); it_x != X.end(); it_x++)
		{
			if ((*it_x) >= x_left_boundary && (*it_x) <= x_right_boundary && (*it_y) == y_fixed)
			{
				fout << (*it_x) << "   " << (*it_z) << endl;
			}
			it_y++;
			it_z++;
		}

		PP_it_res++;
		X.clear();
		Y.clear();
		Z.clear();
		fout.close();
	}
}

void Copy_Data_From_Big_Array3D_To_Small_Array3D(double*** tmpe_for_transform_new, double*** tmpe0, int& Nx_heat, int& Ny_heat, int& Nz_heat, int& left_boundary_x_heat, int& right_boundary_x_heat, int& left_boundary_y_heat, int& right_boundary_y_heat, int& down_boundary_z_heat)
{
	int i_transf = 0, j_transf = 0, k_transf = 0;
	for (int i = 0; i < Nx_heat; i++)
	{
		j_transf = 0;
		for (int j = 0; j < Ny_heat; j++)
		{
			k_transf = 0;
			for (int k = 0; k < Nz_heat; k++)
			{
				if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{
					tmpe_for_transform_new[i_transf][j_transf][k_transf] = tmpe0[i][j][k];
					k_transf++;
				}
			}
			if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
			{
				j_transf++;
			}
		}
		if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
		{
			i_transf++;
		}
	}
}

void Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(double*** tmpe_for_transform_new, double*** tmpe0, int Nx_heat, int Ny_heat, int Nz_heat)
{
	int i_transf = 1, j_transf = 1, k_transf = 1;
	for (int i = 1; i < Nx_heat - 1; i++)
	{
		//j_transf = 1;
		for (int j = 1; j < Ny_heat - 1; j++)
		{
			//k_transf = 1;
			for (int k = 0; k < Nz_heat - 1; k++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{
					tmpe_for_transform_new[i][j][k] = tmpe0[i][j][k];
					//k_transf++;
				}
			}
			//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
			{
				//j_transf++;
			}
		}
		//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
		{
			//i_transf++;
		}
	}
}

void Copy_Data_From_Small_Array3D_To_Big_Array3D(double*** tmpe_for_transform_new, double*** tmpe0, int& left_boundary_x_heat, int& right_boundary_x_heat, int& left_boundary_y_heat, int& right_boundary_y_heat, int& down_boundary_z_heat)
{
	int i_transf, j_transf, k_transf;
	i_transf = 1;
	for (int i = left_boundary_x_heat + 1; i < right_boundary_x_heat; i++)
	{
		j_transf = 1;
		for (int j = left_boundary_y_heat + 1; j < right_boundary_y_heat; j++)
		{
			k_transf = 1;
			for (int k = 1; k < down_boundary_z_heat; k++)
			{
				tmpe0[i][j][k] = tmpe_for_transform_new[i_transf][j_transf][k_transf]; // делаем для того чтобы сшить границы для следующего момета времени
				k_transf++;
			}
			j_transf++;
		}
		i_transf++;
	}
}

void Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(double*** tmpe_for_transform_new, double*** tmpe0, int& left_boundary_x_heat, int& right_boundary_x_heat, int& left_boundary_y_heat, int& right_boundary_y_heat, int& down_boundary_z_heat)
{
	int i_transf, j_transf, k_transf;
	i_transf = 0;
	for (int i = left_boundary_x_heat; i <= right_boundary_x_heat; i++)
	{
		j_transf = 0;
		for (int j = left_boundary_y_heat; j <= right_boundary_y_heat; j++)
		{
			k_transf = 0;
			for (int k = 0; k <= down_boundary_z_heat; k++)
			{
				tmpe0[i][j][k] = tmpe_for_transform_new[i_transf][j_transf][k_transf]; // делаем для того чтобы сшить границы для следующего момета времени
				k_transf++;
			}
			j_transf++;
		}
		i_transf++;
	}
}

void CreateGif2DColor(int numberplot, int time_begin, int time_end, int dt, string firstPartnamefileAnim, string titleplotframe, string xlabel, string ylabel, bool select, double x_left_boundary_for_select, double x_right_boundary_for_select, double y_left_boundary_for_select, double y_right_boundary_for_select, string Time_unit)
{
	vector<int> moments_fix_time_anim;
	vector<string> P_x, P_x_result;
	vector<string>::iterator it_p;
	string time_begin_str, time_end_str, numberplot_str;
	time_begin_str = ConvertNumToString(time_begin);
	time_end_str = ConvertNumToString(time_end);
	numberplot_str = ConvertNumToString(numberplot);
	for (int i = time_begin; i <= time_end; i += dt /*= 10*/) // моменты времени для отрисовки и анимации
	{// аналог n
		moments_fix_time_anim.push_back(i);
		string str, number_frame_Gif;
		number_frame_Gif = ConvertNumToString(i);
		//str = "Data_Gif23 + " + number_frame_Gif + ".txt";
		str = "Data_Gif" + numberplot_str + " + " + number_frame_Gif + ".txt";
		//str = "e(xy) z(index) = " + numberplot_str + " + " + number_frame_Gif + ".txt";
		P_x.push_back(str);
		str.clear();
	}
	//string namefileAnim = firstPartnamefileAnim + " from " + time_begin_str + " fs to " + time_end_str + " fs ";
		string namefileAnim = firstPartnamefileAnim + " from " + time_begin_str + " "+ Time_unit +" to " + time_end_str + " "+ Time_unit +" ";
	if (select)
	{
		SelectDataForPlotColor(P_x, P_x_result, x_left_boundary_for_select, x_right_boundary_for_select, y_left_boundary_for_select, y_right_boundary_for_select);
		GnuPlot Animate2D2(P_x_result);
		Animate2D2.CreateGifOnPlotColor(0, moments_fix_time_anim.size(), moments_fix_time_anim, titleplotframe, xlabel, ylabel, x_right_boundary_for_select, y_right_boundary_for_select, Time_unit, namefileAnim);
	}
	else
	{
		it_p = P_x.begin();
		ifstream fin((*it_p));

		VecD X;
		VecD Y;
		VecD Z;
		VecD::iterator it_x, it_y;

		if (fin.is_open())
		{
			double x, y, z;
			while (fin >> x >> y >> z)
			{
				/*	new_points.push_back(Point{ x, y });*/
				X.push_back(x);
				Y.push_back(y);
				Z.push_back(z);
				//cout << x << "   " << y << endl;
			}
		}
		fin.close();

		it_x = X.end() - 1;
		it_y = Y.end() - 1;

		GnuPlot Animate2D2(P_x); // P_X - вектор названий файлов где лежат данные
		Animate2D2.CreateGifOnPlotColor(0, moments_fix_time_anim.size(), moments_fix_time_anim, titleplotframe, xlabel, ylabel, (*it_x), (*it_y), Time_unit, namefileAnim);
	}
}

void CreateGif_Y_x(int numberplot, int time_begin, int time_end, int dt, string firstPartnamefileAnim, string titleplotframe, string xlabel, string ylabel, vector <string> Legend, bool select, double x_left_boundary_for_select, double x_right_boundary_for_select, string Time_unit)
{
	vector<int> moments_fix_time_anim;
	vector<string> P_x, P_x_result;
	string time_begin_str, time_end_str, numberplot_str;
	time_begin_str = ConvertNumToString(time_begin);
	time_end_str = ConvertNumToString(time_end);
	numberplot_str = ConvertNumToString(numberplot);
	for (int i = time_begin; i <= time_end; i += dt /*= 10*/) // моменты времени для отрисовки и анимации
	{// аналог n
		moments_fix_time_anim.push_back(i);
		string str, number_frame_Gif;
		number_frame_Gif = ConvertNumToString(i);
		//str = "Data_Gif23 + " + number_frame_Gif + ".txt";
		str = "Data_Gif" + numberplot_str + " + " + number_frame_Gif + ".txt";
		P_x.push_back(str);
		str.clear();
	}
	//string namefileAnim = firstPartnamefileAnim + " from " + time_begin_str + " fs to " + time_end_str + " fs ";
	string namefileAnim = firstPartnamefileAnim + " from " + time_begin_str + " " + Time_unit+ " to " + time_end_str + " " + Time_unit + " ";
	if (select)
	{
		SelectDataForPlot2D(P_x, P_x_result, x_left_boundary_for_select, x_right_boundary_for_select);
		GnuPlot Animate2D2(P_x_result);
		Animate2D2.CreateGifOnPlot2D(0, Legend.size(), 3, moments_fix_time_anim.size(), moments_fix_time_anim, titleplotframe, xlabel, ylabel, Legend, namefileAnim, Time_unit);
	}
	else
	{
		GnuPlot Animate2D2(P_x); // P_X - вектор названий файлов где лежат данные
		Animate2D2.CreateGifOnPlot2D(0, Legend.size(), 3, moments_fix_time_anim.size(), moments_fix_time_anim, titleplotframe, xlabel, ylabel, Legend, namefileAnim, Time_unit);
	}
}

/*+++++++++++++*/
void Transform_2D_GridFrom_AcToHeat_xy(int Nx_old, int Ny_old, int Nz_old, double& dx_old, double& dy_old, int Nx_new, int Ny_new, int Nz_new, double& dx_new, double& dy_new, double*** Initial_Data, double*** Final_Data,/* double** Te_old_x_new_y,*/ Splayn spl, int coeff_big_z,
	int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat)
{
	// инетрпол вдоль Y

	VecD X, Y;
	vector<Splayn> vec_spl;
	vector<Splayn>::iterator it_spl;

	for (int i = 0; i < Nx_old; i++) // ac
	{
		for (size_t k = 0; k < Ny_old; k++)
		{
			X.push_back(k * dy_old);
			Y.push_back(Initial_Data[i][k][Nz_old - 1 - coeff_big_z]);
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		/*for (size_t k = 0; k < Ny_new; k++)
		{
			Te_old_x_new_y[i][k] = spl.GetY(k * dy_new);
		}*/
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль X
	int ii = left_boundary_y_heat + 1;
	for (int i = 1; i < Ny_new - 1; i++) // heat
	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < Nx_old; k++) // ac
		{
			X.push_back(k * dx_old);
			Y.push_back((*it_spl).GetY(i * dy_new));//Te_old_x_new_y[k][i]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);
		int kk = 1;
		//cout << endl;
		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_x_heat; k++) // Final = tmpe2
		//for (size_t k = 0; k < Nx_new; k++)
		{
			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
			Final_Data[k][ii][Nz_new - 2] = spl.GetY(kk * dx_new);
			//cout << k << "   " << ii << "   " << Nz_new - 2 <<"   "<< spl.GetY(kk * dx_new) * 300. << endl;
			kk++;
		}
		//system("pause");
		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());
	}
	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());
}

/*+++++++++++++*/
void Transform_2D_GridFrom_AcToHeat_xz(int Nx_old, int Ny_old, int Nz_old, double& dx_old, double& dy_old, double& dz_old, int Nx_new, int left_boundary_y_heat1, int Ny_new, int Nz_new, double& dx_new, double& dz_new, double*** Initial_Data, double*** Final_Data,/* double** Te_old_x_new_z,*/ Splayn spl, int coeff_big_y, int coeff_big_x,
	int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat)
{
	// инетрпол вдоль Z

	VecD X, Y;
	vector<Splayn> vec_spl;
	vector<Splayn>::iterator it_spl;

	/*xz*/

	for (int i = 0; i < Nx_old; i++) // ac
	{
		for (size_t k = 0; k < Nz_old; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back(Initial_Data[i][Ny_old - 1 - coeff_big_y][k]);
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		/*for (size_t k = 0; k < Nz_new; k++)
		{
			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
		}*/
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль X
	int ii = 1;
	for (int i = 1; i < Nz_new - 1; i++) // heat
	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < Nx_old; k++) // ac
		{
			X.push_back(k * dx_old);
			Y.push_back((*it_spl).GetY(i * dz_new)); ;// Y.push_back(Te_old_x_new_z[k][i]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);
		int kk = 1;
		//cout << endl;
		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_x_heat; k++) // Final = tmpe2
		//for (size_t k = 0; k < Nx_new; k++)
		{
			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
			Final_Data[k][Ny_new - 1][ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
			//	cout << k << "   " << Ny_new - 1 << "   " << ii <<"   "<< spl.GetY(kk * dx_new) * 300. << endl;
			kk++;
		}
		//system("pause");
		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());
	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());

	for (int i = 0; i < Nx_old; i++) // ac
	{
		for (size_t k = 0; k < Nz_old; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back(Initial_Data[i][0 + coeff_big_y][k]);
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		/*for (size_t k = 0; k < Nz_new; k++)
		{
			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
		}*/
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль X
	ii = 1;
	for (int i = 1; i < Nz_new - 1; i++) // heat
	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < Nx_old; k++) // ac
		{
			X.push_back(k * dx_old);
			Y.push_back((*it_spl).GetY(i * dz_new)); // Y.push_back(Te_old_x_new_z[k][i]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);
		int kk = 1;
		//cout << endl;
		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_x_heat; k++) // Final = tmpe2
		//for (size_t k = 0; k < Nx_new; k++)
		{
			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
			Final_Data[k][left_boundary_y_heat1 + 1][ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
			//cout << k << "   " << left_boundary_y_heat1 + 1 << "   " << ii << "   " << spl.GetY(kk * dx_new) * 300. << endl;
			kk++;
		}
		//system("pause");
		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());
	}



	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());


	/*yz*/

	for (int i = 0; i < Ny_old; i++) // ac
	{
		for (size_t k = 0; k < Nz_old; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back(Initial_Data[Nx_old - 1 - coeff_big_x][i][k]/*[i][Ny_old - 1 - coeff_big_y][k]*/);
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		/*for (size_t k = 0; k < Nz_new; k++)
		{
			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
		}*/
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль X
	ii = 1;
	for (int i = 1; i < Nz_new - 1; i++) // heat
	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < Ny_old; k++) // ac
		{
			X.push_back(k * dy_old);
			Y.push_back((*it_spl).GetY(i * dz_new)); //Y.push_back(Te_old_x_new_z[k][i]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);
		int kk = 1;
		//cout << endl;
		for (size_t k = left_boundary_y_heat + 1; k < right_boundary_y_heat; k++) // Final = tmpe2
		//for (size_t k = 0; k < Nx_new; k++)
		{
			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
			Final_Data[Nx_new - 1][k]/*[k][Ny_new - 1]*/[ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
			//cout << Nx_new - 1 << "   " << k << "   " << ii << "   " << spl.GetY(kk * dx_new) * 300. << endl;
			kk++;
		}
		//system("pause");
		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());
	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());

	for (int i = 0; i < Ny_old; i++) // ac
	{
		for (size_t k = 0; k < Nz_old; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back(Initial_Data[0 + coeff_big_x][i][k]);
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		/*for (size_t k = 0; k < Nz_new; k++)
		{
			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
		}*/
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль X
	ii = 1;
	for (int i = 1; i < Nz_new - 1; i++) // heat
	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < Ny_old; k++) // ac
		{
			X.push_back(k * dy_old);
			Y.push_back((*it_spl).GetY(i * dz_new)); //Y.push_back(Te_old_x_new_z[k][i]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);
		int kk = 1;
		//cout << endl;
		for (size_t k = left_boundary_y_heat + 1; k < right_boundary_y_heat; k++) // Final = tmpe2
		//for (size_t k = 0; k < Nx_new; k++)
		{
			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
			Final_Data[left_boundary_x_heat + 1][k][ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
			//cout << left_boundary_x_heat + 1 << "   " << k << "   " << ii << "   " << spl.GetY(kk * dx_new) * 300. << endl;
			kk++;
		}
		//system("pause");
		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());
	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());
}

/*+++++++++++++*/
void TransformGrid_2D_xy(int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, double& dx_old, double& dy_old, int Nx_new, int Ny_new, int down_boundary_z_ac, double& dx_new, double& dy_new, double*** Initial_Data, double*** Final_Data, /*double** Te_old_y_new_x,*/ Splayn spl)
{
	VecD X, Y;
	vector<Splayn> vec_spl;
	vector<Splayn>::iterator it_spl;

	// инетрпол вдоль X
	int ii = 0;
	for (int i = left_boundary_y_heat; i <= right_boundary_y_heat; i++)  // Init_data = tmpe2
	{
		int kx = 0;
		for (size_t k = left_boundary_x_heat; k <= right_boundary_x_heat; k++)
		{
			X.push_back(kx * dx_old);
			Y.push_back(Initial_Data[k][i][down_boundary_z_heat]);
			kx++;
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		//for (int k = 0; k < Nx_new; k++) // Nx_new кол-во узлов в tmpe2_ac
		//{
		//	Te_old_y_new_x[k][ii] = spl.GetY(k * dx_new);
		//}

		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль Y

	for (int i = 0; i < Nx_new; i++)
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < ii/*count_y_ht*/; k++)
		{
			X.push_back(k * dy_old);
			Y.push_back((*it_spl).GetY(i * dx_new)); //Y.push_back(Te_old_y_new_x[i][k]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);

		for (int k = 0; k < Ny_new; k++) // tmpe2_ac = final
		{
			Final_Data[i][k][down_boundary_z_ac - 1] = spl.GetY(k * dy_new); // ВСТАВИЛИ ТУДА КУДА НАДО
		}
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());
	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());
}

/*+++++++++++++*/
void TransformGrid_2D_xz(int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, double& dx_old, double& dy_old, double& dz_old, int Nx_new, int Ny_new, int Nz_new, int Nx_part_ac, int Ny_part_ac, double& dx_new, double& dy_new, double& dz_new, double*** Initial_Data, double*** Final_Data, /*double** Te_old_z_new_x,*/ Splayn spl)
{
	VecD X, Y;
	vector<Splayn> vec_spl;
	vector<Splayn>::iterator it_spl;

	// инетрпол вдоль X
	int ii = 0;
	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
	{
		int kx = 0;
		for (size_t k = left_boundary_x_heat; k <= right_boundary_x_heat; k++)
		{
			X.push_back(kx * dx_old);
			Y.push_back(Initial_Data[k][right_boundary_y_heat][i]);
			kx++;
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		//for (int k = 0; k < Nx_new; k++) // Nx_new кол-во узлов в tmpe2_ac
		//{
		//	Te_old_z_new_x[k][ii] = spl.GetY(k * dx_new);
		//}

		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль Z

	for (int i = 0; i < Nx_new; i++)
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < ii; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back((*it_spl).GetY(i * dx_new)); //Y.push_back(Te_old_z_new_x[i][k]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);

		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
		{
			Final_Data[i][Ny_part_ac - 1][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
		}
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());

	// инетрпол вдоль X
	ii = 0;
	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
	{
		int kx = 0;
		for (size_t k = left_boundary_x_heat; k <= right_boundary_x_heat; k++)
		{
			X.push_back(kx * dx_old);
			Y.push_back(Initial_Data[k][left_boundary_y_heat][i]);
			kx++;
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		//for (int k = 0; k < Nx_new; k++) // Nx_new кол-во узлов в tmpe2_ac
		//{
		//	Te_old_z_new_x[k][ii] = spl.GetY(k * dx_new);
		//}

		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль Z

	for (int i = 0; i < Nx_new; i++)
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < ii; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back((*it_spl).GetY(i * dx_new)); //Y.push_back(Te_old_z_new_x[i][k]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);

		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
		{
			Final_Data[i][0][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
		}
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());

	// инетрпол вдоль Y
	ii = 0;
	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
	{
		int kx = 0;
		for (size_t k = left_boundary_y_heat; k <= right_boundary_y_heat; k++)
		{
			X.push_back(kx * dy_old);
			Y.push_back(Initial_Data[right_boundary_x_heat][k][i]);
			kx++;
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		//for (int k = 0; k < Ny_new; k++) // Nx_new кол-во узлов в tmpe2_ac
		//{
		//	Te_old_z_new_x[k][ii] = spl.GetY(k * dy_new);
		//}

		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль Z

	for (int i = 0; i < Ny_new; i++)
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < ii; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back((*it_spl).GetY(i * dy_new)); //Y.push_back(Te_old_z_new_x[i][k]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);

		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
		{
			Final_Data[Nx_part_ac - 1][i][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
		}
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());

	// инетрпол вдоль Y
	ii = 0;
	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
	{
		int kx = 0;
		for (size_t k = left_boundary_y_heat; k <= right_boundary_y_heat; k++)
		{
			X.push_back(kx * dy_old);
			Y.push_back(Initial_Data[left_boundary_x_heat][k][i]);
			kx++;
		}

		spl.InterpolateFast1D(X, Y);
		vec_spl.push_back(spl);

		//for (int k = 0; k < Ny_new; k++) // Nx_new кол-во узлов в tmpe2_ac
		//{
		//	Te_old_z_new_x[k][ii] = spl.GetY(k * dy_new);
		//}

		ii++;
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	// инетрпол вдоль Z

	for (int i = 0; i < Ny_new; i++)
	{
		it_spl = vec_spl.begin();
		for (size_t k = 0; k < ii; k++)
		{
			X.push_back(k * dz_old);
			Y.push_back((*it_spl).GetY(i * dy_new)); //Y.push_back(Te_old_z_new_x[i][k]);
			it_spl++;
		}

		spl.InterpolateFast1D(X, Y);

		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
		{
			Final_Data[0][i][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
		}
		X.clear();
		X.erase(X.begin(), X.end());
		Y.clear();
		Y.erase(Y.begin(), Y.end());

	}

	vec_spl.clear();
	vec_spl.erase(vec_spl.begin(), vec_spl.end());
}

//void Transform_2D_GridFrom_AcToHeat_xy(int Nx_old, int Ny_old, int Nz_old, double& dx_old, double& dy_old, int Nx_new, int Ny_new, int Nz_new, double& dx_new, double& dy_new, double*** Initial_Data, double*** Final_Data, double** Te_old_x_new_y, Splayn spl, int coeff_big_z,
//	int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat)
//{
//	// инетрпол вдоль Y
//
//	VecD X, Y;
//
//	for (int i = 0; i < Nx_old; i++) // ac
//	{
//		for (size_t k = 0; k < Ny_old; k++)
//		{
//			X.push_back(k * dy_old);
//			Y.push_back(Initial_Data[i][k][Nz_old - 1 - coeff_big_z]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (size_t k = 0; k < Ny_new; k++)
//		{
//			Te_old_x_new_y[i][k] = spl.GetY(k * dy_new);
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль X
//	int ii = left_boundary_y_heat + 1;
//	for (int i = 1; i < Ny_new - 1; i++) // heat
//	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
//	{
//		for (size_t k = 0; k < Nx_old; k++) // ac
//		{
//			X.push_back(k * dx_old);
//			Y.push_back(Te_old_x_new_y[k][i]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//		int kk = 1;
//		//cout << endl;
//		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_x_heat; k++) // Final = tmpe2
//		//for (size_t k = 0; k < Nx_new; k++)
//		{
//			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
//			Final_Data[k][ii][Nz_new - 2] = spl.GetY(kk * dx_new);
//			//cout << k << "   " << ii << "   " << Nz_new - 2 <<"   "<< spl.GetY(kk * dx_new) * 300. << endl;
//			kk++;
//		}
//		//system("pause");
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//	}
//}
//
//void Transform_2D_GridFrom_AcToHeat_xz(int Nx_old, int Ny_old, int Nz_old, double& dx_old, double& dy_old, double& dz_old, int Nx_new, int left_boundary_y_heat1, int Ny_new, int Nz_new, double& dx_new, double& dz_new, double*** Initial_Data, double*** Final_Data, double** Te_old_x_new_z, Splayn spl, int coeff_big_y, int coeff_big_x,
//	int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat)
//{
//	// инетрпол вдоль Z
//
//	VecD X, Y;
//
//	/*xz*/
//
//	for (int i = 0; i < Nx_old; i++) // ac
//	{
//		for (size_t k = 0; k < Nz_old; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Initial_Data[i][Ny_old - 1 - coeff_big_y][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (size_t k = 0; k < Nz_new; k++)
//		{
//			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль X
//	int ii = 1;
//	for (int i = 1; i < Nz_new - 1; i++) // heat
//	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
//	{
//		for (size_t k = 0; k < Nx_old; k++) // ac
//		{
//			X.push_back(k * dx_old);
//			Y.push_back(Te_old_x_new_z[k][i]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//		int kk = 1;
//		//cout << endl;
//		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_x_heat; k++) // Final = tmpe2
//		//for (size_t k = 0; k < Nx_new; k++)
//		{
//			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
//			Final_Data[k][Ny_new - 1][ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
//			//	cout << k << "   " << Ny_new - 1 << "   " << ii <<"   "<< spl.GetY(kk * dx_new) * 300. << endl;
//			kk++;
//		}
//		//system("pause");
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//	}
//
//
//
//	for (int i = 0; i < Nx_old; i++) // ac
//	{
//		for (size_t k = 0; k < Nz_old; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Initial_Data[i][0 + coeff_big_y][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (size_t k = 0; k < Nz_new; k++)
//		{
//			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль X
//	ii = 1;
//	for (int i = 1; i < Nz_new - 1; i++) // heat
//	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
//	{
//		for (size_t k = 0; k < Nx_old; k++) // ac
//		{
//			X.push_back(k * dx_old);
//			Y.push_back(Te_old_x_new_z[k][i]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//		int kk = 1;
//		//cout << endl;
//		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_x_heat; k++) // Final = tmpe2
//		//for (size_t k = 0; k < Nx_new; k++)
//		{
//			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
//			Final_Data[k][left_boundary_y_heat1 + 1][ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
//			//cout << k << "   " << left_boundary_y_heat1 + 1 << "   " << ii << "   " << spl.GetY(kk * dx_new) * 300. << endl;
//			kk++;
//		}
//		//system("pause");
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//	}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//	/*yz*/
//
//	for (int i = 0; i < Ny_old; i++) // ac
//	{
//		for (size_t k = 0; k < Nz_old; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Initial_Data[Nx_old - 1 - coeff_big_x][i][k]/*[i][Ny_old - 1 - coeff_big_y][k]*/);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (size_t k = 0; k < Nz_new; k++)
//		{
//			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль X
//	ii = 1;
//	for (int i = 1; i < Nz_new - 1; i++) // heat
//	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
//	{
//		for (size_t k = 0; k < Ny_old; k++) // ac
//		{
//			X.push_back(k * dy_old);
//			Y.push_back(Te_old_x_new_z[k][i]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//		int kk = 1;
//		//cout << endl;
//		for (size_t k = left_boundary_y_heat + 1; k < right_boundary_y_heat; k++) // Final = tmpe2
//		//for (size_t k = 0; k < Nx_new; k++)
//		{
//			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
//			Final_Data[Nx_new - 1][k]/*[k][Ny_new - 1]*/[ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
//			//cout << Nx_new - 1 << "   " << k << "   " << ii << "   " << spl.GetY(kk * dx_new) * 300. << endl;
//			kk++;
//		}
//		//system("pause");
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//	}
//
//
//
//	for (int i = 0; i < Ny_old; i++) // ac
//	{
//		for (size_t k = 0; k < Nz_old; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Initial_Data[0 + coeff_big_x][i][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (size_t k = 0; k < Nz_new; k++)
//		{
//			Te_old_x_new_z[i][k] = spl.GetY(k * dz_new);
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль X
//	ii = 1;
//	for (int i = 1; i < Nz_new - 1; i++) // heat
//	//for (int i = left_boundary_y_heat + 1; i < Ny_new; i++) // heat
//	{
//		for (size_t k = 0; k < Ny_old; k++) // ac
//		{
//			X.push_back(k * dy_old);
//			Y.push_back(Te_old_x_new_z[k][i]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//		int kk = 1;
//		//cout << endl;
//		for (size_t k = left_boundary_x_heat + 1; k < right_boundary_y_heat; k++) // Final = tmpe2
//		//for (size_t k = 0; k < Nx_new; k++)
//		{
//			//Final_Data[k][i][Nz_new - 2] = spl.GetY(k * dx_new);
//			Final_Data[left_boundary_x_heat + 1][k][ii]/*[ii][Nz_new - 2]*/ = spl.GetY(kk * dx_new);
//			//cout << left_boundary_x_heat + 1 << "   " << k << "   " << ii << "   " << spl.GetY(kk * dx_new) * 300. << endl;
//			kk++;
//		}
//		//system("pause");
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//	}
//}
//
//void TransformGrid_2D_xy(int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, double& dx_old, double& dy_old, int Nx_new, int Ny_new, int down_boundary_z_ac, double& dx_new, double& dy_new, double*** Initial_Data, double*** Final_Data, double** Te_old_y_new_x, Splayn spl)
//{
//	VecD X, Y;
//
//	// инетрпол вдоль X
//	int ii = 0;
//	for (int i = left_boundary_y_heat; i <= right_boundary_y_heat; i++)  // Init_data = tmpe2
//	{
//		int kx = 0;
//		for (size_t k = left_boundary_x_heat; k <= right_boundary_x_heat; k++)
//		{
//			X.push_back(kx * dx_old);
//			Y.push_back(Initial_Data[k][i][down_boundary_z_heat]);
//			kx++;
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nx_new; k++) // Nx_new кол-во узлов в tmpe2_ac
//		{
//			Te_old_y_new_x[k][ii] = spl.GetY(k * dx_new);
//		}
//
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль Y
//
//	for (int i = 0; i < Nx_new; i++)
//	{
//		for (size_t k = 0; k < ii/*count_y_ht*/; k++)
//		{
//			X.push_back(k * dy_old);
//			Y.push_back(Te_old_y_new_x[i][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Ny_new; k++) // tmpe2_ac = final
//		{
//			Final_Data[i][k][down_boundary_z_ac - 1] = spl.GetY(k * dy_new); // ВСТАВИЛИ ТУДА КУДА НАДО
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//}
//
//void TransformGrid_2D_xz(int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, double& dx_old, double& dy_old, double& dz_old, int Nx_new, int Ny_new, int Nz_new, int Nx_part_ac, int Ny_part_ac, double& dx_new, double& dy_new, double& dz_new, double*** Initial_Data, double*** Final_Data, double** Te_old_z_new_x, Splayn spl)
//{
//	VecD X, Y;
//
//	// инетрпол вдоль X
//	int ii = 0;
//	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
//	{
//		int kx = 0;
//		for (size_t k = left_boundary_x_heat; k <= right_boundary_x_heat; k++)
//		{
//			X.push_back(kx * dx_old);
//			Y.push_back(Initial_Data[k][right_boundary_y_heat][i]);
//			kx++;
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nx_new; k++) // Nx_new кол-во узлов в tmpe2_ac
//		{
//			Te_old_z_new_x[k][ii] = spl.GetY(k * dx_new);
//		}
//
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль Z
//
//	for (int i = 0; i < Nx_new; i++)
//	{
//		//for (int j = 0; j < Nz_new; j++)
//
//		for (size_t k = 0; k < ii; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Te_old_z_new_x[i][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
//		{
//			Final_Data[i][Ny_part_ac - 1][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//
//
//
//	// инетрпол вдоль X
//	ii = 0;
//	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
//	{
//		int kx = 0;
//		for (size_t k = left_boundary_x_heat; k <= right_boundary_x_heat; k++)
//		{
//			X.push_back(kx * dx_old);
//			Y.push_back(Initial_Data[k][left_boundary_y_heat][i]);
//			kx++;
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nx_new; k++) // Nx_new кол-во узлов в tmpe2_ac
//		{
//			Te_old_z_new_x[k][ii] = spl.GetY(k * dx_new);
//		}
//
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль Z
//
//	for (int i = 0; i < Nx_new; i++)
//	{
//		//for (int j = 0; j < Nz_new; j++)
//
//		for (size_t k = 0; k < ii; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Te_old_z_new_x[i][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
//		{
//			Final_Data[i][0][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//
//
//
//
//
//
//
//
//
//
//
//	// инетрпол вдоль Y
//	ii = 0;
//	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
//	{
//		int kx = 0;
//		for (size_t k = left_boundary_y_heat; k <= right_boundary_y_heat; k++)
//		{
//			X.push_back(kx * dy_old);
//			Y.push_back(Initial_Data[right_boundary_x_heat][k][i]);
//			kx++;
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Ny_new; k++) // Nx_new кол-во узлов в tmpe2_ac
//		{
//			Te_old_z_new_x[k][ii] = spl.GetY(k * dy_new);
//		}
//
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль Z
//
//	for (int i = 0; i < Ny_new; i++)
//	{
//		//for (int j = 0; j < Nz_new; j++)
//
//		for (size_t k = 0; k < ii; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Te_old_z_new_x[i][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
//		{
//			Final_Data[Nx_part_ac - 1][i][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//
//	// инетрпол вдоль Y
//	ii = 0;
//	for (int i = 0; i <= down_boundary_z_heat; i++)  // Init_data = tmpe2
//	{
//		int kx = 0;
//		for (size_t k = left_boundary_y_heat; k <= right_boundary_y_heat; k++)
//		{
//			X.push_back(kx * dy_old);
//			Y.push_back(Initial_Data[left_boundary_x_heat][k][i]);
//			kx++;
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Ny_new; k++) // Nx_new кол-во узлов в tmpe2_ac
//		{
//			Te_old_z_new_x[k][ii] = spl.GetY(k * dy_new);
//		}
//
//		ii++;
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//
//	// инетрпол вдоль Z
//
//	for (int i = 0; i < Ny_new; i++)
//	{
//		//for (int j = 0; j < Nz_new; j++)
//
//		for (size_t k = 0; k < ii; k++)
//		{
//			X.push_back(k * dz_old);
//			Y.push_back(Te_old_z_new_x[i][k]);
//		}
//
//		spl.InterpolateFast1D(X, Y);
//
//		for (int k = 0; k < Nz_new; k++) // tmpe2_ac = final
//		{
//			Final_Data[0][i][k]/*[k][down_boundary_z_ac - 1] */ = spl.GetY(k * dz_new); // ВСТАВИЛИ ТУДА КУДА НАДО
//		}
//		X.clear();
//		X.erase(X.begin(), X.end());
//		Y.clear();
//		Y.erase(Y.begin(), Y.end());
//
//	}
//}

void TransformGridFrom_AcToHeat(int Nx_old, int Ny_old, int Nz_old, double& dx_old, double& dy_old, double& dz_old, int Nx_new, int Ny_new, int Nz_new, double& dx_new, double& dy_new, double& dz_new, double*** Initial_Data, double*** Final_Data, double*** Te_old_xy_new_z, double*** Te_old_y_new_xz, Splayn spl)
{
	// инетрпол вдоль Y

	VecD X, Y;

	for (int i = 0; i < Nx_old; i++) // ac
	{
		for (int j = 0; j < Nz_old; j++) // ac
		{
			for (size_t k = 0; k < Ny_old; k++)
			{
				X.push_back(k * dy_old);
				Y.push_back(Initial_Data[i][k][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (size_t k = 0; k < Ny_new; k++)
			{
				Te_old_y_new_xz[i][k][j] = spl.GetY(k * dy_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	// инетрпол вдоль X

	for (int i = 0; i < Ny_new; i++) // heat
	{
		for (int j = 0; j < Nz_old; j++) // ac
		{
			for (size_t k = 0; k < Nx_old; k++)
			{
				X.push_back(k * dx_old);
				Y.push_back(Te_old_y_new_xz[k][i][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (size_t k = 0; k < Nx_new; k++)
			{
				Te_old_xy_new_z[k][i][j] = spl.GetY(k * dx_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}


	// инетрпол вдоль Z

	for (int i = 0; i < Nx_new; i++)
	{
		for (int j = 0; j < Ny_new; j++)
		{
			for (size_t k = 0; k < Nz_old; k++)
			{
				X.push_back(k * dz_old);
				Y.push_back(Te_old_xy_new_z[i][j][k]);
			}

			spl.InterpolateFast1D(X, Y);
			for (int k = 0; k < Nz_new; k++)
			{
				Final_Data[i][j][k] = spl.GetY(k * dz_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}
}

void TransformGrid(int Nx_old, int Ny_old, int Nz_old, double& dx_old, double& dy_old, double& dz_old, int Nx_new, int Ny_new, int Nz_new, double& dx_new, double& dy_new, double& dz_new, double*** Initial_Data, double*** Final_Data, double*** Te_old_xy_new_z, double*** Te_old_y_new_xz, Splayn spl)
{
	// инетрпол вдоль Z
	VecD X, Y;

	for (size_t i = 0; i < Nx_old; i++)
	{
		for (size_t j = 0; j < Ny_old; j++)
		{
			for (size_t k = 0; k < Nz_old; k++)
			{
				X.push_back(k * dz_old);
				Y.push_back(Initial_Data[i][j][k]);
			}

			spl.InterpolateFast1D(X, Y);

			for (size_t k = 0; k < Nz_new; k++)
			{
				Te_old_xy_new_z[i][j][k] = spl.GetY(k * dz_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	// инетрпол вдоль X

	for (int i = 0; i < Ny_old; i++)
	{
		for (int j = 0; j < Nz_new; j++)
		{
			for (size_t k = 0; k < Nx_old; k++)
			{
				X.push_back(k * dx_old);
				Y.push_back(Te_old_xy_new_z[k][i][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (int k = 0; k < Nx_new; k++)
			{
				Te_old_y_new_xz[k][i][j] = spl.GetY(k * dx_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}

	// инетрпол вдоль Y

	for (int i = 0; i < Nx_new; i++)
	{
		for (int j = 0; j < Nz_new; j++)
		{
			for (size_t k = 0; k < Ny_old; k++)
			{
				X.push_back(k * dy_old);
				Y.push_back(Te_old_y_new_xz[i][k][j]);
			}

			spl.InterpolateFast1D(X, Y);

			for (int k = 0; k < Ny_new; k++)
			{
				Final_Data[i][k][j] = spl.GetY(k * dy_new);
			}
			X.clear();
			X.erase(X.begin(), X.end());
			Y.clear();
			Y.erase(Y.begin(), Y.end());
		}
	}
}

bool operator ==(const Point3D & p1, const Point3D & p2)
{
	return p1.index_x == p2.index_x && p1.index_y == p2.index_y && p1.index_z == p2.index_z;
}

bool comp_x(const Point3D & pt1, const Point3D & pt2)
{
	return (pt1.index_x < pt2.index_x);
}

bool comp_y(const Point3D & pt1, const Point3D & pt2)
{
	return (pt1.index_y < pt2.index_y);
}

bool comp_z(const Point3D & pt1, const Point3D & pt2)
{
	return (pt1.index_z < pt2.index_z);
}

void MySort_Point3D_y(vector<Point3D> & v)
{
	vector<Point3D> copy_v, v_for_partial_sort;
	vector<Point3D>::iterator it;
	copy_v = v;
	v.erase(v.begin(), v.end());
	v.clear();
	vector<Point3D>::iterator begin_pos, end_pos, temp;
	begin_pos = copy_v.begin();
	int count = 0;
	for (it = copy_v.begin(); it != copy_v.end(); it++)// цикл по точкам с которым сравниваем
	{
		if (((*it).index_x == (*(it + 1)).index_x) && (it != copy_v.end()))
		{
			count++;
		}
		else
		{
			end_pos = it;
			Point3D tmp(0, 0, 0);
			for (temp = begin_pos; temp <= end_pos; temp++)
			{
				v_for_partial_sort.push_back(tmp);
			}
			copy(begin_pos, end_pos + 1, v_for_partial_sort.begin());
			sort(v_for_partial_sort.begin(), v_for_partial_sort.end(), comp_y);
			copy(v_for_partial_sort.begin(), v_for_partial_sort.end(), v.end());

			for (temp = v_for_partial_sort.begin(); temp != v_for_partial_sort.end(); temp++)
			{
				v.push_back((*temp));
			}

			begin_pos = (it + 1);
			count = 0;
			v_for_partial_sort.clear();
			v_for_partial_sort.erase(v_for_partial_sort.begin(), v_for_partial_sort.end());
		}
	}
}

void MySort_Point3D_z(vector<Point3D> & v)
{
	vector<Point3D> copy_v, v_for_partial_sort;
	vector<Point3D>::iterator it;
	copy_v = v;
	v.clear();
	v.erase(v.begin(), v.end());
	vector<Point3D>::iterator begin_pos, end_pos, temp;
	begin_pos = copy_v.begin();
	int count = 0;
	for (it = copy_v.begin(); it != copy_v.end(); it++)// цикл по точкам с которым сравниваем
	{
		if (((*it).index_x == (*(it + 1)).index_x) && ((*it).index_y == (*(it + 1)).index_y) && (it != copy_v.end()))
		{
			count++;
		}
		else
		{
			end_pos = it;
			Point3D tmp(0, 0, 0);
			for (temp = begin_pos; temp <= end_pos; temp++)
			{
				v_for_partial_sort.push_back(tmp);
			}
			copy(begin_pos, end_pos + 1, v_for_partial_sort.begin());
			sort(v_for_partial_sort.begin(), v_for_partial_sort.end(), comp_z);
			copy(v_for_partial_sort.begin(), v_for_partial_sort.end(), v.end());

			for (temp = v_for_partial_sort.begin(); temp != v_for_partial_sort.end(); temp++)
			{
				v.push_back((*temp));
			}

			begin_pos = (it + 1);
			count = 0;
			v_for_partial_sort.clear();
			v_for_partial_sort.erase(v_for_partial_sort.begin(), v_for_partial_sort.end());
		}
	}
}

void MySort_Point3D_x(vector<Point3D> & v)
{
	vector<Point3D> copy_v, v_for_partial_sort;
	vector<Point3D>::iterator it;
	copy_v = v;
	v.clear();
	v.erase(v.begin(), v.end());
	vector<Point3D>::iterator begin_pos, end_pos, temp;
	begin_pos = copy_v.begin();
	int count = 0;
	for (it = copy_v.begin(); it != copy_v.end(); it++)// цикл по точкам с которым сравниваем
	{
		if (((*it).index_z == (*(it + 1)).index_z) && ((*it).index_y == (*(it + 1)).index_y) && (it != copy_v.end()))
		{
			count++;
		}
		else
		{
			end_pos = it;
			Point3D tmp(0, 0, 0);
			for (temp = begin_pos; temp <= end_pos; temp++)
			{
				v_for_partial_sort.push_back(tmp);
			}
			copy(begin_pos, end_pos + 1, v_for_partial_sort.begin());
			sort(v_for_partial_sort.begin(), v_for_partial_sort.end(), comp_x);
			copy(v_for_partial_sort.begin(), v_for_partial_sort.end(), v.end());

			for (temp = v_for_partial_sort.begin(); temp != v_for_partial_sort.end(); temp++)
			{
				v.push_back((*temp));
			}

			begin_pos = (it + 1);
			count = 0;
			v_for_partial_sort.clear();
			v_for_partial_sort.erase(v_for_partial_sort.begin(), v_for_partial_sort.end());
		}
	}
}

void My_unique(vector<Point3D> & v)
{
	vector<Point3D>::iterator it1_comp = v.begin();
	vector<Point3D>::iterator it2;

	for (it1_comp = v.begin(); it1_comp != v.end(); it1_comp++)// цикл по точкам с которым сравниваем
	{
		for (it2 = it1_comp + 1; it2 != v.end(); it2++)// перебор по всем точкам остальным
		{
			if ((*it1_comp) == (*it2))
			{
				v.erase(it2);
				it1_comp--;
				it2--;
			}
		}
	}
}

double My_function(double x) // от 0 до п/2
{
	return (pow(x, 4) * exp(x)) / (pow(exp(x) - 1, 2));
}

double delta_function(double Tl, double Tm, double delta)
{
	return (1 / (sqrt(2 * pi) * delta)) * exp(-pow(Tl - Tm, 2) / (2 * delta * delta));
}

double Dependence_k_e_on_T(Metal mt, double T_e) // температура в К
{
	double k_b = 1.38e-23; // постоянная Больцмана
	int i;
	//   Au, Al, Cu, Ni
	//double n_a[4] = { 5.90e+28, 18.10e+28, 8.45e+28, 9.13e+28 };// 18.10e+28;//18.10e+28;//8.45e+28; // концентрация атомов
	double T_D[4] = { 165, 428, 340, 375 };//428;//340; // температура Дебая (К)
	double T_l = 300;
	double v_F[4] = { 1.40e+6, 1.98e+6, 1.57e+6, 2.04e+6 }; // скорость Ферми (м/с)
	double E_F[4] = { 5.53 * 1.6e-19, 11.70 * 1.6e-19, 7.03 * 1.6e-19, 11.67 * 1.6e-19 }; //5.53*1.6e-19; //11.70 * 1.6e-19; // энергия Ферми (Дж)
	double Y[4] = { 68, 91.2, 97, 46 };// 1065;//91.2; // (Дж/м3 К2)
	double A[4] = { 1.18e+7, 0.376e+7, 1.28e+7, 0.59e+7 };//0.376e+7; // (1/(с К2))
	double B[4] = { 1.25e+11, 3.9e+11, 1.23e+11, 1.4e+11 }; // 3.9e+11; // (1/(с К))

	switch (mt)
	{
	case Au:
		i = 0;
		break;
	case Al:
		i = 1;
		break;
	case Cu:
		i = 2;
		break;
	case Ni:
		i = 3;
		break;
	}

	double bb = (B[i] * k_b) / (A[i] * E_F[i]);
	double K = (k_b * pow(v_F[i], 2) * Y[i] / 3) / (0.147 * A[i] * E_F[i]);
	double v_i = k_b * T_l / E_F[i];
	double v_e = k_b * T_e / E_F[i];
	double k_e = K * v_e * (pow(v_e * v_e + 0.16, 5 / 4) * (v_e * v_e + 0.44)) / (pow(v_e * v_e + 0.092, 1 / 2) * (v_e * v_e + bb * v_i));

	return k_e;
}

double Dependence_C_l_on_T(Metal mt, double T_l)
{
	int i;
	//   Au, Al, Cu, Ni
	//double n_a[4] = { 5.90e+28, 18.10e+28, 8.45e+28, 9.13e+28 };// 18.10e+28;//18.10e+28;//8.45e+28; // концентрация атомов
	double n_a[4] = { 5.90e+28, 6.02e+28, 8.45e+28, 9.13e+28 };// 18.10e+28;//18.10e+28;//8.45e+28; // концентрация атомов
	double T_D[4] = { 165, 428, 340, 375 };//428;//340; // температура Дебая (К)
	double k_b = 1.38e-23; // постоянная Больцмана

	switch (mt)
	{
	case Au:
		i = 0;
		break;
	case Al:
		i = 1;
		break;
	case Cu:
		i = 2;
		break;
	case Ni:
		i = 3;
		break;
	}

	double a = 0;// 1e-5;
	double b = T_D[i] / T_l;

	Integral C_l_fun(My_function, a, b);

	return 9 * n_a[i] * k_b* pow(T_l / T_D[i], 3)* C_l_fun.GaussLegendre_IntervalVariety(5);
}

void GRU1(double*** e2, int Nx, int Ny)
{ // Boundary condition of
  // GRU 1

	for (int k = 0; k < Nx; k++)
	{
		for (int kk = 0; kk < Ny; kk++)
		{
			e2[k][kk][0] = 1e-16;
		}
	}
}

void GRU2(double*** b2x, double*** b2y, double*** b2z, int& Nx, int& Ny, int& Nz)
{// Boundary condition of
 // GRU 2

	for (int k = 0; k < Nx; k++)
	{
		for (int kk = 0; kk < Ny; kk++)
		{
			b2x[k][kk][Nz - 1] = 1e-16;
			b2y[k][kk][Nz - 1] = 1e-16;
			b2z[k][kk][Nz - 1] = 1e-16;
		}
	}

	for (int kk = 0; kk < Ny; kk++)
	{
		for (int kkk = 0; kkk < Nz; kkk++)
		{
			b2x[0][kk][kkk] = 1e-16;
			b2y[0][kk][kkk] = 1e-16;
			b2z[0][kk][kkk] = 1e-16;
			b2x[Nx - 1][kk][kkk] = 1e-16;
			b2y[Nx - 1][kk][kkk] = 1e-16;
			b2z[Nx - 1][kk][kkk] = 1e-16;
		}
	}

	for (int k = 0; k < Nx; k++)
	{
		for (int kkk = 0; kkk < Nz; kkk++)
		{
			b2x[k][0][kkk] = 1e-16;
			b2y[k][0][kkk] = 1e-16;
			b2z[k][0][kkk] = 1e-16;
			b2x[k][Ny - 1][kkk] = 1e-16;
			b2y[k][Ny - 1][kkk] = 1e-16;
			b2z[k][Ny - 1][kkk] = 1e-16;
		}
	}
}

double GaussBeam(int& i, int& j, int& k, int& Nx, int& Ny, double& dx, double& dy, double& dz, double& dt, int& n, double& beta)
{
	return exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);
}

double VortexBeam(int& i, int& j, int& k, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, int& n, double& beta, int m)
{
	double h11 = 0, h12 = 0, hf1 = 0;

	h11 = dx * (i - Nx / 2);
	h12 = dy * (j - Ny / 2);
	hf1 = atan(h12 / h11);
	if (h11 == 0 && h12 == 0)
	{
		hf1 = 0;
	}

	if (h11 == 0)
	{
		hf1 = 0;
	}

	return (pow(sqrt(pow(h11 * h11 + h12 * h12, abs(m))) * cos(m * hf1) * (exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2))), 2) +
		pow(sqrt(pow(h11 * h11 + h12 * h12, abs(m))) * sin(m * hf1) * (exp(-dx * dx * (i - Nx / 2) * (i - Nx / 2)) * exp(-dy * dy * (j - Ny / 2) * (j - Ny / 2))), 2)) * exp(-dz * k) * dt * n * exp(-beta * dt * n);
}

void Calculation00(Metal mt, Parametrs param, Splayn spl_C_e_on_T, Splayn spl_G_e_on_T, TypeBeam tbeam, double*** tmpe0, double*** tmpe1, double*** tmpi0, double*** tmpi1, int& Nx, int& Ny, int& Nz, double& dx, double& dy, double& dz, double& dt, double& A2, int& n, double& beta)
{// Calculation of heat conduction by explicit schem

	double A111 = 0, B111 = 0, A1i111 = 0, CC1111 = 0, CC2111 = 0, FFF = 0;

	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int k = 1; k < Nz - 1; k++)
			{
				A111 = (Dependence_k_e_on_T(mt, param.T00 * tmpe0[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe0[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
				B111 = param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe0[i][j][k]) / (100. * 100. * 100.)) * param.T00);
				A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi0[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
				if (tbeam == Gauss)
				{
					FFF = GaussBeam(i, j, k, Nx, Ny, dx, dy, dz, dt, n, beta);
				}

				if (tbeam == Vortex)
				{
					FFF = VortexBeam(i, j, k, Nx, Ny, Nz, dx, dy, dz, dt, n, beta, 1);
				}

				tmpe1[i][j][k] = tmpe0[i][j][k] +
					dt * A111 * ((tmpe0[i + 1][j][k] + tmpe0[i - 1][j][k] - 2 * tmpe0[i][j][k]) / (dx * dx) + (tmpe0[i][j + 1][k] + tmpe0[i][j - 1][k] - 2 * tmpe0[i][j][k]) / (dy * dy) + A2 * (tmpe0[i][j][k - 1] + tmpe0[i][j][k + 1] - 2 * tmpe0[i][j][k]) / (dz * dz))
					+ dt * B111 * FFF;

				tmpi1[i][j][k] = tmpi0[i][j][k] +
					dt * A1i111 * ((tmpi0[i + 1][j][k] + tmpi0[i - 1][j][k] - 2 * tmpi0[i][j][k]) / (dx * dx) + (tmpi0[i][j + 1][k] + tmpi0[i][j - 1][k] - 2 * tmpi0[i][j][k]) / (dy * dy) + A2 * (tmpi0[i][j][k - 1] + tmpi0[i][j][k + 1] - 2 * tmpi0[i][j][k]) / (dz * dz));

			}
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			tmpe1[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpi1[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpe1[i][j][0] = tmpe1[i][j][1];
			tmpi1[i][j][0] = tmpi1[i][j][1];
		}
	}

	for (int j = 0; j < Ny; j++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe1[0][j][k] = tmpe1[1][j][k];
			tmpi1[0][j][k] = tmpi1[1][j][k];
			tmpe1[Nx - 1][j][k] = tmpe1[Nx - 2][j][k];
			tmpi1[Nx - 1][j][k] = tmpi1[Nx - 2][j][k];
		}
	}

	for (int i = 0; i < Nx; i++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe1[i][0][k] = tmpe1[i][1][k];
			tmpi1[i][0][k] = tmpi1[i][1][k];
			tmpe1[i][Ny - 1][k] = tmpe1[i][Ny - 2][k];
			tmpi1[i][Ny - 1][k] = tmpi1[i][Ny - 2][k];
		}
	}
}

void fun2_Tei_for_mixed_grids(int Nx, int Ny, int Nz, int Nx_ac, int Ny_ac, int Nz_ac, double& dx, double& dy, double& dz, double& dx_ac, double& dy_ac, double& dz_ac, double dt, int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double*** tmpe_ac0, double*** tmpe_ac1, double*** tmpe_ac2, double*** tmpi_ac0, double*** tmpi_ac1, double*** tmpi_ac2, double*** tmpe_for_transform_new, vector<Point3D> & points_rupture_ac, double A2, Parametrs & param, vector<Melting> & Melt_metal, Metal mt, /*double*** Ci, double*** gamma,*/ Splayn & spl_G_e_on_T, Splayn spl_C_e_on_T, double beta, int n, Splayn spl_C_l_on_T, TypeBeam tbeam,
	int Nx_part_ht, int Ny_part_ht, int Nx_part_ac, int Ny_part_ac, int down_boundary_z_acoustic, double*** Tei_old_xy_new_z, double*** Tei_old_y_new_xz, Splayn spl_Te, double*** a2x_ac_new)
{
	vector<Point3D>::iterator it;
	it = points_rupture_ac.begin();
	double delta = 1.;

	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				tmpe2[i][j][k] = GaussBeam(i, j, k, Nx, Ny, dx, dy, dz, dt, n, beta);
			}
		}
	}

	// tmpe2 , a2x_ac_new - в данном случае играют роль исочника импульса
	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, tmpe2, Nx, Ny, Nz, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, Nx_part_ac, Ny_part_ac, down_boundary_z_acoustic, dx_ac, dy_ac, dz_ac, tmpe_for_transform_new, a2x_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	double A111 = 0, B111 = 0, A1i111 = 0, CC1111 = 0, CC2111 = 0, FFF = 0;
	// расчет в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int k = 1; k < Nz - 1; k++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && k < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
					// (i <= 30 && i>= 70 && j <= 30 && j >= 70 && k >= 25) // dx = dy = 10 мкм dz = 4 нм
				{
					//cout << " Calcul TeHT " << endl;
					A1i111 = ((Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
					A111 = ((Dependence_k_e_on_T(mt, param.T00 * tmpe1[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
					B111 = (param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.T00));
					CC1111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.))));
					CC2111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.))));
					if (tbeam == Gauss)
					{
						FFF = GaussBeam(i, j, k, Nx, Ny, dx, dy, dz, dt, n, beta);
					}

					if (tbeam == Vortex)
					{
						FFF = VortexBeam(i, j, k, Nx, Ny, Nz, dx, dy, dz, dt, n, beta, 1);
					}

					tmpe2[i][j][k] = tmpe0[i][j][k] * (1 - 2 * dt * A111 / (dx * dx) - 2 * dt * A111 / (dy * dy) - 2 * dt * A111 * A2 / (dz * dz))
						+ 2 * dt * A111 * ((tmpe1[i + 1][j][k] + tmpe1[i - 1][j][k]) / (dx * dx) + (tmpe1[i][j + 1][k] + tmpe1[i][j - 1][k]) / (dy * dy) + A2 * (tmpe1[i][j][k - 1] + tmpe1[i][j][k + 1]) / (dz * dz))
						+ 2 * dt * B111 * FFF - 2 * dt * CC1111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

					tmpe2[i][j][k] = tmpe2[i][j][k] / (1 + 2 * dt * A111 / (dx * dx) + 2 * dt * A111 / (dy * dy) + 2 * dt * A111 * A2 / (dz * dz));

					double delta = 1.;
					if ((tmpi1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
					{

						tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
							+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
							+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

						tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));

					}

					if ((tmpi1[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi1[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
					{
						//melting = true;
						// старый способ расплава
						/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)));*/

						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						// новый способ расплава
						A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) +*/ spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + */ spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])));

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
							+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
							+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

						tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));

					}

					if ((tmpi1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
					{
						A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)));

						/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])));*/

						tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
							+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
							+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

						tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));
					}

				}
			}
		}
	}


	// // расчет в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int k = 1; k < Nz_ac - 1; k++)
			{
				//	int left_bounday_x_heat, int right_bounday_x_heat, int down_bounday_z_heat, int left_bounday_x_acoustic, int right_bounday_x_acoustic, int down_bounday_z_acoustic
				//if (i >= left_bounday_x_acoustic && i <= right_bounday_x_acoustic && j >= left_bounday_y_acoustic && j <= right_bounday_y_acoustic && k <= down_bounday_z_acoustic) // dx = dy = 5 мкм dz = 2 нм
				{
					//cout << " Calcul TeAC " << endl;
					/*A1i111 = ((Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
					A111 = ((Dependence_k_e_on_T(mt, param.T00 * tmpe_ac1[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
					B111 = (param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.T00));
					CC1111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.))));
					CC2111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1[i][j][k]) / (100. * 100. * 100.))));*/
					if (tbeam == Gauss)
					{
						FFF = a2x_ac_new[i][j][k]; //GaussBeam(i, j, k, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, dt, n, beta);
					}

					if (tbeam == Vortex)
					{
						FFF = VortexBeam(i, j, k, Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, dt, n, beta, 1);
					}

					if (!points_rupture_ac.empty())
					{
						if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == k)
						{
							tmpe_ac2[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
							tmpi_ac2[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
							it++;
						}
						else
						{
							if (tmpe_ac1[i - 1][j][k] == 1e-16 || tmpe_ac1[i + 1][j][k] == 1e-16 || tmpe_ac1[i][j - 1][k] == 1e-16 || tmpe_ac1[i][j + 1][k] == 1e-16 || tmpe_ac1[i][j][k - 1] == 1e-16 || tmpe_ac1[i][j][k + 1] == 1e-16 ||
								tmpi_ac1[i - 1][j][k] == 1e-16 || tmpi_ac1[i + 1][j][k] == 1e-16 || tmpi_ac1[i][j - 1][k] == 1e-16 || tmpi_ac1[i][j + 1][k] == 1e-16 || tmpi_ac1[i][j][k - 1] == 1e-16 || tmpi_ac1[i][j][k + 1] == 1e-16 ||
								tmpe_ac1[i][j][k] == 1e-16 || tmpi_ac1[i][j][k] == 1e-16 || tmpe_ac0[i][j][k] == 1e-16 || tmpi_ac0[i][j][k] == 1e-16)
							{
							}
							else
							{
								A1i111 = ((Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
								A111 = ((Dependence_k_e_on_T(mt, param.T00 * tmpe_ac1[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
								B111 = (param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.T00));
								CC1111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.))));
								CC2111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1[i][j][k]) / (100. * 100. * 100.))));

								tmpe_ac2[i][j][k] = tmpe_ac0[i][j][k] * (1 - 2 * dt * A111 / (dx_ac * dx_ac) - 2 * dt * A111 / (dy_ac * dy_ac) - 2 * dt * A111 * A2 / (dz_ac * dz_ac))
									+ 2 * dt * A111 * ((tmpe_ac1[i + 1][j][k] + tmpe_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpe_ac1[i][j + 1][k] + tmpe_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpe_ac1[i][j][k - 1] + tmpe_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
									+ 2 * dt * B111 * FFF - 2 * dt * CC1111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

								tmpe_ac2[i][j][k] = tmpe_ac2[i][j][k] / (1 + 2 * dt * A111 / (dx_ac * dx_ac) + 2 * dt * A111 / (dy_ac * dy_ac) + 2 * dt * A111 * A2 / (dz_ac * dz_ac));

								if ((tmpi_ac1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
								{
									tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
										+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
										+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

									tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));

								}

								if ((tmpi_ac1[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi_ac1[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
								{
									//melting = true;
								// старый способ расплава
								/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)) * param.r0 * param.r0);
								CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)));*/

								/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

								// новый способ расплава
									A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) +*/ spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])) * param.r0 * param.r0);
									CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + */ spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])));

									////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
									tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
										+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
										+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

									tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));
								}

								if ((tmpi_ac1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
								{
									/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
									CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)));*/

									A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])) * param.r0 * param.r0);
									CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])));

									tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
										+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
										+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

									tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));

								}
							}
						}
					}
					else
					{
						A1i111 = ((Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
						A111 = ((Dependence_k_e_on_T(mt, param.T00 * tmpe_ac1[i][j][k]) / (100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0));
						B111 = (param.kabs * param.P0 * param.t0 * beta / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.T00));
						CC1111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.))));
						CC2111 = ((spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1[i][j][k]) / (100. * 100. * 100.))));

						tmpe_ac2[i][j][k] = tmpe_ac0[i][j][k] * (1 - 2 * dt * A111 / (dx_ac * dx_ac) - 2 * dt * A111 / (dy_ac * dy_ac) - 2 * dt * A111 * A2 / (dz_ac * dz_ac))
							+ 2 * dt * A111 * ((tmpe_ac1[i + 1][j][k] + tmpe_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpe_ac1[i][j + 1][k] + tmpe_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpe_ac1[i][j][k - 1] + tmpe_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
							+ 2 * dt * B111 * FFF - 2 * dt * CC1111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

						tmpe_ac2[i][j][k] = tmpe_ac2[i][j][k] / (1 + 2 * dt * A111 / (dx_ac * dx_ac) + 2 * dt * A111 / (dy_ac * dy_ac) + 2 * dt * A111 * A2 / (dz_ac * dz_ac));

						if ((tmpi_ac1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
						{
							tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
								+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
								+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

							tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));

						}

						if ((tmpi_ac1[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi_ac1[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
						{
							//melting = true;
						// старый способ расплава
						/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)));*/

						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						// новый способ расплава
							A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) +*/ spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])) * param.r0 * param.r0);
							CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + */ spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])));

							////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
							tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
								+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
								+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

							tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));
						}

						if ((tmpi_ac1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
						{
							/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
							CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)));*/

							A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])) * param.r0 * param.r0);
							CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])));

							tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
								+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
								+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

							tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));

						}

					}
					/*if ((tmpi_ac1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
					{
						if (!points_rupture_ac.empty()) {
							if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == k)
							{
								tmpi_ac2[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
								it++;
							}
						}
						else
						{
							tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
								+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi_ac1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
								+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

							tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));
						}
					}*/

					//дальше не подправлял схемы относиельно акст сетки
					/*if ((tmpi_ac1[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi_ac1[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
					{
						//melting = true;
						// старый способ расплава
						/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)));

						/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

						// новый способ расплава
						A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) + spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe_ac1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((/*Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.) +  spl_C_l_on_T.GetY(param.T00 * tmpi_ac1[i][j][k])));

						////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
						if (!points_rupture_ac.empty()) {
							if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == k)
							{
								tmpi_ac2[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
								it++;
							}
						}
						else
						{
							tmpi_ac2[i][j][k] = tmpi_ac0[i][j][k] * (1 - 2 * dt * A1i111 / (dx_ac * dx_ac) - 2 * dt * A1i111 / (dy_ac * dy_ac) - 2 * dt * A1i111 * A2 / (dz_ac * dz_ac))
								+ 2 * dt * A1i111 * ((tmpi_ac1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx_ac * dx_ac) + (tmpi_ac1[i][j + 1][k] + tmpi_ac1[i][j - 1][k]) / (dy_ac * dy_ac) + A2 * (tmpi_ac1[i][j][k - 1] + tmpi_ac1[i][j][k + 1]) / (dz_ac * dz_ac))
								+ 2 * dt * CC2111 * (tmpe_ac1[i][j][k] - tmpi_ac1[i][j][k]);

							tmpi_ac2[i][j][k] = tmpi_ac2[i][j][k] / (1 + 2 * dt * A1i111 / (dx_ac * dx_ac) + 2 * dt * A1i111 / (dy_ac * dy_ac) + 2 * dt * A1i111 * A2 / (dz_ac * dz_ac));
						}
					}*/

					//дальше не подправлял схемы относиельно акст сетки
					/*if ((tmpi1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
					{
						A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)));

						/*A1i111 = (Dependence_k_e_on_T(mt, 300.) / (100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])) * param.r0 * param.r0);
						CC2111 = (spl_G_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) * param.t0 / ((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])));

						tmpi2[i][j][k] = tmpi0[i][j][k] * (1 - 2 * dt * A1i111 / (dx * dx) - 2 * dt * A1i111 / (dy * dy) - 2 * dt * A1i111 * A2 / (dz * dz))
							+ 2 * dt * A1i111 * ((tmpi1[i + 1][j][k] + tmpi1[i - 1][j][k]) / (dx * dx) + (tmpi1[i][j + 1][k] + tmpi1[i][j - 1][k]) / (dy * dy) + A2 * (tmpi1[i][j][k - 1] + tmpi1[i][j][k + 1]) / (dz * dz))
							+ 2 * dt * CC2111 * (tmpe1[i][j][k] - tmpi1[i][j][k]);

						tmpi2[i][j][k] = tmpi2[i][j][k] / (1 + 2 * dt * A1i111 / (dx * dx) + 2 * dt * A1i111 / (dy * dy) + 2 * dt * A1i111 * A2 / (dz * dz));
					}*/
				}
			}
		}
	}
}

void fun21(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi1, double*** tmpi2, double& dz, Parametrs & param, Metal mt)
{
	double a4, a5, T_air = 300, eps = 0.1, sigma = 5.67e-12; // Вт/см2 К4 //???????????????
	double h_c = 0.0015; // Вт/cm2 K
	double f1;
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			tmpe2[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpi2[i][j][Nz - 1] = 1.0;// 1e-16;
			tmpe2[i][j][0] = tmpe2[i][j][1];
			tmpi2[i][j][0] = tmpi2[i][j][1];
			//эта строчка заменяются на следующие
		   //a4 = h_c / (param.kabs * (Dependence_k_e_on_T(mt, 300.) / (100. * 100.))); // kl = 1% от ke (при комнат темп)
		   //a5 = eps * sigma * pow(param.T00, 3) / (param.kabs * (Dependence_k_e_on_T(mt, 300.) / (100. * 100.))); //  kl = 1% от ke (при комнат темп)
		   //f1 = a4 * (tmpi1[i][j][0] - 1.) + a5 * (pow(tmpi1[i][j][0], 4) - 1.);// ГУ-2-го рода
		   //fict_masive[i][j] = tmpi2[i][j][1] - 2 * dz * f1; // С использованием центральной разностной аппроксимации градиента на поверхности z=0 находим значение температуры во внешнем узле Q(i,j,0,2).
		   //tmpi2[i][j][0] = (fict_masive[i][j] + tmpi2[i][j][1]) / 2; // Значение температуры на поверхности раздела сред (узел k=1) находим как среднее между внешним узлом Q(i,j,0,2) и внутренним узлом Q(i,j,2,2

		}
	}
}

void fun21_for_ac(int Nx, int Ny, int Nz, double*** tmpe2, double*** tmpi1, double*** tmpi2, double& dz, Parametrs & param, Metal mt, vector<Point3D> v)
{
	double a4, a5, T_air = 300, eps = 0.1, sigma = 5.67e-12; // Вт/см2 К4 //???????????????
	double h_c = 0.0015; // Вт/cm2 K
	double f1;
	vector<Point3D>::iterator it;
	it = v.begin();
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			tmpe2[i][j][0] = tmpe2[i][j][1];
			tmpi2[i][j][0] = tmpi2[i][j][1];

			//эта строчка заменяются на следующие
		   //a4 = h_c / (param.kabs * (Dependence_k_e_on_T(mt, 300.) / (100. * 100.))); // kl = 1% от ke (при комнат темп)
		   //a5 = eps * sigma * pow(param.T00, 3) / (param.kabs * (Dependence_k_e_on_T(mt, 300.) / (100. * 100.))); //  kl = 1% от ke (при комнат темп)
		   //f1 = a4 * (tmpi1[i][j][0] - 1.) + a5 * (pow(tmpi1[i][j][0], 4) - 1.);// ГУ-2-го рода
		   //fict_masive[i][j] = tmpi2[i][j][1] - 2 * dz * f1; // С использованием центральной разностной аппроксимации градиента на поверхности z=0 находим значение температуры во внешнем узле Q(i,j,0,2).
		   //tmpi2[i][j][0] = (fict_masive[i][j] + tmpi2[i][j][1]) / 2; // Значение температуры на поверхности раздела сред (узел k=1) находим как среднее между внешним узлом Q(i,j,0,2) и внутренним узлом Q(i,j,2,2

		}
	}
}

void fun22(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi2)
{
	for (int j = 0; j < Ny; j++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe2[0][j][k] = tmpe2[1][j][k];
			tmpi2[0][j][k] = tmpi2[1][j][k];
			tmpe2[Nx - 1][j][k] = tmpe2[Nx - 2][j][k];
			tmpi2[Nx - 1][j][k] = tmpi2[Nx - 2][j][k];
		}
	}
}

void fun23(int& Nx, int& Ny, int& Nz, double*** tmpe2, double*** tmpi2)
{
	for (int i = 0; i < Nx; i++)
	{
		for (int k = 0; k < Nz; k++)
		{
			tmpe2[i][0][k] = tmpe2[i][1][k];
			tmpi2[i][0][k] = tmpi2[i][1][k];
			tmpe2[i][Ny - 1][k] = tmpe2[i][Ny - 2][k];
			tmpi2[i][Ny - 1][k] = tmpi2[i][Ny - 2][k];
		}
	}
}

void Calculation0_for_mixed_grids(Metal mt, Parametrs & param, Splayn & spl_C_e_on_T, Splayn & spl_G_e_on_T, TypeBeam tbeam, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double*** tmpe_ac0, double*** tmpe_ac1, double*** tmpe_ac2, double*** tmpi_ac0, double*** tmpi_ac1, double*** tmpi_ac2, double*** tmpe_for_transform_new, double& A2, double& dx, double& dy, double& dz, double& dx_ac, double& dy_ac, double& dz_ac, double& dt, int& Nx, int& Ny, int& Nz, int Nx_ac, int Ny_ac, int Nz_ac, int& n, double& beta, vector<Point3D> & points_rupture_ac, vector<Interval> & new_interv, vector<Melting> & Melt_metal, Splayn spl_C_l_on_T,
	int left_bounday_x_heat, int right_bounday_x_heat, int left_bounday_y_heat, int right_bounday_y_heat, int down_bounday_z_heat, int left_bounday_x_acoustic, int right_bounday_x_acoustic, int left_bounday_y_acoustic, int right_bounday_y_acoustic, int down_bounday_z_acoustic, int Nx_part_ac, int Ny_part_ac, int Nx_part_ht, int Ny_part_ht, double*** Tei_old_xy_new_z, double*** Tei_old_y_new_xz, Splayn spl_Te, double*** a2x_ac_new)
{
	///fun2_for_mixed_grids(ref(Nx), ref(Ny), ref(Nz), ref(Nx_ac), ref(Ny_ac), ref(Nz_ac), ref(dx), ref(dy), ref(dz), ref(dx_ac), ref(dy_ac), ref(dz_ac), ref(dt), tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, tmpe_ac0, tmpe_ac1, tmpe_ac2, tmpi_ac0, tmpi_ac1, tmpi_ac2, points_rupture_ac, ref(A2), ref(new_interv[0].begin), ref(new_interv[0].end), ref(param), ref(Melt_metal), mt, ref(spl_G_e_on_T), spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	///fun2_for_mixed_grids(Nx, Ny, Nz, Nx_ac, Ny_ac, Nz_ac, dx, dy, dz, dx_ac, dy_ac, dz_ac, dt, left_bounday_x_heat, right_bounday_x_heat, left_bounday_y_heat, right_bounday_y_heat, down_bounday_z_heat, left_bounday_x_acoustic, right_bounday_x_acoustic, left_bounday_y_acoustic, right_bounday_y_acoustic, down_bounday_z_acoustic, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, tmpe_ac0, tmpe_ac1, tmpe_ac2, tmpi_ac0, tmpi_ac1, tmpi_ac2, points_rupture_ac, A2,/*double*** A1, double*** A1i, double& A2, double*** B1, double*** CC1, double*** CC2,*/ /*double*** F,*/ /*double*** copy_tmpe1, double*** copy_tmpi1,*/ 0, 0, param, Melt_metal, mt, spl_G_e_on_T, spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam);
	fun2_Tei_for_mixed_grids(Nx, Ny, Nz, Nx_part_ac, Ny_part_ac, down_bounday_z_acoustic, dx, dy, dz, dx_ac, dy_ac, dz_ac, dt, left_bounday_x_heat, right_bounday_x_heat, left_bounday_y_heat, right_bounday_y_heat, down_bounday_z_heat,
		tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, tmpe_ac0, tmpe_ac1, tmpe_ac2, tmpi_ac0, tmpi_ac1, tmpi_ac2, tmpe_for_transform_new,
		points_rupture_ac, A2, param, Melt_metal, mt, spl_G_e_on_T, spl_C_e_on_T, beta, n, spl_C_l_on_T, tbeam, Nx_part_ht, Ny_part_ht, Nx_part_ac, Ny_part_ac, down_bounday_z_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te, a2x_ac_new);

	// ПОДПРАВИТ ГРАНИЧНЫЕ УСЛОВИЯ
	thread t51(fun21, ref(Nx), ref(Ny), ref(Nz), tmpe2, tmpi1, tmpi2,/* fict_masive,*/ ref(dz), ref(param), mt);
	thread t52(fun22, ref(Nx), ref(Ny), ref(Nz), tmpe2, tmpi2);
	thread t53(fun23, ref(Nx), ref(Ny), ref(Nz), tmpe2, tmpi2);
	thread t54(fun21_for_ac, ref(Nx_ac), ref(Ny_ac), ref(Nz_ac), tmpe_ac2, tmpi_ac1, tmpi_ac2,/* fict_masive,*/ ref(dz_ac), ref(param), mt, points_rupture_ac);
	//thread t55(fun22, ref(Nx_ac), ref(Ny_ac), ref(Nz_ac), tmpe_ac2, tmpi_ac2);
	//thread t56(fun23, ref(Nx_ac), ref(Ny_ac), ref(Nz_ac), tmpe_ac2, tmpi_ac2);

	t51.join();
	t52.join();
	t53.join(); ////////////////////////////////////
	t54.join();
	//t55.join();
	//t56.join(); ////////////////////////////////////
}

void Calculation1(double*** a1x, /*double*** a2x*/ double*** b1x, double*** b2x, double*** a1z, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** a1y, double*** e1, double*** b1x_ac, double*** b1y_ac, double*** b1z_ac, double*** b2x_ac, double*** b2y_ac, double*** b2z_ac,
	double*** a1x_ac, double*** a1y_ac, double*** a1z_ac, /*double*** a2x_ac*/ double*** a2y_ac, double*** a2z_ac, double*** e1_ac, double& dx, double& dy, double& dz, double& dx_ac, double& dy_ac, double& dz_ac, double& dt, double& CC0, int& Nx, int& Ny, int& Nz, int left_boundary_x_heat, int right_boundary_x_heat,
	int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, int Nx_ac, int Ny_ac, int Nz_ac, vector<Interval> & new_interv_z, vector<Point3D> & points_rupture_ac)
{
	// Calculation 1
	// solving the equation of motion

	vector<Point3D>::iterator it_ac;

	// расчет b2x в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					b2x[i][j][l] = b1x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * (-(a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) + (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz); //	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz); //a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
				}
			}
		}
	}

	it_ac = points_rupture_ac.begin();
	// расчет b2x в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == l)
					{
						b2x_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						b2x_ac[i][j][l] = b1x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l]) * (-(a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) + (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l]) * ((a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac); //	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
						b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1]) * ((a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac); //a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
					}
				}
				else
				{
					b2x_ac[i][j][l] = b1x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l]) * (-(a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) + (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l]) * ((a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac); //	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
					b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1]) * ((a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac); //a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
				}
			}
		}
	}


	// расчет b2y в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					b2y[i][j][l] = b1y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}

	it_ac = points_rupture_ac.begin();
	// расчет b2y в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == l)
					{
						b2y_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						b2y_ac[i][j][l] = b1y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l]) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l]) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1]) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					}
				}
				else
				{
					b2y_ac[i][j][l] = b1y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l]) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l]) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1]) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
				}
			}
		}
	}



	// расчет b2z(граница) в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			//for (int l = 1; l < Nz - 1; l++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					b2z[i][j][0] = b1z[i][j][0] + dt * (e1[i][j][0] - e1[i - 1][j][0]) * ((a1x[i][j][1] - a1x[i][j][0]) * (a1y[i][j + 1][0] - a1y[i][j][0]) - (a1x[i][j + 1][0] - a1x[i][j][0]) * (a1y[i][j][1] - a1y[i][j][0])) / (dx * dy * dz);
					b2z[i][j][0] = b2z[i][j][0] + dt * (e1[i][j][0] - e1[i][j - 1][0]) * ((a1x[i + 1][j][0] - a1x[i][j][0]) * (a1y[i][j][1] - a1y[i][j][0]) - (a1x[i][j][1] - a1x[i][j][0]) * (a1y[i + 1][j][0] - a1y[i][j][0])) / (dx * dy * dz);
					b2z[i][j][0] = b2z[i][j][0] + dt * (e1[i][j][1] - e1[i][j][0]) * ((a1x[i][j + 1][0] - a1x[i][j][0]) * (a1y[i + 1][j][0] - a1y[i][j][0]) - (a1x[i + 1][j][0] - a1x[i][j][0]) * (a1y[i][j + 1][0] - a1y[i][j][0])) / (dx * dy * dz);
				}
			}
		}
	}

	it_ac = points_rupture_ac.begin();
	// расчет b2z(граница) в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			//for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == 0)
					{
						b2z_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						b2z_ac[i][j][0] = b1z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i - 1][j][0]) * ((a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0]) - (a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
						b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i][j - 1][0]) * ((a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0]) - (a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
						b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][1] - e1_ac[i][j][0]) * ((a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0]) - (a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
					}
				}
				else
				{
					b2z_ac[i][j][0] = b1z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i - 1][j][0]) * ((a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0]) - (a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
					b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i][j - 1][0]) * ((a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0]) - (a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
					b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][1] - e1_ac[i][j][0]) * ((a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0]) - (a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
				}
			}
		}
	}


	// расчет b2z в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l]) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l]) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1]) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}

	it_ac = points_rupture_ac.begin();
	// расчет b2z в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == l)
					{
						b2z_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						b2z_ac[i][j][l] = b1z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i - 1][j][l]) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j - 1][l]) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j][l - 1]) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					}
				}
				else
				{
					b2z_ac[i][j][l] = b1z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i - 1][j][l]) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j - 1][l]) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j][l - 1]) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
				}
			}
		}
	}

}

void Calculation10(double*** a1x,/* double*** a2x*/ double*** b1x, double*** b2x, double*** a1z, double*** b1y, double*** b2y, double*** b1z, double*** b2z,
	double*** a1y, double*** e1, double*** V1, double*** V2,
	double*** a1x_ac, /*double*** a2x_ac,*/ double*** b1x_ac, double*** b2x_ac, double*** a1z_ac, double*** b1y_ac, double*** b2y_ac, double*** b1z_ac,
	double*** b2z_ac, double*** a1y_ac, double*** e1_ac, double*** V1_ac, double*** V2_ac, double& dx, double& dy, double& dz, double& dx_ac, double& dy_ac, double& dz_ac,
	double& dt, double& CC0, int& Nx, int& Ny, int& Nz, int left_boundary_x_heat, int right_boundary_x_heat,
	int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, int Nx_ac, int Ny_ac, int Nz_ac, vector<Interval> & new_interv_z, vector<Point3D> & points_rupture_ac)
{
	// Calculation 1
	// solving the equation of motion

	double a = 1.5; //к-т в искусственной вязкости
	double slag1 = 0;
	double slag2 = 0;
	double q = 0;
	vector<Point3D>::iterator it;

	/*for (int j = 0; j < Ny; j++)//????????????????????????????????????????????????????????????????????????
	{// Или тут рассчитываются только внутренние узлы???????????????????????????????????И почему тут часть кода рабоатет???????????
		for (int l = 0; l < Nz; l++)
		{
			if ((j - 1) < 0 || (j + 1) > (Ny - 1) || (l - 1) < 0 || (l + 1) > (Nz - 1))
			{
				//cout << j << "  " << l << endl;
				//b2x[0][j][l] = b1x[0][j][l] + dt*CC0*(e1[1][j][l] - e1[0][j][l])*(-(tmp - a1y[0][j][l])*(tmp - a1z[0][j][l]) + (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((a1y[1][j][l] - a1y[0][j][l])*(tmp - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
				//b2x[0][j][l] = b2x[0][j][l] + dt*CC0*(e1[0][j][l] - tmp)*((tmp - a1y[0][j][l])*(a1z[1][j][l] - a1z[0][j][l]) - (tmp - a1y[0][j][l])*(tmp - a1z[0][j][l])) / (dx*dy*dz);
			}
			else {

				// 1-е слагаемое

				if ((b1x[1][j][l] - b1x[0][j][l]) < 0)
				{
					slag1 = pow(b1x[1][j][l] - b1x[0][j][l], 2) / (V2[1][j][l] + V1[0][j][l]);
				}
				else
				{
					slag1 = 0;
				}

				/*if ((b1x[i][j][l] - b1x[i - 1][j][l]) < 0) //???
				{
					slag2 = pow(b1x[i][j][l] - b1x[i - 1][j][l], 2) / (V2[i - 1][j - 1][l - 1] + V1[i - 1][j - 1][l - 1]);
				}
				else
				{
					slag2 = 0;
				}

				q = 2 * pow(a, 2) * (slag1 - slag2);
				b2x[0][j][l] = b1x[0][j][l] + dt * CC0 * (e1[1][j][l] - e1[0][j][l]) * (-(a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) + (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz);



				// 2-е слагаемое

				if ((b1x[0][j + 1][l] - b1x[0][j][l]) < 0)
				{
					slag1 = pow(b1x[0][j + 1][l] - b1x[0][j][l], 2) / (V2[0][j][l] + V1[0][j][l]);
				}
				else
				{
					slag1 = 0;
				}

				if ((b1x[0][j][l] - b1x[0][j - 1][l]) < 0)  // ??
				{
					slag2 = pow(b1x[0][j][l] - b1x[0][j - 1][l], 2) / (V2[0][j - 1][l] + V1[0][j - 1][l]);
				}
				else
				{
					slag2 = 0;
				}

				q = 2 * pow(a, 2) * (slag1 - slag2);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j - 1][l]) * ((a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j][l + 1] - a1z[0][j][l]) - (a1y[0][j][l + 1] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l])) / (dx * dy * dz); // было a1z[0][j][l + 1] - a1z[0][j][l])



				// 3-е слагаемое

				if ((b1x[0][j][l + 1] - b1x[0][j][l]) < 0)
				{
					slag1 = pow(b1x[0][j][l + 1] - b1x[0][j][l], 2) / (V2[0][j][l] + V1[0][j][l]);
				}
				else
				{
					slag1 = 0;
				}

				if ((b1x[0][j][l] - b1x[0][j][l - 1]) < 0)//???
				{
					slag2 = pow(b1x[0][j][l] - b1x[0][j][l - 1], 2) / (V2[0][j][l - 1] + V1[0][j][l - 1]);
				}
				else
				{
					slag2 = 0;
				}

				q = 2 * pow(a, 2) * (slag1 - slag2);
				b2x[0][j][l] = b2x[0][j][l] + dt * CC0 * (e1[0][j][l] - e1[0][j][l - 1]) * ((a1y[0][j + 1][l] - a1y[0][j][l]) * (a1z[1][j][l] - a1z[0][j][l]) - (a1y[1][j][l] - a1y[0][j][l]) * (a1z[0][j + 1][l] - a1z[0][j][l])) / (dx * dy * dz); // было a1y[0][j+1][l] - a1y[0][j][l]

			}
		}
	}*/

	// расчет b2x в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					// 1-е слагаемое

					if ((b1x[i + 1][j][l] - b1x[i][j][l]) < 0)
					{
						slag1 = pow(b1x[i + 1][j][l] - b1x[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x[i][j][l] - b1x[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1x[i][j][l] - b1x[i - 1][j][l], 2) / (V2[i - 1][j][l] + V1[i - 1][j][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2x[i][j][l] = b1x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l] + q) * (-(a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) + (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);


					// 2-е слагаемое

					if ((b1x[i][j + 1][l] - b1x[i][j][l]) < 0)
					{
						slag1 = pow(b1x[i][j + 1][l] - b1x[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x[i][j][l] - b1x[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1x[i][j][l] - b1x[i][j - 1][l], 2) / (V2[i][j - 1][l] + V1[i][j - 1][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);																																					//	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l] + q) * ((a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1y[i][j][l + 1] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);


					// 3-е слагаемое

					if ((b1x[i][j][l + 1] - b1x[i][j][l]) < 0)
					{
						slag1 = pow(b1x[i][j][l + 1] - b1x[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x[i][j][l] - b1x[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1x[i][j][l] - b1x[i][j][l - 1], 2) / (V2[i][j][l - 1] + V1[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);																														//a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
					b2x[i][j][l] = b2x[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1] + q) * ((a1y[i][j + 1][l] - a1y[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1y[i + 1][j][l] - a1y[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);

				}
			}
		}
	}

	// расчет b2x в акустической сетке в зоне отрыва
	it = points_rupture_ac.begin();
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == l)
					{
						b2x_ac[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
						it++;
					}
					else
					{
						if (b1x_ac[i + 1][j][l] == 1e-16 || b1x_ac[i - 1][j][l] == 1e-16 || b1x_ac[i][j + 1][l] == 1e-16 || b1x_ac[i][j - 1][l] == 1e-16 || b1x_ac[i][j][l + 1] == 1e-16 || b1x_ac[i][j][l - 1] == 1e-16 ||
							a1y_ac[i + 1][j][l] == 1e-16 || a1y_ac[i - 1][j][l] == 1e-16 || a1y_ac[i][j + 1][l] == 1e-16 || a1y_ac[i][j - 1][l] == 1e-16 || a1y_ac[i][j][l + 1] == 1e-16 || a1y_ac[i][j][l - 1] == 1e-16 ||
							a1z_ac[i + 1][j][l] == 1e-16 || a1z_ac[i - 1][j][l] == 1e-16 || a1z_ac[i][j + 1][l] == 1e-16 || a1z_ac[i][j - 1][l] == 1e-16 || a1z_ac[i][j][l + 1] == 1e-16 || a1z_ac[i][j][l - 1] == 1e-16)
						{
						}
						else
						{
							// 1-е слагаемое
							if ((b1x_ac[i + 1][j][l] - b1x_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1x_ac[i + 1][j][l] - b1x_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1x_ac[i][j][l] - b1x_ac[i - 1][j][l]) < 0)
							{
								slag2 = pow(b1x_ac[i][j][l] - b1x_ac[i - 1][j][l], 2) / (V2_ac[i - 1][j][l] + V1_ac[i - 1][j][l]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2x_ac[i][j][l] = b1x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l] + q) * (-(a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) + (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

							// 2-е слагаемое

							if ((b1x_ac[i][j + 1][l] - b1x_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1x_ac[i][j + 1][l] - b1x_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1x_ac[i][j][l] - b1x_ac[i][j - 1][l]) < 0)
							{
								slag2 = pow(b1x_ac[i][j][l] - b1x_ac[i][j - 1][l], 2) / (V2_ac[i][j - 1][l] + V1_ac[i][j - 1][l]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);																																					//	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
							b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l] + q) * ((a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

							// 3-е слагаемое

							if ((b1x_ac[i][j][l + 1] - b1x_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1x_ac[i][j][l + 1] - b1x_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1x_ac[i][j][l] - b1x_ac[i][j][l - 1]) < 0)
							{
								slag2 = pow(b1x_ac[i][j][l] - b1x_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);																														//a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
							b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1] + q) * ((a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						}
					}
				}
				else
				{
					if ((b1x_ac[i + 1][j][l] - b1x_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1x_ac[i + 1][j][l] - b1x_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x_ac[i][j][l] - b1x_ac[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1x_ac[i][j][l] - b1x_ac[i - 1][j][l], 2) / (V2_ac[i - 1][j][l] + V1_ac[i - 1][j][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2x_ac[i][j][l] = b1x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l] + q) * (-(a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) + (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

					// 2-е слагаемое

					if ((b1x_ac[i][j + 1][l] - b1x_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1x_ac[i][j + 1][l] - b1x_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x_ac[i][j][l] - b1x_ac[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1x_ac[i][j][l] - b1x_ac[i][j - 1][l], 2) / (V2_ac[i][j - 1][l] + V1_ac[i][j - 1][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);																																					//	a1z[i + 1][j][l] - a1z[i][j][l] ???????????? а было a1z[i][j][l + 1] - a1z[i][j][l])
					b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l] + q) * ((a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

					// 3-е слагаемое

					if ((b1x_ac[i][j][l + 1] - b1x_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1x_ac[i][j][l + 1] - b1x_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1x_ac[i][j][l] - b1x_ac[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1x_ac[i][j][l] - b1x_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);																														//a1y[i + 1][j][l] - a1y[i][j][l] ??????????? а было  a1y[i][j + 1][l] - a1y[i][j][l]
					b2x_ac[i][j][l] = b2x_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1] + q) * ((a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

				}
			}
		}
	}


	// расчет b2y в тепловой сетке не включая область отрыва/расплава 
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					// 1-е слагаемое

					if ((b1y[i + 1][j][l] - b1y[i][j][l]) < 0)
					{
						slag1 = pow(b1y[i + 1][j][l] - b1y[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y[i][j][l] - b1y[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1y[i][j][l] - b1y[i - 1][j][l], 2) / (V2[i - 1][j][l] + V1[i - 1][j][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y[i][j][l] = b1y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i - 1][j][l] + q) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l])) / (dx * dy * dz);


					// 2-е слагаемое

					if ((b1y[i][j + 1][l] - b1y[i][j][l]) < 0)
					{
						slag1 = pow(b1y[i][j + 1][l] - b1y[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y[i][j][l] - b1y[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1y[i][j][l] - b1y[i][j - 1][l], 2) / (V2[i][j - 1][l] + V1[i][j - 1][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j - 1][l] + q) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j][l + 1] - a1z[i][j][l])) / (dx * dy * dz);


					// 3-е слагаемое

					if ((b1y[i][j][l + 1] - b1y[i][j][l]) < 0)
					{
						slag1 = pow(b1y[i][j][l + 1] - b1y[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y[i][j][l] - b1y[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1y[i][j][l] - b1y[i][j][l - 1], 2) / (V2[i][j][l - 1] + V1[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y[i][j][l] = b2y[i][j][l] + dt * CC0 * (e1[i][j][l] - e1[i][j][l - 1] + q) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1z[i][j + 1][l] - a1z[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1z[i + 1][j][l] - a1z[i][j][l])) / (dx * dy * dz);

				}
			}
		}
	}

	// расчет b2y в акустической сетке в зоне отрыва
	it = points_rupture_ac.begin();
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == l)
					{
						b2y_ac[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
						it++;
					}
					else
					{
						if (b1y_ac[i + 1][j][l] == 1e-16 || b1y_ac[i - 1][j][l] == 1e-16 || b1y_ac[i][j + 1][l] == 1e-16 || b1y_ac[i][j - 1][l] == 1e-16 || b1y_ac[i][j][l + 1] == 1e-16 || b1y_ac[i][j][l - 1] == 1e-16 ||
							a1x_ac[i + 1][j][l] == 1e-16 || a1x_ac[i - 1][j][l] == 1e-16 || a1x_ac[i][j + 1][l] == 1e-16 || a1x_ac[i][j - 1][l] == 1e-16 || a1x_ac[i][j][l + 1] == 1e-16 || a1x_ac[i][j][l - 1] == 1e-16 ||
							a1z_ac[i + 1][j][l] == 1e-16 || a1z_ac[i - 1][j][l] == 1e-16 || a1z_ac[i][j + 1][l] == 1e-16 || a1z_ac[i][j - 1][l] == 1e-16 || a1z_ac[i][j][l + 1] == 1e-16 || a1z_ac[i][j][l - 1] == 1e-16)
						{
						}
						else
						{
							// 1-е слагаемое
							if ((b1y_ac[i + 1][j][l] - b1y_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1y_ac[i + 1][j][l] - b1y_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1y_ac[i][j][l] - b1y_ac[i - 1][j][l]) < 0)
							{
								slag2 = pow(b1y_ac[i][j][l] - b1y_ac[i - 1][j][l], 2) / (V2_ac[i - 1][j][l] + V1_ac[i - 1][j][l]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2y_ac[i][j][l] = b1y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l] + q) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

							// 2-е слагаемое

							if ((b1y_ac[i][j + 1][l] - b1y_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1y_ac[i][j + 1][l] - b1y_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1y_ac[i][j][l] - b1y_ac[i][j - 1][l]) < 0)
							{
								slag2 = pow(b1y_ac[i][j][l] - b1y_ac[i][j - 1][l], 2) / (V2_ac[i][j - 1][l] + V1_ac[i][j - 1][l]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l] + q) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

							// 3-е слагаемое

							if ((b1y_ac[i][j][l + 1] - b1y_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1y_ac[i][j][l + 1] - b1y_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1y_ac[i][j][l] - b1y_ac[i][j][l - 1]) < 0)
							{
								slag2 = pow(b1y_ac[i][j][l] - b1y_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1] + q) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						}
					}
				}
				else
				{
					// 1-е слагаемое

					if ((b1y_ac[i + 1][j][l] - b1y_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1y_ac[i + 1][j][l] - b1y_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y_ac[i][j][l] - b1y_ac[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1y_ac[i][j][l] - b1y_ac[i - 1][j][l], 2) / (V2_ac[i - 1][j][l] + V1_ac[i - 1][j][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y_ac[i][j][l] = b1y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i - 1][j][l] + q) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

					// 2-е слагаемое

					if ((b1y_ac[i][j + 1][l] - b1y_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1y_ac[i][j + 1][l] - b1y_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y_ac[i][j][l] - b1y_ac[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1y_ac[i][j][l] - b1y_ac[i][j - 1][l], 2) / (V2_ac[i][j - 1][l] + V1_ac[i][j - 1][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j - 1][l] + q) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j][l + 1] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

					// 3-е слагаемое

					if ((b1y_ac[i][j][l + 1] - b1y_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1y_ac[i][j][l + 1] - b1y_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1y_ac[i][j][l] - b1y_ac[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1y_ac[i][j][l] - b1y_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2y_ac[i][j][l] = b2y_ac[i][j][l] + dt * CC0 * (e1_ac[i][j][l] - e1_ac[i][j][l - 1] + q) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1z_ac[i][j + 1][l] - a1z_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1z_ac[i + 1][j][l] - a1z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

				}
			}
		}
	}


	// расчет b2z(граница) в тепловой сетке не включая область отрыва/расплава ()
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			//for (int l = 1; l < Nz - 1; l++)
			{
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					// 1-е слагаемое

					if ((b1z[i + 1][j][0] - b1z[i][j][0]) < 0)
					{
						slag1 = pow(b1z[i + 1][j][0] - b1z[i][j][0], 2) / (V2[i][j][0] + V1[i][j][0]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][0] - b1z[i - 1][j][0]) < 0)
					{
						slag2 = pow(b1z[i][j][0] - b1z[i - 1][j][0], 2) / (V2[i - 1][j][0] + V1[i - 1][j][0]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][0] = b1z[i][j][0] + dt * (e1[i][j][0] - e1[i - 1][j][0] + q) * ((a1x[i][j][1] - a1x[i][j][0]) * (a1y[i][j + 1][0] - a1y[i][j][0]) - (a1x[i][j + 1][0] - a1x[i][j][0]) * (a1y[i][j][1] - a1y[i][j][0])) / (dx * dy * dz);

					// 2-е слагаемое

					if ((b1z[i][j + 1][0] - b1z[i][j][0]) < 0)
					{
						slag1 = pow(b1z[i][j + 1][0] - b1z[i][j][0], 2) / (V2[i][j][0] + V1[i][j][0]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][0] - b1z[i][j - 1][0]) < 0)
					{
						slag2 = pow(b1z[i][j][0] - b1z[i][j - 1][0], 2) / (V2[i][j - 1][0] + V1[i][j - 1][0]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][0] = b2z[i][j][0] + dt * (e1[i][j][0] - e1[i][j - 1][0] + q) * ((a1x[i + 1][j][0] - a1x[i][j][0]) * (a1y[i][j][1] - a1y[i][j][0]) - (a1x[i][j][1] - a1x[i][j][0]) * (a1y[i + 1][j][0] - a1y[i][j][0])) / (dx * dy * dz);

					// 3-е слагаемое

					if ((b1z[i][j][1] - b1z[i][j][0]) < 0) // l = 0
					{
						slag1 = pow(b1z[i][j][1] - b1z[i][j][0], 2) / (V2[i][j][0] + V1[i][j][0]);
					}
					else
					{
						slag1 = 0;
					}

					/*if ((b1z[i][j][l] - b1z[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i][j][l - 1], 2) / (V2[i][j][l - 1] + V1[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}*/

					slag2 = 0;

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][0] = b2z[i][j][0] + dt * (e1[i][j][1] - e1[i][j][0] + q) * ((a1x[i][j + 1][0] - a1x[i][j][0]) * (a1y[i + 1][j][0] - a1y[i][j][0]) - (a1x[i + 1][j][0] - a1x[i][j][0]) * (a1y[i][j + 1][0] - a1y[i][j][0])) / (dx * dy * dz);

				}
			}
		}
	}

	// расчет b2z(граница) в акустической сетке в зоне отрыва 
	it = points_rupture_ac.begin();
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			//for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == 0)
					{
						b2z_ac[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
						it++;
					}
					else
					{
						if (b1z_ac[i + 1][j][0] == 1e-16 || b1z_ac[i - 1][j][0] == 1e-16 || b1z_ac[i][j + 1][0] == 1e-16 || b1z_ac[i][j - 1][0] == 1e-16 || b1z_ac[i][j][1] == 1e-16 ||
							a1x_ac[i + 1][j][0] == 1e-16 || a1x_ac[i - 1][j][0] == 1e-16 || a1x_ac[i][j + 1][0] == 1e-16 || a1x_ac[i][j - 1][0] == 1e-16 || a1x_ac[i][j][1] == 1e-16 ||
							a1y_ac[i + 1][j][0] == 1e-16 || a1y_ac[i - 1][j][0] == 1e-16 || a1y_ac[i][j + 1][0] == 1e-16 || a1y_ac[i][j - 1][0] == 1e-16 || a1y_ac[i][j][1] == 1e-16)
						{
						}
						else
						{
							// 1-е слагаемое
							if ((b1z_ac[i + 1][j][0] - b1z_ac[i][j][0]) < 0)
							{
								slag1 = pow(b1z_ac[i + 1][j][0] - b1z_ac[i][j][0], 2) / (V2_ac[i][j][0] + V1_ac[i][j][0]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1z_ac[i][j][0] - b1z_ac[i - 1][j][0]) < 0)
							{
								slag2 = pow(b1z_ac[i][j][0] - b1z_ac[i - 1][j][0], 2) / (V2_ac[i - 1][j][0] + V1_ac[i - 1][j][0]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2z_ac[i][j][0] = b1z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i - 1][j][0] + q) * ((a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0]) - (a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);

							// 2-е слагаемое

							if ((b1z_ac[i][j + 1][0] - b1z_ac[i][j][0]) < 0)
							{
								slag1 = pow(b1z_ac[i][j + 1][0] - b1z_ac[i][j][0], 2) / (V2_ac[i][j][0] + V1_ac[i][j][0]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1z_ac[i][j][0] - b1z_ac[i][j - 1][0]) < 0)
							{
								slag2 = pow(b1z_ac[i][j][0] - b1z_ac[i][j - 1][0], 2) / (V2_ac[i][j - 1][0] + V1_ac[i][j - 1][0]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i][j - 1][0] + q) * ((a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0]) - (a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);

							// 3-е слагаемое

							if ((b1z_ac[i][j][1] - b1z_ac[i][j][0]) < 0)
							{
								slag1 = pow(b1z_ac[i][j][1] - b1z_ac[i][j][0], 2) / (V2_ac[i][j][0] + V1_ac[i][j][0]);
							}
							else
							{
								slag1 = 0;
							}

							/*if ((b1z_ac[i][j][l] - b1z_ac[i][j][l - 1]) < 0)
							{
								slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
							}
							else
							{
								slag2 = 0;
							}*/

							slag2 = 0;

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][1] - e1_ac[i][j][0] + q) * ((a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0]) - (a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);
						}
					}
				}
				else
				{
					// 1-е слагаемое

					if ((b1z_ac[i + 1][j][0] - b1z_ac[i][j][0]) < 0)
					{
						slag1 = pow(b1z_ac[i + 1][j][0] - b1z_ac[i][j][0], 2) / (V2_ac[i][j][0] + V1_ac[i][j][0]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z_ac[i][j][0] - b1z_ac[i - 1][j][0]) < 0)
					{
						slag2 = pow(b1z_ac[i][j][0] - b1z_ac[i - 1][j][0], 2) / (V2_ac[i - 1][j][0] + V1_ac[i - 1][j][0]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z_ac[i][j][0] = b1z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i - 1][j][0] + q) * ((a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0]) - (a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);

					// 2-е слагаемое

					if ((b1z_ac[i][j + 1][0] - b1z_ac[i][j][0]) < 0)
					{
						slag1 = pow(b1z_ac[i][j + 1][0] - b1z_ac[i][j][0], 2) / (V2_ac[i][j][0] + V1_ac[i][j][0]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z_ac[i][j][0] - b1z_ac[i][j - 1][0]) < 0)
					{
						slag2 = pow(b1z_ac[i][j][0] - b1z_ac[i][j - 1][0], 2) / (V2_ac[i][j - 1][0] + V1_ac[i][j - 1][0]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][0] - e1_ac[i][j - 1][0] + q) * ((a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j][1] - a1y_ac[i][j][0]) - (a1x_ac[i][j][1] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);

					// 3-е слагаемое

					if ((b1z_ac[i][j][1] - b1z_ac[i][j][0]) < 0)
					{
						slag1 = pow(b1z_ac[i][j][1] - b1z_ac[i][j][0], 2) / (V2_ac[i][j][0] + V1_ac[i][j][0]);
					}
					else
					{
						slag1 = 0;
					}

					/*if ((b1z_ac[i][j][l] - b1z_ac[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}*/

					slag2 = 0;

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z_ac[i][j][0] = b2z_ac[i][j][0] + dt * (e1_ac[i][j][1] - e1_ac[i][j][0] + q) * ((a1x_ac[i][j + 1][0] - a1x_ac[i][j][0]) * (a1y_ac[i + 1][j][0] - a1y_ac[i][j][0]) - (a1x_ac[i + 1][j][0] - a1x_ac[i][j][0]) * (a1y_ac[i][j + 1][0] - a1y_ac[i][j][0])) / (dx_ac * dy_ac * dz_ac);

				}
			}
		}
	}


	// расчет b2z в тепловой сетке не включая область отрыва/расплава ()
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					// 1-е слагаемое

					if ((b1z[i + 1][j][l] - b1z[i][j][l]) < 0)
					{
						slag1 = pow(b1z[i + 1][j][l] - b1z[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][l] - b1z[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i - 1][j][l], 2) / (V2[i - 1][j][l] + V1[i - 1][j][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][l] = b1z[i][j][l] + dt * (e1[i][j][l] - e1[i - 1][j][l] + q) * ((a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l]) - (a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l])) / (dx * dy * dz);



					// 2-е слагаемое

					if ((b1z[i][j + 1][l] - b1z[i][j][l]) < 0)
					{
						slag1 = pow(b1z[i][j + 1][l] - b1z[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][l] - b1z[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i][j - 1][l], 2) / (V2[i][j - 1][l] + V1[i][j - 1][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j - 1][l] + q) * ((a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j][l + 1] - a1y[i][j][l]) - (a1x[i][j][l + 1] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l])) / (dx * dy * dz);


					// 3-е слагаемое

					if ((b1z[i][j][l + 1] - b1z[i][j][l]) < 0)
					{
						slag1 = pow(b1z[i][j][l + 1] - b1z[i][j][l], 2) / (V2[i][j][l] + V1[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z[i][j][l] - b1z[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1z[i][j][l] - b1z[i][j][l - 1], 2) / (V2[i][j][l - 1] + V1[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z[i][j][l] = b2z[i][j][l] + dt * (e1[i][j][l] - e1[i][j][l - 1] + q) * ((a1x[i][j + 1][l] - a1x[i][j][l]) * (a1y[i + 1][j][l] - a1y[i][j][l]) - (a1x[i + 1][j][l] - a1x[i][j][l]) * (a1y[i][j + 1][l] - a1y[i][j][l])) / (dx * dy * dz);

				}
			}
		}
	}

	// расчет b2z в акустической сетке в зоне отрыва 
	it = points_rupture_ac.begin();
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it).index_x == i && (*it).index_y == j && (*it).index_z == l)
					{
						b2z_ac[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
						it++;
					}
					else
					{
						if (b1z_ac[i + 1][j][l] == 1e-16 || b1z_ac[i - 1][j][l] == 1e-16 || b1z_ac[i][j + 1][l] == 1e-16 || b1z_ac[i][j - 1][l] == 1e-16 || b1z_ac[i][j][l + 1] == 1e-16 || b1z_ac[i][j][l - 1] == 1e-16 ||
							a1x_ac[i + 1][j][l] == 1e-16 || a1x_ac[i - 1][j][l] == 1e-16 || a1x_ac[i][j + 1][l] == 1e-16 || a1x_ac[i][j - 1][l] == 1e-16 || a1x_ac[i][j][l + 1] == 1e-16 || a1x_ac[i][j][l - 1] == 1e-16 ||
							a1y_ac[i + 1][j][l] == 1e-16 || a1y_ac[i - 1][j][l] == 1e-16 || a1y_ac[i][j + 1][l] == 1e-16 || a1y_ac[i][j - 1][l] == 1e-16 || a1y_ac[i][j][l + 1] == 1e-16 || a1y_ac[i][j][l - 1] == 1e-16)
						{
						}
						else
						{
							// 1-е слагаемое
							if ((b1z_ac[i + 1][j][l] - b1z_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1z_ac[i + 1][j][l] - b1z_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1z_ac[i][j][l] - b1z_ac[i - 1][j][l]) < 0)
							{
								slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i - 1][j][l], 2) / (V2_ac[i - 1][j][l] + V1_ac[i - 1][j][l]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2z_ac[i][j][l] = b1z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i - 1][j][l] + q) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

							// 2-е слагаемое

							if ((b1z_ac[i][j + 1][l] - b1z_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1z_ac[i][j + 1][l] - b1z_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1z_ac[i][j][l] - b1z_ac[i][j - 1][l]) < 0)
							{
								slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i][j - 1][l], 2) / (V2_ac[i][j - 1][l] + V1_ac[i][j - 1][l]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j - 1][l] + q) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

							// 3-е слагаемое

							if ((b1z_ac[i][j][l + 1] - b1z_ac[i][j][l]) < 0)
							{
								slag1 = pow(b1z_ac[i][j][l + 1] - b1z_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
							}
							else
							{
								slag1 = 0;
							}

							if ((b1z_ac[i][j][l] - b1z_ac[i][j][l - 1]) < 0)
							{
								slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
							}
							else
							{
								slag2 = 0;
							}

							q = 2 * pow(a, 2) * (slag1 - slag2);
							b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j][l - 1] + q) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						}
					}
				}
				else
				{
					// 1-е слагаемое

					if ((b1z_ac[i + 1][j][l] - b1z_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1z_ac[i + 1][j][l] - b1z_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z_ac[i][j][l] - b1z_ac[i - 1][j][l]) < 0)
					{
						slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i - 1][j][l], 2) / (V2_ac[i - 1][j][l] + V1_ac[i - 1][j][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z_ac[i][j][l] = b1z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i - 1][j][l] + q) * ((a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l]) - (a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

					// 2-е слагаемое

					if ((b1z_ac[i][j + 1][l] - b1z_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1z_ac[i][j + 1][l] - b1z_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z_ac[i][j][l] - b1z_ac[i][j - 1][l]) < 0)
					{
						slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i][j - 1][l], 2) / (V2_ac[i][j - 1][l] + V1_ac[i][j - 1][l]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j - 1][l] + q) * ((a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j][l + 1] - a1y_ac[i][j][l]) - (a1x_ac[i][j][l + 1] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

					// 3-е слагаемое

					if ((b1z_ac[i][j][l + 1] - b1z_ac[i][j][l]) < 0)
					{
						slag1 = pow(b1z_ac[i][j][l + 1] - b1z_ac[i][j][l], 2) / (V2_ac[i][j][l] + V1_ac[i][j][l]);
					}
					else
					{
						slag1 = 0;
					}

					if ((b1z_ac[i][j][l] - b1z_ac[i][j][l - 1]) < 0)
					{
						slag2 = pow(b1z_ac[i][j][l] - b1z_ac[i][j][l - 1], 2) / (V2_ac[i][j][l - 1] + V1_ac[i][j][l - 1]);
					}
					else
					{
						slag2 = 0;
					}

					q = 2 * pow(a, 2) * (slag1 - slag2);
					b2z_ac[i][j][l] = b2z_ac[i][j][l] + dt * (e1_ac[i][j][l] - e1_ac[i][j][l - 1] + q) * ((a1x_ac[i][j + 1][l] - a1x_ac[i][j][l]) * (a1y_ac[i + 1][j][l] - a1y_ac[i][j][l]) - (a1x_ac[i + 1][j][l] - a1x_ac[i][j][l]) * (a1y_ac[i][j + 1][l] - a1y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);

				}
			}
		}
	}

}

void Calculation2(double*** a1x, double*** a1y, double*** a1z, double*** a2x, double*** a2y, double*** a2z, double*** b2x, double*** b2y, double*** b2z, double*** V2, double*** V1, double*** V00, double*** e2,
	double*** tmpe1, double*** tmpi1, double*** a1x_ac, double*** a2x_ac, double*** a1y_ac, double*** a2y_ac, double*** a1z_ac, double*** a2z_ac, double*** b2x_ac, double*** b2y_ac, double*** b2z_ac,
	double*** V2_ac, double*** V1_ac, double*** V00_ac, double*** tmpe_ac1_new, double*** tmpi_ac1_new, double*** e2_ac_new, double*** e1_ac_new, double& dx, double& dy, double& dz, int Nx_part_ht, int  Ny_part_ht, double& dx_ac,
	double& dy_ac, double& dz_ac, double& dt, double& CC0,
	int& Nx, int& Ny, int& Nz, int left_boundary_x_heat, int right_boundary_x_heat,
	int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat, int Nx_ac, int Ny_ac, int Nz_ac, double& V0, P0V & p0v, /*double*** tmpe_for_transform_new,*/
	/*double*** Tei_old_xy_new_z, double*** Tei_old_y_new_xz, double*** Te_acoustic, double*** Tei_old_xy_new_z_for_total_ac, double*** Tei_old_y_new_xz_for_total_ac,*/ /*double** Te_old_x_new_y, double** Te_old_y_new_x, double** Te_old_x_new_z, double** Te_old_z_new_x,*/ int coeff_big_x, int coeff_big_y, int coeff_big_z,
	vector<Melting> & Melt_metal, Metal mt, Parametrs & param,
	Splayn & spl_C_e_on_T, Splayn & spl_Te, Splayn & spl_Ti, Splayn & spl_e, /*double*** Tei_old_xy_new_z, double*** Tei_old_y_new_xz,*/ int& current_count_frame,
	vector<Point3D> & points_rupture_ac, vector<Point3D> & points_rupture_heat, vector<Interval> & new_interv_z, Splayn spl_C_l_on_T)
{
	// Calculation of 
 // Calculatuion 2

	//for (int i = 0; i < Nx; i++)
	//{
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		for (int l = 0; l < Nz; l++) // расчет в тепловой задаче
	//		{
	//			//calculation of the change in the Euler coordinates of Lagrangian particles
	//			a2x[i][j][l] = a1x[i][j][l] + dt * b2x[i][j][l] / CC0;
	//			a2y[i][j][l] = a1y[i][j][l] + dt * b2y[i][j][l] / CC0;
	//			a2z[i][j][l] = a1z[i][j][l] + dt * b2z[i][j][l];
	//		}
	//	}
	//}
	vector<Point3D>::iterator it_ac, it_heat;
	vector<int> for_close_and_open = { 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43,44,45,46,47,48, 49 };

	// расчет Эйлеровых координат в тепловой сетке не включая область отрыва/расплава
	for (int i = 0; i < Nx; i++)
	{
		for (int j = 0; j < Ny; j++)
		{
			for (int k = 0; k < Nz; k++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && k < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					a2x[i][j][k] = a1x[i][j][k] + dt * b2x[i][j][k] / CC0;
					a2y[i][j][k] = a1y[i][j][k] + dt * b2y[i][j][k] / CC0;
					a2z[i][j][k] = a1z[i][j][k] + dt * b2z[i][j][k];
				}
			}
		}
	}

	it_ac = points_rupture_ac.begin();

	// // расчет Эйлеровых координат в акустической сетке в зоне отрыва
	for (int i = 0; i < Nx_ac; i++)
	{
		for (int j = 0; j < Ny_ac; j++)
		{
			for (int k = 0; k < Nz_ac; k++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == k)
					{
						a2x_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						a2y_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						a2z_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						a2x_ac[i][j][k] = a1x_ac[i][j][k] + dt * b2x_ac[i][j][k] / CC0;
						a2y_ac[i][j][k] = a1y_ac[i][j][k] + dt * b2y_ac[i][j][k] / CC0;
						a2z_ac[i][j][k] = a1z_ac[i][j][k] + dt * b2z_ac[i][j][k];
					}
				}
				else
				{
					a2x_ac[i][j][k] = a1x_ac[i][j][k] + dt * b2x_ac[i][j][k] / CC0;
					a2y_ac[i][j][k] = a1y_ac[i][j][k] + dt * b2y_ac[i][j][k] / CC0;
					a2z_ac[i][j][k] = a1z_ac[i][j][k] + dt * b2z_ac[i][j][k];
				}
			}
		}
	}

	//GnuPlot plt_TE_ht(50);
	//GnuPlot plt_TE_ac(50);

	//cout << " Calc 2 " << endl;
	//plt_TE_ht.SetParametrsOnPlotColor(31, "a1x", "z,mkm", "x,mkm", (1e+4 * dz * Nz / 5e+5), (1e+4 * 1e-2 * dx * Nx));
	//plt_TE_ht.SetDataOnPlotColor3D(31, 100, Nx, Ny, Nz, 1e+4 * dx * 1e-2, 1e+4 * dy * 1e-2, 1e+4 * dz / 5e+5, a1x, 1., Ny / 2, xz);
	//plt_TE_ht.ShowDataOnPlotColor(31, "file a1x in ht (x,y div 2,z) ", true);
	//Sleep(1000);
	//plt_TE_ht.SetParametrsOnPlotColor(32, "a2x", "z,mkm", "x,mkm", (1e+4* dz* Nz / 5e+5), (1e+4* 1e-2* dx* Nx));
	//	plt_TE_ht.SetDataOnPlotColor3D(32, 100, Nx, Ny, Nz, 1e+4* dx* 1e-2, 1e+4* dy* 1e-2, 1e+4* dz / 5e+5, a2x, 1., Ny / 2, xz);
	//	plt_TE_ht.ShowDataOnPlotColor(32, "file a2x in ht (x,y div 2,z) ", true);
	//	Sleep(1000);
	//	plt_TE_ac.SetParametrsOnPlotColor(33, "a1x", "z,mkm", "x,mkm", (1e+4 * dz_ac * Nz_ac / 5e+5), (1e+4 * 1e-2 * dx_ac * Nx_ac));
	//	plt_TE_ac.SetDataOnPlotColor3D(33, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4 * dx_ac * 1e-2, 1e+4 * dy_ac * 1e-2, 1e+4 * dz_ac / 5e+5, a1x_ac, 1., Ny_ac / 2, xz);
	//	plt_TE_ac.ShowDataOnPlotColor(33, "file a1x_ac part in ac (x,y div 2,z) ", true);
	//	Sleep(1000);
	//	plt_TE_ac.SetParametrsOnPlotColor(34, "a2x", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	//	plt_TE_ac.SetDataOnPlotColor3D(34, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, a2x_ac, 1., Ny_ac / 2, xz);
	//	plt_TE_ac.ShowDataOnPlotColor(34, "file a2x_ac part in ac (x,y div 2,z) ", true);
	//	Sleep(1000);
	//	cout << " Stop a2x " << endl;
	//	system("pause");

	/* сшивка границ координат (1-й этап) (2-я версия)*/
	/*a2x*/
	Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, a2x_ac, a2x, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, a2x_ac, a2x, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// from ht to ac /* сшивка границ координат (2-й этап) */
	TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
		dx_ac, dy_ac, a2x, a2x_ac, /*Te_old_y_new_x,*/ spl_Te);

	TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, a2x, a2x_ac, /*Te_old_z_new_x,*/ spl_Te);


	/*a2y*/
	/* сшивка границ координат (1-й этап) (2-я версия)*/
	Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, a2y_ac, a2y, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, a2y_ac, a2y, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// from ht to ac /* сшивка границ координат (2-й этап) */
	TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
		dx_ac, dy_ac, a2y, a2y_ac, /*Te_old_y_new_x,*/ spl_Te);

	TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, a2y, a2y_ac, /*Te_old_z_new_x,*/ spl_Te);


	/*a2z*/
	/* сшивка границ координат (1-й этап) (2-я версия)*/
	Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, a2z_ac, a2z, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, a2z_ac, a2z, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// from ht to ac /* сшивка границ координат (2-й этап) */
	TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
		dx_ac, dy_ac, a2z, a2z_ac, /*Te_old_y_new_x,*/ spl_Te);

	TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, a2z, a2z_ac, /*Te_old_z_new_x,*/ spl_Te);

	//сшивка координат (внутренняя граница) (старая версия)
	/*TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, a2x_ac, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, a2x, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, a2y_ac, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, a2y, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, a2z_ac, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, a2z, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);*/

	//plt_TE_ht.SetParametrsOnPlotColor(38, "a2x", "z,mkm", "x,mkm", (1e+4* dz* Nz / 5e+5), (1e+4* 1e-2* dx* Nx));
	//	plt_TE_ht.SetDataOnPlotColor3D(38, 100, Nx, Ny, Nz, 1e+4* dx* 1e-2, 1e+4* dy* 1e-2, 1e+4* dz / 5e+5, a2x, 1., Ny / 2, xz);
	//	plt_TE_ht.ShowDataOnPlotColor(38, "file a2x in HT", true);
	//	Sleep(1000);
	//	cout << " STOP STOP STOP STOP STOP" << endl;
	//	system("pause");

	// расчет ур-я непрерывности (объема) в тепловой сетке не включая область отрыва/расплава
	for (int i = 1; i < Nx - 1; i++)
	{
		for (int j = 1; j < Ny - 1; j++)
		{
			for (int l = 1; l < Nz - 1; l++)
			{
				//if (i >= left_boundary_x_heat && i <= (right_boundary_x_heat /*+ 1*/) && j >= left_boundary_y_heat && j <= (right_boundary_y_heat /*+ 1*/) && k <= (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && l < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					V2[i][j][l] = (a2x[i + 1][j][l] - a2x[i][j][l]) * ((a2y[i][j + 1][l] - a2y[i][j][l]) * (a2z[i][j][l + 1] - a2z[i][j][l]) - (a2y[i][j][l + 1] - a2y[i][j][l]) * (a2z[i][j + 1][l] - a2z[i][j][l])) / (dx * dy * dz);
					V2[i][j][l] = V2[i][j][l] - (a2y[i + 1][j][l] - a2y[i][j][l]) * ((a2x[i][j + 1][l] - a2x[i][j][l]) * (a2z[i][j][l + 1] - a2z[i][j][l]) - (a2x[i][j][l + 1] - a2x[i][j][l]) * (a2z[i][j + 1][l] - a2z[i][j][l])) / (dx * dy * dz);
					V2[i][j][l] = V2[i][j][l] + (a2z[i + 1][j][l] - a2z[i][j][l]) * ((a2x[i][j + 1][l] - a2x[i][j][l]) * (a2y[i][j][l + 1] - a2y[i][j][l]) - (a2x[i][j][l + 1] - a2x[i][j][l]) * (a2y[i][j + 1][l] - a2y[i][j][l])) / (dx * dy * dz);
				}
			}
		}
	}

	//plt_TE_ht.SetParametrsOnPlotColor(35, "V2", "z,mkm", "x,mkm", (1e+4* dz* Nz / 5e+5), (1e+4* 1e-2* dx* Nx));
	//plt_TE_ht.SetDataOnPlotColor3D(35, 100, Nx, Ny, Nz, 1e+4* dx* 1e-2, 1e+4* dy* 1e-2, 1e+4* dz / 5e+5, V2, 1., Ny / 2, xz);
	//plt_TE_ht.ShowDataOnPlotColor(35, "file V2 in ht (x,y div 2,z) ", true);
	//Sleep(1000);
	//system("pause");

	/*CopyDataArray3D(a1x_ac, a2x_ac, Nx_ac, Ny_ac, Nz_ac);
	CopyDataArray3D(a1y_ac, a2y_ac, Nx_ac, Ny_ac, Nz_ac);
	CopyDataArray3D(a1z_ac, a2z_ac, Nx_ac, Ny_ac, Nz_ac);*/

	// в 37-ом границы (общие) 32-го  
	///TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, a2x, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);

	//plt_TE_ac.SetParametrsOnPlotColor(37, "a2x", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_acoustic_total / 5e+5), (1e+4* 1e-2* dx_ac* Nx_acoustic_total));
	//plt_TE_ac.SetDataOnPlotColor3D(37, 100, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, Te_acoustic, 1., Ny_ac / 2, xz);
	//plt_TE_ac.ShowDataOnPlotColor(37, "file a2x (for shiv 1) in ac total (x,y div 2,z) ", true);
	//Sleep(1000);
	//system("pause");

	///Copy_Data_From_Big_Array3D_To_Small_Array3D(a2x_ac, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);

	// в 39-ом главнео что границы сохранены
	//plt_TE_ht.SetParametrsOnPlotColor(39, "a2x", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	//plt_TE_ht.SetDataOnPlotColor3D(39, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, a2x_ac, 1., Ny_ac / 2, xz);
	//plt_TE_ht.ShowDataOnPlotColor(39, "file a2x in Ac bez vnutr tolko graniz (for shiv 2) ", true);
	//Sleep(1000);
	//cout << " STOP STOP STOP STOP STOP" << endl;
	//system("pause");
	//
	///Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(a2x_ac, a1x_ac, Nx_ac, Ny_ac, Nz_ac);

	//plt_TE_ht.SetParametrsOnPlotColor(40, "a2x", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	//	plt_TE_ht.SetDataOnPlotColor3D(40, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, a2x_ac, 1., Ny_ac / 2, xz);
	//	plt_TE_ht.ShowDataOnPlotColor(40, "file a2x in Ac vnutr with boundary", true);
	//	Sleep(1000);
	//	cout << " STOP STOP STOP STOP STOP a2x" << endl;
	//	system("pause");

	// ниже 4 действия (1-е выше там копирование) делаем для того, чтобы  в акуст массиве дать значения смежных границы
	// (старая версия сшивки)
	//TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, a2y, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
	//Copy_Data_From_Big_Array3D_To_Small_Array3D(a2y_ac, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);
	//Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(a2y_ac, a1y_ac, Nx_ac, Ny_ac, Nz_ac);

	//TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, a2z, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
	//Copy_Data_From_Big_Array3D_To_Small_Array3D(a2z_ac, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);
	//Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(a2z_ac, a1z_ac, Nx_ac, Ny_ac, Nz_ac);

	it_ac = points_rupture_ac.begin();
	// // расчет ур-я непрерывности (объема) в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int l = 1; l < Nz_ac - 1; l++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == l)
					{
						///V2_ac[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						if (a2x_ac[i + 1][j][l] == 1e-16 || a2x_ac[i - 1][j][l] == 1e-16 || a2x_ac[i][j + 1][l] == 1e-16 || a2x_ac[i][j - 1][l] == 1e-16 || a2x_ac[i][j][l + 1] == 1e-16 || a2x_ac[i][j][l - 1] == 1e-16 ||
							a2y_ac[i + 1][j][l] == 1e-16 || a2y_ac[i - 1][j][l] == 1e-16 || a2y_ac[i][j + 1][l] == 1e-16 || a2y_ac[i][j - 1][l] == 1e-16 || a2y_ac[i][j][l + 1] == 1e-16 || a2y_ac[i][j][l - 1] == 1e-16 ||
							a2z_ac[i + 1][j][l] == 1e-16 || a2z_ac[i - 1][j][l] == 1e-16 || a2z_ac[i][j + 1][l] == 1e-16 || a2z_ac[i][j - 1][l] == 1e-16 || a2z_ac[i][j][l + 1] == 1e-16 || a2z_ac[i][j][l - 1] == 1e-16)
						{
						}
						else
						{
							V2_ac[i][j][l] = (a2x_ac[i + 1][j][l] - a2x_ac[i][j][l]) * ((a2y_ac[i][j + 1][l] - a2y_ac[i][j][l]) * (a2z_ac[i][j][l + 1] - a2z_ac[i][j][l]) - (a2y_ac[i][j][l + 1] - a2y_ac[i][j][l]) * (a2z_ac[i][j + 1][l] - a2z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
							V2_ac[i][j][l] = V2_ac[i][j][l] - (a2y_ac[i + 1][j][l] - a2y_ac[i][j][l]) * ((a2x_ac[i][j + 1][l] - a2x_ac[i][j][l]) * (a2z_ac[i][j][l + 1] - a2z_ac[i][j][l]) - (a2x_ac[i][j][l + 1] - a2x_ac[i][j][l]) * (a2z_ac[i][j + 1][l] - a2z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
							V2_ac[i][j][l] = V2_ac[i][j][l] + (a2z_ac[i + 1][j][l] - a2z_ac[i][j][l]) * ((a2x_ac[i][j + 1][l] - a2x_ac[i][j][l]) * (a2y_ac[i][j][l + 1] - a2y_ac[i][j][l]) - (a2x_ac[i][j][l + 1] - a2x_ac[i][j][l]) * (a2y_ac[i][j + 1][l] - a2y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
						}
					}
				}
				else
				{
					V2_ac[i][j][l] = (a2x_ac[i + 1][j][l] - a2x_ac[i][j][l]) * ((a2y_ac[i][j + 1][l] - a2y_ac[i][j][l]) * (a2z_ac[i][j][l + 1] - a2z_ac[i][j][l]) - (a2y_ac[i][j][l + 1] - a2y_ac[i][j][l]) * (a2z_ac[i][j + 1][l] - a2z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					V2_ac[i][j][l] = V2_ac[i][j][l] - (a2y_ac[i + 1][j][l] - a2y_ac[i][j][l]) * ((a2x_ac[i][j + 1][l] - a2x_ac[i][j][l]) * (a2z_ac[i][j][l + 1] - a2z_ac[i][j][l]) - (a2x_ac[i][j][l + 1] - a2x_ac[i][j][l]) * (a2z_ac[i][j + 1][l] - a2z_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
					V2_ac[i][j][l] = V2_ac[i][j][l] + (a2z_ac[i + 1][j][l] - a2z_ac[i][j][l]) * ((a2x_ac[i][j + 1][l] - a2x_ac[i][j][l]) * (a2y_ac[i][j][l + 1] - a2y_ac[i][j][l]) - (a2x_ac[i][j][l + 1] - a2x_ac[i][j][l]) * (a2y_ac[i][j + 1][l] - a2y_ac[i][j][l])) / (dx_ac * dy_ac * dz_ac);
				}
			}
		}
	}

	/*plt_TE_ht.SetParametrsOnPlotColor(41, "V2", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	plt_TE_ht.SetDataOnPlotColor3D(41, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, V2_ac, 1., Ny_ac / 2, xz);
	plt_TE_ht.ShowDataOnPlotColor(41, "file V2 part in ac (x,y div 2,z)", true);
	Sleep(1000);
	cout << " STOP STOP STOP STOP STOP V2" << endl;
	system("pause");*/
	/**/

	//V2 (1-й этап)
	Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, V2_ac, V2, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, V2_ac, V2, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// from ht to ac /* сшивка границ объема (2-й этап) */
	TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
		dx_ac, dy_ac, V2, V2_ac, /*Te_old_y_new_x,*/ spl_Te);

	TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, V2, V2_ac, /*Te_old_z_new_x,*/ spl_Te);

	/*V1 (1-й этап) */
	/*//Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, V1_ac, V1, Te_old_x_new_y, spl_Te, coeff_big_z,
	//	left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	//Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, V1_ac, V1, Te_old_x_new_z, spl_Te, coeff_big_y, coeff_big_x,
	//	left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	//// from ht to ac /* сшивка границ объема (2-й этап) //
	//TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
	//	dx_ac, dy_ac, V1, V1_ac, Te_old_y_new_x, spl_Te);

	//TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, V1, V1_ac, Te_old_z_new_x, spl_Te);

	///*V00 (1-й этап) //
	//Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, V00_ac, V00, Te_old_x_new_y, spl_Te, coeff_big_z,
	//	left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	//Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, V00_ac, V00, Te_old_x_new_z, spl_Te, coeff_big_y, coeff_big_x,
	//	left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	//// from ht to ac /* сшивка границ объема (2-й этап) //

	//TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
	//	dx_ac, dy_ac, V00, V00_ac, Te_old_y_new_x, spl_Te);

	//TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, V00, V00_ac, Te_old_z_new_x, spl_Te);
	*/

	// сшивка объема V2 (внутр границы для крупной сетки) (старая версия)
	/*TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, V2_ac, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, V2, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// сшивка объема V2 (внеш границы для мелкой сетки)
	CopyDataArray3D(a1x_ac, V2_ac, Nx_ac, Ny_ac, Nz_ac);
	TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, V2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
	Copy_Data_From_Big_Array3D_To_Small_Array3D(V2_ac, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);
	Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(V2_ac, a1x_ac, Nx_ac, Ny_ac, Nz_ac);

	// сшивка объема V1 (внутр границы для крупной сетки)
	TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, V1_ac, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, V1, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// сшивка объема V1 (внеш границы для мелкой сетки)
	CopyDataArray3D(a1x_ac, V1_ac, Nx_ac, Ny_ac, Nz_ac);
	TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, V1, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
	Copy_Data_From_Big_Array3D_To_Small_Array3D(V1_ac, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);
	Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(V1_ac, a1x_ac, Nx_ac, Ny_ac, Nz_ac);

	// сшивка объема V00 (внутр границы для крупной сетки)
	TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, V00_ac, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, V00, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// сшивка объема V00 (внеш границы для мелкой сетки)
	CopyDataArray3D(a1x_ac, V00_ac, Nx_ac, Ny_ac, Nz_ac);
	TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, V00, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
	Copy_Data_From_Big_Array3D_To_Small_Array3D(V00_ac, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);
	Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(V00_ac, a1x_ac, Nx_ac, Ny_ac, Nz_ac);*/

	//plt_TE_ht.SetParametrsOnPlotColor(42, "V2", "z,mkm", "x,mkm", (1e+4* dz* Nz / 5e+5), (1e+4* 1e-2* dx* Nx));
	//plt_TE_ht.SetDataOnPlotColor3D(42, 100, Nx, Ny, Nz, 1e+4* dx* 1e-2, 1e+4* dy* 1e-2, 1e+4* dz / 5e+5, V2, 1., Ny / 2, xz);
	//plt_TE_ht.ShowDataOnPlotColor(42, "file V2 in HT total ", true);
	//Sleep(1000);
	//cout << " STOP STOP STOP STOP STOP V2 V2 V2" << endl;
	//system("pause");
	//system("pause");

	for (int j = 0; j < Nx; j++)
	{
		for (int l = 0; l < Ny; l++)
		{
			V2[j][l][Nz - 1] = V0;
		}
	}

	/*for (int j = 0; j < Ny; j++)
	{
		for (int l = 0; l < Nz; l++)
		{
			e2[0][j][l] = 1e-16;
		}
	}*/

	double delta = 1.;

	// способы расчета дваления для случая когда еще не учитывали отрыв
	/*if ((Ti_acoustic[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta) )
	{
		e2[i][j][k] = (1. - V2[i][j][k]) +
			(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
			(Te_acoustic[i][j][k] - 1.);
	}

	if ((Ti_acoustic[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (Ti_acoustic[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
	{
		//Ci[J/cm3 K] = ((Dependence_C_l_on_T(mt, param.T00 * tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))
		// spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])
		// spl_Ti.GetY(k * dz) = tmpi2

		/*e2[i][j][k] = (1. - V2[i][j][k]) +
			(Melt_metal[mt].gi * param.T00 * (((((Dependence_C_l_on_T(mt, param.T00* tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
			(Te_acoustic[i][j][k] - 1.);// *//*/

		e2[i][j][k] = (1. - V2[i][j][k]) +
			(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
			(Te_acoustic[i][j][k] - 1.);
	}

	if ((Ti_acoustic[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
	{
		e2[i][j][k] = (1. - V2[i][j][k]) +
			(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * Ti_acoustic[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
			(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
			(Te_acoustic[i][j][k] - 1.);
	}*/

	it_ac = points_rupture_ac.begin();
	ofstream fout_point("Proverka_sovpad.txt");

	//  расчет давления в тепловой сетке не включая область отрыва/расплава
	for (int i = 0; i < Nx; i++) // Nz_heat, dz_heat - не нужны
	{
		for (int j = 0; j < Ny; j++)// циклы по узлам акустики
		{
			for (int k = 1; k < Nz; k++)//l=1
			{
				//if (points_rupture_ac.empty()) // вообще ещё нигде не отровало
				if (i > left_boundary_x_heat && i < (right_boundary_x_heat /*+ 1*/) && j > left_boundary_y_heat && j < (right_boundary_y_heat /*+ 1*/) && k < (down_boundary_z_heat /*+ 1*/)) // dx = dy = 5 мкм dz = 2 нм
				{

				}
				else
				{
					if ((tmpi1[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
					{
						e2[i][j][k] = (1. - V2[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
							(tmpi1[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(tmpe1[i][j][k] - 1.);
					}

					if ((tmpi1[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi1[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
					{
						//Ci[J/cm3 K] = ((Dependence_C_l_on_T(mt, param.T00 * tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))
							// spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])
							// spl_Ti.GetY(k * dz) = tmpi2

						/*e2[i][j][k] = (1. - V2[i][j][k]) +
								(Melt_metal[mt].gi * param.T00 * (((((Dependence_C_l_on_T(mt, param.T00* tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
								(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
								(Te_acoustic[i][j][k] - 1.);// */

						e2[i][j][k] = (1. - V2[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * tmpi1[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(tmpi1[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
							(tmpe1[i][j][k] - 1.);
					}

					if ((tmpi1[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
					{
						e2[i][j][k] = (1. - V2[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * tmpi1[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
							(tmpi1[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe1[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
							(tmpe1[i][j][k] - 1.);
					}
				}
			}
		}
	}

	it_ac = points_rupture_ac.begin();
	// // расчет давлеия в акустической сетке в зоне отрыва
	for (int i = 1; i < Nx_ac - 1; i++)
	{
		for (int j = 1; j < Ny_ac - 1; j++)
		{
			for (int k = 1; k < Nz_ac - 1; k++)
			{
				if (!points_rupture_ac.empty())
				{
					if ((*it_ac).index_x == i && (*it_ac).index_y == j && (*it_ac).index_z == k)
					{
						e2_ac_new[(*it_ac).index_x][(*it_ac).index_y][(*it_ac).index_z] = 1e-16;
						it_ac++;
					}
					else
					{
						if ((tmpi_ac1_new[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
						{
							if (tmpi_ac1_new[i][j][k] == 1e-16 || tmpe_ac1_new[i][j][k] == 1e-16)
							{
							}
							else
							{
								e2_ac_new[i][j][k] = (1. - V2_ac[i][j][k]) +
									(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
									(tmpi_ac1_new[i][j][k] - 1.) / V2_ac[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
									(tmpe_ac1_new[i][j][k] - 1.);
							}
						}

						if ((tmpi_ac1_new[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi_ac1_new[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
						{
							//Ci[J/cm3 K] = ((Dependence_C_l_on_T(mt, param.T00 * tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))
								// spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])
								// spl_Ti.GetY(k * dz) = tmpi2

							/*e2[i][j][k] = (1. - V2[i][j][k]) +
									(Melt_metal[mt].gi * param.T00 * (((((Dependence_C_l_on_T(mt, param.T00* tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
									(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
									(Te_acoustic[i][j][k] - 1.);// */

							if (tmpi_ac1_new[i][j][k] == 1e-16 || tmpe_ac1_new[i][j][k] == 1e-16)
							{
							}
							else
							{
								e2_ac_new[i][j][k] = (1. - V2_ac[i][j][k]) +
									(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * tmpi_ac1_new[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
									(tmpi_ac1_new[i][j][k] - 1.) / V2_ac[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
									(tmpe_ac1_new[i][j][k] - 1.);
							}
						}

						if ((tmpi_ac1_new[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
						{
							if (tmpi_ac1_new[i][j][k] == 1e-16 || tmpe_ac1_new[i][j][k] == 1e-16)
							{
							}
							else
							{
								e2_ac_new[i][j][k] = (1. - V2_ac[i][j][k]) +
									(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1_new[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
									(tmpi_ac1_new[i][j][k] - 1.) / V2_ac[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
									(tmpe_ac1_new[i][j][k] - 1.);
							}
						}
					}
				}
				else
				{
					if ((tmpi_ac1_new[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - delta))
					{
						e2_ac_new[i][j][k] = (1. - V2_ac[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6)) *
							(tmpi_ac1_new[i][j][k] - 1.) / V2_ac[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * 1e-6)) *
							(tmpe_ac1_new[i][j][k] - 1.);
					}

					if ((tmpi_ac1_new[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - delta) && (tmpi_ac1_new[i][j][k] * param.T00) <= (Melt_metal[mt].T_melting + delta))
					{
						//Ci[J/cm3 K] = ((Dependence_C_l_on_T(mt, param.T00 * tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))
							// spl_C_l_on_T.GetY(param.T00 * Ti_acoustic[i][j][k])
							// spl_Ti.GetY(k * dz) = tmpi2

						/*e2[i][j][k] = (1. - V2[i][j][k]) +
								(Melt_metal[mt].gi * param.T00 * (((((Dependence_C_l_on_T(mt, param.T00* tmpi2[i][j][k]) / (100. * 100. * 100.) + (Melt_metal[mt].Q_fusion * Melt_metal[mt].Density) / (2 * delta)))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
								(Ti_acoustic[i][j][k] - 1.) / V2[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * Te_acoustic[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
								(Te_acoustic[i][j][k] - 1.);// */

						e2_ac_new[i][j][k] = (1. - V2_ac[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * (((((spl_C_l_on_T.GetY(param.T00 * tmpi_ac1_new[i][j][k])))))) / (pow((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2, 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * V0 * 1e-6)) *
							(tmpi_ac1_new[i][j][k] - 1.) / V2_ac[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * 1e-6)) *
							(tmpe_ac1_new[i][j][k] - 1.);
					}

					if ((tmpi_ac1_new[i][j][k] * param.T00) > (Melt_metal[mt].T_melting + delta))
					{
						e2_ac_new[i][j][k] = (1. - V2_ac[i][j][k]) +
							(Melt_metal[mt].gi * param.T00 * ((Dependence_C_l_on_T(mt, param.T00 * tmpi_ac1_new[i][j][k]) / (100. * 100. * 100.))) / (Melt_metal[mt].DensityLiquid * Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * V0 * 1e-6)) *
							(tmpi_ac1_new[i][j][k] - 1.) / V2_ac[i][j][k] + (Melt_metal[mt].ge * param.T00 * (spl_C_e_on_T.GetY(param.T00 * tmpe_ac1_new[i][j][k]) / (100. * 100. * 100.)) / (Melt_metal[mt].DensityLiquid * pow(Melt_metal[mt].u0_Liquid, 2) * 1e-6)) *
							(tmpe_ac1_new[i][j][k] - 1.);
					}
				}
			}
		}
	}

	//plt_TE_ht.SetParametrsOnPlotColor(43, "e2", "z,mkm", "x,mkm", (1e+4* dz* Nz / 5e+5), (1e+4* 1e-2* dx* Nx));//p0v.p0v_s
	//plt_TE_ht.SetDataOnPlotColor3D(43, 100, Nx, Ny, Nz, 1e+4* dx* 1e-2, 1e+4* dy* 1e-2, 1e+4* dz / 5e+5, e2, p0v.p0v_s, Ny / 2, xz);
	//plt_TE_ht.ShowDataOnPlotColor(43, "file e2 in ht (x,y div 2,z) ", true);
	//Sleep(1000);

	//plt_TE_ac.SetParametrsOnPlotColor(44, "e2", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	//plt_TE_ac.SetDataOnPlotColor3D(44, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, e2_ac_new, p0v.p0v_s, Ny_ac / 2, xz);
	//plt_TE_ac.ShowDataOnPlotColor(44, "file e2 in ac part (x,y div 2,z) ", true);
	//Sleep(1000);


	/*e2 (1-й этап) */
	Transform_2D_GridFrom_AcToHeat_xy(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx, dy, e2_ac_new, e2, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

	Transform_2D_GridFrom_AcToHeat_xz(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx, dz, e2_ac_new, e2, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
		left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

	// from ht to ac /* сшивка границ давления (2-й этап) */
	TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, Nx_ac, Ny_ac, Nz_ac,
		dx_ac, dy_ac, e2, e2_ac_new, /*Te_old_y_new_x,*/ spl_Te);

	TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx, dy, dz, Nx_ac, Ny_ac, Nz_ac, Nx_ac, Ny_ac, dx_ac, dy_ac, dz_ac, e2, e2_ac_new, /*Te_old_z_new_x,*/ spl_Te);


	//// сшивка давления (внутренн границы) (старая версия)
	/*TransformGridFrom_AcToHeat(Nx_ac, Ny_ac, Nz_ac, dx_ac, dy_ac, dz_ac,
		Nx_part_ht, Ny_part_ht, down_boundary_z_heat + 1, dx, dy, dz, e2_ac_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, e2, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);*/

	//plt_TE_ht.SetParametrsOnPlotColor(45, "e2", "z,mkm", "x,mkm", (1e+4* dz* Nz / 5e+5), (1e+4* 1e-2* dx* Nx));
	//plt_TE_ht.SetDataOnPlotColor3D(45, 100, Nx, Ny, Nz, 1e+4* dx* 1e-2, 1e+4* dy* 1e-2, 1e+4* dz / 5e+5, e2, p0v.p0v_s, Ny / 2, xz);
	//plt_TE_ht.ShowDataOnPlotColor(45, "file e2 in HT", true);
	//Sleep(1000);
	//cout << " STOP STOP STOP STOP STOP e P " << endl;
	//system("pause");

	// ниже 4 действия делаем для того, чтобы  в акуст массиве дать значения смежных границы
	/*CopyDataArray3D(a1x_ac, e2_ac_new, Nx_ac, Ny_ac, Nz_ac);
	TransformGrid(Nx, Ny, Nz, dx, dy, dz, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, dx_ac, dy_ac, dz_ac, e2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
	Copy_Data_From_Big_Array3D_To_Small_Array3D(e2_ac_new, Te_acoustic, Nx_acoustic_total, Ny_acoustic_total, Nz_acoustic_total, left_boundary_x_ac, right_boundary_x_ac, left_boundary_y_ac, right_boundary_y_ac, down_boundary_z_ac);
	Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(e2_ac_new, a1x_ac, Nx_ac, Ny_ac, Nz_ac);*/

	/*plt_TE_ac.SetParametrsOnPlotColor(46, "e2", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	plt_TE_ac.SetDataOnPlotColor3D(46, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, e2_ac_new, p0v.p0v_s, Ny_ac / 2, xz);
	plt_TE_ac.ShowDataOnPlotColor(46, "file e2 in ac part with boundary(x,y div 2,z) ", true);
	Sleep(1000);
	system("pause");

	plt_TE_ac.SetParametrsOnPlotColor(47, "a1x", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	plt_TE_ac.SetDataOnPlotColor3D(47, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, a1x_ac, 1., Ny_ac / 2, xz);
	plt_TE_ac.ShowDataOnPlotColor(47, "file a1x in Ac vnutr with boundary pered vyhod", true);
	Sleep(1000);
	cout << " STOP STOP STOP STOP STOP a2x pered vyhod" << endl;
	system("pause");

	plt_TE_ac.SetParametrsOnPlotColor(48, "a2x", "z,mkm", "x,mkm", (1e+4* dz_ac* Nz_ac / 5e+5), (1e+4* 1e-2* dx_ac* Nx_ac));
	plt_TE_ac.SetDataOnPlotColor3D(48, 100, Nx_ac, Ny_ac, Nz_ac, 1e+4* dx_ac* 1e-2, 1e+4* dy_ac* 1e-2, 1e+4* dz_ac / 5e+5, a2x_ac, 1., Ny_ac / 2, xz);
	plt_TE_ac.ShowDataOnPlotColor(48, "file a2x in Ac vnutr with boundary pered vyhod", true);
	Sleep(1000);
	cout << " STOP STOP STOP STOP STOP a2x pered vyhod" << endl;
	system("pause");

	plt_TE_ac.Close_and_open_files_for_replot(for_close_and_open);
	plt_TE_ht.Close_and_open_files_for_replot(for_close_and_open);*/

	//ofstream fout_ac("e2_ac_new.txt");
	for (int i = 0; i < Nx_ac /*- 1*/; i++) // циклы по узлам акустики (определяем где произошел отрыв)
	{
		for (int j = 0; j < Ny_ac /*- 1*/; j++)
		{
			for (int k = 1; k < Nz_ac /*- 1*/; k++)
			{
				if ((e2_ac_new[i][j][k] * p0v.p0v_s) <= -1400)
				{
					for (int ii = i - 1; ii <= i + 1; ii++)
					{
						for (int jj = j - 1; jj <= j + 1; jj++)
						{
							for (int kk = k - 1; kk <= k + 1; kk++)
							{
								e2_ac_new[ii][jj][kk] = 1e-16;// индикаор того, что в точке уже нет вещества// 1e-16;
								Point3D tmp1(ii, jj, kk);
								points_rupture_ac.push_back(tmp1);
							}
						}
					}
				}
			}
		}
	}


	//for (int i = 0; i < Nx_ac /*- 1*/; i++) // циклы по узлам акустики (определяем где произошел отрыв)
	//{
	//	for (int j = 0; j < Ny_ac /*- 1*/; j++)
	//	{
	//		for (int k = 0; k < Nz_ac /*- 1*/; k++)
	//		{
	//			if ((e2_ac_new[i][j][k] * p0v.p0v_s) <= -1400)
	//			{
	//				e2_ac_new[i][j][k] = 0.;// индикаор того, что в точке уже нет вещества// 1e-16;
	//				e2_ac_new[i][j][k + 1] = 1e-16; // устанавливаем ГУ - свободная граница P = 0
	//				Point3D tmp1(i, j, k);
	//				Point3D tmp2(i, j, k + 1);
	//				points_rupture_ac.push_back(tmp1);
	//				points_rupture_ac.push_back(tmp2);
	//				for (int m = 0; m < k; m++)
	//				{
	//					e2_ac_new[i][j][m] = 0.;
	//					Point3D tmp(i, j, m);
	//					points_rupture_ac.push_back(tmp);
	//				}
	//				// старый варинат
	//				//e2_ac_new[i][j][k] = 1e-16;
	//				//Point3D tmp(i, j, k);
	//				//points_rupture_ac.push_back(tmp);
	//			}
	//		}
	//	}
	//}

	if (!points_rupture_ac.empty() /*&& current_count_frame <= 1599*/)
	{
		string current_count_frame_str = ConvertNumToString(current_count_frame);
		ofstream file_animAc("PointsAc" + current_count_frame_str + ".txt");
		//file_animAc << "Ac" << endl;
		for (it_ac = points_rupture_ac.begin(); it_ac != points_rupture_ac.end(); it_ac++)
		{
			file_animAc << "(" << (*it_ac).index_x << "," << (*it_ac).index_y << "," << (*it_ac).index_z << ")" << endl;
		}
	}

}

void MainProcedure(Metal mt, TypeBeam tbeam, double kte_kte, double ro0, double CeCe, double CiCi, double gammagamma, double g_e, double g_i, double u00, double tptp, double P00, double*** V00,
	double*** V1, double*** V2, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** a2x, double*** a1x, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z,
	double*** b2z, double*** e2, double*** e1, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, int Nx_heat, int  Ny_heat, int  Nz_heat,
	int  Nx_acoustic, int  Ny_acoustic, int  Nz_acoustic, Splayn spl_C_l_on_T, /*double*** Te_acoustic*/ /*double*** Tei_old_xy_new_z_for_total_ac, double*** Tei_old_y_new_xz_for_total_ac,*/ string current_namefile)    // davlenie, gorizontal'naja skorost', vertikal'naja skorost', temperatura
{
	// parametrs of laser and area
	// 20 mkm = 0.002      200 mkm = 0.02 cm 
	double r0 = 1e-2;// 1e-2; // radius of light beam, cm
	double kabs = 5e+5; // absorption coeff, cm-1
	double P0 = P00;//1e+8; // input intensity, W/cm2
	double tp = tptp;// 1e-13; // pulse duration, s
	cout << " tp, s = " << tp << endl;
	cout << " tp, fs = " << tp * 1e+15 << endl;
	double xy0 = 1e-1; // transverse size, cm
	double z0 = 1e-4;  // longsize, cm
	cout << " Fluence = " << tp * P0 << endl;

	double T00 = 300;    // char temperature, K
	double u0 = u00;// 3.24e+5; // velocity of sound, cm/s
	double t0 = 1 / (kabs * u0);  // char time, s
	double beta = t0 / tp;  // parameter of pulse
	double ge = g_e;
	double gi = g_i;

	cout << endl;

	double A2 = kabs * kabs * r0 * r0;
	Parametrs param;
	param.beta = beta;
	param.kabs = kabs;
	param.P0 = P0;
	param.r0 = r0;
	param.t0 = t0;
	param.T00 = T00;

	string metall, type_beam;
	vector<Melting> Melt_metal(4);

	switch (mt)
	{
	case Au:
		metall = "Au";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 63.7; // Дж/г
		Melt_metal[mt].T_melting = 1338.; // K
		Melt_metal[mt].Density = 19.3; // g/cm3
		Melt_metal[mt].DensityLiquid = 17.; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		Melt_metal[mt].u0_Liquid = 2.567e+5;
		break;
	case Al:
		metall = "Al";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 397.; // Дж/г
		Melt_metal[mt].T_melting = 934.; // K
		Melt_metal[mt].Density = 2.7; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		break;
	case Cu:
		metall = "Cu";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 209.; // Дж/г
		Melt_metal[mt].T_melting = 1358.; // K
		Melt_metal[mt].Density = 8.96; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		break;
	case Ni:
		metall = "Ni";
		Melt_metal[mt].mt = mt;
		Melt_metal[mt].Q_fusion = 298.; // Дж/г
		Melt_metal[mt].T_melting = 1726.; // K
		Melt_metal[mt].Density = 8.9; // g/cm3
		Melt_metal[mt].ge = g_e;
		Melt_metal[mt].gi = g_i;
		Melt_metal[mt].u0 = u00;
		break;
	}

	switch (tbeam)
	{
	case Gauss:
		type_beam = "Gauss";
		break;
	case Vortex:
		type_beam = "Vortex";
		break;
	}

	/////////////////////////////////////////

	double V0 = 1;
	double CC0 = 1 / (kabs * r0);

	//double p0v = ro0 * 1000 * pow(u0 * 0.01, 2) * 1e-5;//202.8e+3;   // c00*c00*ro0 - normirovka davlenija bar (r0 *u0* u0) (perevod v SI i v bar)
	double p0v_s = Melt_metal[mt].Density * Melt_metal[mt].u0 * Melt_metal[mt].u0 * (1e-6); // 
	cout << " p0v = " << p0v_s << endl;
	double  p0v_sl = ((Melt_metal[mt].Density + Melt_metal[mt].DensityLiquid) / 2) * pow((Melt_metal[mt].u0 + Melt_metal[mt].u0_Liquid) / 2, 2) * (1e-6);
	cout << " p0v 1 = " << p0v_sl << endl;
	double  p0v_l = Melt_metal[mt].DensityLiquid * Melt_metal[mt].u0_Liquid * Melt_metal[mt].u0_Liquid * (1e-6);
	cout << " p0v = " << p0v_l << endl;

	P0V p0v;
	p0v.p0v_s = p0v_s;
	p0v.p0v_sl = p0v_sl;
	p0v.p0v_l = p0v_l;

	double sigma = 1400; // bar Au - предел прочности

	cout << Nx_heat << "   " << Ny_heat << "   " << Nz_heat << endl;
	cout << Nx_acoustic << "   " << Ny_acoustic << "   " << Nz_acoustic << endl;

	double dxy_heat = xy0 / (r0 * Nx_heat);
	double dx_heat = dxy_heat;
	double dy_heat = dxy_heat;
	double dz_heat = z0 * kabs / Nz_heat;
	double dxy_acoustic = xy0 / (r0 * Nx_acoustic);
	double dx_acoustic = dxy_acoustic;
	double dy_acoustic = dxy_acoustic;
	double dz_acoustic = z0 * kabs / 10000;//10000;//Nz_acoustic; // здесь даем 2000 значит Nz_ac = 400 узлов

	cout << dx_heat << "  " << dy_heat << "  " << dz_heat << endl;
	cout << dx_acoustic << "  " << dy_acoustic << "  " << dz_acoustic << endl;
	cout << " Razmern shag x y (mkm) heat = " << 1e+4 * r0 * dx_heat << endl;
	cout << " Razmern shag z (mkm) heat = " << 1e+4 * dz_heat / kabs << endl;
	cout << " Razmern shag x y (mkm) acoustic = " << 1e+4 * r0 * dx_acoustic << endl;
	cout << " Razmern shag z (mkm) acoustic = " << 1e+4 * dz_acoustic / kabs;
	cout << endl << endl;
	cout << " Razmern shag x y (nm) heat = " << 1e+4 * r0 * dx_heat * 1000 << endl;
	cout << " Razmern shag z (nm) heat = " << 1e+4 * dz_heat / kabs * 1000 << endl;
	cout << " Razmern shag x y (nm) acoustic = " << 1e+4 * r0 * dx_acoustic * 1000 << endl;
	cout << " Razmern shag z (nm) acoustic = " << 1e+4 * dz_acoustic / kabs * 1000;
	cout << endl << endl;

	int left_boundary_x_heat = 300. / (1e+4 * r0 * dx_heat); // 300 мкм
	int right_boundary_x_heat = 700. / (1e+4 * r0 * dx_heat);
	int left_boundary_y_heat = 300. / (1e+4 * r0 * dx_heat); // 300 мкм
	int right_boundary_y_heat = 700. / (1e+4 * r0 * dx_heat);
	int down_boundary_z_heat = 100. / (1e+4 * dz_heat / kabs * 1000); // 100 нм

	int left_boundary_x_acoustic = 300. / (1e+4 * r0 * dx_acoustic);
	int right_boundary_x_acoustic = 700. / (1e+4 * r0 * dx_acoustic);
	int left_boundary_y_acoustic = 300. / (1e+4 * r0 * dx_acoustic);
	int right_boundary_y_acoustic = 700. / (1e+4 * r0 * dx_acoustic);
	int down_boundary_z_acoustic = 100. / (1e+4 * dz_acoustic / kabs * 1000);

	cout << "left_boundary_x_heat  = " << left_boundary_x_heat << endl;
	cout << "right_boundary_x_heat  = " << right_boundary_x_heat << endl;
	cout << "left_boundary_y_heat  = " << left_boundary_y_heat << endl;
	cout << "right_boundary_y_heat  = " << right_boundary_y_heat << endl;
	cout << " down_boundary_z_heat  = " << down_boundary_z_heat << endl;
	cout << endl << endl;

	cout << "left_boundary_x_acoustic  = " << left_boundary_x_acoustic << endl;
	cout << "right_boundary_x_acoustic  = " << right_boundary_x_acoustic << endl;
	cout << "left_boundary_y_acoustic  = " << left_boundary_y_acoustic << endl;
	cout << "right_boundary_y_acoustic  = " << right_boundary_y_acoustic << endl;
	cout << " down_boundary_z_acoustic = " << down_boundary_z_acoustic << endl;

	int Nx_part_ht = right_boundary_x_heat - left_boundary_x_heat;
	int Ny_part_ht = right_boundary_y_heat - left_boundary_y_heat;
	cout << Nx_part_ht << "    " << Ny_part_ht << "   " << down_boundary_z_heat << endl;
	int Nx_part_ac = right_boundary_x_acoustic - left_boundary_x_acoustic;
	int Ny_part_ac = right_boundary_y_acoustic - left_boundary_y_acoustic;
	cout << Nx_part_ac << "    " << Ny_part_ac << "   " << down_boundary_z_acoustic << endl;

	int coeff_big_z = down_boundary_z_acoustic / down_boundary_z_heat;
	int coeff_big_y = Ny_part_ac / Ny_part_ht;
	int coeff_big_x = coeff_big_y;
	cout << " coeff_big_z = " << coeff_big_z << " coeff_big_y = " << coeff_big_y << endl;

	/////////////////////////////

	int Nx_acoustic_show = Nx_part_ac + 1;
	int Ny_acoustic_show = Ny_part_ac + 1;
	int Nz_acoustic_show = down_boundary_z_acoustic + 1;

	/////////////////////////////////

	double dt = 1. / (t0 * 1e+15);// 2e-4;
	cout << " dt = " << dt << endl;
	int n = 1;
	double tt = 0;
	// 32 переменные double

	Splayn spl_G_e_on_T;
	Splayn spl_C_e_on_T;
	spl_C_e_on_T.InterpolateFast1D("Ce_" + metall + "_new.txt");
	spl_G_e_on_T.InterpolateFast1D("Ge_" + metall + "_new.txt");
	Splayn spl_Te, spl_Ti, spl_e;

	double*** V00_ac_new, ***V1_ac_new, ***V2_ac_new, ***a1y_ac_new, ***a2y_ac_new, ***a1z_ac_new, ***a2z_ac_new, ***a2x_ac_new, ***a1x_ac_new, ***b1x_ac_new,
		***b2x_ac_new, ***b1y_ac_new, ***b2y_ac_new, ***b1z_ac_new, ***b2z_ac_new, ***e2_ac_new, ***e1_ac_new, ***Tei_old_xy_new_z, ***Tei_old_y_new_xz;
	double*** tmpe_ac0_new, ***tmpe_ac1_new, ***tmpe_ac2_new, ***tmpi_ac0_new, ***tmpi_ac1_new, ***tmpi_ac2_new;
	double*** tmpe_for_transform_new, **Array3D_for_Plot; // массив куда перебрасываем часть поля (например темепраутры) и эти данные пеерводим из крупной в мелкую сетку
	//double** Te_old_x_new_y, ** Te_old_x_new_z;

	//Te_old_x_new_y = new double* [Nx_part_ac + 1]; 
	//for (int j = 0; j < Nx_part_ac + 1; j++)
	//{
	//	Te_old_x_new_y[j] = new double[Ny_part_ht + 1];
	//}

	//Te_old_x_new_z = new double* [Nx_part_ac + 1];
	//for (int j = 0; j < Nx_part_ac + 1; j++)
	//{
	//	Te_old_x_new_z[j] = new double[down_boundary_z_heat + 1];
	//}

	tmpe_for_transform_new = new double**[Nx_part_ht + 1]; //Nx_part_ht + 2

	for (int i = 0; i < Nx_part_ht + 1; i++)
	{
		tmpe_for_transform_new[i] = new double*[Ny_part_ht + 1];
		for (int j = 0; j < Ny_part_ht + 1; j++)
		{
			tmpe_for_transform_new[i][j] = new double[down_boundary_z_heat + 1];
		}
	}

	Tei_old_xy_new_z = new double**[Nx_part_ht + 1];//1
	Tei_old_y_new_xz = new double**[Nx_part_ac + 1];

	for (int i = 0; i < Nx_part_ht + 1; i++)
	{
		Tei_old_xy_new_z[i] = new double*[Ny_part_ht + 1];
		for (int j = 0; j < Ny_part_ht + 1; j++)
		{
			Tei_old_xy_new_z[i][j] = new double[down_boundary_z_acoustic + 1];
		}
	}

	for (int i = 0; i < Nx_part_ac + 1; i++)
	{
		Tei_old_y_new_xz[i] = new double*[Ny_part_ht + 1];
		for (int j = 0; j < Ny_part_ht + 1; j++)
		{
			Tei_old_y_new_xz[i][j] = new double[down_boundary_z_acoustic + 1];
		}
	}

	tmpe_ac0_new = new double**[Nx_part_ac + 1];
	tmpe_ac1_new = new double**[Nx_part_ac + 1];
	tmpe_ac2_new = new double**[Nx_part_ac + 1];
	tmpi_ac0_new = new double**[Nx_part_ac + 1];
	tmpi_ac1_new = new double**[Nx_part_ac + 1];
	tmpi_ac2_new = new double**[Nx_part_ac + 1];

	V00_ac_new = new double**[Nx_part_ac + 1];
	V1_ac_new = new double**[Nx_part_ac + 1];
	V2_ac_new = new double**[Nx_part_ac + 1];
	a1x_ac_new = new double**[Nx_part_ac + 1];// Euler coordinates
	a2x_ac_new = new double**[Nx_part_ac + 1];// Euler coordinates
	a1y_ac_new = new double**[Nx_part_ac + 1];// Euler coordinates
	a2y_ac_new = new double**[Nx_part_ac + 1];// Euler coordinates
	a1z_ac_new = new double**[Nx_part_ac + 1];// Euler coordinates
	a2z_ac_new = new double**[Nx_part_ac + 1];// Euler coordinates
	b1x_ac_new = new double**[Nx_part_ac + 1]; // particle velocity
	b2x_ac_new = new double**[Nx_part_ac + 1]; // particle velocity
	b1y_ac_new = new double**[Nx_part_ac + 1]; // particle velocity
	b2y_ac_new = new double**[Nx_part_ac + 1]; // particle velocity
	b1z_ac_new = new double**[Nx_part_ac + 1]; // particle velocity
	b2z_ac_new = new double**[Nx_part_ac + 1]; // particle velocity
	e1_ac_new = new double**[Nx_part_ac + 1]; // Pressure
	e2_ac_new = new double**[Nx_part_ac + 1]; // Pressure

	Array3D_for_Plot = new double*[down_boundary_z_acoustic + 1];

	for (int i = 0; i < down_boundary_z_acoustic + 1; i++)
	{
		Array3D_for_Plot[i] = new double[Nx_part_ac + 1];
	}

	for (int i = 0; i < Nx_part_ac + 1; i++)
	{
		tmpe_ac0_new[i] = new double*[Ny_part_ac + 1];
		tmpe_ac1_new[i] = new double*[Ny_part_ac + 1];
		tmpe_ac2_new[i] = new double*[Ny_part_ac + 1];
		tmpi_ac0_new[i] = new double*[Ny_part_ac + 1];
		tmpi_ac1_new[i] = new double*[Ny_part_ac + 1];
		tmpi_ac2_new[i] = new double*[Ny_part_ac + 1];

		V00_ac_new[i] = new double*[Ny_part_ac + 1];
		V1_ac_new[i] = new double*[Ny_part_ac + 1];
		V2_ac_new[i] = new double*[Ny_part_ac + 1];
		a1x_ac_new[i] = new double*[Ny_part_ac + 1];// Euler coordinates
		a2x_ac_new[i] = new double*[Ny_part_ac + 1];// Euler coordinates
		a1y_ac_new[i] = new double*[Ny_part_ac + 1];// Euler coordinates
		a2y_ac_new[i] = new double*[Ny_part_ac + 1];// Euler coordinates
		a1z_ac_new[i] = new double*[Ny_part_ac + 1];// Euler coordinates
		a2z_ac_new[i] = new double*[Ny_part_ac + 1];// Euler coordinates
		b1x_ac_new[i] = new double*[Ny_part_ac + 1]; // particle velocity
		b2x_ac_new[i] = new double*[Ny_part_ac + 1]; // particle velocity
		b1y_ac_new[i] = new double*[Ny_part_ac + 1]; // particle velocity
		b2y_ac_new[i] = new double*[Ny_part_ac + 1]; // particle velocity
		b1z_ac_new[i] = new double*[Ny_part_ac + 1]; // particle velocity
		b2z_ac_new[i] = new double*[Ny_part_ac + 1]; // particle velocity
		e1_ac_new[i] = new double*[Ny_part_ac + 1]; // Pressure
		e2_ac_new[i] = new double*[Ny_part_ac + 1]; // Pressure

		//Array3D_for_Plot[i] = new double*[Nx_part_ac + 1];

		for (int j = 0; j < Ny_part_ac + 1; j++)
		{
			tmpe_ac0_new[i][j] = new double[down_boundary_z_acoustic + 1];
			tmpe_ac1_new[i][j] = new double[down_boundary_z_acoustic + 1];
			tmpe_ac2_new[i][j] = new double[down_boundary_z_acoustic + 1];
			tmpi_ac0_new[i][j] = new double[down_boundary_z_acoustic + 1];
			tmpi_ac1_new[i][j] = new double[down_boundary_z_acoustic + 1];
			tmpi_ac2_new[i][j] = new double[down_boundary_z_acoustic + 1];

			V00_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];
			V1_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];
			V2_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];
			a1x_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];// Euler coordinates
			a2x_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];// Euler coordinates
			a1y_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];// Euler coordinates
			a2y_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];// Euler coordinates
			a1z_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];// Euler coordinates
			a2z_ac_new[i][j] = new double[down_boundary_z_acoustic + 1];// Euler coordinates
			b1x_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // particle velocity
			b2x_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // particle velocity
			b1y_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // particle velocity
			b2y_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // particle velocity
			b1z_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // particle velocity
			b2z_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // particle velocity
			e1_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // Pressure
			e2_ac_new[i][j] = new double[down_boundary_z_acoustic + 1]; // Pressure

			//Array3D_for_Plot[i][j] = new double[down_boundary_z_acoustic + 1];
		}
	}

	// initial conditions

	for (int i = 0; i < Nx_heat; i++)
	{
		for (int j = 0; j < Ny_heat; j++)
		{
			for (int k = 0; k < Nz_heat; k++)
			{
				tmpe0[i][j][k] = 1.0;// 1e-16;
				tmpe1[i][j][k] = 1e-16;
				tmpe2[i][j][k] = 1e-16;
				tmpi0[i][j][k] = 1.0;// 1e-16;
				tmpi1[i][j][k] = 1.0;// 1e-16;
				tmpi2[i][j][k] = 1e-16;
			}
		}
	}

	for (int i = 0; i < Nx_heat; i++)
	{
		for (int j = 0; j < Ny_heat; j++)
		{
			for (int k = 0; k < Nz_heat; k++)
			{
				e1[i][j][k] = 1e-16;
				e2[i][j][k] = 1e-16;
				a1x[i][j][k] = (i - 1) * dx_heat;
				a2x[i][j][k] = 1e-16;
				a1y[i][j][k] = (j - 1) * dy_heat;
				a2y[i][j][k] = 1e-16;
				a1z[i][j][k] = (k - 1) * dz_heat;
				a2z[i][j][k] = 1e-16;
				b1x[i][j][k] = 1e-16;
				b2x[i][j][k] = 1e-16;
				b1y[i][j][k] = 1e-16;
				b2y[i][j][k] = 1e-16;
				b1z[i][j][k] = 1e-16;
				b2z[i][j][k] = 1e-16;
				V2[i][j][k] = V0;
				V1[i][j][k] = V0;
				V00[i][j][k] = V0;
				///massive_melting_tmp1[k][i] = 1e-16;
				/*massive_melting_tmp2[i][j][k] = 1e-16;*/
			}
		}
	}

	SetDataArray3D(tmpe_ac0_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1.0);
	Null_Array3D(tmpe_ac1_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(tmpe_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);/////
	SetDataArray3D(tmpi_ac0_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1.0);
	SetDataArray3D(tmpi_ac1_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1.0);
	Null_Array3D(tmpi_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

	Null_Array3D(e1_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(e2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	//Null_Array3D(a1x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(a2x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	//Null_Array3D(a1y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(a2y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	//Null_Array3D(a1z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(a2z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(b1x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(b2x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(b1y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(b2y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(b1z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
	Null_Array3D(b2z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

	SetDataArray3D(V00_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, V0);
	SetDataArray3D(V1_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, V0);
	SetDataArray3D(V2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, V0);

	//Null_Array3D(Te_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic);
	Null_Array(Array3D_for_Plot, down_boundary_z_acoustic + 1, Nx_part_ac + 1);
	//Null_Array3D(Ti_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic);

	for (int i = 0; i < Nx_part_ac + 1; i++)
	{
		for (int j = 0; j < Ny_part_ac + 1; j++)
		{
			for (int k = 0; k < down_boundary_z_acoustic + 1; k++)
			{
				a1x_ac_new[i][j][k] = (i - 1) * dx_acoustic;
				a1y_ac_new[i][j][k] = (j - 1) * dy_acoustic;
				a1z_ac_new[i][j][k] = (k - 1) * dz_acoustic;
			}
		}
	}

	cout << " hello" << endl;
	//system("pause");

	//для распаралелл тепловой задачи
	vector<Interval> new_interv = {
		Interval{ 1, Nz_heat / 10 },
		Interval{ Nz_heat / 10, 2 * Nz_heat / 10 },
		Interval{ 2 * Nz_heat / 10, 3 * Nz_heat / 10 },
		Interval{ 3 * Nz_heat / 10, 4 * Nz_heat / 10 },
		Interval{ 4 * Nz_heat / 10, 5 * Nz_heat / 10 },
		Interval{ 5 * Nz_heat / 10, 6 * Nz_heat / 10 },
		Interval{ 6 * Nz_heat / 10, 7 * Nz_heat / 10 },
		Interval{ 7 * Nz_heat / 10, 8 * Nz_heat / 10 },
		Interval{ 8 * Nz_heat / 10, 9 * Nz_heat / 10 },
		Interval{ 9 * Nz_heat / 10,  Nz_heat - 1 }
	};

	vector<Interval> new_interv_acoustic_z = {
		Interval{ 1, Nz_acoustic / 10 },
		Interval{ Nz_acoustic / 10, 2 * Nz_acoustic / 10 },
		Interval{ 2 * Nz_acoustic / 10, 3 * Nz_acoustic / 10 },
		Interval{ 3 * Nz_acoustic / 10, 4 * Nz_acoustic / 10 },
		Interval{ 4 * Nz_acoustic / 10, 5 * Nz_acoustic / 10 },
		Interval{ 5 * Nz_acoustic / 10, 6 * Nz_acoustic / 10 },
		Interval{ 6 * Nz_acoustic / 10, 7 * Nz_acoustic / 10 },
		Interval{ 7 * Nz_acoustic / 10, 8 * Nz_acoustic / 10 },
		Interval{ 8 * Nz_acoustic / 10, 9 * Nz_acoustic / 10 },
		Interval{ 9 * Nz_acoustic / 10,  Nz_acoustic - 1 }
	};

	vector<Point3D> points_rupture_ac, points_rupture_heat;

	Calculation00(mt, param, spl_C_e_on_T, spl_G_e_on_T, tbeam, tmpe0, tmpe1, tmpi0, tmpi1, Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, dt, A2, n, beta);
	///Calculation1(a1x, a2x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, new_interv_acoustic_z);// вязкости еще не будет т.к. поле скоростей берутся из н/у 
	///GRU1(e2, Nx_acoustic, Ny_acoustic);
	///GRU2(b2x, b2y, b2z, Nx_acoustic, Ny_acoustic, Nz_acoustic);
	///Calculation2(a1x, a1y, a1z, a2x, a2y, a2z, b2x, b2y, b2z, V2, e2, tmpe1, tmpi1, dx_acoustic, dy_acoustic, dz_acoustic, dx_heat, dy_heat, dz_heat, dt, CC0, Nx_acoustic, Ny_acoustic, Nz_acoustic, Nx_heat, Ny_heat, Nz_heat, V0, p0v, Melt_metal, mt, param, spl_C_e_on_T, spl_Te, spl_Ti, spl_e, /*e_heat, */Te_acoustic, Ti_acoustic, Tei_old_xy_new_z, Tei_old_y_new_xz, current_count_frame_begin, points_rupture_ac, points_rupture_heat, new_interv_acoustic_z, spl_C_l_on_T);

	tt = n * dt * t0 * 1e+15; // vremja v fs

	int count_of_lines = 9;
	int npxzt = 1;
	int current_number_line_melting = 1;
	int current_time_frame = 5;//8200;// 6650;// 600; //650 // 900; // номер кадра 
	int total_time_frame = 2500000;// 2000; // формальный апарметр теперь (уже)
	vector<double> moments_time; // хранит моменты времени
	vector<int> moments_fix_time_anim;
	for (int i = 1; i <= 28000/*2000*/; i++) // сетка времени
	{
		tt = i * dt * t0 * 1e+15; // vremja v fs
		moments_time.push_back(tt);
	}

	int dt_step = 5;
	for (int i = current_time_frame; i <= 27900; i += dt_step/*i++*/ /*= 10*/) // моменты времени для отрисовки и анимации
	{// аналог n
		moments_fix_time_anim.push_back(i);
	}

	//MyCreateFile(list_namefile.size(), list_namefile);
	vector<double>::iterator it_time, it_last_time; // хранит моменты времени
	vector<int>::iterator it_fix_time_anim;
	it_last_time = moments_time.end() - 1;
	it_fix_time_anim = moments_fix_time_anim.begin();
	cout << endl << endl << " last time = " << (*it_last_time) << endl;
	cout << endl << endl << " first fix time = " << (*(moments_fix_time_anim.begin())) << endl;
	cout << endl << endl << " last fix time = " << (*(moments_fix_time_anim.end() - 1)) << endl;

	cout << endl << endl;
	double part_depth = 1.;
	int number_plots = 29;
	GnuPlot plt(number_plots);
	plt.SetParametrs2D(0, 10, 3, "Ti", "x,mkm", "Ti, (K)");
	plt.SetParametrs2D(1, 2, 3, "Te,Ti", "time,fs", "Te,Ti (K)");
	plt.SetParametrsOnPlotColor(2, "Te", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(3, "Ti", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(4, "P(x,N/2,z)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(5, "P(x,y,1)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(6, "Te", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(7, "Ti", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrs2D(8, count_of_lines, 3, "P,bar", "x,mkm", "P,bar");
	plt.SetParametrs2D(9, count_of_lines, 3, "P, bar", "z,mkm", "P, bar");
	plt.SetParametrs2D(10, 1, 3, "Temperature in central zone", "time,fs", "Ti, K"); // копии 1 и 10 графиков
	plt.SetParametrs2D(11, 11, 3, "Temperature T(z)", "z,mkm", "Ti, K");
	plt.SetParametrs2D(12, 11, 3, "Temperature T(x)", "x,mkm", "Ti, K");
	plt.SetParametrsOnPlotColor(13, "Melting zone (xz)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(14, "Melting zone (xy)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrs2D(15, 2, 3, "Te,Ti", "x,mkm", "Te,Ti (K)");
	plt.SetParametrs2D(16, 2, 3, "Te,Ti", "z,mkm", "Te,Ti (K)");

	/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	plt.SetParametrs2D(17, 10, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(18, 1, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(19, 1, 3, "P", "z,mkm", "P (bar)");

	plt.SetParametrsOnPlotColor(20, "Preasure zone (xz)", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(21, "Preasure zone (xy)", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));

	plt.SetParametrs2D(22, 10, 3, "P", "x,mkm", "P (bar)");
	plt.SetParametrs2D(23, 10, 3, "Ti", "x,mkm", "Ti, (K)");

	plt.SetParametrsOnPlotColor(24, "a2z", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(25, "a2x", "z,mkm", "x,mkm", (1e+4 * dz_acoustic * (Nz_acoustic_show * part_depth) / kabs), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(26, "a2z", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrsOnPlotColor(27, "a2x", "y,mkm", "x,mkm", (1e+4 * r0 * dy_acoustic * Ny_acoustic_show), (1e+4 * r0 * dx_acoustic * Nx_acoustic_show));
	plt.SetParametrs2D(28, 2, 3, "a2z", "z,mkm", "a2z, mkm");

	///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// возможо здесь нужно отобразить в начальный момент времени, но это не обзательно

	plt.SetGridOnPlot3D(8, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, count_of_lines, fun_on_x);
	plt.SetGridOnPlot3D(9, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, count_of_lines, fun_on_z);
	plt.SetGridOnPlot3D(11, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 11, fun_on_z);
	plt.SetGridOnPlot3D(12, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 11, fun_on_x);
	plt.SetGridOnPlot3D(0, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 100, fun_on_x);
	plt.SetGridOnPlot3D(17, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 100, fun_on_x);
	plt.SetGridOnPlot3D(22, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 100, fun_on_x);
	plt.SetGridOnPlot3D(23, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 100, fun_on_x);

	clock_t t;
	double sum_t = 0;

	//vector<double***> TeTi = { tmpe2 , tmpi2 };//????
	//vector<double***> Ti = { tmpi2 };//????
	//vector<double***> TeTi = { tmpe_ac2_new , tmpi_ac2_new };//????
	//vector<double***> Ti = { tmpi_ac2_new };//????
	//vector<double***> P = { e1 };/////
	vector<double***> empty;
	vector<string> TeTiTfs = { "Te","Ti" };
	vector<string> Pfs = { "P" };
	vector<string> TiTfs = { "Ti" };
	vector<string> a2z_e_l = { "Lagrange", "Euler" };
	vector<string> Legenda_melting_tmp;
	vector<string> Legenda;
	vector<string> Legenda_Ti_x;
	vector<string> Legenda_P_x;
	vector<string> FilesForWritting;  // имена файлов куда запишутся поля темепратур и т.д
	vector<string> FilesForReading;
	//В векторе не указываем номера плотов, где отображаем 2 и более линий или строим зависимость величины от времени 
	vector<int> for_close_and_open = { 0, 2, 3, 4, 5, 6, 7, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28};
	vector<int> for_close_and_open_29 = { 29 };
	vector<int> for_close_and_open_ac = { 0, 1, 2 };
	string namefile;
	string fileoflistnamefiles = "List.txt";
	string depth = "Depth.txt";
	vector<Point3D> points_rupture_for_plot_z;
	vector<Point3D> points_rupture_for_plot_x;
	vector<Point3D>::iterator it;
	vector<Point3D>::iterator it_z, it_x;
	//points_rupture_ac.clear();
	//points_rupture_ac.erase(points_rupture_ac.begin(), points_rupture_ac.end());

	vector<string> list_namefile = { "P1.txt", "P2.txt", "P3.txt", "P4.txt", "P5.txt" };//?????

	string depth2 = current_namefile;// current_namefile;//"Depth2.txt"; // I0 = 0.25e+10
	ofstream fout_depth2(depth2);

	vector<string> Pacteti = { "Pac.txt", "PTi.txt", "PTe.txt", "Pacx.txt", "PTix.txt", "PTex.txt" };
	MyCreateFile(Pacteti.size(), Pacteti);
	GnuPlot pac(6, Pacteti);
	pac.SetParametrs2D(0, 1, 3, "Pacoustic", "z,mkm", "P (bar)");
	pac.SetParametrs2D(1, 1, 3, "PTi", "z,mkm", "P (bar)");
	pac.SetParametrs2D(2, 1, 3, "PTe", "z,mkm", "P (bar)");
	pac.SetParametrs2D(3, 1, 3, "Pacoustic", "x,mkm", "P (bar)");
	pac.SetParametrs2D(4, 1, 3, "PTi", "x,mkm", "P (bar)");
	pac.SetParametrs2D(5, 1, 3, "PTe", "x,mkm", "P (bar)");
	ofstream fout_Pac(Pacteti[0]);
	ofstream fout_PTi(Pacteti[1]);
	ofstream fout_PTe(Pacteti[2]);
	ofstream fout_Pacx(Pacteti[3]);
	ofstream fout_PTix(Pacteti[4]);
	ofstream fout_PTex(Pacteti[5]);
	Pacteti.clear();
	Pacteti = { "Points.txt" };
	MyCreateFile(Pacteti.size(), Pacteti);
	ofstream fout_point("Points.txt");

	int count_depth = 1;
	cout << endl << endl;
	//GnuPlot plt_TE_ht(50);
	//GnuPlot plt_TE_ac(50);

	// подготвка массивов для расчета на разных сетках (возможно сюда же подготовку температуры перетянуть)
	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, tmpe0, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, tmpe_ac0_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, tmpe1, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, tmpe_ac1_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, tmpi0, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, tmpi_ac0_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, tmpi1, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, tmpi_ac1_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, a1x, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, a1x_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, a1y, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, a1y_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, a1z, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, a1z_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, b1x, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, b1x_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, b1y, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, b1y_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, b1z, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, b1z_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, e1, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, e1_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, V00, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, V00_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, V1, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
	TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, V1_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

	it_time = moments_time.begin();
	cout << endl << "start" << endl;
	system("pause");

	do // цикл по времени
	{
		if (points_rupture_ac.empty()) {}
		else
		{
			My_unique(points_rupture_ac);
			sort(points_rupture_ac.begin(), points_rupture_ac.end(), comp_x);
			MySort_Point3D_y(points_rupture_ac);
			MySort_Point3D_z(points_rupture_ac);

			for (it = points_rupture_ac.begin(); it != points_rupture_ac.end(); it++)
			{
				tmpe_ac1_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				tmpi_ac1_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				b1x_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				b1y_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				b1z_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				a1x_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				a1y_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				a1z_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				//e1_ac_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
			}
		}

		t = clock();

		//GnuPlot plt_TE_ht(23);
		//GnuPlot plt_TE_ac(23);
		//plt_TE_ht.SetParametrsOnPlotColor(0, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		//plt_TE_ht.SetParametrsOnPlotColor(1, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		//plt_TE_ht.SetParametrsOnPlotColor(2, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		//plt_TE_ht.SetParametrsOnPlotColor(3, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* down_boundary_z_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_part_ht));
		//plt_TE_ht.SetDataOnPlotColor3D(0, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe0, 300., Ny_heat / 2, xz);
		//plt_TE_ht.SetDataOnPlotColor3D(1, 101, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe1, 300., Ny_heat / 2, xz);
		//plt_TE_ht.SetDataOnPlotColor3D(2, 101, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., Ny_heat / 2, xz);
		//plt_TE_ht.ShowDataOnPlotColor(0, "file tmpe0 ht (x,y div 2,z) ", true);
		//Sleep(1000);
		//plt_TE_ht.ShowDataOnPlotColor(1, "file tmpe1 ht (x,y div 2,z) ", true);
		//Sleep(1000);
		//plt_TE_ht.ShowDataOnPlotColor(2, "file tmpe2 ht at start (x,y div 2,z) ", true);
		//Sleep(1000);
		//system("pause");

		//plt_TE_ht.SetDataOnPlotColor3D(3, 100, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe_for_transform_new, 300., Ny_part_ht / 2, xz);
		//plt_TE_ht.ShowDataOnPlotColor(3, "file PART ht tmpe1 (x,y div 2,z) ", true);
		//Sleep(1000);

		//system("pause");

		//plt_TE_ac.SetParametrsOnPlotColor(4, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		//plt_TE_ac.SetDataOnPlotColor3D(4, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, tmpe_ac1_new, 300., Ny_part_ac / 2, xz);
		//plt_TE_ac.ShowDataOnPlotColor(4, "file PART ac tmpe1 (x,y div 2,z) ", true);
		//Sleep(1000);
		//cout << " STOP STOP VERNO " << endl;
		//system("pause");

		// отсюда получаем tmpe2 и  tmpi2  в тепл сетке без зоны отрыва и tmpe2_ac и tmpi2_ac - только зона отрыва																																																							// NX_acoustic, Ny_acoustic, Nz_acoustic
		Calculation0_for_mixed_grids(mt, param, spl_C_e_on_T, spl_G_e_on_T, tbeam, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, tmpe_ac0_new, tmpe_ac1_new, tmpe_ac2_new, tmpi_ac0_new, tmpi_ac1_new, tmpi_ac2_new, tmpe_for_transform_new, A2, dx_heat, dy_heat, dz_heat, dx_acoustic, dy_acoustic, dz_acoustic, dt, Nx_heat, Ny_heat, Nz_heat, Nx_part_ac + 1, Nx_part_ac + 1, down_boundary_z_acoustic + 1, n, beta, points_rupture_ac, new_interv, Melt_metal, spl_C_l_on_T,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic + 1, Nx_part_ac + 1, Ny_part_ac + 1, Nx_part_ht + 1, Ny_part_ht + 1, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te, a2x_ac_new);

		/*plt_TE_ht.SetParametrsOnPlotColor(5, "Te", "x,mkm", "x,mkm", (1e+4* 1e-2* dx_heat* Nx_heat), (1e+4* 1e-2* dy_heat* Ny_heat));
		plt_TE_ht.SetDataOnPlotColor3D(5, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300.,down_boundary_z_heat, xy);
		plt_TE_ht.ShowDataOnPlotColor(5, "file tmpe2 ht (x,y) ", true);*/

		//plt_TE_ht.SetParametrsOnPlotColor(5, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		//plt_TE_ht.SetDataOnPlotColor3D(5, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., Ny_heat / 2, xz);
		//plt_TE_ht.ShowDataOnPlotColor(5, "file tmpe2 ht (x,y div 2,z) ", true); // кусок
		//Sleep(1000);
		//plt_TE_ac.SetParametrsOnPlotColor(6, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		//plt_TE_ac.SetDataOnPlotColor3D(6, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, tmpe_ac2_new, 300., Ny_part_ac / 2, xz);
		//plt_TE_ac.ShowDataOnPlotColor(6, "file tmpe2 ac (x,y div 2,z) ", true); // кусок
		//Sleep(1000);

		//plt_TE_ht.SetParametrsOnPlotColor(5, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		//plt_TE_ht.SetDataOnPlotColor3D(5, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., left_boundary_y_heat + 1/*right_boundary_y_heat - 1*/, xz);
		//plt_TE_ht.ShowDataOnPlotColor(5, "file tmpe2 ht (x,y div 2,z) ", true); // кусок
		//Sleep(1000);
		//plt_TE_ac.SetParametrsOnPlotColor(6, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		//plt_TE_ac.SetDataOnPlotColor3D(6, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, tmpe_ac2_new, 300., 1/*Ny_part_ac / 2*/, xz);
		//plt_TE_ac.ShowDataOnPlotColor(6, "file tmpe2 ac (x,y div 2,z) ", true); // кусок
		//Sleep(1000);

		//Null_Array3D(tmpe_for_transform_new, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1);
		//cout << " STOP STOP STOP " << endl;
		//system("pause");


		/* сшивка границ tmpe (1-й этап) (2-я версия)*/

		Transform_2D_GridFrom_AcToHeat_xy(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, tmpe_ac2_new, tmpe2, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

		Transform_2D_GridFrom_AcToHeat_xz(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx_heat, dz_heat, tmpe_ac2_new, tmpe2,/* Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		fun21(Nx_heat, Ny_heat, Nz_heat, tmpe2, tmpi1, tmpi2, dz_heat, param, mt);

		/*	plt_TE_ht.SetParametrsOnPlotColor(3, "Te", "x,mkm", "y,mkm", (1e+4* 1e-2* dx_heat* Nx_heat), (1e+4* 1e-2* dy_heat* Ny_heat));
			plt_TE_ht.SetDataOnPlotColor3D(3, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., 3, xy);
			plt_TE_ht.ShowDataOnPlotColor(3, "file tmpe2 eith bound ht (x,y div 2,z) ", true);*/


			//plt_TE_ht.SetParametrsOnPlotColor(3, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
			//plt_TE_ht.SetDataOnPlotColor3D(3, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., Ny_heat / 2, xz);
			//plt_TE_ht.ShowDataOnPlotColor(3, "file tmpe2 with bound ht (x,y div 2,z) ", true);
			//Sleep(1000);

			//plt_TE_ht.SetParametrsOnPlotColor(3, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* down_boundary_z_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_part_ht));
			//plt_TE_ht.SetDataOnPlotColor3D(3, 100, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe_for_transform_new, 300., Ny_part_ht / 2, xz);
			//plt_TE_ht.ShowDataOnPlotColor(3, "file tmpe2 ht (x,y div 2,z) ", true);

			/*plt_TE_ht.SetParametrsOnPlotColor(3, "Te", "x,mkm", "y,mkm", (1e+4* 1e-2* dx_heat* Nx_part_ht), (1e+4* 1e-2* dy_heat* Ny_part_ht));
			plt_TE_ht.SetDataOnPlotColor3D(3, 100, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe_for_transform_new, 300., down_boundary_z_heat - 1, xy);
			plt_TE_ht.ShowDataOnPlotColor(3, "file tmpe2 ht (x,y,z_last) ", true);*/
			//cout << "Stop proverka " << endl;
			//system("pause");

			// работа с данными в tmpe2 tmpi2
			// 1-е 2 действия делаем для того, чтобы сшить границы для следующего момета времени
			// 2-е 2 действия ниже - для расчета давления в акуст сетке

			/* сшивка границ tmpe (1-й этап) (старая версия) */
			/*TransformGridFrom_AcToHeat(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic,
				Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, tmpe_ac2_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);*/
				// вставляем без границ
			//	Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, tmpe2, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
				// fun21(Nx_heat, Ny_heat, Nz_heat, tmpe2, tmpi1, tmpi2, dz_heat, param, mt);

				/*plt_TE_ht.SetParametrsOnPlotColor(7, "Te", "x,mkm", "y,mkm", (1e+4* 1e-2* dx_heat* Nx_heat), (1e+4* 1e-2* dy_heat* Ny_heat));
				plt_TE_ht.SetDataOnPlotColor3D(7, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., 3, xy);
				plt_TE_ht.ShowDataOnPlotColor(7, "file tmpe2 ht vo l(x,y div 2,z) ", true);*/

				//plt_TE_ht.SetParametrsOnPlotColor(7, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
				//plt_TE_ht.SetDataOnPlotColor3D(7, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., Ny_heat / 2, xz);
				//plt_TE_ht.ShowDataOnPlotColor(7, "file tmpe2 ht vol (x,y div 2,z) ", true);
				//Sleep(1000);

				// from ht to ac /* сшивка границ tmpe (2-й этап) */
		TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, tmpe2, tmpe_ac2_new, /*Te_old_x_new_y,*/ spl_Te);

		TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1,
			down_boundary_z_acoustic + 1, Nx_part_ac + 1, Ny_part_ac + 1, dx_acoustic, dy_acoustic, dz_acoustic, tmpe2, tmpe_ac2_new, /*Te_old_x_new_z,*/ spl_Te);
		fun21_for_ac(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, tmpe_ac2_new, tmpi_ac1_new, tmpi_ac2_new, dz_acoustic, param, mt, points_rupture_ac);

		//plt_TE_ac.SetParametrsOnPlotColor(7, "Te", "x,mkm", "y,mkm", (1e+4* 1e-2* dx_acoustic* Nx_part_ac), (1e+4* 1e-2* dy_acoustic* Ny_part_ac));
		//plt_TE_ac.SetDataOnPlotColor3D(7, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, tmpe_ac2_new, 300., down_boundary_z_acoustic, xy);
		//plt_TE_ac.ShowDataOnPlotColor(7, "file tmpe2 ac (x,y) ", true);

		//plt_TE_ac.SetParametrsOnPlotColor(7, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		//plt_TE_ac.SetDataOnPlotColor3D(7, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, tmpe_ac2_new, 300., Ny_part_ac / 2, xz);
		//plt_TE_ac.ShowDataOnPlotColor(7, "file tmpe2 ac (x,y div 2,z) ", true); // кусок
		//Sleep(1000);
		//cout << "Stop proverka 3" << endl;
		//system("pause");

		/* сшивка границ tmpe (2-й этап) (старая версия)*/

		//CopyDataArray3D(a2x_ac_new, tmpe_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1); // сохраняем внутренноть
		////plt_TE_ac.SetParametrsOnPlotColor(7, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		////plt_TE_ac.SetDataOnPlotColor3D(7, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, a2x_ac_new, 300., Ny_part_ac / 2, xz);
		////plt_TE_ac.ShowDataOnPlotColor(7, "file tmpe22 ac (x,y div 2,z) ", true);
		////Sleep(1000);
		//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpe2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_ac2_new, Te_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(tmpe_ac2_new, a2x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		//fun21_for_ac(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, tmpe_ac2_new, tmpi_ac1_new, tmpi_ac2_new, dz_acoustic, param, mt, points_rupture_ac);

		//plt_TE_ht.SetParametrsOnPlotColor(8, "Te", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		//plt_TE_ht.SetDataOnPlotColor3D(8, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe2, 300., Ny_heat / 2, xz);
		//plt_TE_ht.ShowDataOnPlotColor(8, "Final in HT Total", true);
		//Sleep(1000);
		//cout << " STOP STOP STOP STOP STOP" << endl;
		//system("pause");

		//plt_TE_ac.SetParametrsOnPlotColor(9, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		//plt_TE_ac.SetDataOnPlotColor3D(9, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, tmpe_ac2_new, 300., Ny_part_ac / 2, xz);
		//plt_TE_ac.ShowDataOnPlotColor(9, "file tmpe2 ac (x,y div 2,z) ", true);
		//Sleep(1000);
		//cout << " STOP STOP STOP STOP STOP tmpe2_ac" << endl;
		//system("pause");

		/* сшивка границ tmpi (1-й этап) (2-я версия)*/
		Transform_2D_GridFrom_AcToHeat_xy(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, tmpi_ac2_new, tmpi2, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

		Transform_2D_GridFrom_AcToHeat_xz(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx_heat, dz_heat, tmpi_ac2_new, tmpi2, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		fun21(Nx_heat, Ny_heat, Nz_heat, tmpe2, tmpi1, tmpi2, dz_heat, param, mt);

		// from ht to ac /* сшивка границ tmpi (2-й этап) */
		TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, tmpi2, tmpi_ac2_new, /*Te_old_x_new_y,*/ spl_Te);

		TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1,
			down_boundary_z_acoustic + 1, Nx_part_ac + 1, Ny_part_ac + 1, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, tmpi_ac2_new, /*Te_old_x_new_z,*/ spl_Te);
		fun21_for_ac(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, tmpe_ac2_new, tmpi_ac1_new, tmpi_ac2_new, dz_acoustic, param, mt, points_rupture_ac);


		/*// сшивка границ tmpi (1-й этап) (старая версия)
		TransformGridFrom_AcToHeat(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic,
			Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, tmpi_ac2_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
		Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, tmpi2, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		fun21(Nx_heat, Ny_heat, Nz_heat, tmpe2, tmpi1, tmpi2, dz_heat, param, mt);

		 //сшивка границ tmpe (2-й этап)
		CopyDataArray3D(a2y_ac_new, tmpi_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1); // сохраняем внутренноть
		TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpi_ac2_new, Te_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(tmpi_ac2_new, a2y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		fun21_for_ac(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, tmpe_ac2_new, tmpi_ac1_new, tmpi_ac2_new, dz_acoustic, param, mt, points_rupture_ac);*/

		//// подготвка массивов для расчета акустики
		//Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, a1x, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		//TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
		//	dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, a1x_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

		/*plt_TE_ht.SetParametrsOnPlotColor(10, "a1x", "z,mkm", "x,mkm", (1e+4* dz_heat* down_boundary_z_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_part_ht));
		plt_TE_ht.SetDataOnPlotColor3D(10, 100, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe_for_transform_new, 1., Ny_part_ht / 2, xz);
		plt_TE_ht.ShowDataOnPlotColor(10, "file part a1x ht 1  (x,y div 2,z) ", true);
		Sleep(1000);

		system("pause");

		plt_TE_ac.SetParametrsOnPlotColor(11, "a1x", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		plt_TE_ac.SetDataOnPlotColor3D(11, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, a1x_ac_new, 1., Ny_part_ac / 2, xz);
		plt_TE_ac.ShowDataOnPlotColor(11, "file part ac 1 (x,y div 2,z) ", true);
		Sleep(1000);
		cout << " STOP STOP VERNO " << endl;
		system("pause");*/

		/*Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, a1y, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, a1y_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

		Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, a1z, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, a1z_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

		Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, b1x, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
		TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, b1x_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);*/

			/*plt_TE_ht.SetParametrsOnPlotColor(10, "b1x", "z,mkm", "x,mkm", (1e+4* dz_heat* down_boundary_z_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_part_ht));
			plt_TE_ht.SetDataOnPlotColor3D(10, 100, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, tmpe_for_transform_new, 1., Ny_part_ht / 2, xz);
			plt_TE_ht.ShowDataOnPlotColor(10, "file b1x ht start (x,y div 2,z) ", true);
			Sleep(1000);

			system("pause");

			plt_TE_ac.SetParametrsOnPlotColor(11, "b1x", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
			plt_TE_ac.SetDataOnPlotColor3D(11, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, b1x_ac_new, 1., Ny_part_ac / 2, xz);
			plt_TE_ac.ShowDataOnPlotColor(11, "file part b1x ac start (x,y div 2,z) ", true);
			Sleep(1000);*/

			/*Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, b1y, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
			TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
				dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, b1y_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

			Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, b1z, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
			TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
				dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, b1z_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

			Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, e1, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
			TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
				dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, e1_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

			Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, V00, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
			TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
				dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, V00_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);

			Copy_Data_From_Big_Array3D_To_Small_Array3D(tmpe_for_transform_new, V1, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);
			TransformGrid(Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
				dx_acoustic, dy_acoustic, dz_acoustic, tmpe_for_transform_new, V1_ac_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);*/

		if (n == 1)
		{
			// расчет скорости   b2x, b2y, b2z
			Calculation1(a1x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, b1x_ac_new, b1y_ac_new, b1z_ac_new, b2x_ac_new, b2y_ac_new, b2z_ac_new,
				a1x_ac_new, a1y_ac_new, a1z_ac_new, a2y_ac_new, a2z_ac_new, e1_ac_new, dx_heat, dy_heat, dz_heat, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_heat, Ny_heat, Nz_heat, left_boundary_x_heat, right_boundary_x_heat,
				left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, new_interv_acoustic_z, points_rupture_ac);
		}
		else
		{
			// расчет скорости b2x, b2y, b2z
			Calculation10(a1x, b1x, b2x, a1z, b1y, b2y, b1z, b2z, a1y, e1, V00, V1, a1x_ac_new, b1x_ac_new, b2x_ac_new, a1z_ac_new, b1y_ac_new, b2y_ac_new, b1z_ac_new,
				b2z_ac_new, a1y_ac_new, e1_ac_new, V00_ac_new, V1_ac_new, dx_heat, dy_heat, dz_heat, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_heat, Ny_heat, Nz_heat,
				left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
				new_interv_acoustic_z, points_rupture_ac);
		}

		/*plt_TE_ht.SetParametrsOnPlotColor(5, "b2x", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		plt_TE_ht.SetDataOnPlotColor3D(5, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, b2x, 1., Ny_heat / 2, xz);
		plt_TE_ht.ShowDataOnPlotColor(5, "file b2x ht (x,y div 2,z) ", true);
		Sleep(1000);
		plt_TE_ac.SetParametrsOnPlotColor(6, "b2x", "z,mkm", "x,mkm", (1e+4* dz_acoustic* down_boundary_z_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_part_ac));
		plt_TE_ac.SetDataOnPlotColor3D(6, 100, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, b2x_ac_new, 1., Ny_part_ac / 2, xz);
		plt_TE_ac.ShowDataOnPlotColor(6, "file b2x ac (x,y div 2,z) ", true);
		Sleep(1000);
		cout << " b2x STOP " << endl;
		system("pause");*/

		/*TransformGridFrom_AcToHeat(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic,
			Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, b2x_ac_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
		Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, b2x, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);*/

		/*plt_TE_ht.SetParametrsOnPlotColor(8, "b2x", "z,mkm", "x,mkm", (1e+4* dz_heat* Nz_heat / 5e+5), (1e+4* 1e-2* dx_heat* Nx_heat));
		plt_TE_ht.SetDataOnPlotColor3D(8, 100, Nx_heat, Ny_heat, Nz_heat, 1e+4* dx_heat* 1e-2, 1e+4* dy_heat* 1e-2, 1e+4* dz_heat / 5e+5, b2x, 1., Ny_heat / 2, xz);
		plt_TE_ht.ShowDataOnPlotColor(8, "Final in HT", true);
		Sleep(1000);
		cout << " STOP STOP STOP STOP STOP" << endl;*/
		//system("pause");

		//cout << " STOP STOP VERNO " << endl;
		//("pause");

		//Null_Array3D(tmpe_for_transform_new, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1);
		//cout << " STOP STOP STOP " << endl;
		//system("pause");


			/* сшивка границ скорости (1-й этап) (2-я версия)*/
		Transform_2D_GridFrom_AcToHeat_xy(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, b2x_ac_new, b2x, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

		Transform_2D_GridFrom_AcToHeat_xz(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx_heat, dz_heat, b2x_ac_new, b2x, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

		// from ht to ac /* сшивка границ скорости (2-й этап) */
		TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, b2x, b2x_ac_new, /*Te_old_x_new_y,*/ spl_Te);

		TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1,
			down_boundary_z_acoustic + 1, Nx_part_ac + 1, Ny_part_ac + 1, dx_acoustic, dy_acoustic, dz_acoustic, b2x, b2x_ac_new, /*Te_old_x_new_z,*/ spl_Te);

		/* сшивка границ скорости (1-й этап) (2-я версия)*/
		Transform_2D_GridFrom_AcToHeat_xy(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, b2y_ac_new, b2y, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

		Transform_2D_GridFrom_AcToHeat_xz(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx_heat, dz_heat, b2y_ac_new, b2y, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

		// from ht to ac /* сшивка границ скорости (2-й этап) */
		TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, b2y, b2y_ac_new, /*Te_old_x_new_y,*/ spl_Te);

		TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1,
			down_boundary_z_acoustic + 1, Nx_part_ac + 1, Ny_part_ac + 1, dx_acoustic, dy_acoustic, dz_acoustic, b2y, b2y_ac_new, /*Te_old_x_new_z,*/ spl_Te);

		/* сшивка границ скорости (1-й этап) (2-я версия)*/
		Transform_2D_GridFrom_AcToHeat_xy(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, b2z_ac_new, b2z, /*Te_old_x_new_y,*/ spl_Te, coeff_big_z,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat);

		Transform_2D_GridFrom_AcToHeat_xz(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat + 1, dx_heat, dz_heat, b2z_ac_new, b2z, /*Te_old_x_new_z,*/ spl_Te, coeff_big_y, coeff_big_x,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

		// from ht to ac /* сшивка границ скорости (2-й этап) */
		TransformGrid_2D_xy(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			dx_acoustic, dy_acoustic, b2z, b2z_ac_new, /*Te_old_x_new_y,*/ spl_Te);

		TransformGrid_2D_xz(left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, dx_heat, dy_heat, dz_heat, Nx_part_ac + 1, Ny_part_ac + 1,
			down_boundary_z_acoustic + 1, Nx_part_ac + 1, Ny_part_ac + 1, dx_acoustic, dy_acoustic, dz_acoustic, b2z, b2z_ac_new, /*Te_old_x_new_z,*/ spl_Te);

		/*//сшивка скорости (1-й этап)
		//TransformGridFrom_AcToHeat(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic,
		//	Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, b2x_ac_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
		//Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, b2x, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

		//TransformGridFrom_AcToHeat(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic,
		//	Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, b2y_ac_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
		//Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, b2y, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

		//TransformGridFrom_AcToHeat(Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1, dx_acoustic, dy_acoustic, dz_acoustic,
		//	Nx_part_ht + 1, Ny_part_ht + 1, down_boundary_z_heat + 1, dx_heat, dy_heat, dz_heat, b2z_ac_new, tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, spl_Te);
		//Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_for_transform_new, b2z, left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat);

		//// сшивка скорости (3 компоненты ниже) (внеш границы для мелкой сетки) (2-й этап)
		//CopyDataArray3D(b1x_ac_new, b2x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, b2x, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D(b2x_ac_new, Te_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(b2x_ac_new, b1x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		//CopyDataArray3D(b1x_ac_new, b2y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, b2y, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D(b2y_ac_new, Te_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(b2y_ac_new, b1x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		//CopyDataArray3D(b1x_ac_new, b2z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, b2z, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D(b2z_ac_new, Te_acoustic, Nx_acoustic, Ny_acoustic, Nz_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		//Copy_Data_From_Big_Array3D_To_Small_Array3D_Withoutboundary(b2z_ac_new, b1x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);*/

		GRU1(e2, Nx_heat, Ny_heat);
		GRU1(e2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1);
		GRU2(b2x, b2y, b2z, Nx_heat, Ny_heat, Nz_heat);

		Calculation2(a1x, a1y, a1z, a2x, a2y, a2z, b2x, b2y, b2z, V2, V1, V00, e2, tmpe1, tmpi1, a1x_ac_new, a2x_ac_new, a1y_ac_new, a2y_ac_new, a1z_ac_new, a2z_ac_new, b2x_ac_new, b2y_ac_new,
			b2z_ac_new, V2_ac_new, V1_ac_new, V00_ac_new, tmpe_ac1_new, tmpi_ac1_new, e2_ac_new, e1_ac_new, dx_heat, dy_heat, dz_heat, Nx_part_ht, Ny_part_ht, dx_acoustic, dy_acoustic, dz_acoustic, dt, CC0, Nx_heat, Ny_heat, Nz_heat,
			left_boundary_x_heat, right_boundary_x_heat, left_boundary_y_heat, right_boundary_y_heat, down_boundary_z_heat, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1,
			V0, p0v, /*tmpe_for_transform_new, Tei_old_xy_new_z, Tei_old_y_new_xz, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac,*//* Te_old_x_new_y, Te_old_x_new_y, Te_old_x_new_z, Te_old_x_new_z, */coeff_big_x, coeff_big_y, coeff_big_z, Melt_metal, mt, param, spl_C_e_on_T, spl_Te, spl_Ti, spl_e, current_time_frame, points_rupture_ac, points_rupture_heat, new_interv_acoustic_z, spl_C_l_on_T
		/*Nx_acoustic, Ny_acoustic, Nz_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic*/);

		//// для отрисоки
		//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, e2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		//Copy_Data_From_Small_Array3D_To_Big_Array3D(e2_ac_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

		//plt_TE_ht.SetParametrsOnPlotColor(29, "P", "z,mkm", "x,mkm", (1e+4* dz_acoustic* Nz_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_acoustic));
		//plt_TE_ht.SetDataOnPlotColor3D(29, 100, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, Te_acoustic, p0v_s, Ny_acoustic / 2, xz);
		//plt_TE_ht.ShowDataOnPlotColor(29, "Final P in AC", true);
		//Sleep(1000);
		//cout << " STOP STOP STOP STOP STOP ACoustic" << endl;
		//cout << " dalshe otrisovka " << endl;
		////system("pause");
		//plt_TE_ht.Close_and_open_files_for_replot(for_close_and_open_29);


		// действия ниже - это вынужденая мера по причине того, что гнуплот не хочет рисовать данные с разных сеток
		// выкинуть Te_acoustic и оставить Te_acousrtic - этот массив будет служить плацдармом для отрисовки
		/*TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpe2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpe_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, Ti_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		Copy_Data_From_Small_Array3D_To_Big_Array3D(tmpi_ac2_new, Ti_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		fun21_for_ac(Nx_acoustic, Ny_acoustic, Nz_acoustic, Te_acoustic, Ti_acoustic, Ti_acoustic, dz_acoustic, param, mt, points_rupture_ac);*/

		/*plt_TE_ac.SetParametrsOnPlotColor(9, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* Nz_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_acoustic));
		plt_TE_ac.SetDataOnPlotColor3D(9, 100, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, Te_acoustic, 300., Ny_acoustic / 2, xz);
		plt_TE_ac.ShowDataOnPlotColor(9, "Final Final Te in AC  Total ", true);
		Sleep(1000);
		plt_TE_ac.SetParametrsOnPlotColor(10, "Te", "z,mkm", "x,mkm", (1e+4* dz_acoustic* Nz_acoustic / 5e+5), (1e+4* 1e-2* dx_acoustic* Nx_acoustic));
		plt_TE_ac.SetDataOnPlotColor3D(10, 100, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4* dx_acoustic* 1e-2, 1e+4* dy_acoustic* 1e-2, 1e+4* dz_acoustic / 5e+5, Ti_acoustic, 300., Ny_acoustic / 2, xz);
		plt_TE_ac.ShowDataOnPlotColor(10, "Final Final  Ti in AC  Total ", true);
		Sleep(1000);
		cout << " STOP STOP STOP STOP STOP" << endl;

		cout << "FINAL FINAL FINAL " << endl;
		system("pause");*/

		/*cout << "It took me %d clicks (%f seconds)." << endl <<
			(std::clock() - t) / (double)CLOCKS_PER_SEC << endl;
		sum_t += (std::clock() - t) / (double)CLOCKS_PER_SEC;
		cout << " General duration calculation = " << sum_t << endl;
		cout << endl;*/

		double level_depth = 0.;

		tt = n * dt * t0 * 1e+15; // vremja v fs
		cout << " time, fs = " << tt << endl;

		if (points_rupture_ac.empty()) {}
		else
		{
			My_unique(points_rupture_ac);
			sort(points_rupture_ac.begin(), points_rupture_ac.end(), comp_x);
			MySort_Point3D_y(points_rupture_ac);
			MySort_Point3D_z(points_rupture_ac);

			//tmpe1_ac - проверить на предмет отрыва
			for (it = points_rupture_ac.begin(); it != points_rupture_ac.end(); it++)
			{
				tmpe_ac2_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
				tmpi_ac2_new[(*it).index_x][(*it).index_y][(*it).index_z] = 1e-16;
			}
		}

		// Отрисовка

		/*if (!points_rupture_ac.empty())
		{
			TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(tmpi_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);
		}*/

		//if ((tmpi1[Nx_heat / 2][Ny_heat / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // 2-й вариант определения макс глубины расплава
		//if ((Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // 2-й вариант определения макс глубины расплава
		if (!points_rupture_ac.empty())
		{
			for (int i = 0; i < Nz_acoustic_show; i++)
			{
				if ((tmpi_ac2_new[Nx_acoustic_show / 2][Ny_acoustic_show / 2][i] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
				{
					level_depth = (1e+4 * dz_acoustic * i / kabs);
				}
			}
			count_depth++;

			fout_depth2 << level_depth << endl;
		}

		//vector<double***> TeTi = { tmpe_ac2_new , tmpi_ac2_new };//????
		//vector<double***> Ti = { tmpi_ac2_new };//????
		//plt.SetDataOnPlot3D(1, 0, NULL, T00, Nx_acoustic / 2, Ny_acoustic / 2, 1, 2, 0, tt, TeTi, fun_on_t);
		//plt.SetDataOnPlot3D(10, 0, NULL, T00, Nx_acoustic / 2, Ny_acoustic / 2, 1, 1, 0, tt, Ti, fun_on_t);

		//if ((tmpi1[Nx / 2][Ny / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Зона плавления
		//if ((Te_acoustic[Nx_acoustic / 2][Ny_acoustic / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) 

		// Зона плавления
		if (/*!points_rupture_ac.empty()*/ ((*it_fix_time_anim) >= 30000000) && (it_fix_time_anim != moments_fix_time_anim.end())) // Зона плавления + явление отрыва
		{
			for (int i = 0; i < Nx_acoustic_show; i++)
			{
				for (int k = 0; k < Nz_acoustic_show; k++)
				{
					if ((tmpi_ac2_new[i][Ny_acoustic_show / 2][k] * param.T00) >= (Melt_metal[mt].T_melting - 1.))
					{
						Array3D_for_Plot[i][k] = (Melt_metal[mt].T_melting + 200.) / param.T00;
					}

					if (((tmpi_ac2_new[i][Ny_acoustic_show / 2][k] * param.T00) < (Melt_metal[mt].T_melting - 1.)) && ((tmpi_ac2_new[i][Ny_acoustic_show / 2][k] * param.T00) >= param.T00))
					{
						Array3D_for_Plot[i][k] = (Melt_metal[mt].T_melting * 0.5) / param.T00;
					}

					if ((tmpi_ac2_new[i][Ny_acoustic_show / 2][k]) == 1e-16)
					{
						Array3D_for_Plot[i][k] = 1e-16;
					}
				}
			}

			plt.SetDataOnPlotColor2D(13, current_time_frame, Nz_acoustic_show, Nx_acoustic_show, 1e+4 * dz_acoustic / kabs, 1e+4 * dx_acoustic * r0, /*Te_acoustic*/ Array3D_for_Plot, param.T00);

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(tmpi_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			///CopyDataArray3D(Array3D_for_Plot, tmpi_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

			//for (int i = 0; i < Nx_acoustic; i++)
			//{
			//	for (int j = 0; j < Ny_acoustic; j++)
			//	{
			//		if ((Te_acoustic[i][j][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // точно расплав (температуру устанавливаем фиктивную для контраста цветов)
			//		{
			//			//Te_acoustic[i][j][1] = (Melt_metal[mt].T_melting + 200.) / param.T00;
			//		}

			//		if (((Te_acoustic[i][j][1] * param.T00) < (Melt_metal[mt].T_melting + 1.)) && ((Te_acoustic[i][j][1] * param.T00) >= param.T00))
			//		{
			//			Te_acoustic[i][j][1] = (Melt_metal[mt].T_melting * 0.5) / param.T00;
			//		}
			//	}
			//}

			//for (int k = 0; k < Nz_acoustic_show; k++)
			//{
			//	for (int i = 0; i < Nx_acoustic_show; i++)
			//	{
			//		for (int j = 0; j < Ny_acoustic_show; j++)
			//		{
			//			if ((/*Te_acoustic*/tmpi_ac2_new[i][j][k] * param.T00) >= (Melt_metal[mt].T_melting - 1.)) // точно расплав (температуру устанавливаем фиктивную для контраста цветов)
			//			{
			//				/*Te_acoustic*/ Array3D_for_Plot[i][j][k] = (Melt_metal[mt].T_melting + 200.) / param.T00;
			//			}

			//			if (((/*Te_acoustic*/tmpi_ac2_new[i][j][k] * param.T00) < (Melt_metal[mt].T_melting - 1.)) && ((/*Te_acoustic*/tmpi_ac2_new[i][j][k] * param.T00) >= param.T00))
			//			{
			//				/*Te_acoustic*/  Array3D_for_Plot[i][j][k] = (Melt_metal[mt].T_melting * 0.5) / param.T00;
			//			}
			//		}
			//	}
			//}

			//int z_depth = 100; // нм
			//int z_max = z_depth / (1e+4 * dz_acoustic / kabs * 1000);
			//for (int index_z = 0; index_z < z_max; index_z++)
			//{
			//	plt.SetDataOnPlotColor3D(14, current_time_frame, Nx_acoustic, Ny_acoustic, Nz_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, Te_acoustic, param.T00, index_z, xy);
			//	namefile.clear();
			//	namefile = "Melting zone (yx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, z = " + ConvertNumToString(z_depth) + " nm";
			//	plt.ShowDataOnPlotColor(14, namefile, true);
			//	Sleep(1500);
			//	plt.Close_and_open_files_for_replot(for_close_and_open);
			//}

			namefile.clear();
			namefile = "Melting zone (zx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(13, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "Melting zone (yx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(14, namefile, true);
			//Sleep(1000);

			/*Null_Array(massive_melting_tmp1, Nz_acoustic, Nx_acoustic);
			Null_Array(massive_melting_tmp2, Nx_acoustic, Ny_acoustic);*/

			plt.Close_and_open_files_for_replot(for_close_and_open);
		}

		// Preasure zone
		if (!points_rupture_ac.empty())/*(e2[Nx_acoustic / 2][Ny_acoustic / 2][1] * p0v) <= -sigma*/
		{
			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, e2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(e2_ac_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

		///	CopyDataArray3D(Array3D_for_Plot, e2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

			/*SetDataArray(massive_melting_tmp1, Nz_acoustic, Nx_acoustic, 2000. / p0v_s);
			SetDataArray(massive_melting_tmp2, Nx_acoustic, Ny_acoustic, 2000. / p0v_s);*/
			///for (int i = 0; i < Nx_acoustic_show; i++)
			{
				///for (int k = 0; k < Nz_acoustic_show; k++)
				{
					//if ((/*Te_acoustic*/e2_ac_new[i][Ny_acoustic_show / 2][k] * p0v_s) == (1e-16 * p0v_s))
					{
						//massive_melting_tmp1[i][Ny_acoustic / 2][k] = tmpi1[i][Ny_acoustic / 2][k];
						///massive_melting_tmp1[k][i] = Te_acoustic[i][Ny_acoustic / 2][k];
					}
					//else
					{
						/*Te_acoustic*/ //Array3D_for_Plot[i][Ny_acoustic_show / 2][k] = 2000. / p0v_s;
					}
				}
			}

			///plt.SetDataOnPlotColor3D(20, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ Array3D_for_Plot, p0v_s, Ny_acoustic_show / 2, xz);

			///CopyDataArray3D(Array3D_for_Plot, e2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, e2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(e2_ac_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			///for (int i = 0; i < Nx_acoustic_show; i++)
			{
				///for (int j = 0; j < Ny_acoustic_show; j++)
				{
					//if ((tmpi1[i][j][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.))
					//if ((/*Te_acoustic*/// e2_ac_new[i][j][1] * p0v_s) == (1e-16 * p0v_s))
					{
						//massive_melting_tmp2[i][j][1] = tmpi1[i][j][1];
						///massive_melting_tmp2[i][j] = Te_acoustic[i][j][1];
					}
					//else
					{
						/*Te_acoustic*/ //Array3D_for_Plot[i][j][1] = 2000. / p0v_s;
					}
				}
			}

			//plt.SetDataOnPlotColor2D(20, current_time_frame, Nz_acoustic, Nx_acoustic, 1e+4 * dz_acoustic / kabs, 1e+4 * dy_acoustic * r0, massive_melting_tmp1, p0v_s);
			//plt.SetDataOnPlotColor2D(21, current_time_frame, Nx_acoustic, Ny_acoustic, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, massive_melting_tmp2, p0v_s);
			///plt.SetDataOnPlotColor3D(21, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ Array3D_for_Plot, p0v_s, 1, xy);

			namefile.clear();
			namefile = "Preasure zone (zx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			///	plt.ShowDataOnPlotColor(20, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "Preasure zone (yx) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			///	plt.ShowDataOnPlotColor(21, namefile, true);
			//Sleep(1000);

			/*Null_Array(massive_melting_tmp1, Nz_acoustic, Nx_acoustic);
			Null_Array(massive_melting_tmp2, Nx_acoustic, Ny_acoustic);*/

			plt.Close_and_open_files_for_replot(for_close_and_open);
		}

		// Профили температур ионной решетки в разные моменты
		if ((tmpi1[Nx_heat / 2][Ny_heat / 2][1] * param.T00) >= (Melt_metal[mt].T_melting + 1.)) // Профили температур ионной решетки в разные моменты
		{
			if ((int(ceil(tt)) % 1000 == 0 || int(trunc(tt)) % 1000 == 0) && current_number_line_melting <= 10)
			{
				//plt.SetDataOnPlot3D(12, 0, tmpi1, param.T00, NULL, Ny_heat / 2, 1, 12, current_number_line_melting, NULL, empty, fun_on_x);
				//plt.SetDataOnPlot3D(11, 0, tmpi1, param.T00, Nx_heat / 2, Ny_heat / 2, NULL, 12, current_number_line_melting, NULL, empty, fun_on_z);
				current_number_line_melting++;
				Legenda_melting_tmp.push_back(ConvertNumToStringdouble(tt) + " fs");
			}
		}

		if ((n == (*it_fix_time_anim))) // блок отрисовки данных и запись данных в файлы
		{
			namefile.clear();
			namefile = "Te,Ti (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			//plt.ShowDataOnPlot2D(1, false, 2, TeTiTfs, namefile, true);
			namefile.clear();
			namefile = "Ti (t,fms) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs, 1";
			//plt.ShowDataOnPlot2D(10, false, 1, TiTfs, namefile, true);

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpe2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(tmpe_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			/*пока скрыл*/
			///plt.SetDataOnPlotColor3D(2, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ tmpe_ac2_new, param.T00, Ny_acoustic_show / 2, xz);
			///plt.SetDataOnPlotColor3D(6, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ tmpe_ac2_new, param.T00, 1, xy);

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(tmpi_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			plt.SetDataOnPlotColor3D(3, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ tmpi_ac2_new, param.T00, Ny_acoustic_show / 2, xz);
			plt.SetDataOnPlotColor3D(7, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ tmpi_ac2_new, param.T00, 1, xy);

			int z_point = 0;// 5;
			Legenda_Ti_x.clear();
			int total_count_z = 30;
			for (int Ti_P_x_profile = 1; Ti_P_x_profile <= total_count_z; Ti_P_x_profile++)// цикл по номеру столбца 
			{
				// Ny_heat / 2 = 50   Ny_ac / 2 = 100
				// если меняю общее кол-во линий, то меняем в методе plt.SetGridOnPlot3D
				//plt.SetDataOnPlot3D(23, current_time_frame, /*Te_acoustic*/ tmpi_ac2_new, param.T00, NULL, Ny_acoustic_show / 2, z_point, total_count_z, Ti_P_x_profile, NULL, empty, fun_on_x);
				Legenda_Ti_x.push_back(ConvertNumToStringdouble(1e+4 * z_point * dz_acoustic / kabs * 1000) + " nm");
				z_point++;// = 5;// номер узла по легенде глубин (по z)
			}

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, e2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(e2_ac_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			plt.SetDataOnPlotColor3D(4, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ e2_ac_new, p0v_s, Ny_acoustic_show / 2, xz);
			plt.SetDataOnPlotColor3D(5, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, /*Te_acoustic*/ e2_ac_new, p0v_s, 1, xy);

			/*plt.SetDataOnPlotColor3D(6, current_time_frame, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpe2, param.T00, 1, xy);
			plt.SetDataOnPlotColor3D(7, current_time_frame, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, tmpi2, param.T00, 1, xy);*/

			plt.SetDataOnPlotColor3D(24, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, a2z_ac_new, pow(kabs, -1), Ny_acoustic_show / 2, xz);
			plt.SetDataOnPlotColor3D(25, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, a2x_ac_new, r0, Ny_acoustic_show / 2, xz);
			plt.SetDataOnPlotColor3D(26, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, a2z_ac_new, pow(kabs, -1), 1, xy);
			plt.SetDataOnPlotColor3D(27, current_time_frame, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, a2x_ac_new, r0, 1, xy);

			z_point = 0;// 5;
			Legenda_P_x.clear();
			total_count_z = 30;
			for (int Ti_P_x_profile = 1; Ti_P_x_profile <= total_count_z; Ti_P_x_profile++)// цикл по номеру столбца 
			{
				// Ny_heat / 2 = 50   Ny_ac / 2 = 100
				// если меняю общее кол-во линий, то меняем в методе plt.SetGridOnPlot3D
				//plt.SetDataOnPlot3D(17, current_time_frame, /*Te_acoustic*/ e2_ac_new, p0v_s, NULL, Ny_acoustic_show / 2, z_point, total_count_z, Ti_P_x_profile, NULL, empty, fun_on_x);
				Legenda_P_x.push_back(ConvertNumToStringdouble(1e+4 * z_point * dz_acoustic / kabs * 1000) + " nm");
				z_point++;// = 5;// номер узла по легенде глубин (по z)
			}

			namefile.clear();
			namefile = "P (x,mkm) Legenda z" + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(17, false, 10, Legenda_Ti_x, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "Ti (x,mkm) Legenda z" + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(23, false, 10, Legenda_P_x, namefile, true);
			//Sleep(1000);

			if (npxzt <= count_of_lines) // если нцжно будет то обертываем в определенные моменты времени
			{
				//plt.SetDataOnPlot3D(8, 0, e2, p0v_s, NULL, Ny_acoustic / 2, 1, count_of_lines, npxzt, NULL, empty, fun_on_x);
				//plt.SetDataOnPlot3D(9, 0, e2, p0v_s, Nx_acoustic / 2, Ny_acoustic / 2, NULL, count_of_lines, npxzt, NULL, empty, fun_on_z);
				npxzt++;
				Legenda.push_back(ConvertNumToStringdouble(tt) + " fs");
			}

			namefile.clear();
			namefile = "file Te(x, N div 2, z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(2, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file Ti(x, N div 2, z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(3, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file P (x,N div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(4, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file P (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(5, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file Te (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(6, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file Ti (x,Ny,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(7, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file a2z (x,Ny div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(24, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file a2x (x,Ny div 2,z) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(25, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file a2z (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(26, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "file a2x (x,y,1) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlotColor(27, namefile, true);
			//Sleep(1000);

			namefile.clear();

			// отрисовка температуры в тепл сетке
			/*plt.SetGridOnPlot3D(15, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 2, fun_on_x);
			plt.SetGridOnPlot3D(16, Nx_heat, Ny_heat, Nz_heat, 1e+4 * dx_heat * r0, 1e+4 * dy_heat * r0, 1e+4 * dz_heat / kabs, 2, fun_on_z);*/

			// отрисовка температуры в акуст сетке

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpe2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(tmpe_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			// отрисовка температуры в акуст сетке
			plt.SetGridOnPlot3D(15, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 2, fun_on_x);
			plt.SetGridOnPlot3D(16, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 2, fun_on_z);
			plt.SetDataOnPlot3D(15, current_time_frame, /*Te_acoustic*/ tmpe_ac2_new, param.T00, NULL, Ny_acoustic_show / 2, 1, 2, 1, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(16, current_time_frame, /*Te_acoustic*/ tmpe_ac2_new, param.T00, Nx_acoustic_show / 2, Ny_acoustic_show / 2, NULL, 2, 1, NULL, empty, fun_on_z);
			plt.SetDataOnPlot3D(15, current_time_frame, /*Te_acoustic*/ tmpi_ac2_new, param.T00, NULL, Ny_acoustic_show / 2, 1, 2, 2, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(16, current_time_frame, /*Te_acoustic*/ tmpi_ac2_new, param.T00, Nx_acoustic_show / 2, Ny_acoustic_show / 2, NULL, 2, 2, NULL, empty, fun_on_z);


			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, tmpi2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
			//Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(tmpi_ac2_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			//TransformGrid(Nx_heat, Ny_heat, Nz_heat, dx_heat, dy_heat, dz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, dx_acoustic, dy_acoustic, dz_acoustic, e2, Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac, spl_Te);
		//	Copy_Data_From_Small_Array3D_To_Big_Array3D_Withboundary(e2_ac_new, Te_acoustic, left_boundary_x_acoustic, right_boundary_x_acoustic, left_boundary_y_acoustic, right_boundary_y_acoustic, down_boundary_z_acoustic);

			plt.SetGridOnPlot3D(18, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_x);
			plt.SetGridOnPlot3D(19, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 1, fun_on_z);
			plt.SetDataOnPlot3D(18, current_time_frame, /*Te_acoustic*/ e2_ac_new, p0v_s, NULL, Ny_acoustic_show / 2, 1, 1, 1, NULL, empty, fun_on_x);
			plt.SetDataOnPlot3D(19, current_time_frame,/* Te_acoustic*/ e2_ac_new, p0v_s, Nx_acoustic_show / 2, Ny_acoustic_show / 2, NULL, 1, 1, NULL, empty, fun_on_z);


			plt.SetGridOnPlot3D(28, Nx_acoustic_show, Ny_acoustic_show, Nz_acoustic_show, 1e+4 * dx_acoustic * r0, 1e+4 * dy_acoustic * r0, 1e+4 * dz_acoustic / kabs, 2, fun_on_z);
			plt.SetDataOnPlot3D(28, current_time_frame,/* Te_acoustic*/ a2z_ac_new, pow(kabs, -1) * 1e+4, Nx_acoustic_show / 2, Ny_acoustic_show / 2, NULL, 2, 2, NULL, empty, fun_on_z);
			//??plt.SetDataOnPlot3D(28, current_time_frame,/* Te_acoustic*/ a2z_ac_new, 1., Nx_acoustic_show / 2, Ny_acoustic_show / 2, NULL, 2, 2, NULL, empty, fun_on_z);

			//plt.SetDataOnPlot1D(28, current_time_frame,/* Te_acoustic*/ a2z_ac_new, 1., Nx_acoustic_show / 2, Ny_acoustic_show / 2, NULL, 2, 2, NULL, empty, fun_on_z);

			// ЗАПИСЬ ОТДЕЛЬНЫХ СЛАГАЕМЫХ УРАВНЕНИЯ СОСТОЯНИЯ МИ-ГРЮНАЙЗЕНА 
			///spl_Te.InterpolateFast(1, tmpe2, Nx_acoustic / 2, Ny_acoustic / 2);
		///	spl_Ti.InterpolateFast(1, tmpi2, Nx_acoustic / 2, Ny_acoustic / 2);
			// может значения Т перезаписать в матрицу под сетку акустики?
			// если есть изменения в акустике не забывать менять там где идет расчет

			/*if (!points_rupture_ac.empty())
			//{
			//	for (it = points_rupture_ac.begin(); it != points_rupture_ac.end(); it++)// цикл по точкам где произошел разрыв
			//	{
			//		if ((*it).index_x == (Nx_acoustic / 2) && (*it).index_y == (Ny_acoustic / 2))
			//		{
			//			points_rupture_for_plot_z.push_back((*it));
			//		}

			//		if ((*it).index_y == (Ny_acoustic / 2) && (*it).index_z == 1)
			//		{
			//			points_rupture_for_plot_x.push_back((*it));
			//		}
			//	}
			//	sort(points_rupture_for_plot_z.begin(), points_rupture_for_plot_z.end(), comp_z);
			//	it_z = points_rupture_for_plot_z.begin();
			//	it_x = points_rupture_for_plot_x.begin();

			//	cout << " Points z for plotting " << endl;
			//	for (it = points_rupture_for_plot_z.begin(); it != points_rupture_for_plot_z.end(); it++)// цикл по точкам где произошел разрыв
			//	{
			//		cout << (*it).index_x << "   " << (*it).index_y << "   " << (*it).index_z << endl;
			//	}

			//	cout << endl << endl;
			//	cout << " Points x for plotting " << endl;

			//	for (it = points_rupture_for_plot_x.begin(); it != points_rupture_for_plot_x.end(); it++)// цикл по точкам где произошел разрыв
			//	{
			//		cout << (*it).index_x << "   " << (*it).index_y << "   " << (*it).index_z << endl;
			//	}

			//	cout << endl;
			//}*/

			//system("pause");

			namefile.clear();
			namefile = "Pacoustic (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(0, 1, Pfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "PTi (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(1, 1, Pfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "PTe (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(2, 1, Pfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "Pacoustic (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(3, 1, Pfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "PTi (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(4, 1, Pfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "PTe (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//pac.ShowDataOnPlot2D(5, 1, Pfs, namefile, true);
			//Sleep(1000);

			fout_Pac.close();
			fout_PTe.close();
			fout_PTi.close();
			fout_Pacx.close();
			fout_PTex.close();
			fout_PTix.close();
			points_rupture_for_plot_z.clear();
			points_rupture_for_plot_z.erase(points_rupture_for_plot_z.begin(), points_rupture_for_plot_z.end());
			points_rupture_for_plot_x.clear();
			points_rupture_for_plot_x.erase(points_rupture_for_plot_x.begin(), points_rupture_for_plot_x.end());
			fout_Pac.open(Pacteti[0]);
			fout_PTi.open(Pacteti[1]);
			fout_PTe.open(Pacteti[2]);
			fout_Pacx.open(Pacteti[3]);
			fout_PTix.open(Pacteti[4]);
			fout_PTex.open(Pacteti[5]);

			// открытие файлов произвожу выше

			namefile.clear();
			namefile = "Te,Ti (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(15, false, 2, TeTiTfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "Te,Ti (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(16, false, 2, TeTiTfs, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "P (x,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(18, false, 1, Pfs, namefile, true);
			//Sleep(1000);
			namefile.clear();
			namefile = "P (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(19, false, 1, Pfs, namefile, true);
			//Sleep(1000);

			namefile.clear();
			namefile = "a2z (z,mkm) " + metall + "," + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
			//plt.ShowDataOnPlot2D(28, false, 2, a2z_e_l, namefile, true);
			//Sleep(1000);

			current_time_frame += dt_step;
			it_fix_time_anim++;
			plt.Close_and_open_files_for_replot(for_close_and_open);
		}
		// next step
		n++;
		cout << " step = " << n << endl;

		// Perehod na sledujuschij shag po vremeni
		/*tmpe0 = tmpe1;
		tmpe1 = tmpe2;
		tmpe_ac0_new = tmpe_ac1_new;
		tmpe_ac1_new = tmpe_ac2_new;*/

		CopyDataArray3D(tmpe0, tmpe1, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(tmpe1, tmpe2, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(tmpe_ac0_new, tmpe_ac1_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(tmpe_ac1_new, tmpe_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		/*tmpi0 = tmpi1;
		tmpi1 = tmpi2;
		tmpi_ac0_new = tmpi_ac1_new;
		tmpi_ac1_new = tmpi_ac2_new;*/

		CopyDataArray3D(tmpi0, tmpi1, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(tmpi1, tmpi2, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(tmpi_ac0_new, tmpi_ac1_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(tmpi_ac1_new, tmpi_ac2_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		/*a1x = a2x;
		a1y = a2y;
		a1z = a2z;
		b1x = b2x;
		b1y = b2y;
		b1z = b2z;*/

		CopyDataArray3D(a1x, a2x, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(a1y, a2y, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(a1z, a2z, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(b1x, b2x, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(b1y, b2y, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(b1z, b2z, Nx_heat, Ny_heat, Nz_heat);

		CopyDataArray3D(a1x_ac_new, a2x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(a1y_ac_new, a2y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(a1z_ac_new, a2z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(b1x_ac_new, b2x_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(b1y_ac_new, b2y_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(b1z_ac_new, b2z_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		/*a1x_ac_new = a2x_ac_new;
		a1y_ac_new = a2y_ac_new;
		a1z_ac_new = a2z_ac_new;
		b1x_ac_new = b2x_ac_new;
		b1y_ac_new = b2y_ac_new;
		b1z_ac_new = b2z_ac_new;*/

		CopyDataArray3D(V00, V1, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(V1, V2, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(V00_ac_new, V1_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);
		CopyDataArray3D(V1_ac_new, V2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		//V00 = V1;
		//V1 = V2;

		CopyDataArray3D(e1, e2, Nx_heat, Ny_heat, Nz_heat);
		CopyDataArray3D(e1_ac_new, e2_ac_new, Nx_part_ac + 1, Ny_part_ac + 1, down_boundary_z_acoustic + 1);

		//e1 = e2;
		//e1_ac_new = e2_ac_new;
		it_time++;

		cout << "It took me %d clicks (%f seconds)." << endl <<
			(std::clock() - t) / (double)CLOCKS_PER_SEC << endl;
		sum_t += (std::clock() - t) / (double)CLOCKS_PER_SEC;
		cout << " General duration calculation = " << sum_t << endl;
		cout << endl;

		cout << endl;
	} while (tt <= (*it_last_time) /*8000*//*nn < 1000*/);

	t = clock() - t;

	cout << endl << endl;
	cout << "It took me %d clicks (%f seconds)." << endl <<
		(int)t << "   " << ((double)t) / CLOCKS_PER_SEC;
	cout << endl;

	// могут приодитьс (если вдруг моментов времени меньше заявленных линий на графике, и как результат мы не запишем данные из матрицы в файл )
	if (npxzt < count_of_lines)
	{
		npxzt = count_of_lines;
		plt.SetDataOnPlot3D(8, 0, e1, p0v_s, NULL, Ny_acoustic / 2, 1, count_of_lines, npxzt, NULL, empty, fun_on_x);
		plt.SetDataOnPlot3D(9, 0, e1, p0v_s, Nx_acoustic / 2, Ny_acoustic / 2, NULL, count_of_lines, npxzt, NULL, empty, fun_on_z);
	}

	//double*** Field_melting;
	//Field_melting = new double** [Nx];
	//for (int i = 0; i < Nx; i++)
	//{
	//	Field_melting[i] = new double* [Ny];
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		Field_melting[i][j] = new double[Nz];
	//	}
	//}

	//for (int i = 0; i < Nx; i++)
	//{
	//	for (int j = 0; j < Ny; j++)
	//	{
	//		for (int k = 0; k < Nz; k++)
	//		{
	//			Field_melting[i][j][k] = Melt_metal[mt].T_melting;// 1e-16;
	//		}
	//	}
	//}

	//if (current_number_line_melting < 15)
	//{
	//	current_number_line_melting = 11;
	//	plt.SetDataOnPlot3D(12, Nx, Ny, Nz, 1e+4 * dx * param.r0, 1e+4 * dy * param.r0, 1e+4 * dz / param.kabs, Field_melting, 1., NULL, Ny / 2, 1, 11, current_number_line_melting, NULL, empty, fun_on_x);
	//	plt.SetDataOnPlot3D(11, Nx, Ny, Nz, 1e+4 * dx * param.r0, 1e+4 * dy * param.r0, 1e+4 * dz / param.kabs, Field_melting, 1., Nx / 2, Ny / 2, NULL, 11, current_number_line_melting, NULL, empty, fun_on_z);
	//}

	//vector<string> Legenda = { "500 fs", "1000 fs", "1500 fs","2000 fs" , "3000 fs" , "4000 fs" , "5000 fs" , "6000 fs", "7000 fs", "8000 fs" };
	namefile.clear();
	namefile = "file P(x), bar " + metall + ", " + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
	//plt.ShowDataOnPlot2D(8, count_of_lines, Legenda, namefile, true);
	namefile.clear();
	namefile = "file P(z), bar " + metall + ", " + type_beam + ", pulse duration = " + ConvertNumToStringdouble(tp * 1e+15) + "fs, moment time = " + ConvertNumToStringdouble(tt) + "fs";
	//plt.ShowDataOnPlot2D(9, count_of_lines, Legenda, namefile, true);
	namefile.clear();
	//vector<string> Legenda_melting = { "2500 fs" , "3000 fs", "3500 fs", "4000 fs", "4500 fs", "5000 fs", "5500 fs", "6000 fs", "6500 fs", "7000 fs","T melting" };
	vector<string> Legenda_melting = Legenda_melting_tmp;
	Legenda_melting.push_back("T melting");
	namefile = "file Profile Teperature T(x)";
	//plt.ShowDataOnPlot2D(12, 11, Legenda_melting, namefile, true);
	namefile.clear();
	namefile = "file Profile Teperature T(z)";
	//plt.ShowDataOnPlot2D(11, 11, Legenda_melting, namefile, true);
	namefile.clear();

	sum_t = clock() - t;

	plt.Close_all_files_and_plots(number_plots);
}

inline int random_ab(int a, int b, std::default_random_engine & g)
{
	//return a + (b - a) * ((double)rand() / (double)RAND_MAX);
	std::uniform_int_distribution<int> distribution(a, b);
	return distribution(g);
}

int main()
{
	//CreateGif2DColor(3, 5, 290/*7640*/, 5, "Anim Ti (z,x) (28.09.2025)", "Ti(z,x)", "z,mkm", "x,mkm", false, NULL, NULL, NULL, NULL, "fs");
	//Sleep(1000);
	//CreateGif2DColor(4, 5, 290/*7640*/, 5, "Anim P (x,N div 2,z) (28.09.2025)", "P(x,z)", "z,mkm", "x,mkm", false, NULL, NULL, NULL, NULL, "fs");
	//Sleep(1000);
	//CreateGif2DColor(5, 5, 290/*7640*/, 5, "Anim P (x,y,1) (28.09.2025)", "P(x,y)", "y,mkm", "x,mkm", false, NULL, NULL, NULL, NULL, "fs");
	//Sleep(1000);
	//CreateGif2DColor(7, 5, 290/*7640*/, 5, "Anim Ti (y,x) (28.09.2025)", "Ti(y,x)", "y,mkm", "x,mkm", false, NULL, NULL, NULL, NULL, "fs");

	vector<string> legP = { "P" };
	vector<string> legT = { "Te", "Ti" };
	/*CreateGif_Y_x(15, 5, 290, 5, "Anim Te,Ti (x) (28.09.2025) Y(x)", "Te,Ti(x)", "x,mkm", "Te, Ti, K", legT, false, NULL, NULL, "fs");
	Sleep(1000);
	CreateGif_Y_x(16, 5, 290, 5, "Anim Te,Ti (z) (28.09.2025) Y(x)", "Te,Ti(z)", "z,mkm", "Te, Ti, K", legT, false, NULL, NULL, "fs");
	Sleep(1000);
	CreateGif_Y_x(18, 5, 290, 5, "Anim P (x) (28.09.2025) Y(x)", "P(x)", "x,mkm", "P, bar", legP, false, NULL, NULL, "fs");
	Sleep(1000);
	CreateGif_Y_x(19, 5, 290, 5, "Anim P (z) (28.09.2025) Y(x)", "P(z)", "z,mkm", "P, bar", legP, false, NULL, NULL, "fs");
	Sleep(1000);*/

	system("pause");

	//setlocale(LC_ALL, "rus");
	TypeBeam tp = Gauss;
	cout << " tp = " << tp << endl;
	cout << endl;

	int Nx_heat = 100;// 100;// dxy; // 100 //100; // chislo uzlov po x // 100
	int Ny_heat = 100;// 100;// dxy;// 100// 100; // chislo uzlov po y //100
	int Nz_heat = 250;// 250;//250; //!!!!!!!!!!!!!!

	int Nx_acoustic = 200;// 500;//200;// 200;//400; //dxy; // 100 //100; // chislo uzlov po x // 100
	int Ny_acoustic = 200;// 500;// 200;// 200;// 400;// dxy;// 100// 100; // chislo uzlov po y //100
	int Nz_acoustic = 2000;// 400;// 400;// 10000; // !!!!!!

	double*** V00, *** V1, *** V2, *** a1y, *** a2y, *** a1z, *** a2z, *** a2x,
		*** a1x, *** b1x, *** b2x, *** b1y, *** b2y, *** b1z, *** b2z, *** e2, *** e1;
	double*** tmpe0, *** tmpe1, *** tmpe2, *** tmpi0, *** tmpi1, *** tmpi2;
	//double*** Te_acoustic, *** Tei_old_xy_new_z_for_total_ac, *** Tei_old_y_new_xz_for_total_ac; // массивы для переопрелеления сетки температуры под акустику

	// 21 массива

	V00 = new double**[Nx_heat];
	V1 = new double**[Nx_heat];
	V2 = new double**[Nx_heat];
	a1y = new double**[Nx_heat]; // Euler coordinates
	a2y = new double**[Nx_heat]; // Euler coordinates
	a1z = new double**[Nx_heat]; // Euler coordinates
	a2z = new double**[Nx_heat]; // Euler coordinates
	a1x = new double**[Nx_heat]; // Euler coordinates
	a2x = new double**[Nx_heat]; // Euler coordinates
	b1x = new double**[Nx_heat]; // particle velocity
	b2x = new double**[Nx_heat]; // particle velocity
	b1y = new double**[Nx_heat]; // particle velocity
	b2y = new double**[Nx_heat]; // particle velocity
	b1z = new double**[Nx_heat]; // particle velocity
	b2z = new double**[Nx_heat]; // particle velocity
	e1 = new double**[Nx_heat]; // Pressure
	e2 = new double**[Nx_heat]; // Pressure

	tmpe0 = new double**[Nx_heat]; // temperature of electrons 
	tmpe1 = new double**[Nx_heat]; // temperature of electrons 
	tmpe2 = new double**[Nx_heat]; // temperature of electrons 
	tmpi0 = new double**[Nx_heat]; // temperature of ions
	tmpi1 = new double**[Nx_heat]; // temperature of ions
	tmpi2 = new double**[Nx_heat]; // temperature of ions

	//Te_acoustic = new double** [Nx_acoustic];
	//Tei_old_xy_new_z_for_total_ac = new double** [Nx_heat];
	//Tei_old_y_new_xz_for_total_ac = new double** [Nx_acoustic];
	//system("pause");

	for (int i = 0; i < Nx_acoustic; i++)
	{
		//Te_acoustic[i] = new double* [Ny_acoustic];
		for (int j = 0; j < Ny_acoustic; j++)
		{
			//	Te_acoustic[i][j] = new double[Nz_acoustic];
		}
	}

	for (int i = 0; i < Nx_heat; i++)
	{
		tmpe0[i] = new double*[Ny_heat];
		tmpe1[i] = new double*[Ny_heat];
		tmpe2[i] = new double*[Ny_heat];
		tmpi0[i] = new double*[Ny_heat];
		tmpi1[i] = new double*[Ny_heat];
		tmpi2[i] = new double*[Ny_heat];

		V00[i] = new double*[Ny_heat];
		V1[i] = new double*[Ny_heat];
		V2[i] = new double*[Ny_heat];
		a1y[i] = new double*[Ny_heat];
		a2y[i] = new double*[Ny_heat];
		a1z[i] = new double*[Ny_heat];
		a2z[i] = new double*[Ny_heat];
		a1x[i] = new double*[Ny_heat];
		a2x[i] = new double*[Ny_heat];
		b1x[i] = new double*[Ny_heat];
		b2x[i] = new double*[Ny_heat];
		b1y[i] = new double*[Ny_heat];
		b2y[i] = new double*[Ny_heat];
		b1z[i] = new double*[Ny_heat];
		b2z[i] = new double*[Ny_heat];
		e1[i] = new double*[Ny_heat];
		e2[i] = new double*[Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			tmpe0[i][j] = new double[Nz_heat];
			tmpe1[i][j] = new double[Nz_heat];
			tmpe2[i][j] = new double[Nz_heat];
			tmpi0[i][j] = new double[Nz_heat];
			tmpi1[i][j] = new double[Nz_heat];
			tmpi2[i][j] = new double[Nz_heat];

			V00[i][j] = new double[Nz_heat];
			V1[i][j] = new double[Nz_heat];
			V2[i][j] = new double[Nz_heat];
			a1y[i][j] = new double[Nz_heat];
			a2y[i][j] = new double[Nz_heat];
			a1z[i][j] = new double[Nz_heat];
			a2z[i][j] = new double[Nz_heat];
			a1x[i][j] = new double[Nz_heat];
			a2x[i][j] = new double[Nz_heat];
			b1x[i][j] = new double[Nz_heat];
			b2x[i][j] = new double[Nz_heat];
			b1y[i][j] = new double[Nz_heat];
			b2y[i][j] = new double[Nz_heat];
			b1z[i][j] = new double[Nz_heat];
			b2z[i][j] = new double[Nz_heat];
			e1[i][j] = new double[Nz_heat];
			e2[i][j] = new double[Nz_heat];
		}
	}

	for (int i = 0; i < Nx_heat; i++)
	{
		//Tei_old_xy_new_z_for_total_ac[i] = new double* [Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			//Tei_old_xy_new_z_for_total_ac[i][j] = new double[Nz_acoustic];
		}
	}

	for (int i = 0; i < Nx_acoustic; i++)
	{
		//Tei_old_y_new_xz_for_total_ac[i] = new double* [Ny_heat];
		for (int j = 0; j < Ny_heat; j++)
		{
			//Tei_old_y_new_xz_for_total_ac[i][j] = new double[Nz_acoustic];
		}
	}

	double TTT = 300.;
	/*  отрисовка графиков */

	Splayn spl_G_e_on_T, spl_C_e_on_T; //= Calculation_Interpolation("Ge_Au_new.txt");
	spl_C_e_on_T.InterpolateFast1D("Ce_Au_new.txt");
	spl_G_e_on_T.InterpolateFast1D("Ge_Au_new.txt");

	cout << "C_e_on_T = " << spl_C_e_on_T.GetY(TTT) / (100. * 100. * 100.) << "   " << 5.472e-2 << endl;
	cout << "G_e_on_T = " << spl_G_e_on_T.GetY(TTT) / (100. * 100. * 100.) << "   " << 2.5e+10 << endl;
	cout << " kte = " << Dependence_k_e_on_T(Au, TTT) / (100.) << "   " << 3.115 << endl;
	cout << " Cl = " << Dependence_C_l_on_T(Al, TTT) / (100. * 100. * 100.) << "   " << 2.550 << endl;

	ofstream fout1("Ce_interp.txt");
	ofstream fout2("Ge_interp.txt");
	ofstream fout3("ke_depend.txt");
	ofstream fout4("Cl_depend.txt");

	for (double i = 300.; i <= 50250.; i++)
	{
		fout1 << i << "   " << spl_C_e_on_T.GetY(i) << endl;
		fout2 << i << "   " << spl_G_e_on_T.GetY(i) << endl;
		fout3 << i << "   " << Dependence_k_e_on_T(Au, i) << endl;
	}

	/*for (double i = 300.; i < 3000.; i++)
	{
		fout4 << i << "   " << Dependence_C_l_on_T(Au, i) << endl;
	}*/
	fout4.precision(10);
	double d = 0.5;//121.;

	for (double i = 1338 - d - 2; i <= 1338. + d + 2; i += 0.5)
	{
		if (i < (1338. - d))
		{
			//////fout4 << i << "   " << (1.11 * 300. * (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) / (19.3 * 19.3 * pow(3.24e+5, 2) * 1e-6)) << endl;
			fout4 << i << "   " << (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) << endl;
		}

		if (i >= (1338. - d) && i <= (1338. + d))
		{
			///	//fout4 << i << "   " << (1.11* 300.* (((Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.) + (63.7 * 19.3) / (2 * 1.)))) / (pow((19.3 + 17.) / 2, 2) * pow((3.24e+5 + 2.567e+5) / 2, 2) * 1e-6)) << endl;
			fout4 << i << "   " << (((Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.) + (63.7 * 19.3) * delta_function(i, 1338., d)))) << endl;
			///	//fout4 << i << "   " << 0. << endl;
		}

		if (i > (1338. + d))
		{
			fout4 << i << "   " << (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) << endl;
			///////fout4 << i << "   " << (1.11 * 300.* (Dependence_C_l_on_T(Au, i) / (100. * 100. * 100.)) / (17. * 17. * pow(2.567e+5, 2) * 1e-6)) << endl;
		}
		//(Melt_metal[mt].gi * param.T00 * (Dependence_C_l_on_T(mt, param.T00 * spl_Ti.GetY(k * dz)) / (100. * 100. * 100.)) / (Melt_metal[mt].Density * Melt_metal[mt].Density * pow(Melt_metal[mt].u0, 2) * V0 * 1e-6))
	}

	//system("pause");
	Splayn spl_C_l_on_T;
	spl_C_l_on_T.InterpolateFast1D("Cl_depend.txt");
	fout4.close();
	fout4.open("Cl_depend.txt");

	for (double i = 1338 - d - 2.; i <= 1338. + d + 2.; i += 0.001)
		//for (double i = 1287.99; i <= 1288.05; i += 0.001)
	{
		fout4 << i << "   " << spl_C_l_on_T.GetY(i) << endl;
	}

	///////////////////////////////////////!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!?////////////////////////////////

	vector<string> vv = { "Ce_interp.txt", "Ge_interp.txt", "ke_depend.txt", "Cl_depend.txt", "Dependence Hmax on Fluence.txt" };
	GnuPlot plt(5, vv); // объект хранит 5 плотиков
	plt.SetParametrs2D(0, 1, 3, "C_e(T_e)", "T_e, K", "C_e, J/m^3K");
	plt.SetParametrs2D(1, 1, 3, "G_e(T_e)", "T_e, K", "G_e, W/m^3K");
	plt.SetParametrs2D(2, 1, 3, "k_e(T_e)", "T_e, K", "k_e, W/mK");
	plt.SetParametrs2D(3, 1, 3, "C_i(T_i)", "T_i, K", "C_i, J/m^3K");
	//plt.SetParametrs2D(4, 1, 3, "h max(Fluence)", "F, J/m^2", "z, nm");

	string  namefile;
	vector<string> Metal = { "1" };
	namefile = "Dependence Ce on Te";
	plt.ShowDataOnPlot2D(0, false, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence Ge on Te";
	plt.ShowDataOnPlot2D(1, false, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence ke on Te";
	plt.ShowDataOnPlot2D(2, false, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence Cl on Tl";
	plt.ShowDataOnPlot2D(3, false, 1, Metal, namefile, true);
	namefile.clear();
	namefile = "Dependence hmax on Fluence";
	//plt.ShowDataOnPlot2D(4, false, 1, Metal, namefile, true);
	namefile.clear();

	//fout1.close();
	//fout2.close();
	//fout3.close();
	fout4.close();

	//system("pause");

	vector<string> vvv = { "100 fs.txt", "1 ps.txt", "3 ps.txt", "5 ps.txt", "10 ps.txt" };
	GnuPlot pltt(vvv);
	pltt.SetParametrs2D(0, 5, 3, "h_max(F)", "F, J/cm^2", "h_max, nm");
	vector<string> Melting1 = { "100 fs", "1 ps", "3 ps", "5 ps", "10 ps" };
	namefile = "Dependence h_max on F";
	pltt.ShowDataOnPlot2D(0, true, 5, Melting1, namefile, true);
	namefile.clear();

	vector<Point3D> v;
	std::default_random_engine g(std::chrono::system_clock::now().time_since_epoch().count());
	for (int i = 0; i < 60; i++)
	{
		srand(static_cast<unsigned int>(time(nullptr)));
		int x = random_ab(47, 53, g);//47 + rand() % 53; // случайное число - формула
		srand(static_cast<unsigned int>(time(nullptr)));
		int y = random_ab(47, 53, g); ;// 47 + rand() % 53; // случайное число - формула
		srand(static_cast<unsigned int>(time(nullptr)));
		int z = random_ab(1, 6, g); ;// 1 + rand() % 6; // случайное число - формула

		Point3D tmp(x, y, z);
		//cout << tmp.index_x << "   " << tmp.index_y << "   " << tmp.index_z << endl;
		v.push_back(tmp);
	}

	if (v.empty())
	{
		cout << " empty " << endl;
	}
	else
	{
		cout << " not empty " << endl;
	}
	//cout << v.empty() << endl;
	//cout << endl << " size " << v.size() << endl;
	//Point3D tmp();
	vector<Point3D>::iterator iter;
	My_unique(v);
	//cout << endl << v.size() << endl;

	//cout << endl << endl << endl << endl;
	//system("pause");

	sort(v.begin(), v.end(), comp_x);
	//My_unique(v);
	for (iter = v.begin(); iter != v.end(); iter++)
	{
		//cout << (*iter).index_x << "   " << (*iter).index_y << "   " << (*iter).index_z << endl;
	}
	//cout << endl << endl;
	//cout << endl << v.size() << endl;

	//system("pause");

	MySort_Point3D_y(v);
	cout << endl;
	for (iter = v.begin(); iter != v.end(); iter++)
	{
		//cout << (*iter).index_x << "   " << (*iter).index_y << "   " << (*iter).index_z << endl;
	}
	//cout << endl << v.size() << endl;

	//system("pause");

	MySort_Point3D_z(v);
	//cout << endl;
	//cout << endl << v.size() << endl;
	for (iter = v.begin(); iter != v.end(); iter++)
	{
		//cout << (*iter).index_x << "   " << (*iter).index_y << "   " << (*iter).index_z << endl;
	}
	//cout << endl << v.size() << endl;

	cout << "the end " << endl;
	//system("pause");

	//vector<string> moments_time_string; // хранит моменты времени
	//map<string, string> momtime_n_t;// задает соответствие времени в фс с номером кадра в аниации
	//map<string, string>::iterator it_map;
	//it_map = momtime_n_t.begin();
	//int count_time_frame = 1;
	//for (double i = 100; i <= 100000; i += 100)
	//{
	//	//moments_time_string.push_back(ConvertNumToStringdouble(i) + "fs");
	//	string tmp_first = ConvertNumToStringdouble(i) + "fs";
	//	//(*it_map).second = tmp_second + "fs";
	//	string tmp_second = ConvertNumToString(count_time_frame);
	//	///cout << tmp_first << "   " << tmp_second << endl;
	//	//momtime_n_t[tmp_first] = tmp_second + "fs";
	//	momtime_n_t.insert(make_pair(tmp_first, tmp_second));
	//	cout << tmp_first << "   " << tmp_second << endl;
	//	//cout << (*it_map).first << "  =  " << (*it_map).second << endl;
	//	count_time_frame++;
	//	it_map++;
	//	//list_namefile.push_back()
	//}

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


	vector<int> moments_fix_time_anim;
	vector<string> Ti_x, P_x, P_X_Result;
	string time_begin_str, time_end_str;
	int time_begin = 2400;
	int time_end = 3000;
	time_begin_str = ConvertNumToString(time_begin);
	time_end_str = ConvertNumToString(time_end);
	for (int i = time_begin; i <= time_end; i++ /*= 10*/) // моменты времени для отрисовки и анимации
	{// аналог n
		moments_fix_time_anim.push_back(i);
		string str, number_frame_Gif;
		number_frame_Gif = ConvertNumToString(i);
		str = "Data_Gif23 + " + number_frame_Gif + ".txt";
		Ti_x.push_back(str);
		str.clear();
		str = "Data_Gif17 + " + number_frame_Gif + ".txt";
		// e(xy) z(index) = 1 + 1378
		//cout << str << endl;
		P_x.push_back(str);
	}

	///GnuPlot Animate2D(Ti_x);
	vector <string> Tix = { "0 nm", "2 nm", "4 nm", "6 nm", "8 nm", "10 nm" };
	namefile = "Animate Ti (x,mkm) Legend z from " + time_begin_str + " fs to " + time_end_str + " fs ";
	//Animate2D.CreateGifOnPlot2D(0, 6, 3, moments_fix_time_anim.size(), moments_fix_time_anim, "Ti", "x,mkm", "Ti (K)", Tix, namefile, true);
	Sleep(1000);
	namefile.clear();
	namefile = "Animate P (x,mkm) Legend z";
	//SelectDataFromPlotColorForPlot2D(P_x, P_X_Result, )
	///GnuPlot Animate2D_P(P_x);
	//Animate2D_P.CreateGifOnPlot2D(0, 2, 3, moments_fix_time_anim.size(), moments_fix_time_anim, "P", "x,mkm", "P (bar)", Tix, namefile, true);

	moments_fix_time_anim.clear();
	P_x.clear();
	Ti_x.clear();

	string fixed_z_index_str, fixed_z_depth_str;
	int fixed_z_index = 5; // меняем номер узла сетки по z
	int dz_depth = 2;
	time_begin = 2400;
	time_end = 3000;
	for (int i = time_begin; i <= time_end; i++ /*= 10*/) // моменты времени для отрисовки и анимации
	{// аналог n
		moments_fix_time_anim.push_back(i);
		string str, number_frame_Gif;
		number_frame_Gif = ConvertNumToString(i);
		fixed_z_index_str = ConvertNumToString(fixed_z_index);
		//str = "e(xy) z(index) = 1 + " + number_frame_Gif + ".txt";
		str = "e(xy) z(index) = " + fixed_z_index_str + " + " + number_frame_Gif + ".txt";
		P_x.push_back(str);
		str.clear();
	}

	namefile.clear();
	time_begin_str = ConvertNumToString(time_begin);
	time_end_str = ConvertNumToString(time_end);
	fixed_z_depth_str = ConvertNumToString(fixed_z_index * dz_depth);
	namefile = "Animate P (x) in fix z = " + fixed_z_depth_str + " nm, fix y = 500 mkm from " + time_begin_str + " fs to " + time_end_str + " fs ";
	vector <string> Pxz = { fixed_z_depth_str + " nm" };
	//SelectDataFromPlotColorForPlot2D(P_x, P_X_Result, 0., 1000., 500.);
	//GnuPlot Animate2D3(P_X_Result);
	//Animate2D3.CreateGifOnPlot2D(0, 1, 3, moments_fix_time_anim.size(), moments_fix_time_anim, "P", "x,mkm", "P (bar)", Pxz, namefile, true);
	namefile.clear();
	namefile = "Animate P (x,y) in fix z = " + fixed_z_depth_str + " nm from " + time_begin_str + " fs to " + time_end_str + " fs ";
	//GnuPlot Animate2D2(P_x);
	//Animate2D2.CreateGifOnPlotColor(0, moments_fix_time_anim.size(), moments_fix_time_anim, "P (x,y) ", "X,mkm", "P, (bar)", 1000., 1000., namefile);



	///////////////////////////////////////////////////////////////////////////////////

	//vector<string> symbols;
	//vector<Point3D> pt;
	//ifstream fin("PointsAc4783.txt");
	//if (fin.is_open())
	//{
	//	double x, y, z;
	//	char symbol1, symbol2, symbol3, symbol4;
	//	while (fin >> symbol1 >> x >> symbol2 >> y >> symbol3 >> z >> symbol4)
	//	{
	//		Point3D tmp(x, y, z);
	//		pt.push_back(tmp);
	//		//cout << symbol1 <<"   "<< x << "   " << y << "   " << z << endl;
	//		/*	new_points.push_back(Point{ x, y });*/
	//		//X.push_back(x);
	//		//Y1.push_back(y1);
	//		//cout << x << "   " << y << endl;
	//	}
	//	sort(pt.begin(), pt.end(), comp_x);
	//	cout << (*(pt.begin())).index_x << "   " << (*(pt.begin())).index_y << "   " << (*(pt.begin())).index_z << endl;
	//	cout << (*(pt.end() - 1)).index_x << "   " << (*(pt.end() - 1)).index_y << "   " << (*(pt.end() - 1)).index_z << endl;
	//	sort(pt.begin(), pt.end(), comp_y);
	//	cout << (*(pt.begin())).index_x << "   " << (*(pt.begin())).index_y << "   " << (*(pt.begin())).index_z << endl;
	//	cout << (*(pt.end() - 1)).index_x << "   " << (*(pt.end() - 1)).index_y << "   " << (*(pt.end() - 1)).index_z << endl;
	//	sort(pt.begin(), pt.end(), comp_z);
	//	cout << (*(pt.end() - 1)).index_x << "   " << (*(pt.end() - 1)).index_y << "   " << (*(pt.end() - 1)).index_z << endl;
	//}
	//fin.close();


	//vector<string> P_x = { "Data_Gif18 + 4.txt", "Data_Gif18 + 9.txt", "Data_Gif18 + 19.txt", "Data_Gif18 + 49.txt", "Data_Gif18 + 59.txt",
	//"Data_Gif18 + 69.txt", "Data_Gif18 + 79.txt", "Data_Gif18 + 119.txt", "Data_Gif18 + 159.txt", "Data_Gif18 + 199.txt" };
	//vector<string> P_z_result;

	//string namefile1;
	//GnuPlot plt_P_X(P_x);
	//plt_P_X.SetParametrs2D(0, 10, 3, "P", "x, mkm", "P (bar)");
	//vector<string> plt_P_X_legend = { "50 fs", "100 fs", "200 fs", "500 fs", "600 fs", "700 fs", "800 fs",  "1200 fs", "1600 fs", "2000 fs" };
	//namefile1 = "Anim P(x) in dif mom time";
	//plt_P_X.ShowDataOnPlot2D(0, 0, 10, plt_P_X_legend, namefile1, true);
	//namefile1.clear();

	//vector<string> P_z = { "Data_Gif19 + 4.txt", "Data_Gif19 + 9.txt", "Data_Gif19 + 19.txt", "Data_Gif19 + 29.txt", "Data_Gif19 + 39.txt",
	//"Data_Gif19 + 59.txt", "Data_Gif19 + 79.txt", "Data_Gif19 + 119.txt", "Data_Gif19 + 159.txt", "Data_Gif19 + 199.txt" };
	//SelectDataForPlot2D(P_z, P_z_result, 0., 0.14);
	//GnuPlot plt_P_Z(P_z_result);
	//plt_P_Z.SetParametrs2D(0, 10, 3, "P", "z, mkm", "P (bar)");
	////"100 fs", "200 fs", "300 fs", "400 fs", "500 fs", "600 fs", "700 fs", "800 fs", "900 fs", "1000 fs"
	//vector<string> plt_P_Z_legend = { "50 fs", "100 fs", "200 fs", "300 fs", "400 fs", "600 fs", "800 fs", "1200 fs", "1600 fs", "2000 fs" };
	//namefile1 = "Anim P(z) in dif mom time";
	//plt_P_Z.ShowDataOnPlot2D(0, 0, 10, plt_P_Z_legend, namefile1, true);
	//namefile1.clear();


	vector<double> data1;
	vector<double> data2;
	vector<double> data3;
	vector<double> ::iterator it_db1, it_db2, it_db3;

	//ifstream fin_data("Data_Gif28 + 670.txt");
	//if (fin_data.is_open())
	//{
	//	double x, y, z;
	//	while (fin_data >> x >> y >> z)
	//	{
	//		//cout << x << "   " << y << "   " << z << endl;
	//		data1.push_back(x);
	//		data2.push_back(y);
	//		data3.push_back(z);
	//	}
	//}

	//it_db3 = data3.begin();
	for (it_db3 = data3.begin(); it_db3 != data3.end(); it_db3++)
	{
		(*it_db3) *= (1e+4);
		//(*it_db3) /= (1e+4);
		//(*it_db1) *= pow(5e+5, -1);
	}

	//ofstream fout_data("Data_draw.txt");
	//cout << endl;
	//it_db1 = data1.begin();
	//it_db3 = data3.begin();
	//for (it_db2 = data2.begin(); it_db2 != data2.end(); it_db2++)
	//{
	//	fout_data << (*it_db1) << "   " << (*it_db2) << "   " << (*it_db3) << endl;
	//	//	cout << (*it_db1) << "   " << (*it_db2) << "   " << (*it_db3) << endl;
	//	it_db1++;
	//	it_db3++;
	//}

	cout << endl;

	vector<string> vv_data = { "Data_draw.txt" };
	///GnuPlot plt_a(vv_data.size(), vv_data);
	///plt_a.SetParametrs2D(0, 2, 3, "a2z(z)", "z, mkm", "a2z, mkm");
	//plt.SetParametrs2D(4, 1, 3, "h max(Fluence)", "F, J/m^2", "z, nm");

	//string  namefile;
	namefile.clear();
	vector<string> Metal_a = { "Lagrange", "Euler" };
	namefile = "Dependence a2z on z";
	//plt_a.ShowDataOnPlot2D(0, false, 2, Metal_a, namefile, true);
	namefile.clear();

	////////////////////////////////////////////////////////////////////////////////

	//system("pause");

	//GnuPlot plt_experim(24,40);
	cout << endl << endl << "Gnuplot" << endl;
	system("pause");
	cout << " STOP " << endl;
	system("pause");
	//type metal, type beam, kte, Ce,Ci,	   gamma,    ge,   gi, us	tp s	
																								//I0 = 1e+7 0.15e+12
	// double*** V2, double*** a1y, double*** a2y, double*** a1z, double*** a2z, double*** a2x, double*** a1x, double*** b1x, double*** b2x, double*** b1y, double*** b2y, double*** b1z, double*** b2z, double*** e2, double*** e1, double*** F, double*** tmpe0, double*** tmpe1, double*** tmpe2, double*** tmpi0, double*** tmpi1, double*** tmpi2, double** Pxt, double** Pzt, int Nx, int Ny, int Nz, Splayn spl_C_l_on_T)
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, 0.75e+11, V00, V1, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, Nx_heat, Ny_heat, Nz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic,
		spl_C_l_on_T,/* Te_acoustic, Tei_old_xy_new_z_for_total_ac, Tei_old_y_new_xz_for_total_ac,*/ "Depth1.txt");
	system("pause");
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-13, 1e+7, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, tmpe0, tmpe1, tmpe2, tmpi0, tmpi1, tmpi2, Nx_heat, Ny_heat, Nz_heat, Nx_acoustic, Ny_acoustic, Nz_acoustic, spl_C_l_on_T, Te_acoustic, Ti_acoustic);
	system("pause");

	/*MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 5e-14, 0.2e+8, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 2e-13, 0.5e+7, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 5e-13, 0.2e+7, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	*///MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Gauss, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	system("pause");

	//type metal, type beam, kte, ro0, 	Ce,			Ci,	   gamma,  ge, gi,		us	  tp s	
//	MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Au, Vortex, 3.115, 19.32, 5.472e-2, 2.550, 2.5e+10, 1.5, 1.11, 3.24e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);

	//system("pause");
	//type metal, type beam, kte,   Ce,			Ci,	   gamma,  ge,  gi,     us	   tp s	
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Gauss, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);

	//type metal, type beam,  kte,   Ce,			Ci,	   gamma,  ge,  gi,     us	   tp s		
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	//MainProcedure(Ni, Vortex, 0.900, 8.96, 2.733e-2, 3.910, 4.05e+12, 1.5, 1.83, 5.63e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);	// продублировать по разным длительностям импульса
	system("pause");


	//type metal, type beam, kte,   Ce,	???		Ci,	   gamma,  ge,  gi,     us	   tp s	
	/*MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Gauss, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);

	//type metal, type beam,  kte,   Ce,??			Ci,	   gamma,  ge,  gi,     us	   tp s
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 5e-14, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 1e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);
	MainProcedure(Cu, Vortex, 3.930, 8.92, 2.733e-2, 3.568, 1.0e+11, 1.5, 1.96, 4.73e+5, 2e-13, V2, a1y, a2y, a1z, a2z, a2x, a1x, b1x, b2x, b1y, b2y, b1z, b2z, e2, e1, F, tmpe0, tmpe1, tmpe2, tmpi1, tmpi2, Pxt, Pzt, Nx, Ny, Nz);	// продублировать по разным длительностям импульса
	*/
	system("pause");
	return 0;
}



