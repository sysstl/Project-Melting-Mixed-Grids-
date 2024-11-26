#pragma once
#include <vector>
#include <string>
#include <fstream>
#include <algorithm>
#include <iostream>
#include "Matrix.h"

using namespace std;

struct DataforColorPlot
{
	DataforColorPlot(double index_x, double index_y, double index_z) : index_x{ index_x }, index_y{ index_y }, index_z{ index_z } {}
	double index_x;
	double index_y;
	double index_z;
};

enum Plane { xy, xz, yz };
enum Profile { fun_on_x, fun_on_y, fun_on_z, fun_on_t };

class GnuPlot {

	int number_of_plots;
	int count_frame = 0;

	vector<FILE*> gnuplotPipe; //���� � ������� ������������ ������� ��� ��������
	vector<ofstream> file; // ������ ��� ������ ������ � ���� (����������� ��������)
	//vector<ifstream> file_in;
	vector<string> filename; // ���� ������� �������� ������ ��� ���������  // ����� ������ ��� ������ ������ � ���� (����������� ��������)

	vector<FILE*> gnuplotPipe_Gif;
	//vector<vector<ofstream>> file_Gif; // ������ ��� ������ ������ � ���� (��������)
	vector<ofstream> file_Gif; // ������ ��� ������ ������ � ���� (��������)
	vector<string> filename_Gif;  // ����� ������ ��� ������ ������ � ���� (�������)
	vector<Matrix> P_xyzt_Gif;

	vector<Matrix> P_xyzt;

	string GetGPPath();
	string GetGPTerminal();
	size_t GetTermnalWidth(); //������ ������ 

public:
	//  �������� ������� ������� � ������� ������� � �������� ���� (����������: fun � SetDataOnPlotColor3D)!!!!!

	// ����������� ��� �������� �������� � ������� ����� ���������� ������ ���������� ����� � ������� � ������� � ��������
	// ������������ ������ ��� ����������� ��������
	GnuPlot(int number_of_plots_); // ��� ������ ������� �� ���������� ������ // ++++++++++++++++++

	// ������������ ��� ����������� �������� � ��������
	GnuPlot(int number_of_plots_, int count_frame);

	// ����������� ��� �������� �������� � ������� ����� ���������� ����������������� ������ ������� � ������ ����������
	GnuPlot(int number_of_plots_, vector<string> filename_);
	// ����������� ��� �������� 1 �����, ������� ������� ������ � ������ ������ (+ ������ ����� ������ �����)
	// ����������� ��� �������� 1 �����, ������� ������� �������� (gif) �� ������ � ������ ������
	GnuPlot(vector<string> filename_);

	// ���������� ����� ��� ����������� �������(��) y(x)
	void SetParametrs2D(int Current_number_plot, int number_of_lines, int width_line, string title_plot, string xlabel, string ylabel); // +++ ��������� � ����� �����  // Current_number_plot=0,1,2, ....
	// ���������� ����� ��� ����������� ������ �������� �������� (������)
	void SetParametrsOnPlotColor(int Current_number_plot, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y); // Current_number_plot=0,1,2, ....

	// ��� ������� ����� ������� ������ 1 ��� ����� SetDataOnPlot2D (��� 2-� �-��� ���� � ����� � SetDataOnPlot1D(2D, 3D))
	// ���������� ���� 2-�� �-���: �-��� SetDataOnPlot3D ���������� � ����� �� ������� + ��������� ��������� ������ ��� ������ � ����
	void SetGridOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, int number_of_lines, Profile prof);
	// ��� ������� ����� ������� ������ 1 ��� ����� SetDataOnPlot3D (��� 2-� �-��� ���� � ����� � SetDataOnPlot1D(2D, 3D))
	// ���������� ���� 2-�� �-���: �-��� SetDataOnPlot3D ���������� � ����� �� ������� + ��������� ��������� ������ ��� ������ � ����
	void SetGridOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, int number_of_lines, Profile prof);

	// 1D, 2D, 3D - ����������� ������
   ///(����������� �������� ������) y(x) (Ce �� Te) (�������� ��� ��������� ������� �� ������������ ����� ������� ��� � ������(�� ���������� � �����))
   // ��� ��������� ������� ��� 1D - ������ ���� (U(x), U(t))  
	void SetDataOnPlot1D(int Current_number_plot, int Nx, double dx, double** fun_dimensionless, double parametr_for_dimension, int moment_of_time, double dt, Profile prof);
	// ��� ��������� ������� ������� (y(x)) ��� 2D - ������ ���� (U(x), U(y), U(t))  
	void SetDataOnPlot2D(int Current_number_plot, double** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double**> vec, Profile prof);// ������� ������ ���� �������, ������ ������������� � �������� y(x), � �� ���������
	// ��� ��������� ������� ������� (y(x) ��� ����� t) ��� 3D - ������ ���� (U(x), U(y), U(z),U(t)) 
	void SetDataOnPlot3D(int Current_number_plot, int current_count_frame, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int fixed_point_on_axis_z, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double***> vec, Profile prof); // linetype � ���������� �� ������ � ������ (�� ������ ����� ����)

	// ��� ��������� ���������� ������� (���� ����������� �������� 2D)
	void SetDataOnPlotColor2D(int Current_number_plot, int current_count_frame, int Nx, int Ny, double dx, double dy, double** fun, double parametr_for_dimension);
	// ��� ��������� ���������� ������� (������� ����������� � ������� (���������) 3-� ������� ������� - ����)
	void SetDataOnPlotColor3D(int Current_number_plot, int current_count_frame, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis, Plane plane);
	void SetDataOnPlotColor3D_For_Mixed_Grids(int Current_number_plot, int current_count_frame, int Nx_heat, int Ny_heat, int Nz_heat, double dx_heat, double dy_heat, double dz_heat, int Nx_ac, int Ny_ac, int Nz_ac, double dx_ac, double dy_ac, double dz_ac, int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat,
		double*** fun_dimensionless_heat, double*** fun_dimensionless_ac, double parametr_for_dimension, int fixed_point_on_axis_heat, int fixed_point_on_axis_ac, Plane plane);

	// ��������� ������ ������� (���� ���� ����� ����� �� �������) (� ������ ����� �������� 1 ���� ������)
	//void ShowDataOnPlot2D(int Current_number_plot, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_);
	// ��������� ������ ������� (���� ���� ����� ����� �� �������) (� ������ ����� ��������� ��������� ������ � �������)
	void ShowDataOnPlot2D(int Current_number_plot, bool DataFromDiffFiles, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_);
	void CreateGifOnPlot2D(int Current_number_plot, int number_of_lines, int width_line, int count_frame, vector<int> moments_fix_time_anim, string title_plot, string xlabel, string ylabel, vector<string> list_name_line, string name_of_file, bool gif_);
	void CreateGifOnPlotColor(int Current_number_plot, int count_frame, vector<int> moments_fix_time_anim, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y, string name_of_file);

	// ��������� ������ �������� �������� (������)
	void ShowDataOnPlotColor(int Current_number_plot, string name_of_file, bool png_);
	// ����� ���� ��������� ���� �� ������� ������� (��� ��������� �������� �������� ������ �������)
	void Close_and_open_files_for_replot(vector<int> Current_number_of_plots_); //Current_number_of_plots_ - ������ ������� ���������� ������ ������ ������� ����� ������� � ����� �������
	void Close_all_files_and_plots(int number_of_plots_); // ++++++++++++++++++

};