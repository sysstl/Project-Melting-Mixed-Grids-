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

	vector<FILE*> gnuplotPipe; //файл в который записываются команды для гнуплота
	vector<ofstream> file; // потоки для вывода данных в файл (статические картинки)
	//vector<ifstream> file_in;
	vector<string> filename; // файл который содержит данные для отрисовки  // имена файлов для вывода данных в файл (статические картинки)

	vector<FILE*> gnuplotPipe_Gif;
	//vector<vector<ofstream>> file_Gif; // потоки для вывода данных в файл (Анимация)
	vector<ofstream> file_Gif; // потоки для вывода данных в файл (Анимация)
	vector<string> filename_Gif;  // имена файлов для вывода данных в файл (Анимция)
	vector<Matrix> P_xyzt_Gif;

	vector<Matrix> P_xyzt;

	string GetGPPath();
	string GetGPTerminal();
	size_t GetTermnalWidth(); //ширина экрана 

public:
	//  передаем размеры области и зачение функций в размером виде (иселючение: fun в SetDataOnPlotColor3D)!!!!!

	// конструктор для создания плотиков в которых будем отображать данные полученные здесь в проекте и лежащие в массивах
	// используется только для статических картинок
	GnuPlot(int number_of_plots_); // наш объект состоит из нескольких плотов // ++++++++++++++++++

	// используется для статических картинок и анимации
	GnuPlot(int number_of_plots_, int count_frame);

	// конструктор для создания плотиков в которых будем отображать экспериментальные данные лежащие в файлах изначально
	GnuPlot(int number_of_plots_, vector<string> filename_);
	// конструктор для создания 1 плота, который покажет данные с разных файлов (+ данные имеют разные сетки)
	// конструктор для создания 1 плота, который покажет анимацию (gif) по данным с разных файлов
	GnuPlot(vector<string> filename_);

	// Подготовка плота для отображения графика(ов) y(x)
	void SetParametrs2D(int Current_number_plot, int number_of_lines, int width_line, string title_plot, string xlabel, string ylabel); // +++ проверить а риеал даных  // Current_number_plot=0,1,2, ....
	// Подготовка плота для отображения данных цветовой палитрой (радуга)
	void SetParametrsOnPlotColor(int Current_number_plot, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y); // Current_number_plot=0,1,2, ....

	// эту функцию нужно вызвать только 1 раз перед SetDataOnPlot2D (эти 2-у ф-ции идут в свзке с SetDataOnPlot1D(2D, 3D))
	// потребость этих 2-ух ф-ций: ф-ция SetDataOnPlot3D вызывается в цикле по времени + сохранить структуру данных при записи в файл
	void SetGridOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, int number_of_lines, Profile prof);
	// эту функцию нужно вызвать только 1 раз перед SetDataOnPlot3D (эти 2-у ф-ции идут в свзке с SetDataOnPlot1D(2D, 3D))
	// потребость этих 2-ух ф-ций: ф-ция SetDataOnPlot3D вызывается в цикле по времени + сохранить структуру данных при записи в файл
	void SetGridOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, int number_of_lines, Profile prof);

	// 1D, 2D, 3D - размерность задачи
   ///(зависимости эксперим данных) y(x) (Ce от Te) (подумать над методамми которые бы отрисовывали даные лежащие уже в файлах(не полученные в проге))
   // Для получения графика для 1D - задача тепл (U(x), U(t))  
	void SetDataOnPlot1D(int Current_number_plot, int Nx, double dx, double** fun_dimensionless, double parametr_for_dimension, int moment_of_time, double dt, Profile prof);
	// Для получения графика функции (y(x)) для 2D - задача тепл (U(x), U(y), U(t))  
	void SetDataOnPlot2D(int Current_number_plot, double** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double**> vec, Profile prof);// аналоги нижних двух функций, только пприменительо к графикам y(x), а не красочным
	// Для получения графика функции (y(x) при серии t) для 3D - задача тепл (U(x), U(y), U(z),U(t)) 
	void SetDataOnPlot3D(int Current_number_plot, int current_count_frame, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int fixed_point_on_axis_z, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double***> vec, Profile prof); // linetype и поглдывать на рафики с задачи (от времеи могут быть)

	// Для получения красочного рисунка (поле температуры пластины 2D)
	void SetDataOnPlotColor2D(int Current_number_plot, int current_count_frame, int Nx, int Ny, double dx, double dy, double** fun, double parametr_for_dimension);
	// Для получения красочного рисунка (профиль температуры в сечении (плоскости) 3-х мерного объекта - куба)
	void SetDataOnPlotColor3D(int Current_number_plot, int current_count_frame, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis, Plane plane);
	void SetDataOnPlotColor3D_For_Mixed_Grids(int Current_number_plot, int current_count_frame, int Nx_heat, int Ny_heat, int Nz_heat, double dx_heat, double dy_heat, double dz_heat, int Nx_ac, int Ny_ac, int Nz_ac, double dx_ac, double dy_ac, double dz_ac, int left_boundary_x_heat, int right_boundary_x_heat, int left_boundary_y_heat, int right_boundary_y_heat, int down_boundary_z_heat,
		double*** fun_dimensionless_heat, double*** fun_dimensionless_ac, double parametr_for_dimension, int fixed_point_on_axis_heat, int fixed_point_on_axis_ac, Plane plane);

	// Отриосвка данных линиями (если есть много линий то легенда) (к одному плоту привязан 1 файл данных)
	//void ShowDataOnPlot2D(int Current_number_plot, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_);
	// Отриосвка данных линиями (если есть много линий то легенда) (к одному плоту привязано несколько файлов с данными)
	void ShowDataOnPlot2D(int Current_number_plot, bool DataFromDiffFiles, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_);
	void CreateGifOnPlot2D(int Current_number_plot, int number_of_lines, int width_line, int count_frame, vector<int> moments_fix_time_anim, string title_plot, string xlabel, string ylabel, vector<string> list_name_line, string name_of_file, bool gif_);
	void CreateGifOnPlotColor(int Current_number_plot, int count_frame, vector<int> moments_fix_time_anim, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y, string name_of_file);

	// Отриосвка данных цветовой палитрой (радуга)
	void ShowDataOnPlotColor(int Current_number_plot, string name_of_file, bool png_);
	// нужно если отрисовка идет во времени кадрами (при повтонром открытии удалются старые даннные)
	void Close_and_open_files_for_replot(vector<int> Current_number_of_plots_); //Current_number_of_plots_ - вектор хранщий конкретные номера плотов которые нужно закрыть и сразу открыть
	void Close_all_files_and_plots(int number_of_plots_); // ++++++++++++++++++

};