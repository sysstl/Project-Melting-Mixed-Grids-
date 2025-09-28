#include "GnuPlot.h"
#include <iostream>
#include <iomanip>
#include <string>
#include <cmath>
#include <vector>

string ConvertNumToStringint(int i)
{
	string tmp_num_line;

	int length = snprintf(NULL, 0, "%i", i + 1);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%i", i + 1);
	for (int j = 0; j < length; j++)
	{
		tmp_num_line.push_back(str[j]); // номер столбца
	}

	return tmp_num_line;
}

string ConvertNumToStringint_1(int i)
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

string ConvertNumToStringdouble_gnuplot(double i)
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

using namespace std;
// ORIGINAL
GnuPlot::GnuPlot(int number_of_plots_)
{
	number_of_plots = number_of_plots_;
	Matrix matr(0, 0);

	string text;
	FILE* gnuplotPipe_tmp;
	FILE* gnuplotPipe_tmp_Gif;
	for (int i = 0; i < number_of_plots; ++i)
	{
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			text.push_back(str[j]);
		}

		string file_name = "gnuplot" + text + ".txt";
		text.clear();
		filename.push_back(file_name);
		file.emplace_back(ofstream{ file_name }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);
		gnuplotPipe_tmp = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
		gnuplotPipe.push_back(gnuplotPipe_tmp); // вектор файлов
		gnuplotPipe[i] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом  (файл системный в который записываются команды для гнуплота)
		P_xyzt.push_back(matr);

		gnuplotPipe_tmp_Gif = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
		gnuplotPipe_Gif.push_back(gnuplotPipe_tmp_Gif); // вектор файлов
		gnuplotPipe_Gif[i] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом  (файл системный в который записываются команды для гнуплота)
		P_xyzt_Gif.push_back(matr);
	}
}

GnuPlot::GnuPlot(int number_of_plots_, int count_frame_)
{
	number_of_plots = number_of_plots_;
	count_frame = count_frame_;
	//file_Gif.reserve(number_of_plots);
	Matrix matr(0, 0);

	string text_Gif;
	string text;
	FILE* gnuplotPipe_tmp_Gif;
	FILE* gnuplotPipe_tmp;
	for (int i = 0; i < number_of_plots; ++i)// цикл для создания системных файлов и файлов для данных для статических картинок
	{
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			text.push_back(str[j]);
		}
		string file_name = "gnuplot" + text + ".txt";
		text.clear();
		filename.push_back(file_name);
		file.emplace_back(ofstream{ file_name }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);

		gnuplotPipe_tmp = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
		gnuplotPipe.push_back(gnuplotPipe_tmp); // вектор файлов
		gnuplotPipe[i] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом  (файл системный в который записываются команды для гнуплота)
		P_xyzt.push_back(matr);

		//!!!!!
		gnuplotPipe_tmp_Gif = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
		gnuplotPipe_Gif.push_back(gnuplotPipe_tmp_Gif); // вектор файлов
		gnuplotPipe_Gif[i] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом  (файл системный в который записываются команды для гнуплота)
		P_xyzt_Gif.push_back(matr);
	}

	//if (count_frame != 0)
	//{
	//	for (int i = 0; i < number_of_plots; ++i)// цикл для создания системных файлов и файлов для данных для Анимации
	//	{
	//		string type_plot_Gif;
	//		string file_name;
	//		int length = snprintf(NULL, 0, "%i", i);
	//		char* str = new char[length + 1];
	//		snprintf(str, length + 1, "%i", i);
	//		for (int j = 0; j < length; j++)
	//		{
	//			type_plot_Gif.push_back(str[j]);
	//		}

	//		int j = 0;
	//		//for (int j = 0; j < count_frame_; ++j)
	//		{
	//			string number_frame_Gif /*= "GIF"*/;

	//			int length = snprintf(NULL, 0, "%i", j);
	//			char* str = new char[length + 1];
	//			snprintf(str, length + 1, "%i", j);
	//			for (int k = 0; k < length; k++)
	//			{
	//				number_frame_Gif.push_back(str[k]);
	//			}

	//			file_name = "Data_Gif" + type_plot_Gif + " + " + number_frame_Gif + ".txt";
	//			number_frame_Gif.clear();
	//			//filename_Gif.push_back(file_name);

	//			//v1[j] = ofstream{ file_name };
	//			//v1.emplace_back(ofstream{ file_name });
	//			///file_Gif.emplace_back(ofstream{ file_name }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);
	//			//file_Gif.push_back(ofstream{ file_name }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);
	//		}

	//		//file_Gif.emplace_back(ofstream{ type_plot_Gif }); // создаие потока вывода в файл и самого файла дл записи туда данных и ео открытие file.open(filename);
	//		type_plot_Gif.clear();

			//gnuplotPipe_tmp_Gif = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
			//gnuplotPipe_Gif.push_back(gnuplotPipe_tmp_Gif); // вектор файлов
			//gnuplotPipe_Gif[i] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом  (файл системный в который записываются команды для гнуплота)
			//P_xyzt_Gif.push_back(matr);
	//	}
	//}
}

GnuPlot::GnuPlot(int number_of_plots_, vector<string> filename_)
{
	number_of_plots = number_of_plots_;
	FILE* gnuplotPipe_tmp;
	for (int i = 0; i < number_of_plots; ++i)
	{
		gnuplotPipe_tmp = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
		gnuplotPipe.push_back(gnuplotPipe_tmp); // вектор файлов
		gnuplotPipe[i] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом (файл системный в который записываются команды для гнуплота)
		filename = filename_;
	}
}

GnuPlot::GnuPlot(vector<string> filename_)
{
	number_of_plots = 1;
	FILE* gnuplotPipe_tmp;
	gnuplotPipe_tmp = new FILE; // создание файла дл гнуплота чтобы записывать туда комады дл гнуплота
	gnuplotPipe.push_back(gnuplotPipe_tmp); // вектор файлов
	gnuplotPipe[0] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом (файл системный в который записываются команды для гнуплота)
	filename = filename_;

	gnuplotPipe_Gif.push_back(gnuplotPipe_tmp); // вектор файлов
	gnuplotPipe_Gif[0] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот и свзываем его с его же файлом (файл системный в который записываются команды для гнуплота)
	filename = filename_;
}

void GnuPlot::SetParametrs2D(int Current_number_plot, int number_of_lines, int width_line, string title_plot, string xlabel, string ylabel) // Current_number_plot=0,1,2, ....
{
	//gnuplotPipe[Current_number_plot] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот
	//fprintf(gnuplotPipe[Current_number_plot], "set xrange [%lf:%lf]\n", 0., 50.); 
	//fprintf(gnuplotPipe[Current_number_plot], "set yrange [%lf:%lf]\n", 0., 50.);

	string tmp_width_line;
	int length = snprintf(NULL, 0, "%i", width_line);
	char* str = new char[length + 1];
	snprintf(str, length + 1, "%i", width_line);
	for (int j = 0; j < length; j++)
	{
		tmp_width_line.push_back(str[j]);
	}
	string color;

	for (int i = 1; i <= number_of_lines; i++) // number_of_lines - число линий на плоте, 
	{
		string tmp_num_line;
		int length = snprintf(NULL, 0, "%i", i);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", i);
		for (int j = 0; j < length; j++)
		{
			tmp_num_line.push_back(str[j]); // номер линии это уже строка
		}

		if (i == 1)
		{
			color = "red";
		}

		if (i == 2)
		{
			color = "blue";
		}

		if (i == 3)
		{
			color = "orange";
		}

		if (i == 4)
		{
			color = "green";
		}

		if (i == 5)
		{
			color = "cyan";
		}

		if (i == 6)
		{
			color = "yellow";
		}

		if (i == 7)
		{
			color = "violet";
		}

		if (i == 8)
		{
			color = "pink";
		}

		if (i == 9)
		{
			color = "dark-red";
		}

		if (i == 10)
		{
			color = "dark-yellow";
		}

		if (i == 11)
		{
			color = "dark-orange";
		}

		if (i == 12)
		{
			color = "dark-green";
		}

		if (i == 13)
		{
			color = "dark-cyan";
		}

		if (i == 14)
		{
			color = "dark-blue";
		}

		if (i == 15)
		{
			color = "dark-violet";
		}

		if (i == 16)
		{
			color = "dark-pink";
		}

		string str_str = "set linetype " + tmp_num_line + " lw " + tmp_width_line + " lc rgb " + "\"" + color + "\"\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str.c_str());
	}

	string title = "set title \"" + title_plot + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str()); //передадим текущее название плота (менть Te Ti)
	title.clear();
	title = "set xlabel \"" + xlabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str()); // передаем текущее название осей (подправить универсалом)
	title.clear();
	title = "set ylabel \"" + ylabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());

	//!!!!!!!!!!!!!!!!!!!!!
	/////////////////////////////////////////
	/*title.clear();
	title = "set label \" (200,28) \" at 200,27 \"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	title.clear();
	title = "set label \" (250,40) \" at 242,42 \"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	title.clear();
	title = "set label \" (300,52) \" at 286,52 \"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());*/
	/////////////////////////////////////////////////

	///fprintf(gnuplotPipe[Current_number_plot], ("set term " + GetGPTerminal() + " position 0,0 size %zu,%zu\n").c_str(), GetTermnalWidth(), 415);
	//fprintf(gnuplotPipe[Current_number_plot], "plot [][0:1] 2\n");
	//printf("Press enter when window will appear");
	fflush(gnuplotPipe[Current_number_plot]); // Пропихиваем данные туда
}

void GnuPlot::SetParametrsOnPlotColor(int Current_number_plot, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y) // Current_number_plot=0,1,2, ....
{
	//gnuplotPipe[Current_number_plot] = _popen(GetGPPath().c_str(), "w"); // открываем текущий плот
	fprintf(gnuplotPipe_Gif[Current_number_plot], "reset\n");
	// вместе с русской кодировкой
	/*string stroka1 = ConvertNumToStringdouble_gnuplot(0.);
	string stroka2 = ConvertNumToStringdouble_gnuplot(right_bondary_x);
	string stroka3 = ConvertNumToStringdouble_gnuplot(0.);
	string stroka4 = ConvertNumToStringdouble_gnuplot(top_bondary_y);

	stroka1.replace(stroka1.find(","), 1,".");
	stroka2.replace(stroka2.find(","), 1, ".");
	stroka3.replace(stroka3.find(","), 1, ".");
	stroka4.replace(stroka4.find(","), 1, ".");

	string total_str = "set xrange ["+ stroka1  + ":" + stroka2 +"]\n";
	fprintf(gnuplotPipe[Current_number_plot], total_str.c_str());
	total_str.clear();
	total_str = "set yrange [" + stroka3 + ":" + stroka4 + "]\n";
	fprintf(gnuplotPipe[Current_number_plot], total_str.c_str());*/

	///fprintf(gnuplotPipe[Current_number_plot], "set yrange [%f:%f]\n", 0., top_bondary_y);
	fprintf(gnuplotPipe[Current_number_plot], "set xrange [%lf:%lf]\n", 0., right_bondary_x);
	fprintf(gnuplotPipe[Current_number_plot], "set yrange [%lf:%lf]\n", 0., top_bondary_y);

	fprintf(gnuplotPipe[Current_number_plot], "set palette defined ( 0 \"dark-violet\", 1 \"blue\", 2 \"cyan\", 3 \"green\", 4 \"yellow\", 5 \"orange\", 6 \"red\")\n");
	if (Current_number_plot == 13 || Current_number_plot == 14)
	{
		fprintf(gnuplotPipe[Current_number_plot], "set palette defined (0 \"black\", 1 \"green\", 1 \"red\")\n");
	}

	if (Current_number_plot == 20 || Current_number_plot == 21)
	{
		fprintf(gnuplotPipe[Current_number_plot], "set palette defined (0 \"blue\", 1 \"red\")\n");
	}

	string title = "set title \"" + title_plot + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	title.clear();
	title = "set xlabel \"" + xlabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	title.clear();
	title = "set ylabel \"" + ylabel + "\"\n";
	fprintf(gnuplotPipe[Current_number_plot], title.c_str());
	///fprintf(gnuplotPipe[Current_number_plot], ("set term " + GetGPTerminal() + " position 0,0 size %zu,%zu\n").c_str(), GetTermnalWidth(), 415);
	//fprintf(gnuplotPipe[Current_number_plot], "plot [][0:1] 2\n");
	//printf("Press enter when window will appear");
	fflush(gnuplotPipe[Current_number_plot]); // Пропихиваем данные туда
}

void GnuPlot::SetGridOnPlot2D(int Current_number_plot, int Nx, int Ny, double dx, double dy, int number_of_lines, Profile prof)
{
	if (prof == fun_on_x)
	{
		//Matrix mtr1(Nx, number_of_lines + 1);
		//P_xyzt[Current_number_plot] = mtr1; // рабочий вариант
		P_xyzt[Current_number_plot].Resize(Nx, number_of_lines + 1);// если out of range  будет то еще на 1 добавим

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dx;
		}
	}

	if (prof == fun_on_y)
	{
		P_xyzt[Current_number_plot].Resize(Ny, number_of_lines + 1);

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dy;
		}
	}
}

void GnuPlot::SetGridOnPlot3D(int Current_number_plot, int Nx, int Ny, int Nz, double dx, double dy, double dz, int number_of_lines, Profile prof)
{
	if (prof == fun_on_x)
	{
		//Matrix mtr1(Nx, number_of_lines + 1);
		//P_xyzt[Current_number_plot] = mtr1; // рабочий вариант
		P_xyzt[Current_number_plot].Resize(Nx, number_of_lines + 1);// если out of range  будет то еще на 1 добавим

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dx;
		}


		//if (count_frame != 0)
		{
			P_xyzt_Gif[Current_number_plot].Resize(Nx, number_of_lines + 1);// если out of range  будет то еще на 1 добавим

			for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
			{
				P_xyzt_Gif[Current_number_plot][i][0] = i * dx;
			}
		}
	}

	if (prof == fun_on_y)
	{
		P_xyzt[Current_number_plot].Resize(Ny, number_of_lines + 1);

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dy;
		}


		//if (count_frame != 0)
		{
			P_xyzt_Gif[Current_number_plot].Resize(Ny, number_of_lines + 1);

			for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
			{
				P_xyzt_Gif[Current_number_plot][i][0] = i * dy;
			}
		}
	}

	if (prof == fun_on_z)
	{
		P_xyzt[Current_number_plot].Resize(Nz, number_of_lines + 1);

		for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			P_xyzt[Current_number_plot][i][0] = i * dz;
		}

		//if (count_frame != 0)
		{
			P_xyzt_Gif[Current_number_plot].Resize(Nz, number_of_lines + 1);

			for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
			{
				P_xyzt_Gif[Current_number_plot][i][0] = i * dz;
			}
		}



		if (Current_number_plot == 28) // сделано под конкретную задачу
		{
			for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
			{
				P_xyzt[Current_number_plot][i][1] = i * dz;
			}

			//if (count_frame != 0)
			{
				P_xyzt_Gif[Current_number_plot].Resize(Nz, number_of_lines + 1);

				for (size_t i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					P_xyzt_Gif[Current_number_plot][i][1] = i * dz;
				}
			}
		}

	}
}

void GnuPlot::SetDataOnPlot1D(int Current_number_plot, int Nx, double dx, double** fun_dimensionless, double parametr_for_dimension, int moment_of_time, double dt, Profile prof)
{
	if (prof == fun_on_x) // current_number_of_line должа начинатс с 1 
	{ // такого плана графики лучше создавать на последих по номеру плотах т.к. мб поадобитс закрытие файлов предыдущих плотов 
		P_xyzt[Current_number_plot].Resize(Nx, moment_of_time + 1);// если out of range  будет то еще на 1 добавим

		for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++)
		{
			P_xyzt[Current_number_plot][j][0] = j * dx;
			for (size_t k = 1; k < moment_of_time; k++) {
				P_xyzt[Current_number_plot][j][k] = parametr_for_dimension * fun_dimensionless[j][k];//????
			}
		}

		for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			for (int j = 0; j < moment_of_time + 1; j++)
			{
				file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
			}
			file[Current_number_plot] << endl;
		}
	}

	if (prof == fun_on_t)
	{
		P_xyzt[Current_number_plot].Resize(moment_of_time + 1, Nx);// если out of range  будет то еще на 1 добавим

		for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++)
		{
			P_xyzt[Current_number_plot][j][0] = j * dt;
			for (size_t k = 1; k < Nx; k++) {
				P_xyzt[Current_number_plot][j][k] = parametr_for_dimension * fun_dimensionless[j][k];
			}
		}

		for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
		{
			for (int j = 0; j < Nx; j++)
			{
				file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
			}
			file[Current_number_plot] << endl;
		}
	}
}

void GnuPlot::SetDataOnPlot2D(int Current_number_plot, double** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double**> vec, Profile prof)
{
	if (prof == fun_on_x) // current_number_of_line должа начинатс с 1 
	{ // такого плана графики лучше создавать на последих по номеру плотах т.к. мб поадобитс закрытие файлов предыдущих плотов 
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[j][fixed_point_on_axis_y];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_y)
	{
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_x][j];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}
	}

	if (prof == fun_on_t)
	{
		//  moment_of_time - размерна величина
		file[Current_number_plot] << moment_of_time << "   ";
		for (int i = 0; i < vec.size(); i++)
		{
			file[Current_number_plot] << parametr_for_dimension * vec[i][fixed_point_on_axis_x][fixed_point_on_axis_y] << "   ";
		}
		file[Current_number_plot] << endl;
	}
}

void GnuPlot::SetDataOnPlot3D(int Current_number_plot, int current_count_frame, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis_x, int fixed_point_on_axis_y, int fixed_point_on_axis_z, int number_of_lines, int current_number_of_line, double moment_of_time, vector<double***> vec, Profile prof)
{//при реализации откртыть файлы (P от x сери врем) (P от z сери врем) (T от времени)

	string type_plot_Gif = ConvertNumToStringint(Current_number_plot - 1);
	string current_count_frame_str = ConvertNumToStringint(current_count_frame - 1);
	ofstream file_anim("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");
	//filename_Gif.push_back("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");//????

	if (prof == fun_on_x) // current_number_of_line должа начинатс с 1 
	{ // такого плана графики лучше создавать на последих по номеру плотах т.к. мб поадобитс закрытие файлов предыдущих плотов 

		//для статики
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[j][fixed_point_on_axis_y][fixed_point_on_axis_z];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}

		//для анимаци
		//if (count_frame != 0)
		{
			if (current_number_of_line <= number_of_lines)
			{
				for (size_t j = 0; j < P_xyzt_Gif[Current_number_plot].GetRows(); j++) {
					P_xyzt_Gif[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[j][fixed_point_on_axis_y][fixed_point_on_axis_z];
				}
				if (current_number_of_line == number_of_lines)
				{
					// sprintf('Data_Gif" + type_plot_Gif + " + " + "%%d" + ".txt', i)
					//file_Gif[Current_number_plot].open("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");
					for (int i = 0; i < P_xyzt_Gif[Current_number_plot].GetRows(); i++)
					{
						for (int j = 0; j < number_of_lines + 1; j++)
						{
							file_anim << P_xyzt_Gif[Current_number_plot][i][j] << "   ";
						}
						file_anim << endl;
					}

				}
			}
		}

		if (Current_number_plot == 0 || Current_number_plot == 22)
		{
			//ofstream file_anim("e(xy) z(index) = " + z_str + " + " + current_count_frame_str + ".txt");
			// продолить запись профилей P(x) легенда z но для от y = Ny/2 до y = Ny/2-40 (для тепловой задачи)
		}

		if (Current_number_plot == 17 || Current_number_plot == 23)
		{
			//ofstream file_anim("e(xy) z(index) = " + z_str + " + " + current_count_frame_str + ".txt");

		}

	}

	if (prof == fun_on_y)
	{
		//для статики
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_x][j][fixed_point_on_axis_z];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}


		//для анимаци
		//if (count_frame != 0)
		{
			if (current_number_of_line <= number_of_lines)
			{
				for (size_t j = 0; j < P_xyzt_Gif[Current_number_plot].GetRows(); j++) {
					P_xyzt_Gif[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_x][j][fixed_point_on_axis_z];
				}
				if (current_number_of_line == number_of_lines)
				{
					for (int i = 0; i < P_xyzt_Gif[Current_number_plot].GetRows(); i++)
					{
						for (int j = 0; j < number_of_lines + 1; j++)
						{
							file_anim << P_xyzt_Gif[Current_number_plot][i][j] << "   ";
						}
						file_anim << endl;
					}
				}
			}
		}
	}

	if (prof == fun_on_z)
	{
		if (current_number_of_line <= number_of_lines)
		{
			for (size_t j = 0; j < P_xyzt[Current_number_plot].GetRows(); j++) {
				P_xyzt[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_x][fixed_point_on_axis_y][j];
			}
			if (current_number_of_line == number_of_lines)
			{
				for (int i = 0; i < P_xyzt[Current_number_plot].GetRows(); i++)
				{
					for (int j = 0; j < number_of_lines + 1; j++)
					{
						file[Current_number_plot] << P_xyzt[Current_number_plot][i][j] << "   ";
					}
					file[Current_number_plot] << endl;
				}
			}
		}

		//if (count_frame != 0)
		{
			if (current_number_of_line <= number_of_lines)
			{
				for (size_t j = 0; j < P_xyzt_Gif[Current_number_plot].GetRows(); j++) {
					P_xyzt_Gif[Current_number_plot][j][current_number_of_line] = parametr_for_dimension * fun_dimensionless[fixed_point_on_axis_x][fixed_point_on_axis_y][j];
				}
				if (current_number_of_line == number_of_lines)
				{
					for (int i = 0; i < P_xyzt_Gif[Current_number_plot].GetRows(); i++)
					{
						for (int j = 0; j < number_of_lines + 1; j++)
						{
							file_anim << P_xyzt_Gif[Current_number_plot][i][j] << "   ";
						}
						file_anim << endl;
					}
				}
			}
		}
	}

	if (prof == fun_on_t)
	{
		//  moment_of_time - размерна величина
		file[Current_number_plot] << moment_of_time << "   ";
		for (int i = 0; i < vec.size(); i++)
		{
			file[Current_number_plot] << parametr_for_dimension * vec[i][fixed_point_on_axis_x][fixed_point_on_axis_y][fixed_point_on_axis_z] << "   ";
		}
		file[Current_number_plot] << endl;


		if (count_frame != 0)
		{
			file_anim << moment_of_time << "   ";
			for (int i = 0; i < vec.size(); i++)
			{
				file_anim << parametr_for_dimension * vec[i][fixed_point_on_axis_x][fixed_point_on_axis_y][fixed_point_on_axis_z] << "   ";
			}
			file_anim << endl;
		}
	}
}

void GnuPlot::SetDataOnPlotColor3D(int Current_number_plot, int current_count_frame, int Nx, int Ny, int Nz, double dx, double dy, double dz, double*** fun_dimensionless, double parametr_for_dimension, int fixed_point_on_axis, Plane plane)
{
	string type_plot_Gif = ConvertNumToStringint(Current_number_plot - 1);
	string current_count_frame_str = ConvertNumToStringint(current_count_frame - 1);
	ofstream file_anim("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");

	//  переаем размеры области и зачение функций в размером виде (иселючение: fun в SetDataOnPlotColor3D)!!!!!
	if (plane == xy)
	{
		/*if (Current_number_plot == 14)
		{
			for (size_t i = 0; i < Nx; i++)
				for (size_t j = 0; j < Ny; j++) {
					file[Current_number_plot] << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun_dimensionless[i][j][fixed_point_on_axis] << endl;
				}

			for (size_t i = 0; i < Nx; i++)
				for (size_t j = 0; j < Ny; j++) {
					file_anim << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun_dimensionless[i][j][fixed_point_on_axis] << endl;
				}
		}
		else*/
		{
			//if (Current_number_plot != 21)
			{
				for (size_t i = 0; i < Nx; i++)
					for (size_t j = 0; j < Ny; j++) {
						file[Current_number_plot] << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun_dimensionless[i][j][fixed_point_on_axis] << endl;
					}

				//if (count_frame != 0)
				{
					//cout << "hello ANIMATE " << endl;
					//cout << "Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt" << endl;
					//file_Gif[Current_number_plot].open("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");
					for (size_t i = 0; i < Nx; i++)
						for (size_t j = 0; j < Ny; j++) {
							file_anim << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun_dimensionless[i][j][fixed_point_on_axis] << endl;
						}
				}
			}
		}


		//if (Current_number_plot == 5)
		//{
		//	string z_str;
		//	for (int k = 1; k <= 20; k++) // цикл по индексам оси z
		//	{
		//		z_str = ConvertNumToStringint(k - 1);
		//		ofstream file_anim1("e(xy) z(index) = " + z_str + " + " + current_count_frame_str + ".txt");

		//		for (size_t i = 0; i < Nx; i++)
		//			for (size_t j = 0; j < Ny; j++) {
		//				file_anim1 << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun_dimensionless[i][j][k] << endl;
		//			}

		//		file_anim1.close();
		//	}
		//}
	}

	if (plane == xz)
	{
		for (size_t i = 0; i < Nx; i++)
			for (size_t j = 0; j < Nz; j++) {
				file[Current_number_plot] << j * dz << "   " << i * dx << "   " << parametr_for_dimension * fun_dimensionless[i][fixed_point_on_axis][j] << endl;
			}

		//if (count_frame != 0)
		{
			//file_Gif[Current_number_plot].open("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");
			for (size_t i = 0; i < Nx; i++)
				for (size_t j = 0; j < Nz; j++) {
					file_anim << j * dz << "   " << i * dx << "   " << parametr_for_dimension * fun_dimensionless[i][fixed_point_on_axis][j] << endl;
				}
		}
	}

	if (plane == yz)
	{
		for (size_t i = 0; i < Ny; i++)
			for (size_t j = 0; j < Nz; j++) {
				file[Current_number_plot] << j * dz << "   " << i * dy << "   " << parametr_for_dimension * fun_dimensionless[fixed_point_on_axis][i][j] << endl;
			}

		//if (count_frame != 0)
		{
			//file_Gif[Current_number_plot].open("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");
			for (size_t i = 0; i < Ny; i++)
				for (size_t j = 0; j < Nz; j++) {
					file_anim << j * dz << "   " << i * dy << "   " << parametr_for_dimension * fun_dimensionless[fixed_point_on_axis][i][j] << endl;
				}
		}
	}
}

void GnuPlot::SetDataOnPlotColor2D(int Current_number_plot, int current_count_frame, int Nx, int Ny, double dx, double dy, double** fun, double parametr_for_dimension)
{
	string type_plot_Gif = ConvertNumToStringint(Current_number_plot - 1);
	string current_count_frame_str = ConvertNumToStringint(current_count_frame - 1);
	ofstream file_anim("Data_Gif" + type_plot_Gif + " + " + current_count_frame_str + ".txt");
	for (size_t i = 0; i < Nx; i++)
		for (size_t j = 0; j < Ny; j++) {
			file[Current_number_plot] << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun[i][j] << endl;
		}

	//if (count_frame != 0)
	{
		for (size_t i = 0; i < Nx; i++)
			for (size_t j = 0; j < Ny; j++) {
				file_anim << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun[i][j] << endl;
			}
	}

	//if (Current_number_plot == 21)
	//{
	//	string z_str;
	//	for (int i = 0; i <= 50; i++) // цикл по индексам оси z
	//	{
	//		z_str = ConvertNumToStringint(i - 1);
	//		ofstream file_anim("e(xy) z(index) = " + z_str + " + " + current_count_frame_str + ".txt");

	//		for (size_t i = 0; i < Nx; i++)
	//			for (size_t j = 0; j < Ny; j++) {
	//				file_anim << i * dx << "   " << j * dy << "   " << parametr_for_dimension * fun[i][j] << endl;
	//			}

	//		file_anim.close();
	//	}
	//}

}

void GnuPlot::ShowDataOnPlot2D(int Current_number_plot, bool DataFromDiffFiles, int number_of_lines, vector<string> list_name_line, string name_of_file, bool png_)
{
	if (png_)
	{
		fprintf(gnuplotPipe[Current_number_plot], "set terminal png\n");
		// name_of_file - имя файла png
		string str_str_str = "set output \"" + name_of_file + ".png\" \n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str());
		str_str_str.clear();

		string str_str = "plot ";
		for (int i = 1; i <= number_of_lines; i++) // number_of_lines - число линий на плоте, 
		{
			string tmp_num_line;

			int length = snprintf(NULL, 0, "%i", i + 1);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", i + 1);
			for (int j = 0; j < length; j++)
			{
				tmp_num_line.push_back(str[j]); // номер столбца
			}
			// string name_line - назваие линии в легенде (например когда рисуем значение давлен от икса в различные моменты времени )
			string name_line = list_name_line[i - 1]; // если не заработает, то ужно строки в сишные переводить .c_str(), а стринг это C++
			if (DataFromDiffFiles)
			{
				str_str += "'" + filename[i - 1] + "' u 1:2 w l title \"" + name_line + "\", ";
			}
			else
			{
				str_str += "'" + filename[Current_number_plot] + "' u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ";
			}
			//str_str += "'" + filename[Current_number_plot] + "' u 1:" + tmp_num_line + " w lp lc rgb 'red' lw 2 pt 5 title \"" + name_line + "\", ";
		}
		str_str.pop_back();
		str_str.pop_back();
		str_str += "\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str.c_str()/*, filename[Current_number_filename].c_str()*/);
	}
	else
	{
		// здесь будет только цикл по числу линий без первых 2-ух строчек выше от него
		string str_str = "plot ";
		for (int i = 1; i <= number_of_lines; i++) // number_of_lines - число линий на плоте, 
		{
			string tmp_num_line;

			int length = snprintf(NULL, 0, "%i", i + 1);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", i + 1);
			for (int j = 0; j < length; j++)
			{
				tmp_num_line.push_back(str[j]); // номер столбца
			}
			// string name_line - назваие линий в легенде 
			string name_line = list_name_line[i - 1];
			if (DataFromDiffFiles)
			{
				str_str += "'" + filename[i - 1] + "' u 1:2 w l title \"" + name_line + "\", ";
			}
			else
			{
				str_str += "'" + filename[Current_number_plot] + "' u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ";
			}
		}
		str_str.pop_back();
		str_str.pop_back();
		str_str += "\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str.c_str()/*, filename[Current_number_filename].c_str()*/);
	}
	fflush(gnuplotPipe[Current_number_plot]);
}

void GnuPlot::CreateGifOnPlot2D(int Current_number_plot, int number_of_lines, int width_line, int count_frame, vector<int> moments_fix_time_anim, string title_plot, string xlabel, string ylabel, vector<string> list_name_line, string name_of_file, string Time_unit)
{
	//if (gif_)
	{
		fprintf(gnuplotPipe_Gif[Current_number_plot], "reset\n");
		fprintf(gnuplotPipe_Gif[Current_number_plot], "set terminal gif font arial 11 animate delay 100 loop 0 size 640,480\n"); // delay 10 = 0.1 сек; loop 0 - бесконечный цикл повтора (1 - 1 раз повтор)
		string str_str_str = "set output \"" + name_of_file + ".gif\" \n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], str_str_str.c_str());
		fprintf(gnuplotPipe_Gif[Current_number_plot], "set macros\n"); // для задания строковой переменной
		str_str_str.clear();

		string tmp_width_line;
		int length = snprintf(NULL, 0, "%i", width_line);
		char* str = new char[length + 1];
		snprintf(str, length + 1, "%i", width_line);
		for (int j = 0; j < length; j++)
		{
			tmp_width_line.push_back(str[j]);
		}
		string color;

		for (int i = 1; i <= number_of_lines; i++) // number_of_lines - число линий на плоте, 
		{
			string tmp_num_line;
			int length = snprintf(NULL, 0, "%i", i);
			char* str = new char[length + 1];
			snprintf(str, length + 1, "%i", i);
			for (int j = 0; j < length; j++)
			{
				tmp_num_line.push_back(str[j]); // номер линии это уже строка
			}

			if (i == 1)
			{
				color = "red";
			}

			if (i == 2)
			{
				color = "blue";
			}

			if (i == 3)
			{
				color = "orange";
			}

			if (i == 4)
			{
				color = "green";
			}

			if (i == 5)
			{
				color = "cyan";
			}

			if (i == 6)
			{
				color = "yellow";
			}

			if (i == 7)
			{
				color = "violet";
			}

			if (i == 8)
			{
				color = "pink";
			}

			if (i == 9)
			{
				color = "dark-red";
			}

			if (i == 10)
			{
				color = "dark-yellow";
			}

			if (i == 11)
			{
				color = "dark-orange";
			}

			if (i == 12)
			{
				color = "dark-green";
			}

			if (i == 13)
			{
				color = "dark-cyan";
			}

			if (i == 14)
			{
				color = "dark-blue";
			}

			if (i == 15)
			{
				color = "dark-violet";
			}

			if (i == 16)
			{
				color = "dark-pink";
			}

			string str_str = "set linetype " + tmp_num_line + " lw " + tmp_width_line + " lc rgb " + "\"" + color + "\"\n";
			fprintf(gnuplotPipe_Gif[Current_number_plot], str_str.c_str());
		}

		string MyArray;
		string count_frame_string = ConvertNumToStringint(count_frame - 2);
		//string type_plot_Gif = ConvertNumToStringint(Current_number_plot - 1);

		MyArray = "aGet(name, i) = value(sprintf('_%%s_%%i', name, i)) \n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
		MyArray.clear();
		MyArray = "aSet(name, i,value) = sprintf('_%%s_%%i = %%.16e', name, i, value) \n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
		MyArray.clear();

		vector<int>::iterator it_fix_time_anim;
		vector<string>::iterator it_str;
		int kkk = 0;
		for (it_fix_time_anim = moments_fix_time_anim.begin(); it_fix_time_anim != moments_fix_time_anim.end(); it_fix_time_anim++)
		{
			string fix_time_anim_string = ConvertNumToStringint_1((*it_fix_time_anim));
			string i_str = ConvertNumToStringint_1(kkk);
			MyArray = "eval aSet(\"A\"," + i_str + "," + fix_time_anim_string + ")\n";
			fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
			kkk++;
			MyArray.clear();
		}

		MyArray = "aSet(name, i,value) = sprintf('_%%s_%%i = %%s', name, i, value) \n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
		MyArray.clear();
		kkk = 0;
		for (it_str = filename.begin(); it_str != filename.end(); it_str++)
		{
			string result_str = "colString = \"" + (*it_str) + "\"\n";
			fprintf(gnuplotPipe_Gif[Current_number_plot], result_str.c_str());
			string i_str = ConvertNumToStringint_1(kkk);
			MyArray = "eval aSet(\"B\"," + i_str + ",\"colString\")\n";
			fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
			kkk++;
			MyArray.clear();
		}

		string title;
		//title = "set title \"" + title_plot + "\"\n";
		//fprintf(gnuplotPipe_Gif[Current_number_plot], title.c_str());
		//title.clear();
		title = "set xlabel \"" + xlabel + "\"\n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], title.c_str()); // передаем текущее название осей (подправить универсалом)
		title.clear();
		title = "set ylabel \"" + ylabel + "\"\n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], title.c_str());
		title.clear();
																									
		/*!!!*///string str_str = "do for [i=0:" + count_frame_string + "] {set title sprintf('" + title_plot + " %%d fs', aGet(\"A\",i)); plot sprintf(aGet(\"B\",i))";//!!!!!! ПРОДОЛЖИТЬ С ЭТОГО МЕСТА
		string str_str = "do for [i=0:" + count_frame_string + "] {set title sprintf('" + title_plot + " %%d " + Time_unit + "', aGet(\"A\",i)); plot sprintf(aGet(\"B\",i))";//!!!!!! ПРОДОЛЖИТЬ С ЭТОГО МЕСТА

		for (int ii = 1; ii <= number_of_lines; ii++) // number_of_lines - число линий на плоте, 
		{
			string tmp_num_line = ConvertNumToStringint(ii);
			// string name_line - назваие линии в легенде (например когда рисуем значение давлен от икса в различные моменты времени )
			string name_line = list_name_line[ii - 1]; // если не заработает, то ужно строки в сишные переводить .c_str(), а стринг это C++
			str_str += " u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ''";
		}

		str_str.pop_back();
		str_str.pop_back();
		str_str.pop_back();
		str_str.pop_back();
		str_str += "}"; // аналог delay (значение устанавливается в секундах)
		str_str += "\n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], str_str.c_str());
	}
	//else
	//{
	//	string count_frame_string = ConvertNumToStringint(count_frame - 1);
	//	string str_str = "do for [i=1:" + count_frame_string + "] { plot sprintf('Data%%d.txt', i)";

	//	for (int ii = 1; ii <= number_of_lines; ii++) // number_of_lines - число линий на плоте, 
	//	{
	//		string tmp_num_line = ConvertNumToStringint(ii);
	//		// string name_line - назваие линии в легенде (например когда рисуем значение давлен от икса в различные моменты времени )
	//		string name_line = list_name_line[ii - 1]; // если не заработает, то ужно строки в сишные переводить .c_str(), а стринг это C++
	//		str_str += " u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ";
	//		//str_str += "'" + filename[k - 1] + "' u 1:" + tmp_num_line + " w l title \"" + name_line + "\", ";
	//		//str_str += "'" + filename[Current_number_plot] + "' u 1:" + tmp_num_line + " w lp lc rgb 'red' lw 2 pt 5 title \"" + name_line + "\", ";
	//	}

	//	str_str.pop_back();
	//	str_str.pop_back();
	//	str_str += "; pause 1 }"; // аналог delay (значение устанавливается в секундах)
	//	str_str += "\n";
	//	fprintf(gnuplotPipe[Current_number_plot], str_str.c_str()/*, filename[Current_number_filename].c_str()*/);
	//}
	fflush(gnuplotPipe_Gif[Current_number_plot]);
	//str_str += "'" + filename[i - 1] + "' u 1:2 w l title \"" + name_line + "\", ";
}

void GnuPlot::ShowDataOnPlotColor(int Current_number_plot, string name_of_file, bool png_)
{
	if (png_)
	{
		fprintf(gnuplotPipe[Current_number_plot], "set terminal png\n");
		string str_str_str = "set output \"" + name_of_file + ".png\" \n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str());
		str_str_str.clear();
		str_str_str = "plot '" + filename[Current_number_plot] + "' u 1:2:3 with image\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str(), filename[Current_number_plot].c_str());
	}
	else
	{
		string str_str_str = "plot '" + filename[Current_number_plot] + "' u 1:2:3 with image\n";
		fprintf(gnuplotPipe[Current_number_plot], str_str_str.c_str(), filename[Current_number_plot].c_str());
	}
	fflush(gnuplotPipe[Current_number_plot]);
}

void GnuPlot::CreateGifOnPlotColor(int Current_number_plot, int count_frame, vector<int> moments_fix_time_anim, string title_plot, string xlabel, string ylabel, long float right_bondary_x, long float top_bondary_y, string Time_unit, string name_of_file)
{
	fprintf(gnuplotPipe_Gif[Current_number_plot], "reset\n");
	fprintf(gnuplotPipe_Gif[Current_number_plot], "set terminal gif font arial 11 animate delay 100 loop 0 size 640,480\n"); // delay 10 = 0.1 сек; loop 0 - бесконечный цикл повтора (1 - 1 раз повтор)
	string str_str_str = "set output \"" + name_of_file + ".gif\" \n";
	fprintf(gnuplotPipe_Gif[Current_number_plot], str_str_str.c_str());
	fprintf(gnuplotPipe_Gif[Current_number_plot], "set macros\n"); // для задания строковой переменной
	str_str_str.clear();

	/*string stroka1 = ConvertNumToStringdouble_gnuplot(0.);
	string stroka2 = ConvertNumToStringdouble_gnuplot(right_bondary_x);
	string stroka3 = ConvertNumToStringdouble_gnuplot(0.);
	string stroka4 = ConvertNumToStringdouble_gnuplot(top_bondary_y);

	stroka1.replace(stroka1.find(","), 1, ".");
	stroka2.replace(stroka2.find(","), 1, ".");
	stroka3.replace(stroka3.find(","), 1, ".");
	stroka4.replace(stroka4.find(","), 1, ".");

	string total_str = "set xrange [" + stroka1 + ":" + stroka2 + "]\n";
	fprintf(gnuplotPipe[Current_number_plot], total_str.c_str());
	total_str.clear();
	total_str = "set yrange [" + stroka3 + ":" + stroka4 + "]\n";
	fprintf(gnuplotPipe[Current_number_plot], total_str.c_str());*/

	fprintf(gnuplotPipe_Gif[Current_number_plot], "set xrange [%lf:%lf]\n", 0., right_bondary_x);
	fprintf(gnuplotPipe_Gif[Current_number_plot], "set yrange [%lf:%lf]\n", 0., top_bondary_y);

	fprintf(gnuplotPipe_Gif[Current_number_plot], "set palette defined ( 0 \"dark-violet\", 1 \"blue\", 2 \"cyan\", 3 \"green\", 4 \"yellow\", 5 \"orange\", 6 \"red\")\n");
	if (Current_number_plot == 13 || Current_number_plot == 14)// для отрисовки зоны плавления в двух плоскостях
	{
		fprintf(gnuplotPipe_Gif[Current_number_plot], "set palette defined (0 \"blue\", 1 \"red\")\n");
	}

	if (Current_number_plot == 20 || Current_number_plot == 21)
	{
		fprintf(gnuplotPipe[Current_number_plot], "set palette defined (0 \"red\", 1 \"blue\")\n");
	}

	string MyArray;
	string count_frame_string = ConvertNumToStringint(count_frame - 2);
	//string type_plot_Gif = ConvertNumToStringint(Current_number_plot - 1);

	MyArray = "aGet(name, i) = value(sprintf('_%%s_%%i', name, i)) \n";
	fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
	MyArray.clear();
	MyArray = "aSet(name, i,value) = sprintf('_%%s_%%i = %%.16e', name, i, value) \n";
	fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
	MyArray.clear();

	vector<int>::iterator it_fix_time_anim;
	vector <string>::iterator it_str;
	int kkk = 0;
	for (it_fix_time_anim = moments_fix_time_anim.begin(); it_fix_time_anim != moments_fix_time_anim.end(); it_fix_time_anim++)
	{
		string fix_time_anim_string = ConvertNumToStringint_1((*it_fix_time_anim));
		string i_str = ConvertNumToStringint_1(kkk);
		kkk++;
		MyArray = "eval aSet(\"A\"," + i_str + "," + fix_time_anim_string + ")\n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
		MyArray.clear();
	}

	MyArray = "aSet(name, i,value) = sprintf('_%%s_%%i = %%s', name, i, value) \n";
	fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
	MyArray.clear();
	kkk = 0;
	for (it_str = filename.begin(); it_str != filename.end(); it_str++)
	{
		string result_str = "colString = \"" + (*it_str) + "\"\n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], result_str.c_str());
		string i_str = ConvertNumToStringint_1(kkk);
		MyArray = "eval aSet(\"B\"," + i_str + ",\"colString\")\n";
		fprintf(gnuplotPipe_Gif[Current_number_plot], MyArray.c_str());
		kkk++;
		MyArray.clear();
	}

	string title;// = "set title \"" + title_plot + "\"\n";
	//fprintf(gnuplotPipe_Gif[Current_number_plot], title.c_str());
	//title.clear();
	title = "set xlabel \"" + xlabel + "\"\n";
	fprintf(gnuplotPipe_Gif[Current_number_plot], title.c_str());
	title.clear();
	title = "set ylabel \"" + ylabel + "\"\n";
	fprintf(gnuplotPipe_Gif[Current_number_plot], title.c_str());
	title.clear();

	string str_str = "do for [i=0:" + count_frame_string + "] {set title sprintf('" + title_plot + " %%d " + Time_unit + "', aGet(\"A\",i)); plot sprintf(aGet(\"B\",i)) u 1:2:3 with image}\n";//!!!!!!

	fprintf(gnuplotPipe_Gif[Current_number_plot], str_str.c_str());
	fflush(gnuplotPipe_Gif[Current_number_plot]);
}

void GnuPlot::Close_and_open_files_for_replot(vector<int> Current_number_of_plots_)
{
	for (int i = 0; i < Current_number_of_plots_.size(); ++i)
	{
		file[Current_number_of_plots_[i]].close();
		file[Current_number_of_plots_[i]].open(filename[Current_number_of_plots_[i]]);
	}
}

void GnuPlot::Close_all_files_and_plots(int number_of_plots_)
{
	number_of_plots = number_of_plots_;
	for (int i = 0; i < number_of_plots; ++i)
	{
		file[i].close();
		fprintf(gnuplotPipe[i], "exit\n");
		_pclose(gnuplotPipe[i]);
	}
}


string GnuPlot::GetGPPath()
{
	return "\"C:\\Program Files\\gnuplot\\bin\\gnuplot\""; // где расположена программа
}

string GnuPlot::GetGPTerminal()
{
	return "wxt";
}

size_t GnuPlot::GetTermnalWidth() //ширина экрана 
{
	return 650;
}

