#pragma once
#include <iostream>
#include <vector>

using namespace std;

class Matrix
{
public:
	Matrix(int n_rows, int n_columns = 1);											//объявление
	//Matrix(int n_el);

	int nRows() const;
	int nColumns() const;
	
	double &operator()(int row, int column = 0);									//возвращается не копия, а ссылка, чтобы его можно было менять
	//double &operator()(int index);



private:
	vector<double> el;
	int n_rows;																	//строки
	int n_columns;																//столбцы
};

Matrix operator+(Matrix m1, Matrix m2);						//сложение матриц

Matrix operator-(Matrix m1, Matrix m2);						//вычитание матриц

Matrix operator*(Matrix m, double value);						//умножение матрицы на константу

Matrix operator*(double value, Matrix m);						//умножение константы на матрицу, для идиотов


Matrix operator*(Matrix m1, Matrix m2);							//умножение матрицы на матрицу


ostream &operator<<(ostream &out, Matrix m);				//вывод матрицы

double det(Matrix m);
Matrix jordan_gauss(Matrix a, Matrix b);
Matrix inv_jg(Matrix a);

double vect_norm1(Matrix v);
double vect_norm2(Matrix v);
double vect_norm3(Matrix v);
double matr_norm1(Matrix m);
double matr_norm2(Matrix m);
Matrix jacobi_iter(Matrix m, Matrix b, double eps);				//метод Якоби

//передаем функцию
using Func = double(double);
using Func2 = double(double, double);
/*
using my_int = int;
using dvect = vector<double>;
*/

using Funcs = double(double, vector<double> );
using FuncsV = double( Matrix );

double newton(Func f, double x, double eps = 0.0001);
double m_hord(Func f, double x, double eps = 0.0000001);
double m_iter(Func f, double x, double eps = 0.0001);


double int_rect(Func f, double a, double b, double eps = 0.0001);
double int_trap(Func f, double a, double b, double eps = 0.0001);
double int_simpson(Func f, double a, double b, double eps = 0.0001);

//double ode_Euler(Func2 f, double x0, double y0, double x, double eps = 0.0001);


//vector<double> vode_Euler(vector<Funcs*> f, double x0, vector<double> y0, double x, double eps);



vector<double> s_euler_solver(vector<Funcs*> f, double x0, vector<double> y0, double x, int n);
vector<double> s_ode_euler(vector<Funcs*> f, double x, vector<double> y, double h);
vector<double> s_ode_euler_koshi(vector<Funcs*> f, double x, vector<double> y, double h);

class S_Ode_Solver
{
public:
	using Method = vector<double>(vector<Funcs*>, double, vector<double>, double);
	S_Ode_Solver(Method *method, int order);

	vector<double> operator()(vector<Funcs*> f, double x0, vector<double> y0, double x, double eps);
private:
	vector<double> solve(vector<Funcs*> f, double x0, vector<double> y0, double x, int n);
	Method *method;
	int order;
};


double diff_1_p1(Func f, double x, double eps = 0.0001);
double diff_1_p2(Func f, double x, double eps = 0.0001);
double diff_2_p2(Func f, double x, double eps = 0.0001);

double minimize_dih(Func f, double a, double b, double eps = 0.0001);
double minimize_gold(Func f, double a, double b, double eps = 0.0001);
