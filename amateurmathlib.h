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
/*
using my_int = int;
using dvect = vector<double>;
*/
double newton(Func f, double x, double eps = 0.0001);