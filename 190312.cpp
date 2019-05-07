#include "pch.h"
#include <iostream>
#include <vector>																//библиотека по работе с векторами
#include <string>
#include "amateurmathlib.h"
using namespace std;



double f1(double x)
{
	return x * x * exp(x) - 4;
	
}

double f2(double x, double y)
{
	return x + y;
}



int main()
{

	//выявление ошибок
	/*
	try
	{
		Matrix m1(2, 2);
		m1(0, 0) = 1;
		m1(0, 1) = 2;
		m1(1, 0) = 3;
		m1(1, 1) = 4;

		Matrix m2(2, 2);
		m2(0, 0) = 5;
		m2(0, 1) = 6;
		m2(1, 0) = 7;
		m2(1, 1) = 8;

		Matrix m3 = m1 + m2;

		cout << m3(1, 0) << endl;
	}
	catch (string s)
	{
		cout << s << endl;
	}

	*/
	Matrix a(3, 3);

	a(0, 0) = 10;
	a(0, 1) = 2;
	a(0, 2) = 3;

	a(1, 0) = 4;
	a(1, 1) = 10;
	a(1, 2) = 2;

	a(2, 0) = 3;
	a(2, 1) = 4;
	a(2, 2) = 10;



	Matrix b(3, 1);

	b(0) = 2;
	b(1) = 0;
	b(2) = 1;


	//Matrix x = jordan_gauss(a, b);


	//Matrix c = inv_jg(a);

	//cout << vect_norm1(b);
	/*
	try
	{
		cout << matr_norm1(a);
	}
	catch (string s)
	{
		cout << s << endl;
	}
	*/
	// cout << jacobi_iter(a, b, 0.001);

	//cout << newton(f1, -2.2) << endl;
	//cout << m_hord(f1, -10) << endl;

	//cout << m_hord(f1, -10) << endl;


	/*
	cout << int_rect(f1, 0, 2, 0.01) << endl;
	cout << int_trap(f1, 0, 2, 0.01) << endl;
	cout << int_simpson(f1, 0, 2, 0.01) << endl;

	cout << endl;

	cout << int_rect(f1, 0, 2, 0.00001) << endl;
	cout << int_trap(f1, 0, 2, 0.00001) << endl;
	cout << int_simpson(f1, 0, 2, 0.00001) << endl;
	*/

	//cout << ode_Euler(f2, 0, 1, 2, 0.001);


	
	vector<Funcs*> funcs(2);

	funcs[0] = [](double x, vector<double> y) -> double
	{
		return y[1];
	};
	funcs[1] = [](double x, vector<double> y) -> double
	{
		return x + 2 * y[0] - y[1];
	};
	double x0 = 0;
	vector<double> y0 = { 3. / 4., 1. / 2. };
	S_Ode_Solver euler(s_ode_euler_koshi, 2);
	vector<double> y = euler(funcs, x0, y0, 2., 0.0001);

	for (auto el: y)
	{
		cout << el << endl;
	}

	
	

	return 0;
}
