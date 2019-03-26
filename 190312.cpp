#include "pch.h"
#include <iostream>
#include <vector>																//библиотека по работе с векторами
#include <string>
#include "amateurmathlib.h"
using namespace std;



double f1(double x)
{
	return x * x - 4;
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
	//починить

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

	cout << newton(f1, -2.2) << endl;

	return 0;
}
