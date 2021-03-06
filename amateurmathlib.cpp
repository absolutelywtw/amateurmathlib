#include "pch.h"
#include "amateurmathlib.h"
#include <string>
#include <iomanip>
using namespace std;

/*
namespace mylib
{
	int a;
}
mylib::a
*/




Matrix::Matrix(int n_rows, int n_columns)
/*:
n_rows(n_rows), n_columns(n_columns),
el(n_rows * n_columns)
*/
{
	if (n_rows < 1 || n_columns < 1)
		throw string("Wrong number of rows and columns");						//выкинуть ошибку
																				//this-> указатель на объект класса
	this->n_rows = n_rows;
	this->n_columns = n_columns;
	el.resize(n_rows * n_columns);
}

/*
Matrix::Matrix(int n_el)
{
	if (n_el < 1)
		throw string("Wrong number of rows and columns");						//выкинуть ошибку
																				//this-> указатель на объект класса
	this->n_rows = n_el;
	this->n_columns = 1;
	el.resize(n_rows * n_columns);
}
*/

int Matrix::nRows() const
{
	return n_rows;
}

int Matrix::nColumns() const
{
	return n_columns;
}

double &Matrix::operator()(int row, int column)					//выражаем матрицу в виде вектора
{
	/*
	234
	567
	108

	(2,3,4,5,6,7,1,0,8)
	*/

	if (row < 0 || row >= n_rows)
		throw string("Wrong row index");
	if (column < 0 || column >= n_columns)
		throw string("Wrong column index");
	return el[n_columns * row + column];
}

/*
double & Matrix::operator()(int index)
{
	return el[index];
}
*/

Matrix add(Matrix m1, Matrix m2, int sign)						//костыль для сложения-вычитания
{
	Matrix result(m1.nRows(), m1.nColumns());
	for (int i = 0; i < m1.nRows(); i++)
	{
		for (int j = 0; j < m1.nColumns(); j++)
		{
			result(i, j) = m1(i, j) + sign * m2(i, j);
		}
	}
	return result;
}

Matrix operator+(Matrix m1, Matrix m2)							//сложение матриц
{
	return add(m1, m2, +1);
}

Matrix operator-(Matrix m1, Matrix m2)							//вычитание матриц
{
	return add(m1, m2, -1);
}

Matrix operator*(Matrix m, double value)						//умножение матрицы на константу
{
	Matrix result(m.nRows(), m.nColumns());
	for (int i = 0; i < m.nRows(); i++)
	{
		for (int j = 0; j < m.nColumns(); j++)
		{
			result(i, j) = value * m(i, j);
		}
	}
	return result;
}

Matrix operator*(double value, Matrix m)						//умножение константы на матрицу, для идиотов
{
	return m * value;
}

Matrix operator*(Matrix m1, Matrix m2)							//умножение матрицы на матрицу
{
	Matrix result(m1.nRows(), m2.nColumns());
	for (int i = 0; i < m1.nRows(); i++)
	{

		for (int j = 0; j < m2.nColumns(); j++)
		{
			double tmp = 0;
			for (int k = 0; k < m2.nRows(); k++)
			{
				tmp += m1(i, k) * m2(k, j);
			}
			result(i, j) = tmp;
		}
	}
	return result;
}

ostream &operator<<(ostream &out, Matrix m)					//вывод матрицы
{
	for (int i = 0; i < m.nRows(); i++)
	{
		for (int j = 0; j < m.nColumns(); j++)
		{
			double tmp;										//ммм, костыли
			tmp = m(i, j);
			if (abs(tmp) < 0.000001)
				tmp = 0.;
			out << setw(9) << setprecision(3) << tmp << " \t";
		}
		out << endl;
	}
	return out;
}

double det(Matrix m)
{
	for (int i = 0; i < m.nRows() - 1; i++)
	{
		for (int j = i + 1; j < m.nRows(); j++)
		{
			double coeff = m(j, i) / m(i, i);
			for (int k = i; k < m.nColumns(); k++)
			{
				m(j, k) = m(j, k) - coeff * m(i, k);
			}

		}


	}

	double result = 1;
	for (int i = 0; i < m.nRows(); i++)
	{
		result *= m(i, i);
	}
	return result;

}

//метод Жордана-Гаусса
//A*X=B

Matrix jordan_gauss(Matrix a, Matrix b)
{
	/*
	1)
	723
	345
	410

	2)
	723
	0##
	0##

	*/



	for (int i = 0; i < a.nRows() - 1; i++)					//ведущая строка
	{
		for (int j = i + 1; j < a.nRows(); j++)				//текущая строка
		{
			double coeff = a(j, i) / a(i, i);				//коэффициент
			for (int k = i; k < a.nColumns(); k++)
			{
				a(j, k) = a(j, k) - coeff * a(i, k);
			}
			b(j, 0) = b(j, 0) - coeff * b(i, 0);
		}


	}

	for (int i = a.nRows() - 1; i >= 1; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double coeff = a(j, i) / a(i, i);
			a(j, i) = a(j, i) - coeff * a(i, i);
			b(j, 0) = b(j, 0) - coeff * b(i, 0);
		}
	}

	for (int i = 0; i < a.nRows(); i++)
	{
		b(i, 0) = b(i, 0) / a(i, i);
		a(i, i) = 1;
	}

	return b;

	//X = jordan_gauss(a, b)
	//для проверки считаем A*X, должно получиться B
}

Matrix inv_jg(Matrix a)
{

	if (a.nRows() != a.nColumns())
		throw string("Matrix is not square");

	if (abs(det(a)) < 0.00000000000001)										//вещественное число
		throw string(" Det cannot be equal to 0");

	Matrix result(a.nRows(), a.nColumns());
	for (int i = 0; i < a.nRows(); i++)
	{
		result(i, i) = 1;													//единичная матрица
	}






	for (int i = 0; i < a.nRows() - 1; i++)					//ведущая строка
	{
		for (int j = i + 1; j < a.nRows(); j++)				//текущая строка
		{
			double coeff = a(j, i) / a(i, i);				//коэффициент
			for (int k = i; k < a.nColumns(); k++)
			{
				a(j, k) = a(j, k) - coeff * a(i, k);
			}
			for (int k = 0; k < a.nColumns(); k++)
			{
				result(j, k) = result(j, k) - coeff * result(i, k);
			}


		}


	}

	for (int i = a.nRows() - 1; i >= 1; i--)
	{
		for (int j = i - 1; j >= 0; j--)
		{
			double coeff = a(j, i) / a(i, i);
			a(j, i) = a(j, i) - coeff * a(i, i);
			for (int k = 0; k < a.nColumns(); k++)
			{
				result(j, k) = result(j, k) - coeff * result(i, k);
			}

		}
	}

	for (int i = 0; i < a.nRows(); i++)
	{
		for (int k = 0; k < a.nColumns(); k++)
		{
			result(i, k) /= a(i, i);
		}
		a(i, i) = 1;
	}

	return result;

}

double vect_norm1(Matrix v)
{
	double result = 0;
	for (int i = 0; i < v.nRows(); i++)
	{
		if ( abs( v(i) ) > result )
		{
			result = abs( v(i) );
		}
	}

	return result;
}

double vect_norm2(Matrix v)
{
	double result = 0;
	for (int i = 0; i < v.nRows(); i++)
	{
		result += abs( v(i) );
	}
	return result;
}

double vect_norm3(Matrix v)
{
	double result = 0;
	for (int i = 0; i < v.nRows(); i++)
	{
		result += pow(v(i), 2);
	}
	return sqrt(result);
}

double matr_norm1(Matrix m)
{
	double result = 0;
	for (int i = 0; i < m.nRows(); i++)
	{
		double sum = 0;
		for (int j = 0; j < m.nColumns(); j++)
		{
			sum += abs( m(i,j) );
		}
		if (result < sum)
		{
			result = sum;
		}
	}
	return result;
}

double matr_norm2(Matrix m)
{
	double result = 0;
	for (int j = 0; j < m.nColumns(); j++)
	{
		double sum = 0;
		for (int i = 0; i < m.nRows(); i++)
		{
			sum += abs(m(i, j));
		}
		if (result < sum)
		{
			result = sum;
		}
	}
	return result;
}

Matrix jacobi_iter(Matrix m, Matrix b, double eps)
{
	Matrix result = b;
	for (int k = 0; k < m.nColumns(); k++)
	{
		b(k) /= m(k, k);
	}

	for (int i = 0; i < m.nRows(); i++)
	{
		for (int k = 0; k < m.nColumns(); k++)
		{
			if (k == i)
			{
				continue;
			}

			m(i, k) /= -m(i, i);
		}

		m(i, i) = 0;

	}


	while (42)
	{
		Matrix last = result;
		result = m * last + b;							//метод простых итерация
		if (vect_norm1(result - last) < eps)
		{
			break;
		}
	}


	return result;
}


double newton(Func f, double x, double eps)
{
	//x = x - df(x) / f;
	

	while (42)
	{
		//int n = 0;
		double last = x;
		double df = (f(x + eps) - f(x)) / eps;
		x = last - f(x) / df;
		//n++;
		if ( abs(x-last) < eps )
		{
			break;
		}
	}
	//cout << n << endl;
	return x;
}


double m_hord(Func f, double x, double eps)
{
//x_k+1 x_k-f(x_k) / ( f(a) - (f(x_k)*(a-x_k));

	double a = x - eps; 
	while (42)
	{
		//int n = 0;
		double last = x;
		//double df = (f(x + eps) - f(x)) / eps;
		x = last - f(last) * (a - last) / ( f(a) - f(last) );
		//n++;
		if (abs(x - last) < eps)
		{
			break;
		}
	}
	//cout << n << endl;
	return x;
}

double m_iter(Func f, double x, double eps)
{

	//f(x) = 0; -> x = g(x);
	double a = x - eps;
	while (42)
	{

		double last = x;

		x = last + f(x);

		if (abs(x - last) < eps)
		{
			break;
		}
	}
	//cout << n << endl;
	return x;
}





double helper_int_rect(Func f, double a, double b, int n)
{
	double h = (b - a) / n;
	double sum = 0;
	for (double x = a; x < b; x += h)
	{
		sum += f(x) * h;
	}

	return sum;

}


double int_rect(Func f, double a, double b, double eps)		//интегрирование методом прямоугольников
{
	int n = 2;
	double I = helper_int_rect(f, a, b, n);
	while (42)
	{
		
		n *= 2;
		double I_half = helper_int_rect(f, a, b, n);
		double R = (I - I_half) / 1;
		if (abs(R) < eps)
		{
			return I_half;
		}
		I = I_half;



	}


}

double helper_int_trap(Func f, double a, double b, int n)
{
	//I ~=  h/2 * (y_0 + 2 * sum(y_k,  k = 1..n-1) + y_n ) 
	// = h * ( sum(y_k,  k = 1..n-1) + (y_0 + y_n) / 2 )

	double h = (b - a) / n;
	double sum = ( f(a) + f(b) ) / 2.;
	//int n = ceil((b - a) / h); //округление вверх
	for (double x = a+h; x < b-h; x += h)
	{
		sum += f(x);
	}

	return sum * h;

}


double int_trap(Func f, double a, double b, double eps)		
{
	int n = 2;
	double I = helper_int_trap(f, a, b, n);
	while (42)
	{

		n *= 2;
		double I_half = helper_int_trap(f, a, b, n);
		double R = (I_half - I) / 3.;
		if (abs(R) < eps)
		{
			return I_half + R;
		}
		I = I_half;



	}


}


double helper_int_simpson(Func f, double a, double b, double n)
{
	//I ~=  h/3 * ( (y_0 + y_n) + 4*sum(y_2k-1, k=1..n) + 2*sum(y_2k, k=1..n-1))
	// = h/3 * (result1 + 4*result2 + 2*result3)


	double h = (b - a) / n;		//на какие интервалы делится отрезок
	double double_h = 2 * h;

	double result_1 = f(a) + f(b);

	double result_2 = 0.;
	for (double x = a + h; x < b; x += double_h)
	{
		result_2 += f(x);
	}

	double result_3 = 0.;
	for (double x = a + double_h; x < b; x += double_h)
	{
		result_3 += f(x);
	}

	return h * (result_1 + 4 * result_2 + 2 * result_3) / 3.;

}

double int_simpson(Func f, double a, double b, double eps)
{
	int n = 2;
	double I = helper_int_simpson(f, a, b, n);
	while (42)
	{

		n *= 2;
		double I_half = helper_int_simpson(f, a, b, n);
		double R = (I_half - I) / 15.;
		if (abs(R) < eps)
		{
			return I_half + R;
		}
		I = I_half;



	}


}




/*
double helper_ode_Euler(Func2 f, double x0, double y0, double x, int n)
{

	double y_i = y0;
	double h = (x - x0) / n;


	for (double x_i = x0; x_i <= x; x_i += h)
	{
		y_i = y_i + h * f(x_i, y_i);

	}

	return y_i;

}

double ode_Euler(Func2 f, double x0, double y0, double x, double eps)
{
	int n = 2;
	double I = helper_ode_Euler(f, x0, y0, x, n);
	while (42)
	{

		n *= 2;
		double I_half = helper_ode_Euler(f, x0, y0, x, n);
		double R = (I_half - I) / 1.;
		if (abs(R) < eps)
		{
			return I_half;
		}
		I = I_half;



	}


}
*/

vector<double> s_ode_euler(vector<Funcs*> f, double x, vector<double> y, double h)
{
	vector<double> result( y.size() );
	for (int k = 0; k < y.size(); k++)
	{
		result[k] = y[k] + h * f[k](x, y);
	}

	return result;
}

vector<double> s_ode_euler_koshi(vector<Funcs*> f, double x, vector<double> y, double h)
{
	vector<double> result(y.size());
	for (int k = 0; k < y.size(); k++)
	{
		result[k] = y[k] + h * f[k](x, y);
	}
	double h_half = h / 2;
	for (int k = 0; k < y.size(); k++)
	{
		result[k] = y[k] + h_half * ( f[k](x,y) + f[k](x+h, result));
	}
	return result;
}

//vector<double> s_ode_rk(vector<Funcs*> f, double x, vector<double> y, double h)


/*vector<double> s_euler_solver(vector<Funcs*> f, double x0, vector<double> y0, double x, int n)
{

	vector<double> y_i = y0;
	double h = (x - x0) / n;


	for (double x_i = x0; x_i <= x; x_i += h)
	{
		y_i = method(f, x_i, y_i, h);
	}

	return y_i;

}
*/




/*
vector<double> vode_Euler(vector<Funcs*> f, double x0, vector<double> y0, double x, double eps)
{
	int n = 2;
	vector<double> y = helper_vode_Euler(f, x0, y0, x, n);
	while (42)
	{

		n *= 2;
		vector<double> y_half = helper_vode_Euler(f, x0, y0, x, n);
			   
		double max_R = 0;
		for (int k = 0; k < y.size(); k++)
		{
			double R = abs(y_half[k] - y[k]);
			if (R > max_R)
			{
				max_R = R;
			}
		}

		if (max_R < eps)
		{
			for (int i = 0; i < y.size(); i++)
			{
				y_half[i] += ( (y_half[i] - y[i]) / (1.) );
			}
			return y_half;
		}

		y = y_half;

	}


}
*/




S_Ode_Solver::S_Ode_Solver(Method method, int order): method(method), order(order)			//конструктор
{
}

vector<double> S_Ode_Solver::operator()(vector<Funcs*> f, double x0, vector<double> y0, double x, double eps)
{

	int n = 2;
	vector<double> y = solve(f, x0, y0, x, n);
	const double magic_denom = pow(2, order) - 1;
	while (42)
	{

		n *= 2;
		vector<double> y_half = solve(f, x0, y0, x, n);

		double max_R = 0;
		for (int k = 0; k < y.size(); k++)
		{
			double R = abs(y_half[k] - y[k])/magic_denom;
			if (R > max_R)
			{
				max_R = R;
			}
		}
		//cout << max_R << endl;
		if (max_R < eps)
		{
			for (int i = 0; i < y.size(); i++)
			{
				y_half[i] += ((y_half[i] - y[i]) / (magic_denom));
			}
			return y_half;
		}

		y = y_half;

	}
}

vector<double> S_Ode_Solver::solve(vector<Funcs*> f, double x0, vector<double> y0, double x, int n)
{

	vector<double> y_i = y0;
	double h = (x - x0) / n;


	for (double x_i = x0; x_i <= x; x_i += h)
	{
		y_i = method(f, x_i, y_i, h);
	}

	return y_i;

}

double helper_diff_1_p1(Func f, double x, double h)
{
	return ( f(x+h) - f( x ) ) / h;
}

double diff_1_p1(Func f, double x, double eps)
{
	double h = eps;
	double df = helper_diff_1_p1(f, x, h);
	while (42)
	{

		h /= 2;
		double df_half = helper_diff_1_p1(f, x, h);
		double R = (df_half - df) / 1.;
		if (abs(R) < eps)
		{
			return df_half + R;
		}
		df = df_half;



	}


}

double helper_diff_1_p2(Func f, double x, double h)
{
	return (f(x + h) - f(x - h)) / (h*2);
}

double diff_1_p2(Func f, double x, double eps)
{
	double h = eps;
	double df = helper_diff_1_p2(f, x, h);
	const double magic_denom = (pow(2, 2) - 1);
	while (42)
	{

		h /= 2;
		double df_half = helper_diff_1_p2(f, x, h);
		double R = (df_half - df) / magic_denom;
		if (abs(R) < eps)
		{
			return df_half + R;
		}
		df = df_half;



	}


}

double helper_diff_2_p2(Func f, double x, double h)
{
	return (f(x - h) - f(x)*2 + f(x + h)) / (h * h);
}

double diff_2_p2(Func f, double x, double eps)
{
	double h = eps;
	double df = helper_diff_2_p2(f, x, h);
	const double magic_denom = (pow(2, 2) - 1);
	while (42)
	{

		h /= 2;
		double df_half = helper_diff_2_p2(f, x, h);
		double R = (df_half - df) / magic_denom;
		if (abs(R) < eps)
		{
			return df_half + R;
		}
		df = df_half;



	}


}


double minimize_dih(Func f, double a, double b, double eps)
{
	while (42)
	{
		double delta = eps / 2;
		double half = (a + b) / 2;
		double hme = half - delta;
		double hpe = half + delta;

		if (f(hme) < f(hpe))
		{
			b = hpe;
		}
		else
		{
			a = hme;
		}
		
		if ((b - a)/2 < eps)
		{
			return (a + b) / 2;
		}
	}

}

double minimize_gold(Func f, double a, double b, double eps)
{
	const double magic_l = ( 3 - sqrt(5) ) / 2;
	const double magic_r = ( sqrt(5) - 1 ) / 2;

	double x1 = a + (b - a) * magic_l;
	double x2 = a + (b - a) * magic_r;

	double f1 = f(x1);
	double f2 = f(x2);

	while (42)
	{
		if ((b - a) / 2 < eps)
		{
			return (a + b) / 2;
		}


		if (f1 <= f2)
		{
			b = x2;
			x2 = x1;
			x1 = a + b - x2;

			f2 = f1;
			f1 = f(x1);
		}
		else
		{
			a = x1;
			x1 = x2;
			x2 = (a + b) - x1;

			f1 = f2;
			f2 = f(x2);
		}


	}

}

/*
Matrix grad( FuncsV f, Matrix x )
{
	Matrix result( x.nRows() );
	result[0] = ();
}
*/
