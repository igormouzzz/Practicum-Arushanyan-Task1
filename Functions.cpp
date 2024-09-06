#include "Header.h"
#include "Vector.h"

double f1(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return x1 * x1 - x2 * x2;
}
double f2(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return x1 * x1 + x2 * x2 * x2 * x2;
}
double f3(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return (x1 - 1) * (x1 - 1) + (x2 + 1) * (x2 + 1);
}
double f4(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return x1 * x1 - 2 * x1 * x2 + x2 * x2 + 100;
}
double f5(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1], x3 = x.v[2];
	return x1 * x1 + x2 * x2 + x3 * x3;
}
double f6(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return 0;
}
double f7(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return 0;
}
double f8(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return 0;
}
double f9(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return 0;
}
double f10(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return 0;
}
double f11(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return 0;
}
double f12(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1], x3 = x.v[2];
	return x1*x1 + 5*x2*x2 + 3*x3*x3 + 4*x1*x2 - 2*x1*x3 - 2*x2*x3;
}

double g(Vector x)
{
	double x1 = x.v[0], x2 = x.v[1];
	return -cos(x1 * x1 + x2 * x2);
}