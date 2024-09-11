#include "Header.h"
#include "Matrix.h"

namespace Task1		// 1 1
{
	double StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk);
	double NewtonMethod(double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon);
	double NewtonMethod2(double s, double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon);
	void FirstMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x);
	void SecondMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x);
	void main_func();
}
