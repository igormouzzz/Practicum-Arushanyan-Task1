#include "Header.h"
#include "Matrix.h"

namespace Task1		// 1 1
{
	double Derivative(const std::function<double(double)>& phi, double t, double tau);

	double StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk);
	double NewtonMethod(double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon);
	void FirstMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x);
	void SecondMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x);
	void main_func();
}
