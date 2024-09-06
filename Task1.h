#include "Header.h"
#include "Vector.h"

namespace Task1		// 1 1
{
	double StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk);
	void RunCalculation(double (*f)(Vector x), const double epsilon, Vector& x0);
	void main_func();
}
