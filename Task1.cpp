#include "Task1.h"

double Task1::StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk)
{
	auto phi = [&](double t) { return f(xk + hk * t); };
	double lambda = 0.5, mu = 2;
	double alpha, beta;
	double alpha_k;

	double step_max = 2;

	beta = 0.1;
	alpha = beta;

	if (f(xk + hk * alpha) <= f(xk))
	{
		alpha = mu * beta;
		while (!(phi(alpha) <= phi(beta)) && alpha <= step_max)
		{
			alpha = mu * beta;
			cout << 1 << endl;
			cout << phi(alpha) << "\t" << phi(beta) << endl;
		}
		alpha_k = alpha;
	}
	else
	{
		cout << f(xk + hk * alpha) << "\t" << f(xk) << endl;
		alpha = lambda * beta;
		while (!(f(xk + hk * alpha) < f(xk)))
		{
			alpha = lambda * beta;
			cout << 2 << endl;
			cout << f(xk + hk*alpha) << "\t" << f(xk) << endl;
		}
		alpha_k = alpha;
	}

	return alpha_k;
}

void Task1::RunCalculation(double(*f)(Vector x), const double epsilon, Vector& x0)
{
	size_t n = x0.GetSize();
	Vector x(n);
	Vector xk_1(n), xk(n);
	xk_1 = x0;
	vector<double> df(n);
	Vector Df(n);
	Vector hk_1(n);
	Vector diff(n);
	double alpha_k_1;
	double weak_conv;
	double (*norm_sqr)(Vector& x);
	norm_sqr = l2_norm_square;

	const double tau = 0.1 * sqrt(epsilon);

	//vector<double> d1({ tau,0 }), d2({ 0,tau });
	vector<vector<double>> d(n, vector<double>(n));
	for (int i = 0; i < n; i++) d[i][i] = tau;

	do
	{
		for (int i = 0; i < n; i++)
		{
			df[i] = (f(xk_1 + d[i]) - f(xk_1 - d[i])) / (2 * tau);
		}
		//df[0] = (f(xk_1 + d1) - f(xk_1 - d1)) / (2 * tau);
		//df[1] = (f(xk_1 + d2) - f(xk_1 - d2)) / (2 * tau);

		Df = df;
		hk_1 = Df * (-1);
		cout << "xk_1 = " << xk_1;
		cout << "hk_1 = "<<hk_1;
		alpha_k_1 = StepSplittingMethod(f, xk_1, hk_1);
		cout << "alpha_k_1 = " << alpha_k_1 << endl;
		xk = xk_1 + hk_1 * alpha_k_1;
		diff = xk - xk_1;
		cout << "diff = " << diff;
		weak_conv = std::abs(f(xk) - f(xk_1));
		cout << "weak_conv = " << weak_conv << endl << endl;

		xk_1 = xk;
	}
	while (norm_sqr(diff) > epsilon || weak_conv > epsilon || norm_sqr(Df) > epsilon);

	x = xk;
	cout << "Answer: x = ";
	cout << x << endl;
}

void Task1::main_func()
{
	double (*f)(Vector x);
	f = f3;
	double epsilon = 1e-12;
	vector<double> x00(2);
	cout << "Enter the initial point: " << endl;
	for (int i = 0; i < x00.size(); i++) cin >> x00[i];

	Vector x0(x00);

	RunCalculation(f, epsilon, x0);

	cout << "Program has been completed successfully!" << endl;
}