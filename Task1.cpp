#include "Task1.h"

double Task1::Derivative(const std::function<double(double)>& phi, double t, double tau)
{
	return (phi(t + tau) - phi(t - tau)) / (2 * tau);
}

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

double Task1::NewtonMethod(double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon)
{
	auto phi = [&](double t) { return f(xk + hk * t); };
	double tau = 0.1 * sqrt(epsilon);
	double a = -10, b = 10;
	double d_phi_a = Derivative(phi, a, tau), d_phi_b = Derivative(phi, b, tau);
	if (std::fabs(d_phi_a) < epsilon) return a;
	if (std::fabs(d_phi_b) < epsilon) return b;
	if (d_phi_a > 0 && d_phi_b > 0) return a;
	if (d_phi_a < 0 && d_phi_b < 0) return b;

	double an = a, bn = b;
	double d_phi_an = Derivative(phi, an, tau), d_phi_bn = Derivative(phi, bn, tau);
	double cn = (phi(an) - phi(bn) + d_phi_bn * bn - d_phi_an * an) / (d_phi_bn - d_phi_an);

	double d_phi_cn = Derivative(phi, cn, tau);

	if (std::fabs(d_phi_cn) < epsilon) return cn;
	double an1, bn1;

	while (!(std::fabs(d_phi_cn) < epsilon))
	{
		if (d_phi_cn < 0)
		{
			an1 = cn; bn1 = bn; cout << "<0" << endl;
		}
		else
		{
			an1 = an; bn1 = cn; cout << ">0" << endl;
		}

		an = an1; bn = bn1;
		cout << an << "\t" << bn << endl;
		d_phi_an = Derivative(phi, an, tau), d_phi_bn = Derivative(phi, bn, tau);
		if (std::fabs(d_phi_an) < epsilon) return an;
		if (std::fabs(d_phi_bn) < epsilon) return bn;
		if (d_phi_an > 0 && d_phi_bn > 0) return an;
		if (d_phi_an < 0 && d_phi_bn < 0) return bn;

		cn = (phi(an) - phi(bn) + d_phi_bn * bn - d_phi_an * an) / (d_phi_bn - d_phi_an);
		cout << "cn = " << cn << endl;

		d_phi_cn = Derivative(phi, cn, tau);
		cout << "d_phi_cn = " << d_phi_cn << endl << endl;
	}
	return cn;
}

void Task1::FirstMethod(double(*f)(Vector x), const double epsilon, Vector& x0, Vector& x)
{
	size_t n = x0.GetSize();
	Vector xk_1(n), xk(n);
	xk_1 = x0;
	vector<double> df(n);
	Vector Df(n);
	Vector hk_1(n);
	Vector diff_iter(n);
	double alpha_k_1;
	double diff_f;
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
		diff_iter = xk - xk_1;
		cout << "diff_iter = " << diff_iter;
		diff_f = std::abs(f(xk) - f(xk_1));
		cout << "diff_f = " << diff_f << endl << endl;

		xk_1 = xk;
	}
	while (norm_sqr(diff_iter) > epsilon || diff_f * diff_f > epsilon || norm_sqr(Df) > epsilon);

	x = xk;
	cout << "First method result: x = ";
	cout << x << endl;
}

void Task1::SecondMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x)
{
	size_t n = x0.GetSize();
	Vector xk_1(n), xk(n);
	xk_1 = x0;
	vector<double> df(n);
	Vector Df(n);
	Vector hk_1(n);
	Vector diff_iter(n);
	double alpha_k_1;
	double diff_f;
	double (*norm_sqr)(Vector& x);
	norm_sqr = l2_norm_square;
	const double tau = 0.1 * sqrt(epsilon);

	vector<vector<double>> vectorized_hessian(2); for (int i = 0; i < 2; i++) vectorized_hessian[i].resize(2);
	Matrix He(2, 2);
	//double fy, fy4;
	//double fyx4, fy4x;

	//vector<double> d1({ tau,0 }), d2({ 0,tau });
	vector<vector<double>> d(n, vector<double>(n));
	for (int i = 0; i < n; i++) d[i][i] = tau;

	auto fy = [&](Vector p) { return (f(p + d[1]) - f(p)) / tau; };
	auto fy4 = [&](Vector p) { return (f(p) - f(p - d[1])) / tau; };
	double fyx4, fy4x;

	vector<double> d2f(3);
	do
	{
		for (int i = 0; i < n; i++)
		{
			df[i] = (f(xk_1 + d[i]) - f(xk_1 - d[i])) / (2 * tau);
		}

		Df = df;
		for (int i = 0; i < 2; i++)
		{
			d2f[i] = (f(xk_1 + d[i]) - 2 * f(xk_1) + f(xk_1 - d[i])) / (tau * tau);
		}

		//fy = ( f(xk_1 +d[1]) - f(xk_1) ) / tau;
		//fy4 = ( f(xk_1) - f(xk_1 -d[1]) ) / tau;
		fyx4 = (fy(xk_1) - fy(xk_1 - d[0])) / tau;
		fy4x = (fy4(xk_1 + d[0]) - fy4(xk_1)) / tau;

		d2f[2] = 0.5 * (fyx4 + fy4x);

		vectorized_hessian = { { d2f[0], d2f[2]}, {d2f[2], d2f[1]} };
		He = Matrix(vectorized_hessian);
		cout << He << endl;
		hk_1 = He.Inv2() * Df * (-1);

		cout << "xk_1 = " << xk_1;
		cout << "hk_1 = " << hk_1;
		alpha_k_1 = NewtonMethod(f, xk_1, hk_1, epsilon);
		cout << "alpha_k_1 = " << alpha_k_1 << endl;
		xk = xk_1 + hk_1 * alpha_k_1;
		diff_iter = xk - xk_1;
		cout << "diff_iter = " << diff_iter;
		diff_f = std::abs(f(xk) - f(xk_1));
		cout << "diff_f = " << diff_f << endl << endl;

		xk_1 = xk;
	} while (norm_sqr(diff_iter) > epsilon*epsilon || diff_f > epsilon || norm_sqr(Df) > epsilon*epsilon);

	x = xk;
	cout << "Answer: x = ";
	cout << x << endl;
}

void Task1::main_func()
{
	double (*f)(Vector x);
	f = f4;
	double epsilon = 1e-8;
	vector<double> x00(2);
	cout << "Enter the initial point: " << endl;
	for (int i = 0; i < x00.size(); i++) cin >> x00[i];

	Vector x0(x00); Vector x(x0.GetSize());

	FirstMethod(f, epsilon, x0, x);
	SecondMethod(f, epsilon, x, x);

	cout << "Program has been completed successfully!" << endl;
}