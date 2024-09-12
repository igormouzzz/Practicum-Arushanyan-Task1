#include "Task1.h"


double Task1::StepSplittingMethod(double (*f)(Vector x), Vector& xk, Vector& hk)
{
	auto phi = [&](double t) { return f(xk + hk * t); };
	double lambda = 0.5, mu = 2;
	double alpha, beta;
	double alpha_k;

	double alpha_max = 10;

	beta = 0.1;
	alpha = beta;

	if (f(xk + hk * alpha) < f(xk))
	{
		//cout << phi(alpha) << "\t" << phi(0) << endl;
		while (phi(alpha) < phi(beta) && alpha <= alpha_max);
		{
			alpha_k = alpha;
			alpha *= mu;
			cout << 1 << endl;
		} 
	}
	else
	{
		do
		{
			alpha_k = alpha;
			alpha *= lambda;
			cout << 2 << endl;
		} 
		while (!(f(xk + hk * alpha) < f(xk)));
	}
	return alpha_k;
}

double Task1::NewtonMethod(double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon)
{
	auto phi = [&](double t) { return f(xk + hk * t); };
	const double tau = 0.1 * sqrt(epsilon);
	double a = -10, b = 10;
	auto Derivative = [&](double t) -> double { return (phi(t + tau) - phi(t - tau)) / (2 * tau); };
	double d_phi_a = Derivative(a), d_phi_b = Derivative(b);
	if (std::fabs(d_phi_a) < epsilon) return a;
	if (std::fabs(d_phi_b) < epsilon) return b;
	if (d_phi_a > 0 && d_phi_b > 0) return a;
	if (d_phi_a < 0 && d_phi_b < 0) return b;
	cout << 100 << endl;
	double an = a, bn = b;
	double d_phi_an = Derivative(an), d_phi_bn = Derivative(bn);
	double cn = (phi(an) - phi(bn) + d_phi_bn * bn - d_phi_an * an) / (d_phi_bn - d_phi_an);
	cout << "cn = " << cn << endl;
	double d_phi_cn = Derivative(cn);
	cout << "d_phi_cn = " << d_phi_cn << endl;
	if (std::fabs(d_phi_cn) < epsilon) return cn;
	double an1, bn1;
	cout << 200 << endl;
	do
	{
		if (d_phi_cn < 0)
		{
			an1 = cn; bn1 = bn;
		}
		else
		{
			an1 = an; bn1 = cn;
		}

		an = an1; bn = bn1;
		cout << an << "\t" << bn << endl;
		d_phi_an = Derivative(an), d_phi_bn = Derivative(bn);
		if (std::fabs(d_phi_an) < epsilon) return an;
		if (std::fabs(d_phi_bn) < epsilon) return bn;
		if (d_phi_an > 0 && d_phi_bn > 0) return an;
		if (d_phi_an < 0 && d_phi_bn < 0) return bn;
		cout << 300 << endl;
		cn = (phi(an) - phi(bn) + d_phi_bn * bn - d_phi_an * an) / (d_phi_bn - d_phi_an);
		cout << "cn = " << cn << endl;

		d_phi_cn = Derivative(cn);
		cout << "d_phi_cn = " << d_phi_cn << endl << endl;
	}
	while (!(std::fabs(d_phi_cn) < epsilon));
	return cn;
}

double Task1::NewtonMethod2(double s, double (*f)(Vector x), Vector& xk, Vector& hk, double epsilon)
{
	auto phi = [&](double t) { return f(xk + hk * t); };
	const double tau = 0.1 * sqrt(epsilon);
	auto d_phi = [&](double t) { return (phi(t + tau) - phi(t - tau)) / (2 * tau); };

	while (fabs(phi(s)) > epsilon)
	{
		cout << "Po4erevin" << endl;
		s = s - phi(s) / d_phi(s);
	}
	return s;
}

void Task1::FirstMethod(double(*f)(Vector x), const double epsilon, Vector& x0, Vector& x)
{
	size_t n = x0.GetSize();
	Vector xk(n), xk1(n);
	xk = x0;
	vector<double> df(n);
	Vector Df(n);
	Vector hk(n);
	Vector diff_iter(n);
	double alpha_k;
	double diff_f;
	double (*norm_sqr)(Vector& x);
	norm_sqr = l2_norm_square;
	double difference_norm_sqr, grad_norm_sqr;

	const double tau = 0.1 * sqrt(epsilon);

	vector<vector<double>> d(n, vector<double>(n));
	for (int i = 0; i < n; i++) d[i][i] = tau;

	int iterations = 0;
	do
	{
		for (int i = 0; i < n; i++)
		{
			df[i] = (f(xk + d[i]) - f(xk - d[i])) / (2 * tau);
		}

		Df = df;
		hk = Df * (-1);
		alpha_k = StepSplittingMethod(f, xk, hk);
		xk1 = xk + hk * alpha_k;
		diff_iter = xk1 - xk;
		difference_norm_sqr = norm_sqr(diff_iter);
		diff_f = std::fabs(f(xk1) - f(xk));
		grad_norm_sqr = norm_sqr(Df);

		cout << "Df = " << Df;
		cout << "alpha_k = " << alpha_k << endl;
		cout << "xk = " << xk;
		cout << "xk+1 = " << xk1;
		cout << difference_norm_sqr << "\t" << diff_f << "\t" << grad_norm_sqr << endl << endl;

		xk = xk1;
		iterations++;
	}
	while (difference_norm_sqr > epsilon || diff_f * diff_f > epsilon || grad_norm_sqr > epsilon);

	x = xk1;
	cout << "Iterations: " << iterations << endl << endl;
	cout << "First method result: x = ";
	cout << x << endl;
}

void Task1::SecondMethod(double (*f)(Vector x), const double epsilon, Vector& x0, Vector& x)
{
	size_t n = x0.GetSize();
	Vector xk(n), xk1(n);
	xk = x0;
	vector<double> grad_f(n);
	Vector Df(n);
	Vector hk(n);
	Vector diff_iter(n);
	double alpha_k;
	double diff_f;
	double (*norm_sqr)(Vector& x);
	norm_sqr = l2_norm_square;
	double difference_norm_sqr, grad_norm_sqr;

	const double tau = 0.1 * sqrt(epsilon);

	Matrix He(n, n);

	vector<vector<double>> d(n, vector<double>(n));
	for (int i = 0; i < n; i++) d[i][i] = tau;

	vector<vector<double>> der2f(n, vector<double>(n));

	vector<std::function<double(Vector)>> df(n);
	vector<std::function<double(Vector)>> df4(n);
	double fij4, fi4j;
	
	int iterations = 0;
	do
	{
		for (int i = 0; i < n; i++)
		{
			grad_f[i] = (f(xk + d[i]) - f(xk - d[i])) / (2 * tau);
		}

		Df = grad_f;
		
		for (int i = 0; i < n; i++)
		{
			df[i] = [&](Vector p) { return (f(p + d[i]) - f(p)) / tau; };
			df4[i] = [&](Vector p) { return (f(p) - f(p - d[i])) / tau; };
			der2f[i][i] = (f(xk + d[i]) - 2 * f(xk) + f(xk - d[i])) / (tau * tau);
			for (int j = i + 1; j < n; j++)
			{
				fij4 = (df[i](xk) - df[i](xk - d[j])) / tau;
				fi4j = (df4[i](xk + d[j]) - df4[i](xk)) / tau;
				der2f[i][j] = der2f[j][i] = 0.5 * (fij4 + fi4j);
			}
		}
		He.Move(der2f);
		cout << He;
		hk = He.Gauss(Df) * (-1);
		//alpha_k = NewtonMethod(f, xk, hk, epsilon);
		alpha_k = NewtonMethod2(0.5, f, xk, hk, epsilon);
		xk1 = xk + hk * alpha_k;
		diff_iter = xk1 - xk;
		difference_norm_sqr = norm_sqr(diff_iter);
		diff_f = std::abs(f(xk1) - f(xk));
		grad_norm_sqr = norm_sqr(Df);

		cout << "xk = " << xk << endl;
		cout << "alpha_k = " << alpha_k << endl;
		cout << "Df = " << Df << endl;

		xk = xk1;
		iterations++;
	} while (difference_norm_sqr > epsilon*epsilon || diff_f > epsilon || grad_norm_sqr > epsilon*epsilon);

	x = xk1;
	cout << "Iterations: " << iterations << endl << endl;
	cout << "Answer: x = ";
	cout << x << endl;
}

void Task1::main_func()
{
	double (*f)(Vector x);
	f = f4;
	double epsilon = 1e-4;
	vector<double> x00(2);
	cout << "Enter the initial point: " << endl;
	for (int i = 0; i < x00.size(); i++) cin >> x00[i];

	Vector x0(x00); Vector x(x0.GetSize());

	FirstMethod(f, epsilon, x0, x);
	system("pause");
	SecondMethod(f, epsilon, x, x);

	cout << "Program has been completed successfully!" << endl;
}