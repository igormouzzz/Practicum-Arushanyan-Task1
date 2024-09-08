#include "Header.h"
#include "Vector.h"

using row = vector<double>;

class Matrix
{
protected:
	vector<row> a;
	int N, M;
public:
	Matrix();
	Matrix(int N, int M);
	Matrix(vector<row> rows);
	Matrix(const Matrix& b);
	size_t GetN();
	size_t GetM();
	void SetValue(int i, int j, double value);
	Matrix& operator=(const Matrix& b);
	Matrix operator+(const Matrix& b);
	Matrix operator-(const Matrix& b);
	Matrix operator*(double k);
	Matrix operator*(const Matrix& b);
	vector<double> operator*(vector<double>& b);
	Vector operator*(Vector& b);
	Matrix T();
	Matrix Inv2();
	Matrix Inversed(Matrix& Inv);
	
	friend class Vector;
	friend ostream& operator<<(ostream& cout, const Matrix& b);
};
