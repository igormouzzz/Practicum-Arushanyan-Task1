#pragma once
#include <stdio.h>
#include <math.h>
#include <numeric>
#include <algorithm>
#include <cmath> 

#include <thread>
#include <omp.h>
#include <typeinfo>
#include <chrono>

#include <iostream>
#include <fstream>

#include <vector>
#include <algorithm>

using namespace std;

class Vector;

double l1_norm(Vector& x);
double l2_norm(Vector& x);
double l2_norm_square(Vector& x);
double l_inf_norm(Vector& x);

double f1(Vector x);
double f2(Vector x);
double f3(Vector x);
double f4(Vector x);
double f5(Vector x);
double f6(Vector x);
double f7(Vector x);
double f8(Vector x);
double f9(Vector x);
double f10(Vector x);
double f11(Vector x);
double f12(Vector x);

double g(Vector x);