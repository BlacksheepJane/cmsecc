#pragma once
#include <iostream>
#include <vector>
#include <gmp.h>
#include <gmpxx.h>

using namespace std;
#define mpz mpz_class
#define mpq mpq_class
#define mpf mpf_class
typedef vector<mpq_class> Vectorq;
typedef vector<mpf_class> Vectorf;
typedef vector<mpz_class> Vector;
