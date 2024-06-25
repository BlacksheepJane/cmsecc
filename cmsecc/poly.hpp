#pragma once
#include "options.h"

class poly {
private:

public:
    vector<vector<mpz> > coef;
    //mpz e;
    int m, t;
    poly(int m = 0, int t = 0, mpz e = 0);

    mpz getCoef(int i, int j) const;

    void setCoef(int i, int j, mpz value);

    poly operator+(const poly& other) const;
    poly operator-(const poly& other) const;
    poly operator*(const poly& other) const;

    //poly mulf(mpz A, mpz B, mpz C, mpz X, mpz Y);
    poly powf(int exp);// ¿ìËÙÃÝ

    poly powf1(int exp, mpz A, mpz B, mpz C, mpz X, mpz Y);
    poly mulc(mpz c);
    poly xshift(int i, int k, mpz X, mpz e);

    poly yshift(int j, int k, mpz Y, mpz e);
    void print();
    //poly powy(int exp, mpz Y);
};