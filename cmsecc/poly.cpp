#include "poly.hpp"



poly::poly(int m, int t, mpz e) : m(m), t(t) {
        coef.resize(m + 1, vector<mpz>(m + t + 1, 0));

    }

mpz poly::getCoef(int i, int j) const {
        if (i >= 0 && i <= m && j >= 0 && j <= m + t) {
            return coef[i][j];
        }
        return 0;
    }

void poly::setCoef(int i, int j, mpz value) {
        if (i >= 0 && i <= m && j >= 0 && j <= m + t) {
            coef[i][j] = value;
        }
    }

poly poly::operator+(const poly& other) const {
        poly result(m, t);

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                mpz sum = getCoef(i, j) + other.getCoef(i, j);
                result.setCoef(i, j, sum);
            }
        }

        return result;
    }

poly poly::operator-(const poly& other) const {
        poly result(m, t);

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                mpz diff = getCoef(i, j) - other.getCoef(i, j);
                result.setCoef(i, j, diff);
            }
        }

        return result;
    }
poly poly::operator*(const poly& other) const {
        poly result(m, t);

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                for (int k = 0; k <= m; ++k) {
                    for (int l = 0; l <= m + t; ++l) {
                        mpz product = getCoef(i, j) * other.getCoef(k, l);
                        result.setCoef(i + k, j + l, result.getCoef(i + k, j + l) + product);
                    }
                }
            }
        }
        return result;
    }



    //poly poly::mulf(mpz A, mpz B, mpz C, mpz X, mpz Y) {
    //    poly p1(m, t), p2(m, t), p3(m, t), p4(m, t);

    //    for (int i = 0; i <= m; ++i) {
    //        for (int j = 0; j <= m + t; ++j) {
    //            p1.setCoef(i + 1, j + 1, X * Y * coef[i][j]);
    //        }
    //    }

    //    for (int i = 0; i <= m; ++i) {
    //        for (int j = 0; j <= m + t; ++j) {
    //            p2.setCoef(i + 1, j, A * X * coef[i][j]);
    //        }
    //    }

    //    for (int i = 0; i <= m; ++i) {
    //        for (int j = 0; j <= m + t; ++j) {
    //            p3.setCoef(i, j + 1, B * Y * coef[i][j]);
    //        }
    //    }

    //    for (int i = 0; i <= m; ++i) {
    //        for (int j = 0; j <= m + t; ++j) {
    //            p4.setCoef(i, j, C * coef[i][j]);
    //        }
    //    }

    //    return p1 + p2 + p3 + p4;
    //}
poly poly::powf(int exp) {//¿ìËÙÃÝ
        poly result = *this;
        exp--;
        poly base = *this;
        while (exp > 0) {
            if (exp % 2 == 1) {
                result = result * base;
            }
            base = base * base;
            exp /= 2;
        }
        return result;
    }

poly poly::powf1(int exp, mpz A, mpz B, mpz C, mpz X, mpz Y) {//this*f^(exp-1)
    poly result = *this;
    for (int i = 1; i < exp; i++) {
        poly p1(m, t), p2(m, t), p3(m, t), p4(m, t);

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                p1.setCoef(i + 1, j + 1, X * Y * result.coef[i][j]);
            }
        }

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                p2.setCoef(i + 1, j, A * X * result.coef[i][j]);
            }
        }

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                p3.setCoef(i, j + 1, B * Y * result.coef[i][j]);
            }
        }

        for (int i = 0; i <= m; ++i) {
            for (int j = 0; j <= m + t; ++j) {
                p4.setCoef(i, j, C * result.coef[i][j]);
            }
        }
        result = p1 + p2 + p3 + p4;
    }
    return result;
}
poly poly::mulc(mpz c) {
        poly result(m, t);
        for (int i = 0; i <= m; i++)
            for (int j = 0; j <= m + t; j++)
                result.setCoef(i, j, c * coef[i][j]);
        return result;
    }
poly poly::xshift(int i, int k, mpz X, mpz e) {
        poly result(m, t);

        poly fp = powf(k);//or powf1

        mpz ep;
        mpz_pow_ui(ep.get_mpz_t(), e.get_mpz_t(), m - k);

        mpz xp;
        mpz_pow_ui(xp.get_mpz_t(), X.get_mpz_t(), i);

        for (int ii = 0; ii <= m; ++ii) {
            for (int jj = 0; jj <= m + t; ++jj) {
                result.setCoef(ii + i, jj, ep * xp * fp.coef[ii][jj]);
            }
        }
        return result;
    }

poly poly::yshift(int j, int k, mpz Y, mpz e) {
    poly result(m, t);

    poly fp = powf(k);//or powf1

    mpz ep;
    mpz_pow_ui(ep.get_mpz_t(), e.get_mpz_t(), m - k);

    mpz yp;
    mpz_pow_ui(yp.get_mpz_t(), Y.get_mpz_t(), j);

    for (int ii = 0; ii <= m; ++ii) {
        for (int jj = 0; jj <= m + t; ++jj) {
            result.setCoef(ii, jj + j, ep * yp * fp.coef[ii][jj]);
        }
    }
    return result;
}
void poly::print() {
    for (int i = 0; i <= m; ++i) {
        for (int j = 0; j <= m + t; ++j) {
            cout << coef[i][j] << ' ';
        }
        cout << endl;
    }
    cout << endl;
}
    //poly poly::powy(int exp, mpz Y) {
    //    mpz yp;
    //    mpz_pow_ui(yp.get_mpz_t(), Y.get_mpz_t(), exp);
    //    poly result(m, t);
    //    for (int i = 0; i <= m; ++i) {
    //        for (int j = 0; j <= m + t; ++j) {
    //            if (i + exp <= m + 1) {
    //                result.setCoef(i, j + exp, yp * coef[i][j]);
    //            }
    //        }
    //    }
    //    return result;
    //}