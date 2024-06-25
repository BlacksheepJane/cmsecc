//#pragma once
//#include"LLL.hpp"
//
//// 内积函数
//mpq_class dot(const Vectorq& v1, const Vectorq& v2) {
//    mpq_class result = 0;
//    for (size_t i = 0; i < v1.size(); ++i) {
//        result += v1[i] * v2[i];
//    }
//    return result;
//}
//
//mpq_class dot(const Vector& v1, const Vectorq& v2) {
//    mpq_class result = 0;
//    for (size_t i = 0; i < v1.size(); ++i) {
//        result += v1[i] * v2[i];
//    }
//    return result;
//}
//
////// 计算向量的范数
////mpq_class norm(const Vectorq& v) {
////    mpq_class norm_squared = dot(v, v);
////    mpf_class result = sqrt(norm_squared.get_d());
////    return mpq_class(result);
////}
//
////// 分数表示函数
////string f_to_tex(const mpq_class& value) {
////    ostringstream out;
////    out << value;
////    return out.str();
////}
//Vectorq initq(const Vector v) {
//    int n = v.size();
//    Vectorq result;
//    for (int i = 0; i < n; i++) {
//        mpq q(v[i], 1);
//        result.push_back(q);
//    }
//    return result;
//}
//// Gram-Schmidt 正交化过程
//vector<Vectorq> gramSchmidt(const vector<Vector>& B) {
//    size_t n = B.size();
//    vector<Vectorq> B_star(n, Vectorq(B[0].size()));
//    for (size_t i = 0; i < n; ++i) {
//        //B_star[i] = B[i];
//        B_star[i] = initq(B[i]);
//        for (size_t j = 0; j < i; ++j) {
//            mpq_class mu = dot(B[i], B_star[j]) / dot(B_star[j], B_star[j]);
//            for (size_t k = 0; k < B[i].size(); ++k) {
//                B_star[i][k] -= mu * B_star[j][k];
//            }
//        }
//    }
//    return B_star;
//}
//
//mpq_class mu(Vector B, Vectorq B_star) {
//    if (dot(B_star, B_star) != 0)
//        return (dot(B, B_star) / dot(B_star, B_star));
//    return 0;
//}
//mpq_class mu(Vectorq B, Vectorq B_star) {
//    if (dot(B_star, B_star) != 0)
//        return (dot(B, B_star) / dot(B_star, B_star));
//    return 0;
//}
//// LLL算法实现
//void LLL(vector<Vector>& B, mpz delta) {//B中向量必须线性无关
//    size_t n = B.size();
//    vector<Vectorq> B_star = gramSchmidt(B);
//
//    size_t k = 1;
//    while (k < n) {
//        for (int j = k - 1; j >= 0; --j) {
//            mpq_class mu_ij = mu(B[k], B_star[j]);
//            if (abs(mu_ij.get_d()) > 0.5) {
//                mpz_class r = round(mu_ij.get_d());
//                for (size_t i = 0; i < B[k].size(); ++i) {
//                    B[k][i] -= r * B[j][i];
//                }
//                B_star = gramSchmidt(B);
//            }
//        }
//
//        if (dot(B_star[k], B_star[k]) >= (delta - mu(B_star[k], B_star[k - 1]) * mu(B_star[k], B_star[k - 1])) * dot(B_star[k - 1], B_star[k - 1])) {
//            ++k;
//        }
//        else {
//            swap(B[k], B[k - 1]);
//            B_star = gramSchmidt(B);
//            k = max(k - 1, size_t(1));
//        }
//    }
//}
//
//// 打印向量
//void printVector(const Vector& v) {
//    for (size_t i = 0; i < v.size(); ++i) {
//        cout << v[i];
//        if (i < v.size() - 1)
//            cout << ", ";
//    }
//    cout << endl;
//}
//
//// 打印矩阵
//void printMat(const vector<Vector>& B) {
//    for (const Vector& v : B) {
//        printVector(v);
//    }
//}

#pragma once
#include"LLL.hpp"

// 内积函数
mpq_class dot(Vectorq& v1, Vectorq& v2) {
    mpq_class result = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

mpq_class dot(Vector& v1, Vectorq& v2) {
    mpq_class result = 0;
    for (size_t i = 0; i < v1.size(); ++i) {
        result += v1[i] * v2[i];
    }
    return result;
}

Vectorq initq(const Vector v) {
    int n = v.size();
    Vectorq result;
    for (int i = 0; i < n; i++) {
        mpq q(v[i], 1);
        result.push_back(q);
    }
    return result;
}
mpq_class mu(Vector& B, Vectorq& B_star) {
    if (dot(B_star, B_star) != 0)
        return (dot(B, B_star) / dot(B_star, B_star));
    return 0;
}
mpq_class mu(Vectorq& B, Vectorq& B_star) {
    if (dot(B_star, B_star) != 0)
        return (dot(B, B_star) / dot(B_star, B_star));
    return 0;
}
// Gram-Schmidt 正交化过程
vector<Vectorq> gramSchmidt(vector<Vector>& B) {
    size_t n = B.size();
    vector<Vectorq> B_star(n, Vectorq(B[0].size()));
    for (size_t i = 0; i < n; ++i) {
        //B_star[i] = B[i];
        B_star[i] = initq(B[i]);
        for (size_t j = 0; j < i; ++j) {
            mpq_class mu_ij = mu(B[i], B_star[j]);
            for (size_t k = 0; k < B[i].size(); ++k) {
                B_star[i][k] -= mu_ij * B_star[j][k];
            }
        }
    }
    return B_star;
}

// 打印向量
void printVectorq(const Vectorq& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        cout << v[i];
        if (i < v.size() - 1)
            cout << ", ";
    }
    cout << endl;
}

// 打印矩阵
void printMatq(const vector<Vectorq>& B) {
    for (const Vectorq& v : B) {
        printVectorq(v);
    }
    cout << endl;
}
// LLL算法实现
void LLL(vector<Vector>& B, mpq delta) {//B中向量必须线性无关
    size_t n = B.size();
    cout << n << endl;
    vector<Vectorq> B_star = gramSchmidt(B);
    //printMat(B);
    size_t k = 1;
    int flag = 1;
    while (k < n) {
        for (int j = k - 1; j >= 0; --j) {
            mpq_class mu_ij = mu(B[k], B_star[j]);
            //cout << mu_ij << endl << endl;
            if (abs(mu_ij.get_d()) > 0.5) {
                /*cout << mu_ij << endl;*/
                mpz_class r = round(mu_ij.get_d());
                //cout << r << endl;
                for (size_t i = 0; i < B[k].size(); ++i) {
                    B[k][i] -= r * B[j][i];
                }
                //printMat(B);
                B_star = gramSchmidt(B);
                //printMatq(B_star);
            }
        }
        
        //printMatq(B_star);
        //cout << delta << endl;
        //cout << delta - mu(B_star[k], B_star[k - 1]) << endl;
        //cout << dot(B_star[k - 1], B_star[k - 1]) << endl;
        //cout << endl;
        if (dot(B_star[k], B_star[k]) >= (delta - mu(B_star[k], B_star[k - 1]) * mu(B_star[k], B_star[k - 1])) * dot(B_star[k - 1], B_star[k - 1])) {
            //cout << dot(B_star[k], B_star[k]) << endl;
            //cout << (delta - mu(B_star[k], B_star[k - 1]) * mu(B_star[k], B_star[k - 1])) * dot(B_star[k - 1], B_star[k - 1]) << endl;
            //cout << endl;
            ++k;
            //cout << "t" << endl;
        }
        else {
            //cout << dot(B_star[k], B_star[k]) << endl;
            //cout << (delta - mu(B_star[k], B_star[k - 1]) * mu(B_star[k], B_star[k - 1])) * dot(B_star[k - 1], B_star[k - 1]) << endl;
            //cout << endl;
            swap(B[k], B[k - 1]);
            B_star = gramSchmidt(B);
            k = max(k - 1, size_t(1));
            //printMat(B);
            //cout << "f" << endl;
        }
        if (k == flag) {
            cout << k << endl;
            flag++;
        }
    }
}

// 打印向量
void printVector(const Vector& v) {
    for (size_t i = 0; i < v.size(); ++i) {
        cout << mpf_class(v[i].get_str()) << ' ';
        if (i < v.size() - 1)
            cout << ", ";
    }
    cout << endl;
}

// 打印矩阵
void printMat(const vector<Vector>& B) {
    for (const Vector& v : B) {
        printVector(v);
    }
    cout << endl;
}

