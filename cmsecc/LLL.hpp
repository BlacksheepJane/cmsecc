#pragma once
#include"options.h"
#include <algorithm>
#include <sstream>

// 内积函数
mpq_class dot(const Vectorq& v1, const Vectorq& v2);

mpq_class dot(const Vector& v1, const Vectorq& v2);

//// 计算向量的范数
//mpq_class norm(const Vectorq& v);

//// 分数表示函数
//string f_to_tex(const mpq_class& value);
Vectorq initq(const Vector v);
// Gram-Schmidt 正交化过程
vector<Vectorq> gramSchmidt(const vector<Vector>& B);

// LLL算法实现
void LLL(vector<Vector>& B, mpq delta = 0.75);//B必须线性无关

// 打印向量
void printVector(const Vector& v);

// 打印矩阵
void printMat(const vector<Vector>& B);