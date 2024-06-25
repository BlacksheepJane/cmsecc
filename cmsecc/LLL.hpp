#pragma once
#include"options.h"
#include <algorithm>
#include <sstream>

// �ڻ�����
mpq_class dot(const Vectorq& v1, const Vectorq& v2);

mpq_class dot(const Vector& v1, const Vectorq& v2);

//// ���������ķ���
//mpq_class norm(const Vectorq& v);

//// ������ʾ����
//string f_to_tex(const mpq_class& value);
Vectorq initq(const Vector v);
// Gram-Schmidt ����������
vector<Vectorq> gramSchmidt(const vector<Vector>& B);

// LLL�㷨ʵ��
void LLL(vector<Vector>& B, mpq delta = 0.75);//B���������޹�

// ��ӡ����
void printVector(const Vector& v);

// ��ӡ����
void printMat(const vector<Vector>& B);