#include "T5.hpp"

mpf vectorMagnitude(const Vector& v) {
	mpf_class sum = 0.0;
	for (const auto& value : v) {
		sum += value * value;
	}
	return sqrt(sum);
}

mpf findminm(const vector<Vector>& Mat) {
	if (Mat.empty() || Mat[0].empty()) {
		return 0.0; // 处理空矩阵或空向量的情况
	}

	mpf minMagnitude = vectorMagnitude(Mat[0]); // 初始化为第一个向量的模长

	for (const auto& vec : Mat) {
		mpf mag = vectorMagnitude(vec);
		if (mag < minMagnitude) {
			minMagnitude = mag;
		}
	}

	return minMagnitude;
}

Vector mattovec(vector<Vector> Mat) {
	Vector result;
	int m = Mat.size(), n = Mat[0].size();
	for (int i = 0; i < m; i++)
		for (int j = 0; j < n; j++)
			result.push_back(Mat[i][j]);
	return result;
}
void solve5() {
	mpz e, dl, dm, N;
	e.set_str("124c552642ef2467aaecde51b0f3e1bee2ebe87bae39a956ad56cf7eec669cdc7b9664ea435b4c3492b8e610e0a182e1a76c7af443ca2962672b4e703c4f359cf8d88a67db77be2491b74bcdae58691b69e6ea06d067815b26fc0d669d8c06f11a728154dc8cdf983a056633fecadc417df4304625c3e6f91ec3d655a91a29e9", 16);
	dl.set_str("2b26d177dc20ceea15de6e3c5a03207fb326a42d53a9", 16);
	dm.set_str("acfad4bbb97a99b6bbc82c8b44a5260bcfe9c4a0acf437186ff4d5d1594cc5c1", 16);
	N.set_str("94eab94581f4931a5ea6aabcfe0598600fa3e0a06573887aed69e274f14484472dc3feaf50d4ef384e502f747f5605c1d2a4c8172b6ef134b7e96d6c383a9cb967ccbbd8b3647848d34928982a274999c2df00bd7dd11bf25acd61411e3395637e85dd84ecf785ff1027eed91f3976c8186e2e940edcb5fed8d759a5028b47a1", 16);
	mpz B, X, Y, K, E;
	B = e * (dl + dm << 256) - 1;
	E = e << 176;
	mpz_root(Y.get_mpz_t(), N.get_mpz_t(), 2);
	mpz_root(X.get_mpz_t(), N.get_mpz_t(), 4);
	K = B / (N - 2 * Y + 1);
	poly f(2, 1);
	f.setCoef(0, 0, K * (N + 1) - B);
	f.setCoef(1, 0, X * (N + 1));
	f.setCoef(0, 1, K * Y);
	f.setCoef(1, 1, 1);

	vector<poly> p(8,poly(2,1));
	p[0].setCoef(0, 0, E * E);
	p[1].setCoef(1, 0, X * E * E);
	p[2].setCoef(0, 1, Y * E * E);
	p[3] = f.mulc(E);
	p[4].setCoef(2, 0, E * E);
	p[5] = f.xshift(1, 1, X, E);
	p[6] = f * f;
	p[7] = f.yshift(1, 1, Y, E);

	vector<Vector> Mat(8);
	for(int i=0;i<8;i++)
		Mat[i] = mattovec(p[i].coef);


	//printMat(Mat);
	//cout << endl;
	//cout << findminm(Mat) << endl;
	//cout << endl;

	LLL(Mat);

	printMat(Mat);
	cout << endl;
	cout << findminm(Mat) << endl;
	mpf E2 = mpf_class(E.get_str()) * mpf_class(E.get_str());
	cout << E2;
}
