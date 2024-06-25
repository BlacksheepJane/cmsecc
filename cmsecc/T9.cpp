#include "T9.hpp"
tuple<mpz, mpz, mpz> extendedEuclid(mpz a, mpz b) {//拓展欧几里得
	if (b == 0)
		return make_tuple(a, 1, 0);

	mpz gcd, x, y;
	tie(gcd, x, y) = extendedEuclid(b, a % b);
	return make_tuple(gcd, y, x - (a / b) * y);
}


mpz modInverse(mpz b, mpz a) {// 求 b 在模 a 情况下的逆
	mpz gcd, x, y;
	tie(gcd, x, y) = extendedEuclid(b, a);

	if (gcd != 1) {
		throw invalid_argument("Inverse does not exist.");
	}
	else {
		// x 可能为负数，将其转换为正数
		x = (x % a + a) % a;
		return x;
	}
}


void decompose(mpz k, mpz& a, mpz& b, mpz& p) {//分解k为a*b，其中b为奇数，p=2^a
	a = 0, p = 1;
	while (k % 2 == 0) {
		k = k / 2;
		a = a + 1;
		p = p << 1;
	}
	b = k;
}

mpz solvet(mpz s, mpz m) {//解最后的方程
	if (s % 8 != 1) return -1;
	mpz k = 3;
	mpz t, d, a, p = 1;//p=pow(2,k-3)
	t = 1;
	d = (1 - s) >> 3;
	while (k != m + 1) {
		a = abs(d % 2);
		d = (d + a * t) / 2 + a * p;
		t = t + 4 * p * a;
		k++;
		p = p << 1;
	}
	return t;
}

pair<mpz, mpz> getpq() {
	mpz e = 65537;
	mpz dl, N;
	dl.set_str("20142ae2802b877eb4dfa8a462e7d017c4d348181c367fd1a661ec9b6bbcca9dcb6601ccb6c10416b7f3c20129527346bbc136ee60f9945125cba03a9bba3720f7411", 16);
	N.set_str("cc5b706f373a79c680cec9527aac573fd435129cf16c23334085bf97832e5a6c78b633c2f244b12a62f87ec5295dd89fcf3c808c39e45a9afdbda2f8d2d0b50d61b685c0fe9eb41a7018a40f98892f96d738e2a4e740d4e507bcbd07f68c1ecb2ca10bd780ce65265a7e4da00f1031a5db9d038878a29a5ffefcaf2119720005", 16);
	for (mpz k = 1; k < e; k++) {
		mpz a, b, po;
		decompose(k, a, b, po);//分解k为po*b，其中b为奇数，po=2^a
		if ((e * dl - 1) % po != 0) continue;
		//cout << "a:" << a << "\n" << "b:" << b << "\n" << "k:" << k << "\n" << "pow:" << po << "\n";
		mpz p530 = pow(2, 530);
		mpz bi = modInverse(b, p530);//ok
		mpz r = (N + 1 - (e * dl - 1) / po * bi) % (p530 / po);
		r = (r + p530 / po) % (p530 / po);
		//cout << r % 4 << endl;
		if (r % 4 != 2) continue;
		//cout << k << endl;
		mpz p512 = pow(2, 512);
		mpz m = 512, s = (r * r / 4 - N) / 4;
		mpz p1 = (2 * solvet(s, m) + r / 2 + p512) % p512;
		mpz p2 = (-2 * solvet(s, m) + r / 2 + p512) % p512;
		if (N % p1 == 0) {
			cout << p1 << endl << N / p1;
			return { p1,N / p1 };
		}
		else if (N % p2 == 0) {
			cout << p2 << endl << N / p2;
			return { p2,N / p2 };
		}
	}
	cout << -1;
	return { -1,-1 };
}

//ans
//	mpz p("10846327614507406655792564994667714933899253952298425269758486277699020260863878841336945944423557227322075142932155647161674834513419649086027797728283207");
//  mpz q("13230698921743059551824416344443899403552566005548844066885308262194945624798691552085027494542401017580805838924364362875395655984648652258879913619766611");
