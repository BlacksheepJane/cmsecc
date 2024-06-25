#pragma once
#include"options.h"
#include <string>
#include <sstream>

int hexCharToInt(char ch) {
    if (ch >= '0' && ch <= '9')
        return ch - '0';
    else if (ch >= 'A' && ch <= 'F')
        return ch - 'A' + 10;
    else if (ch >= 'a' && ch <= 'f')
        return ch - 'a' + 10;
    else
        return 0; // 如果输入不是合法的十六进制字符，返回 0 或者其他适当的错误处理
}

// 将十六进制字符串转换为 ASCII 字符串
string hexToAscii(const string& hexStr) {
    string asciiStr;
    stringstream ss;

    // 每两个十六进制字符对应一个 ASCII 字符
    for (size_t i = 0; i < hexStr.length(); i += 2) {
        // 将两个十六进制字符转换为整数
        int highNibble = hexCharToInt(hexStr[i]);
        int lowNibble = hexCharToInt(hexStr[i + 1]);

        // 计算对应的 ASCII 字符并添加到结果字符串中
        char asciiChar = static_cast<char>((highNibble << 4) | lowNibble);
        asciiStr += asciiChar;
    }
    reverse(asciiStr.begin(), asciiStr.end());
    return asciiStr;
}
string decode(mpz e ,mpz c,mpz p,mpz q) {
    //// 计算欧拉函数 phi(N) = (p-1)(q-1)
    mpz_class N = p * q;
    mpz_class phi_N = N + 1 - p - q;
    
    // 计算私钥 d，满足 d * e ≡ 1 (mod phi(N))
    mpz_class d;
    mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi_N.get_mpz_t());

    // 使用私钥 d 解密密文 c
    mpz_class decrypted_message;
    mpz_powm(decrypted_message.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());

    string str = "2e636974656d6874697261206e69206c756665737520646e6120746e6174726f706d692074736f6d2065687420666f20656e6f206562206f74206e776f6e6b2073692073726f7463616620656d697270207269656874206f746e6920737265626d756e20657469736f706d6f632065687420676e69766c6f736552";
    // 转换为 ASCII 字符串
    string asciiString = hexToAscii(str);
    // 输出结果
    return asciiString;
}