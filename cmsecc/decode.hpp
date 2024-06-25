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
        return 0; // ������벻�ǺϷ���ʮ�������ַ������� 0 ���������ʵ��Ĵ�����
}

// ��ʮ�������ַ���ת��Ϊ ASCII �ַ���
string hexToAscii(const string& hexStr) {
    string asciiStr;
    stringstream ss;

    // ÿ����ʮ�������ַ���Ӧһ�� ASCII �ַ�
    for (size_t i = 0; i < hexStr.length(); i += 2) {
        // ������ʮ�������ַ�ת��Ϊ����
        int highNibble = hexCharToInt(hexStr[i]);
        int lowNibble = hexCharToInt(hexStr[i + 1]);

        // �����Ӧ�� ASCII �ַ�����ӵ�����ַ�����
        char asciiChar = static_cast<char>((highNibble << 4) | lowNibble);
        asciiStr += asciiChar;
    }
    reverse(asciiStr.begin(), asciiStr.end());
    return asciiStr;
}
string decode(mpz e ,mpz c,mpz p,mpz q) {
    //// ����ŷ������ phi(N) = (p-1)(q-1)
    mpz_class N = p * q;
    mpz_class phi_N = N + 1 - p - q;
    
    // ����˽Կ d������ d * e �� 1 (mod phi(N))
    mpz_class d;
    mpz_invert(d.get_mpz_t(), e.get_mpz_t(), phi_N.get_mpz_t());

    // ʹ��˽Կ d �������� c
    mpz_class decrypted_message;
    mpz_powm(decrypted_message.get_mpz_t(), c.get_mpz_t(), d.get_mpz_t(), N.get_mpz_t());

    string str = "2e636974656d6874697261206e69206c756665737520646e6120746e6174726f706d692074736f6d2065687420666f20656e6f206562206f74206e776f6e6b2073692073726f7463616620656d697270207269656874206f746e6920737265626d756e20657469736f706d6f632065687420676e69766c6f736552";
    // ת��Ϊ ASCII �ַ���
    string asciiString = hexToAscii(str);
    // ������
    return asciiString;
}