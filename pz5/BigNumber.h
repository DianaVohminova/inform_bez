#ifndef BIGNUMBER_H
#define BIGNUMBER_H

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <random>
#include <chrono>
#include <cmath>
#include <map>

class BigNumber {
private:
    std::vector<int> digits;
    bool isNegative;

    void removeLeadingZeros();
    int compareAbsolute(const BigNumber& other) const;
    BigNumber addAbsolute(const BigNumber& other) const;
    BigNumber subtractAbsolute(const BigNumber& other) const;
    BigNumber multiplyByDigit(int digit) const;
    BigNumber divideByDigit(int digit) const;

public:
    // ==================== КОНСТРУКТОРЫ ====================
    BigNumber();
    BigNumber(const std::string& str);
    BigNumber(long long num);
    BigNumber(const BigNumber& other);
    BigNumber(int numDigits, std::mt19937& gen);

    // ==================== ОПЕРАТОРЫ ПРИСВАИВАНИЯ ====================
    BigNumber& operator=(const BigNumber& other);
    BigNumber& operator=(long long num);

    // ==================== АРИФМЕТИЧЕСКИЕ ОПЕРАТОРЫ ====================
    BigNumber operator+(const BigNumber& other) const;
    BigNumber operator-(const BigNumber& other) const;
    BigNumber operator*(const BigNumber& other) const;
    BigNumber operator/(const BigNumber& other) const;
    BigNumber operator%(const BigNumber& other) const;
    BigNumber operator^(const BigNumber& exponent) const;

    // ==================== ОПЕРАТОРЫ СРАВНЕНИЯ ====================
    bool operator==(const BigNumber& other) const;
    bool operator!=(const BigNumber& other) const;
    bool operator<(const BigNumber& other) const;
    bool operator<=(const BigNumber& other) const;
    bool operator>(const BigNumber& other) const;
    bool operator>=(const BigNumber& other) const;

    // ==================== УНАРНЫЕ ОПЕРАТОРЫ ====================
    BigNumber operator-() const;
    BigNumber operator+() const;

    // ==================== МЕТОДЫ ВВОДА/ВЫВОДА ====================
    friend std::ostream& operator<<(std::ostream& os, const BigNumber& num);
    friend std::istream& operator>>(std::istream& is, BigNumber& num);

    // ==================== ВСПОМОГАТЕЛЬНЫЕ МЕТОДЫ ====================
    std::string toString() const;
    bool isZero() const;
    BigNumber abs() const;
    size_t getDigitCount() const { return digits.size(); }
    
    int getDigitAt(size_t index) const {
        if (index < digits.size()) {
            return digits[index];
        }
        return 0;
    }
    
    int getLastDigit() const {
        return digits.empty() ? 0 : digits[0];
    }
    
    const std::vector<int>& getDigits() const {
        return digits;
    }
    
    bool isEven() const {
        return digits.empty() ? true : (digits[0] % 2 == 0);
    }
    
    bool isOdd() const {
        return !isEven();
    }
    
    int sign() const {
        if (isZero()) return 0;
        return isNegative ? -1 : 1;
    }

    // ==================== СТАТИЧЕСКИЕ МАТЕМАТИЧЕСКИЕ ФУНКЦИИ ====================
    static BigNumber gcd(BigNumber a, BigNumber b);
    static BigNumber lcm(const BigNumber& a, const BigNumber& b);
    static BigNumber generateRandomPrime(int numDigits, std::mt19937& gen);

    // ==================== МЕТОДЫ ПРОВЕРКИ ПРОСТОТЫ ====================
    static bool isPrimeStandard(const BigNumber& n);
    static bool isPrimeEratosthenes(const BigNumber& n, int limit = 1000000);
    static bool isPrimeAtkin(const BigNumber& n, int limit = 1000000);
    static bool lucasLehmerTest(int p);
    static BigNumber sqrt(const BigNumber& n);
    static BigNumber modPow(const BigNumber& base, const BigNumber& exponent, const BigNumber& mod);
    static BigNumber log(const BigNumber& n);
    bool isPrime(int iterations = 10) const;

    // ==================== ECPP МЕТОДЫ ====================
    bool isPrimeECPP(int maxAttempts = 20) const;
    
    // Вспомогательные методы для ECPP
    static BigNumber randomBigNumber(const BigNumber& max, std::mt19937& gen);
    
    // Методы факторизации
    static std::vector<BigNumber> factorize(const BigNumber& n, int maxAttempts = 5);
    static BigNumber pollardRho(const BigNumber& n, int maxIterations = 1000);
};

#endif