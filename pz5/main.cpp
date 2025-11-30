/*
g++ -std=c++11 -O2 -o ecpp_test BigNumber.cpp main.cpp
./ecpp_test
*/

#include "BigNumber.h"
#include <iostream>
#include <chrono>
#include <vector>
#include <iomanip>

using namespace std;
using namespace std::chrono;

template<typename Func>
long long measureTime(Func func, bool& result) {
    auto start = high_resolution_clock::now();
    result = func();
    auto end = high_resolution_clock::now();
    return duration_cast<milliseconds>(end - start).count();
}

void testECPPPrimes() {
    cout << "=== ECPP ПРОВЕРКА ПРОСТЫХ ЧИСЕЛ ===" << endl << endl;
    
    vector<pair<string, string>> testCases = {
        {"1009", "Малое простое"},
        {"10007", "Простое среднего размера"}, 
        {"104729", "10000-е простое"},
        {"1000003", "Простое > 10^6"},
        {"10000019", "Большое простое"}
    };
    
    for (const auto& testCase : testCases) {
        string numStr = testCase.first;
        string description = testCase.second;
        BigNumber num(numStr);
        
        cout << "Тест: " << description << endl;
        cout << "Число: " << num << " (" << numStr.length() << " цифр)" << endl;
        
        bool result;
        long long ecppTime = measureTime([&]() {
            return num.isPrimeECPP(15);
        }, result);
        
        cout << "ECPP результат: " << (result ? "простое" : "составное") 
             << " (время: " << ecppTime << " мс)" << endl;
        cout << string(50, '-') << endl;
    }
}

void testECPPComposites() {
    cout << "=== ECPP ПРОВЕРКА СОСТАВНЫХ ЧИСЕЛ ===" << endl << endl;
    
    vector<pair<string, string>> composites = {
        {"1001", "7 × 11 × 13"},
        {"10001", "73 × 137"}, 
        {"100001", "11 × 9091"},
        {"999999", "3³ × 7 × 11 × 13 × 37"}
    };
    
    for (const auto& test : composites) {
        BigNumber num(test.first);
        cout << "Число: " << num << endl;
        
        bool result;
        long long time = measureTime([&]() {
            return num.isPrimeECPP(10);
        }, result);
        
        cout << "ECPP результат: " << (result ? "простое" : "составное") 
             << " (время: " << time << " мс)" << endl;
        
        if (!result) {
            auto factors = BigNumber::factorize(num, 3);
            cout << "Факторизация: ";
            for (size_t i = 0; i < factors.size(); ++i) {
                cout << factors[i];
                if (i < factors.size() - 1) cout << " × ";
            }
            cout << " (ожидается: " << test.second << ")" << endl;
        }
        cout << endl;
    }
}

void testMersenneWithECPP() {
    cout << "=== ЧИСЛА МЕРСЕННА С ECPP ===" << endl << endl;
    
    vector<int> exponents = {2, 3, 5, 7, 13, 17};
    
    for (int exp : exponents) {
        cout << "2^" << exp << "-1:" << endl;
        
        BigNumber mersenne = (BigNumber(2) ^ BigNumber(exp)) - BigNumber(1);
        cout << "Число: " << mersenne.toString().length() << " цифр" << endl;
        
        bool result;
        long long ecppTime = measureTime([&]() {
            return mersenne.isPrimeECPP(10);
        }, result);
        
        cout << "ECPP: " << (result ? "простое" : "составное")
             << " (" << ecppTime << " мс)" << endl << endl;
    }
}

void testECPPPerformance() {
    cout << "=== ПРОИЗВОДИТЕЛЬНОСТЬ ECPP ===" << endl << endl;
    
    vector<string> testNumbers = {
        "1009", "7919", "10007", "104729", "1000003", "10000019", "1000000007", "2147483647"
    };
    
    cout << setw(12) << "Число" << setw(8) << "\t\tЦифр" << setw(12) << "ECPP(мс)" << endl;
    cout << string(35, '-') << endl;
    
    for (const auto& numStr : testNumbers) {
        BigNumber num(numStr);
        bool result;
        
        long long ecppTime = measureTime([&]() { 
            return num.isPrimeECPP(10); 
        }, result);
        
        cout << setw(12) << numStr << setw(8) << numStr.length()
             << setw(12) << ecppTime << endl;
    }
}

int main() {
    cout << "ПРАКТИЧЕСКАЯ РАБОТА №5 - АЛГОРИТМ ECPP" << endl;
    cout << "Аткин-Морейн на эллиптических кривых" << endl;
    cout << "=============================================" << endl << endl;
    
    testECPPPrimes();
    testECPPComposites(); 
    testMersenneWithECPP();
    testECPPPerformance();
    
    return 0;
}