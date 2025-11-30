
// gcc -o big_primes main.c BigNumber.c Tests.c -lm
// ./big_primes
#include "BigNumber.h"
#include <stdio.h>
#include <stdlib.h>
#include <locale.h>
#include <time.h>
#include <string.h>

int main() {
    
    // Простые числа < 2^32
    BigNum small_prime = MakeBigNum("999983"); 
    // BigNum small_prime = MakeBigNum("2147483647"); // 2^31 - 1 


    // Простые числа > 2^64
    
    BigNum large_prime = MakeBigNum("162259276829213363391578010288127"); // 2^107 - 1
    // BigNum large_prime = MakeBigNum("170141183460469231731687303715884105727"); // 2^127 - 1


    // // Составные числа
    // BigNum small_prime = MakeBigNum("536870911"); // 2^29 - 1 
    // BigNum large_prime = MakeBigNum("147573952589676412927");// 2^67 - 1


    // ==================== ТЕСТ AKS ====================
    
    // Тест 1: Малое простое число
    TestResult result1 = measure_time_with_result(AKSIsPrime, small_prime);
    
    char num1[20];
    get_number_preview(small_prime, num1, sizeof(num1));
    
    print_aks_result_table("AKS Тест 1", num1, "< 2^32", 
                          result1.is_prime ? "Простое" : "Составное", 
                          result1.time);

    // Тест 2: Большое простое число  
    TestResult result2 = measure_time_with_result(AKSIsPrime, large_prime);
    
    char num2[40];
    get_number_preview(large_prime, num2, sizeof(num2));
    
    print_aks_result_table("AKS Тест 2", num2, "> 2^64", 
                          result2.is_prime ? "Простое" : "Составное", 
                          result2.time);

    // ==================== ТЕСТ МИЛЛЕРА ====================
    
    // Тест 1: Малое простое число
    TestResult result1_m = measure_time_with_result(MillerTest, small_prime);
    
    char num1_m[20];
    get_number_preview(small_prime, num1_m, sizeof(num1_m));
    
    print_miller_result_table("Миллер Тест 1", num1_m, "< 2^32", 
                             result1_m.is_prime ? "Простое" : "Составное", 
                             result1_m.time);

    // Тест 2: Большое простое число
    TestResult result2_m = measure_time_with_result(MillerTest, large_prime);
    
    char num2_m[40];
    get_number_preview(large_prime, num2_m, sizeof(num2_m));
    
    print_miller_result_table("Миллер Тест 2", num2_m, "> 2^64", 
                             result2_m.is_prime ? "Простое" : "Составное", 
                             result2_m.time);

    // ==================== СРАВНЕНИЕ ====================
   
    double aks_ratio = result2.time / result1.time;
    double miller_ratio = result2_m.time / result1_m.time;
    
    printf("AKS малое число:     %-17.4f\n", result1.time);
    printf("AKS большое число:   %-17.4f\n", result2.time);
    printf("Отношение AKS:       %-17.2f\n", aks_ratio);
    printf("Миллер малое число:  %-17.4f\n", result1_m.time);
    printf("Миллер большое число:%-17.4f\n", result2_m.time);
    printf("Отношение Миллер:    %-17.2f\n", miller_ratio);

    // Освобождение памяти
    FreeBigNum(&small_prime);
    FreeBigNum(&large_prime);
    
    return 0;
}