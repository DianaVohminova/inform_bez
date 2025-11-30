#ifndef BIGNUMBER_H
#define BIGNUMBER_H

#define M_PI 3.14159265358979323846
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>
#include <wchar.h>
#include <locale.h>

// Используется для реализации быстрого преобразования Фурье (FFT)
#include <complex.h>

// Константы для FFT умножения
#define BASE_DIGS 4
#define FFT_BASE 10000

// Большое число
typedef struct {
    int* digits; // динамический массив для хранения цифр
    int length;             // длина числа
    int sign;              // отрицательное число
                            // 1-полож. -1 - отриц.
} BigNum;

// ==================== ОСНОВНЫЕ ФУНКЦИИ ====================

// Создание BigNum заданной длины
BigNum CreateBigNum(int length);

// Освобождение памяти, занятой BigNum
void FreeBigNum(BigNum *num);

// Копирование BigNum
BigNum CopyBigNum(BigNum a);

// Удаление ведущих нулей
void TrimBigNum(BigNum *num);

// Генерация случайного BigNum заданной длины
void GenerateRandomBigNum(BigNum *num, int length);

// Преобразование строки в BigNum
BigNum MakeBigNum(const char *str);

// Печать BigNum в файл
void PrintBigNum(const char *filename, BigNum a);

// Чтение BigNum из файла
BigNum ReadBigNum(const char *filename);

// ==================== АРИФМЕТИЧЕСКИЕ ОПЕРАЦИИ ====================

// Сложение по модулю (без учета знака)
BigNum AddAbs(BigNum a, BigNum b);

// Вычитание по модулю (a >= b)
BigNum SubAbs(BigNum a, BigNum b);

// Сравнение по модулю (без учета знака)
int CompareAbs(BigNum a, BigNum b);

// Сложение с учетом знака
BigNum Add(BigNum a, BigNum b);

// Вычитание с учетом знака
BigNum Sub(BigNum a, BigNum b);

// Простое умножение (для маленьких чисел)
BigNum SimpleMul(BigNum a, BigNum b);

// Основная функция умножения
BigNum Mul(BigNum a, BigNum b);

// Сравнение с учетом знака
int Compare(BigNum a, BigNum b);

// Деление в столбик
BigNum Div(BigNum a, BigNum b);

// Сдвиг влево (умножение на 10^k)
void ShiftLeftT(BigNum *a, int k);

// Умножение BigNum на цифру (0..9)
BigNum MulByDigit(BigNum a, int d);

// ==================== НОД И НОК ====================

// Наибольший общий делитель (НОД)
BigNum NOD(BigNum a, BigNum b);

// Наименьшее общее кратное (НОК)
BigNum NOK(BigNum a, BigNum b);

// ==================== ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ====================

// Проверка на ноль
int IsZero(BigNum a);

// Проверка на чётность
int IsEven(BigNum a);

// Оптимизированный сдвиг вправо (деление на 2)
void ShiftRight(BigNum *a);

// Оптимизированный сдвиг влево (умножение на 2)
void ShiftLeft(BigNum *a);

// Сдвиг вправо на k позиций (деление на 10^k)
void ShiftRightT(BigNum *a, int k);

// ==================== FFT УМНОЖЕНИЕ ====================

// Быстрое преобразование Фурье
static void fft(long double complex *a, int n, int invert);

// Свертка с помощью FFT
static long long * convolution(const int *a, int na, const int *b, int nb, int *out_len);

// Упаковка BigNum в блоки для FFT
static int * pack_blocks(const BigNum *a, int *blocks_len);

// Распаковка блоков обратно в BigNum
static BigNum unpack_blocks_to_BigNum(long long *blocks, int blocks_len);

// Умножение больших чисел с использованием FFT
BigNum FFTMul(BigNum a, BigNum b);

// ==================== ДОПОЛНИТЕЛЬНЫЕ ФУНКЦИИ ====================

// Создание BigNum из int
BigNum CreateFromInt(int n);

// Возведение в степень
BigNum Pow(BigNum a, int n);

// Деление BigNum на цифру (0..9)
BigNum DivByDigit(BigNum a, int d);

// Оптимизированное возведение в степень по модулю
BigNum PowMod(BigNum a, BigNum b, BigNum m);
// ==================== ТЕСТЫ ====================


// Возведение в степень по модулю (a^b mod m)
BigNum PowMod(BigNum a, BigNum b, BigNum m);
// Вычисление остатка от деления (a mod m)
BigNum Mod(BigNum a, BigNum m);

// Проверка, является ли число точной степенью
int IsPerfectPower(BigNum n);
// Поиск наименьшего r такого, что ord_r(n) > log^2(n)
BigNum FindSmallestR(BigNum n);
// Проверка малых простых делителей
int HasSmallPrimeFactors(BigNum n);
// Основная функция AKS теста
int AKSIsPrime(BigNum n);
// Упрощенная версия AKS для тестирования (использует вероятностные проверки для больших чисел)
// int AKSIsPrimeSimple(BigNum n);
// Функция для измерения времени выполнения
double measure_time(int (*test_func)(BigNum), BigNum n);
// Приближенное вычисление кубического корня
BigNum CbrtApprox(BigNum n);
// Приближенное вычисление квадратного корня
BigNum SqrtApprox(BigNum n);



// Упрощенная версия теста Миллера для очень больших чисел
int MillerTestSimple(BigNum n);
// Детерминированный тест Миллера для проверки простоты
int MillerTest(BigNum n);

// Функция для вывода BigNum в консоль (для отладки)
void PrintBigNumToConsole(BigNum num);

// Функция для получения строкового представления BigNum (первые 10 цифр)
void get_number_preview(BigNum n, char* buffer, int buffer_size);



// Функция для измерения времени с возвратом результата
typedef struct {
    int is_prime;
    double time;
} TestResult;

TestResult measure_time_with_result(int (*test_func)(BigNum), BigNum n);

// Функции для красивого вывода таблиц
void print_aks_header();
void print_miller_header();

void print_table_row(const char* test_name, const char* number, const char* size_info, 
                    const char* result, double time, int is_last);

void print_footer();





// Функции для вывода отдельных таблиц для каждого числа
void print_test_header(const char* algorithm, const char* number_type);

void print_number_details(const char* label, const char* value);
void print_test_result(const char* result, double time);

void print_test_footer();

// Функция для получения полного представления BigNum как строки
void get_full_number_string(BigNum n, char* buffer, int buffer_size);

// Функция для получения информации о размере числа
const char* get_size_info(BigNum n);


// Функция для подсчета видимой ширины строки с UTF-8 символами
int utf8_strwidth(const char *str);

// Функция для выравнивания UTF-8 строки
void print_utf8_aligned(const char *str, int total_width);



// Функция для получения строкового представления BigNum (первые 10 цифр)
void get_number_preview(BigNum n, char* buffer, int buffer_size);
// Функции для вывода таблиц 
void print_aks_result_table(const char* test_name, const char* number, const char* size_info, 
                           const char* result, double time);

void print_miller_result_table(const char* test_name, const char* number, const char* size_info, 
                              const char* result, double time);



#endif // BIGNUMBER_H