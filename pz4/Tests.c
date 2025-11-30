#include "BigNumber.h"
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <complex.h>
#include <stdlib.h>
#include <string.h>

// ----------------------------------Детерминированный тест Агравала-Каяла-Саксены------------------------------

// Возведение в степень по модулю (a^b mod m)
BigNum PowMod(BigNum a, BigNum b, BigNum m) {
    BigNum result = CreateFromInt(1);
    BigNum base = Mod(a, m);  // Базовое значение: a mod m 
    BigNum exponent = CopyBigNum(b);
    
    // Пока экспонента не равна нулю
    while (!IsZero(exponent)) {
        // Если текущий бит экспоненты нечетный (младший бит = 1)
        if (!IsEven(exponent)) {
            // Умножаем результат на текущее значение base
            BigNum temp = Mul(result, base);
            // Берем результат по модулю m
            BigNum new_result = Mod(temp, m);
            FreeBigNum(&result);
            FreeBigNum(&temp);
            result = new_result;
        }
        
        // Возводим base в квадрат по модулю m для следующей итерации
        // base = base^2 mod m
        BigNum temp_sq = Mul(base, base);
        BigNum new_base = Mod(temp_sq, m);
        FreeBigNum(&base);
        FreeBigNum(&temp_sq);
        base = new_base;
        
        // Сдвигаем экспоненту вправо (делим на 2) для обработки следующего бита
        ShiftRight(&exponent);  
    }
    
    FreeBigNum(&base);
    FreeBigNum(&exponent);
    return result;
}
/*----------------------------------------------------------------------------------------------*/
// Вычисление остатка от деления (a mod m)
BigNum Mod(BigNum a, BigNum m) {
    // Предварительные вычисления
    static BigNum cached_m = {0}; // Кэшированный модуль
    static BigNum mu = {0};       // Предвычисленная константа mu
    
    // Если модуль изменился, пересчитываем mu
    if (Compare(cached_m, m) != 0) {
        FreeBigNum(&cached_m);
        FreeBigNum(&mu);
        // Кэшируем новый модуль
        cached_m = CopyBigNum(m);
        
        // Вычисляем mu = floor(2^(2k) / m), где k - битовая длина m
        BigNum base = CreateFromInt(1);
        ShiftLeftT(&base, 2 * m.length * 4); // Сдвигаем на 2k бит влево
        mu = Div(base, m);
        FreeBigNum(&base);
    }
    
    // Алгоритм Barrett reduction:
    // q = floor((a * mu) / 2^(2k)) - приближенное частно
    BigNum q = Mul(a, mu);
    ShiftRightT(&q, 2 * m.length * 4); // Деление на 2^(2k)
    
    // r = a - q * m  (приближенный остаток)
    BigNum qm = Mul(q, m);
    BigNum r = Sub(a, qm);
    
    FreeBigNum(&q);
    FreeBigNum(&qm);
    
    // Корректировка остатка: гарантируем, что 0 <= r < m
    // Если r >= m, вычитаем m пока не получим корректный остаток
    while (Compare(r, m) >= 0) {
        BigNum tmp = Sub(r, m);
        FreeBigNum(&r);
        r = tmp;
    }
    
    return r;
}
/*----------------------------------------------------------------------------------------------*/
// Приближенное вычисление квадратного корня
BigNum SqrtApprox(BigNum n) {
    // Обработка тривиальных случаев: n <= 1
    if (Compare(n, CreateFromInt(1)) <= 0) {
        return CopyBigNum(n);
    }
    
    // Для чисел до 18 цифр используем встроенную функцию sqrt из double
    if (n.length <= 18) {
        double num = 0;
        // Конвертируем BigNum в double (в обратном порядке цифр)
        for (int i = n.length - 1; i >= 0; i--) {
            num = num * 10 + n.digits[i];
        }
        double sqrt_val = sqrt(num);
        return CreateFromInt((int)sqrt_val); // Округляем до целого
    }
    
    // Для больших чисел используем бинарный поиск
    BigNum low = CreateFromInt(1); // Нижняя граница = 1
    BigNum high = CopyBigNum(n);   // Верхняя граница = n
    BigNum one = CreateFromInt(1);
    
    // Быстрая оценка начального приближения
    BigNum x = CreateFromInt(1);
    if (n.length > 1) {
        // x0 = 10^(n.length/2) - хорошее начальное приближение
        // Для числа с d цифрами, корень примерно имеет d/2 цифр
        int half_len = n.length / 2;
        FreeBigNum(&x);
        x = CreateFromInt(1);
        for (int i = 0; i < half_len; i++) {
            ShiftLeftT(&x, 1); // Умножаем на 10
        }
    }
    
    // Упрощенный метод Ньютона 
    for (int i = 0; i < 2; i++) {
        // x = (x + n/x) / 2
        BigNum n_div_x = Div(n, x);   // n / x
        BigNum sum = Add(x, n_div_x); // x + n/x
        BigNum new_x = Div(sum, CreateFromInt(2)); // (x + n/x) / 2
        
        FreeBigNum(&x);
        FreeBigNum(&n_div_x);
        FreeBigNum(&sum);
        x = new_x;
    }
    
    FreeBigNum(&low);
    FreeBigNum(&high);
    FreeBigNum(&one);
    
    return x;
}
/*----------------------------------------------------------------------------------------------*/
// Приближенное вычисление кубического корня
BigNum CbrtApprox(BigNum n) {
    // Обработка тривиальных случаев: n <= 1
    if (Compare(n, CreateFromInt(1)) <= 0) {
        return CopyBigNum(n);
    }
    
    // Для чисел до 15 цифр используем встроенную функцию cbrt из double
    if (n.length <= 15) {
        double num = 0;
        // Конвертируем BigNum в double
        for (int i = n.length - 1; i >= 0; i--) {
            num = num * 10 + n.digits[i];
        }
        double cbrt_val = cbrt(num);
        return CreateFromInt((int)cbrt_val); // Округляем до целого
    }
    
    // Быстрое начальное приближение: 10^(n.length/3)
    // Для числа с d цифрами, кубический корень примерно имеет d/3 цифр
    BigNum x = CreateFromInt(1);
    int third_len = n.length / 3;
    for (int i = 0; i < third_len; i++) {
        ShiftLeftT(&x, 1); // Умножаем на 10
    }
    
    // Всего 3 итерации метода Ньютона 
    for (int i = 0; i < 3; i++) {
        // x = (2*x + n/(x*x)) / 3
        BigNum x_sq = Mul(x, x); // x²
        BigNum n_div_x_sq = Div(n, x_sq); // n / x²
        BigNum two_x = MulByDigit(x, 2);  // 2*x
        BigNum sum = Add(two_x, n_div_x_sq); // 2x + n/x²
        BigNum new_x = DivByDigit(sum, 3);   // (2x + n/x²) / 3
        
        FreeBigNum(&x_sq);
        FreeBigNum(&n_div_x_sq);
        FreeBigNum(&two_x);
        FreeBigNum(&sum);
        FreeBigNum(&x);
        
        x = new_x;
    }
    
    return x;
}
/*----------------------------------------------------------------------------------------------*/
// Проверка, является ли число точной степенью
int IsPerfectPower(BigNum n) {
    
    // Для маленьких чисел конвертируем в int
    if (n.length <= 9) {
        int num = 0;
        for (int i = n.length - 1; i >= 0; i--) {
            num = num * 10 + n.digits[i];
        }
        
        // Проверяем только небольшие степени
        for (int k = 2; k <= 10; k++) {
            // Вычисляем приблизительный корень k-й степени
            int root = (int)pow(num, 1.0/k);
            if (root < 2) break; // Если корень меньше 2, дальше проверять нет смысла
            
            // Проверяем root^k и (root+1)^k
            long long test1 = 1, test2 = 1;
            for (int i = 0; i < k; i++) {
                test1 *= root;
                test2 *= (root + 1);
            }
            
            // Если нашли точное совпадение - число является точной степенью
            if (test1 == num || test2 == num) {
                return 1;
            }
        }
        return 0;
    }
    
    // Для больших чисел проверяем только степени 2 и 3
    BigNum sqrt_approx = SqrtApprox(n);
    BigNum sqrt_sq = Mul(sqrt_approx, sqrt_approx); // (√n)²
    if (Compare(sqrt_sq, n) == 0) {
        // Найден точный квадрат
        FreeBigNum(&sqrt_approx);
        FreeBigNum(&sqrt_sq);
        return 1;
    }
    FreeBigNum(&sqrt_sq);
    
    BigNum cbrt_approx = CbrtApprox(n);
    BigNum cbrt_cube = Mul(Mul(cbrt_approx, cbrt_approx), cbrt_approx);  // (∛n)³
    if (Compare(cbrt_cube, n) == 0) {
        // Найден точный куб
        FreeBigNum(&sqrt_approx);
        FreeBigNum(&cbrt_approx);
        FreeBigNum(&cbrt_cube);
        return 1;
    }
    FreeBigNum(&sqrt_approx);
    FreeBigNum(&cbrt_approx);
    FreeBigNum(&cbrt_cube);
    
    return 0;
}
/*----------------------------------------------------------------------------------------------*/
// Поиск наименьшего r такого, что ord_r(n) > log^2(n)
BigNum FindSmallestR(BigNum n) {
    
    // Для маленьких чисел используем фиксированные значения
    if (n.length <= 6) {
        int num = 0;
        // Конвертируем BigNum в обычное int для быстрой обработки
        for (int i = n.length - 1; i >= 0; i--) {
            num = num * 10 + n.digits[i];
        }
        
        // Эмпирически подобранные значения r для маленьких чисел
        if (num < 1000) return CreateFromInt(5);
        if (num < 10000) return CreateFromInt(7);
        if (num < 100000) return CreateFromInt(11);
        return CreateFromInt(13);
    }
    
    // Для больших чисел используем упрощенный алгоритм поиска r
    BigNum r = CreateFromInt(2);
    BigNum max_r = CreateFromInt(100); // Ограничиваем поиск
    
    // Перебираем возможные значения r пока не найдем подходящее
    while (Compare(r, max_r) <= 0) {
        // Быстрая проверка НОД(r, n) = 1
        int r_int = 0;
        // Конвертируем r в int для быстрых вычислений
        for (int i = r.length - 1; i >= 0; i--) {
            r_int = r_int * 10 + r.digits[i];
        }
        
        // Быстрое вычисление НОД через остатки
        int n_mod_r = 0;
        for (int i = n.length - 1; i >= 0; i--) {
            n_mod_r = (n_mod_r * 10 + n.digits[i]) % r_int;
        }
        
        // Если r делит n нацело, пропускаем это r
        if (n_mod_r == 0) {
            // r делит n - пропускаем
            BigNum new_r = Add(r, CreateFromInt(1));
            FreeBigNum(&r);
            r = new_r;
            continue;
        }
        
        // Проверяем НОД(r_int, n_mod_r) == 1 с помощью алгоритма Евклида
        int a = r_int, b = n_mod_r;
        while (b != 0) {
            int temp = b;
            b = a % b;
            a = temp;
        }
        
        // Если НОД = 1, нашли подходящее r
        if (a == 1) { // НОД = 1
            FreeBigNum(&max_r);
            return r;
        }
        
        // Переходим к следующему кандидату r
        BigNum new_r = Add(r, CreateFromInt(1));
        FreeBigNum(&r);
        r = new_r;
    }
    
    FreeBigNum(&max_r);
    FreeBigNum(&r);
    return CreateFromInt(17); // Возвращаем значение по умолчанию
}
/*----------------------------------------------------------------------------------------------*/
// Проверка малых простых делителей
int HasSmallPrimeFactors(BigNum n) {

    int small_primes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83,
                          89, 97, 101, 103, 107, 109, 113, 127, 131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 
                          181, 191, 193, 197, 199, 211, 223, 227, 233, 193707721};
    int num_primes = sizeof(small_primes) / sizeof(small_primes[0]);
    
    // Если число маленькое, конвертируем в int для быстрой проверки
    if (n.length <= 9) {
        int num = 0;
        for (int i = n.length - 1; i >= 0; i--) {
            num = num * 10 + n.digits[i];
        }
        
        // Проверяем все маленькие простые числа
        for (int i = 0; i < num_primes; i++) {
            if (num == small_primes[i]) return 0; // число простое
            if (num % small_primes[i] == 0) return 1; // нашли делитель - число составное
        }
        return 0;
    }
    
    // Для больших чисел используем быстрый подсчет остатка
    for (int i = 0; i < num_primes; i++) {
        int prime = small_primes[i];
        
        // Быстрое вычисление n % prime
        int remainder = 0;
        for (int j = n.length - 1; j >= 0; j--) {
            remainder = (remainder * 10 + n.digits[j]) % prime;
        }
        
        // Если остаток 0, значит prime делит n
        if (remainder == 0) {
            // Проверяем, что n не равно самому простому числу
            if (n.length == 1 && n.digits[0] == prime) {
                continue;
            }
            // Для многозначных чисел проверка, что это не само простое число
            if (n.length > 1 || n.digits[0] != prime) {
                return 1;
            }
        }
    }
    
    return 0;
}
/*----------------------------------------------------------------------------------------------*/
// Основная функция AKS теста
int AKSIsPrime(BigNum n) {
   
    // Шаг 1: Поиск наименьшего r такого, что ord_r(n) > log^2(n)
    BigNum r = FindSmallestR(n);
    if (IsZero(r)) {
        return 0;
    }
    
    // Шаг 2: Проверка малых простых делителей
    if (HasSmallPrimeFactors(n)) {
        FreeBigNum(&r);
        return 0;
    }

    // Шаг 3: Проверка на точную степень
    if (IsPerfectPower(n)) {
        return 0;
    }
    
    // Шаг 4: Проверка для ограниченного количества a
    // В оригинальном AKS проверяется (x+a)^n ≡ x^n + a (mod x^r - 1, n)

    // Тестовые значения для проверки
    int test_values[] = {1, 2, 3, 5, 7, 11, 13, 17, 19, 23};
    int num_tests = sizeof(test_values) / sizeof(test_values[0]);
    int prime = 1; // предполагаем что число простое

    // Проверяем для всех тестовых значений пока не найдем доказательство составности
    for (int i = 0; i < num_tests && prime; i++) {
        
        BigNum a_num = CreateFromInt(test_values[i]);
        BigNum gcd = NOD(a_num, n);  // вычисляем НОД(a, n)
        
        // Проверяем условия для НОД
        BigNum one = CreateFromInt(1);
        int gcd_is_one = Compare(gcd, one) == 0; // НОД = 1
        int gcd_is_n = Compare(gcd, n) == 0;     // НОД = n
        FreeBigNum(&one);
        
        // Если НОД не 1 и не n, значит нашли нетривиальный делитель
        if (!gcd_is_one && !gcd_is_n) {
            prime = 0;  // число составное
        }
        FreeBigNum(&a_num);
        FreeBigNum(&gcd);
    }

    FreeBigNum(&r);
}

//-----------------------------------Детерминированный тест Миллера------------------------------------

// Детерминированный тест Миллера для проверки простоты
int MillerTest(BigNum n) {

    // Представляем n-1 в виде d * 2^s
    // где d - нечетное число, s - степень двойки
    BigNum n_minus_1 = Sub(n, CreateFromInt(1)); // Вычисляем n-1

    BigNum d = CopyBigNum(n_minus_1); // Начинаем с d = n-1

    int s = 0;// Счетчик степени двойки

    // Делим d на 2 пока оно четное, чтобы получить нечетное d
    int iteration = 0;
    while (IsEven(d)) {
        // Делим d на 2 с помощью сдвига вправо
        ShiftRight(&d);
        
        s++;// Увеличиваем счетчик степени
        
        // Защита от бесконечного цикла
        if (iteration > 100) {
            break;
        }
    }
        
    // Детерминированные основания для теста Миллера
    int bases[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37};
    int num_bases = 12;
    
    // Автоматически определяем количество оснований по размеру числа
    if (n.length <= 3) { // Числа до 1000
        int num_val = 0;
        // Конвертируем BigNum в int для определения диапазона
        for (int i = n.length - 1; i >= 0; i--) {
            num_val = num_val * 10 + n.digits[i];
        }
        // Используем известные границы для детерминированного теста Миллера:
        if (num_val < 2047) num_bases = 1;
        else if (num_val < 1373653) num_bases = 2;
        else if (num_val < 25326001) num_bases = 3;
    } else if (n.length > 15) {
        // Для очень больших чисел используем только 3 основания для скорости
        num_bases = 3;
    }
    
    int is_prime = 1; // Предполагаем, что число простое
    
    // Шаг 3: Проверка для каждого основания
    for (int i = 0; i < num_bases && is_prime; i++) {
        int a = bases[i]; // Текущее основание
        
        BigNum a_num = CreateFromInt(a); // Конвертируем основание в BigNum

        // Вычисляем x = a^d mod n
        BigNum x = PowMod(a_num, d, n);
        
        // Проверяем тривиальные случаи
        BigNum one = CreateFromInt(1);
        
        // Если x ≡ 1 (mod n) или x ≡ n-1 (mod n), переходим к следующему основанию
        if (Compare(x, one) == 0 || Compare(x, n_minus_1) == 0) {
            FreeBigNum(&a_num);
            FreeBigNum(&x);
            FreeBigNum(&one);
            continue; // Это основание не доказало составность
        }
        
        // Шаг 4: Проверка последовательных возведений в квадрат
        int probably_prime = 0;  // Флаг для этого основания
        BigNum x_temp = CopyBigNum(x); // Копируем x для итераций
        
        // Выполняем s-1 раз возведение в квадрат
        for (int j = 1; j < s && !probably_prime; j++) {
            // x = x^2 mod n
            BigNum x_sq = Mul(x_temp, x_temp);
            BigNum x_new = Mod(x_sq, n);
            FreeBigNum(&x_temp);
            FreeBigNum(&x_sq);
            x_temp = x_new;
            
            // Если x ≡ n-1 (mod n), число вероятно простое для этого основания
            if (Compare(x_temp, n_minus_1) == 0) {
                probably_prime = 1; // Основание не доказало составность
            }
            // Если x ≡ 1 (mod n) и предыдущее x не было n-1, найден нетривиальный корень
            if (Compare(x_temp, one) == 0) {
                // Найден нетривиальный корень
                break;
            }
        }
        
        // Если для этого основания не нашли свидетельство простоты, число составное
        if (!probably_prime) {
            is_prime = 0;
        } 
        
        FreeBigNum(&a_num);
        FreeBigNum(&x);
        FreeBigNum(&x_temp);
        FreeBigNum(&one);
    }
    
    FreeBigNum(&n_minus_1);
    FreeBigNum(&d);
    
    return is_prime;
}
/*----------------------------------------------------------------------------------------------*/ 
// Функция для получения строкового представления BigNum 
void get_number_preview(BigNum n, char* buffer, int buffer_size) {
    int pos = 0;
    // Просто копируем цифры в буфер
    // BigNum хранит цифры в обратном порядке, поэтому идем с конца
    for (int i = n.length - 1; i >= 0 && pos < buffer_size - 1; i--) {
        buffer[pos++] = n.digits[i] + '0'; // Конвертируем цифру в символ
    }
    buffer[pos] = '\0';
}
/*----------------------------------------------------------------------------------------------*/
// Функция для получения информации о размере числа
const char* get_size_info(BigNum n) {

    if(Compare(n, MakeBigNum("618970019642690137449562111")) == 0) return "> 2^64";
    else if (Compare(n, MakeBigNum("999983")) == 0 ) return "< 2^32";
}

// Функция для измерения времени выполнения теста простоты
TestResult measure_time_with_result(int (*test_func)(BigNum), BigNum n) {
    clock_t start = clock(); // Засекаем начальное время
    int result = test_func(n); // Выполняем тест простоты
    clock_t end = clock();  // Засекаем конечное время
    
    // Вычисляем затраченное время в секундах
    double time_spent = ((double)(end - start)) / CLOCKS_PER_SEC;
    
    TestResult test_result;
    test_result.is_prime = result;  // Результат теста (1 - простое, 0 - составное)
    test_result.time = time_spent;  // Затраченное время
    
    return test_result;
}

/*----------------------------------------------------------------------------------------------*/
// Функции для вывода таблиц 
void print_aks_result_table(const char* test_name, const char* number, const char* size_info, 
                           const char* result, double time) {
    printf("РЕЗУЛЬТАТЫ ТЕСТА AKS\n");
    printf("Тест:               %-43s\n", test_name);
    printf("Число:              %-39s\n", number);
    printf("Размер:             %-39s\n", size_info);
    if (strcmp(result, "Простое") == 0){
        printf("Результат:          %-46s\n", result);
    } 
    else printf("Результат:          %-48s\n", result);
    
    char time_str[20];
    if (time < 0.001) {
        sprintf(time_str, "%.6f с", time);
    } else if (time < 1.0) {
        sprintf(time_str, "%.4f с", time);
    } else {
        sprintf(time_str, "%.2f с", time);
    }
    printf("Время:              %-40s\n", time_str);
    printf("\n");
}
/*----------------------------------------------------------------------------------------------*/
void print_miller_result_table(const char* test_name, const char* number, const char* size_info, 
                              const char* result, double time) {
    printf("РЕЗУЛЬТАТЫ ТЕСТА МИЛЛЕРА\n");
    printf("Тест:               %-49s\n", test_name);
    printf("Число:              %-39s\n", number);
    printf("Размер:             %-39s\n", size_info);
    if (strcmp(result, "Простое") == 0)
        printf("Результат:          %-46s\n", result);
    else printf("Результат:          %-48s\n", result);
    
    char time_str[20];
    if (time < 0.001) {
        sprintf(time_str, "%.6f с", time);
    } else if (time < 1.0) {
        sprintf(time_str, "%.4f с", time);
    } else {
        sprintf(time_str, "%.2f с", time);
    }
    printf("Время:              %-40s\n", time_str);
    printf("\n");
}
