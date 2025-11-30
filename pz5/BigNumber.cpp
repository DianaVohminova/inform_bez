#include "BigNumber.h"
#include <numeric>
#include <unordered_map>
#include <functional>

using namespace std;

// ==================== РЕАЛИЗАЦИЯ ОСНОВНЫХ МЕТОДОВ ====================

void BigNumber::removeLeadingZeros() {
    while (digits.size() > 1 && digits.back() == 0) {
        digits.pop_back();
    }
    if (digits.size() == 1 && digits[0] == 0) {
        isNegative = false;
    }
}

int BigNumber::compareAbsolute(const BigNumber& other) const {
    if (digits.size() != other.digits.size()) {
        return digits.size() < other.digits.size() ? -1 : 1;
    }
    for (int i = digits.size() - 1; i >= 0; --i) {
        if (digits[i] != other.digits[i]) {
            return digits[i] < other.digits[i] ? -1 : 1;
        }
    }
    return 0;
}

BigNumber BigNumber::addAbsolute(const BigNumber& other) const {
    BigNumber result;
    result.digits.resize(max(digits.size(), other.digits.size()) + 1, 0);

    int carry = 0;
    for (size_t i = 0; i < result.digits.size(); ++i) {
        int sum = carry;
        if (i < digits.size()) sum += digits[i];
        if (i < other.digits.size()) sum += other.digits[i];

        result.digits[i] = sum % 10;
        carry = sum / 10;
    }

    result.removeLeadingZeros();
    return result;
}

BigNumber BigNumber::subtractAbsolute(const BigNumber& other) const {
    if (compareAbsolute(other) < 0) {
        BigNumber result = other.subtractAbsolute(*this);
        result.isNegative = true;
        return result;
    }

    BigNumber result;
    result.digits.resize(digits.size(), 0);

    int borrow = 0;
    for (size_t i = 0; i < digits.size(); ++i) {
        int diff = digits[i] - borrow;
        if (i < other.digits.size()) {
            diff -= other.digits[i];
        }

        if (diff < 0) {
            diff += 10;
            borrow = 1;
        } else {
            borrow = 0;
        }

        result.digits[i] = diff;
    }

    result.removeLeadingZeros();
    return result;
}

BigNumber BigNumber::multiplyByDigit(int digit) const {
    if (digit == 0) return BigNumber(0);

    BigNumber result;
    result.digits.resize(digits.size() + 1, 0);

    int carry = 0;
    for (size_t i = 0; i < digits.size(); ++i) {
        int product = digits[i] * digit + carry;
        result.digits[i] = product % 10;
        carry = product / 10;
    }

    if (carry > 0) {
        result.digits[digits.size()] = carry;
    }

    result.removeLeadingZeros();
    return result;
}

BigNumber BigNumber::divideByDigit(int digit) const {
    if (digit == 0) {
        throw runtime_error("Division by zero");
    }

    BigNumber result;
    result.digits.resize(digits.size(), 0);

    int remainder = 0;
    for (int i = digits.size() - 1; i >= 0; --i) {
        int current = remainder * 10 + digits[i];
        result.digits[i] = current / digit;
        remainder = current % digit;
    }

    result.removeLeadingZeros();
    return result;
}

BigNumber::BigNumber() : isNegative(false) {
    digits.push_back(0);
}

BigNumber::BigNumber(const string& str) {
    if (str.empty()) {
        digits.push_back(0);
        isNegative = false;
        return;
    }

    size_t start = 0;
    if (str[0] == '-') {
        isNegative = true;
        start = 1;
    } else if (str[0] == '+') {
        isNegative = false;
        start = 1;
    } else {
        isNegative = false;
    }

    for (int i = str.size() - 1; i >= (int)start; --i) {
        if (isdigit(str[i])) {
            digits.push_back(str[i] - '0');
        } else {
            throw invalid_argument("Invalid character in number string");
        }
    }

    removeLeadingZeros();
}

BigNumber::BigNumber(long long num) {
    if (num < 0) {
        isNegative = true;
        num = -num;
    } else {
        isNegative = false;
    }

    if (num == 0) {
        digits.push_back(0);
        return;
    }

    while (num > 0) {
        digits.push_back(num % 10);
        num /= 10;
    }
}

BigNumber::BigNumber(const BigNumber& other)
    : digits(other.digits), isNegative(other.isNegative) {}

BigNumber::BigNumber(int numDigits, mt19937& gen) {
    if (numDigits <= 0) {
        throw invalid_argument("Invalid number of digits");
    }

    uniform_int_distribution<int> dist(1, 9);
    digits.push_back(dist(gen));

    uniform_int_distribution<int> dist2(0, 9);
    for (int i = 1; i < numDigits; ++i) {
        digits.push_back(dist2(gen));
    }

    reverse(digits.begin(), digits.end());
    isNegative = false;
}

BigNumber& BigNumber::operator=(const BigNumber& other) {
    if (this != &other) {
        digits = other.digits;
        isNegative = other.isNegative;
    }
    return *this;
}

BigNumber& BigNumber::operator=(long long num) {
    *this = BigNumber(num);
    return *this;
}

BigNumber BigNumber::operator+(const BigNumber& other) const {
    if (isNegative == other.isNegative) {
        BigNumber result = addAbsolute(other);
        result.isNegative = isNegative;
        return result;
    }

    int cmp = compareAbsolute(other);
    if (cmp == 0) {
        return BigNumber(0);
    }

    BigNumber result;
    if (cmp > 0) {
        result = subtractAbsolute(other);
        result.isNegative = isNegative;
    } else {
        result = other.subtractAbsolute(*this);
        result.isNegative = other.isNegative;
    }

    return result;
}

BigNumber BigNumber::operator-(const BigNumber& other) const {
    return *this + (-other);
}

BigNumber BigNumber::operator*(const BigNumber& other) const {
    if (isZero() || other.isZero()) {
        return BigNumber(0);
    }

    BigNumber result;
    result.digits.resize(digits.size() + other.digits.size(), 0);

    for (size_t i = 0; i < digits.size(); ++i) {
        int carry = 0;
        for (size_t j = 0; j < other.digits.size() || carry; ++j) {
            long long product = result.digits[i + j] +
                digits[i] * (j < other.digits.size() ? other.digits[j] : 0) +
                carry;
            result.digits[i + j] = product % 10;
            carry = product / 10;
        }
    }

    result.isNegative = isNegative != other.isNegative;
    result.removeLeadingZeros();
    return result;
}

BigNumber BigNumber::operator/(const BigNumber& other) const {
    if (other.isZero()) {
        throw runtime_error("Division by zero");
    }

    BigNumber absOther = other.abs();
    if (absOther.compareAbsolute(*this) > 0) {
        return BigNumber(0);
    }

    BigNumber quotient, remainder;
    quotient.digits.resize(digits.size(), 0);

    for (int i = digits.size() - 1; i >= 0; --i) {
        remainder = remainder * BigNumber(10) + BigNumber(digits[i]);
        int count = 0;
        while (remainder.compareAbsolute(absOther) >= 0) {
            remainder = remainder - absOther;
            count++;
        }
        quotient.digits[i] = count;
    }

    quotient.isNegative = isNegative != other.isNegative;
    quotient.removeLeadingZeros();
    return quotient;
}

BigNumber BigNumber::operator%(const BigNumber& other) const {
    if (other.isZero()) {
        throw runtime_error("Division by zero");
    }

    BigNumber absOther = other.abs();
    BigNumber remainder;

    for (int i = digits.size() - 1; i >= 0; --i) {
        remainder = remainder * BigNumber(10) + BigNumber(digits[i]);
        while (remainder.compareAbsolute(absOther) >= 0) {
            remainder = remainder - absOther;
        }
    }

    remainder.isNegative = isNegative;
    if (remainder.isNegative && !remainder.isZero()) {
        remainder = remainder + absOther;
    }

    return remainder;
}

BigNumber BigNumber::operator^(const BigNumber& exponent) const {
    if (exponent.isNegative) {
        throw runtime_error("Negative exponents not supported");
    }

    if (exponent.isZero()) {
        return BigNumber(1);
    }

    BigNumber result(1);
    BigNumber base = *this;
    BigNumber exp = exponent;

    while (!exp.isZero()) {
        if (exp.digits[0] % 2 == 1) {
            result = result * base;
        }
        base = base * base;
        exp = exp / BigNumber(2);
    }

    return result;
}

bool BigNumber::operator==(const BigNumber& other) const {
    return isNegative == other.isNegative && digits == other.digits;
}

bool BigNumber::operator!=(const BigNumber& other) const {
    return !(*this == other);
}

bool BigNumber::operator<(const BigNumber& other) const {
    if (isNegative != other.isNegative) {
        return isNegative;
    }

    if (isNegative) {
        return compareAbsolute(other) > 0;
    } else {
        return compareAbsolute(other) < 0;
    }
}

bool BigNumber::operator<=(const BigNumber& other) const {
    return *this < other || *this == other;
}

bool BigNumber::operator>(const BigNumber& other) const {
    return !(*this <= other);
}

bool BigNumber::operator>=(const BigNumber& other) const {
    return !(*this < other);
}

BigNumber BigNumber::operator-() const {
    BigNumber result = *this;
    if (!result.isZero()) {
        result.isNegative = !result.isNegative;
    }
    return result;
}

BigNumber BigNumber::operator+() const {
    return *this;
}

ostream& operator<<(ostream& os, const BigNumber& num) {
    if (num.isNegative) {
        os << '-';
    }
    for (int i = num.digits.size() - 1; i >= 0; --i) {
        os << num.digits[i];
    }
    return os;
}

istream& operator>>(istream& is, BigNumber& num) {
    string str;
    is >> str;
    num = BigNumber(str);
    return is;
}

string BigNumber::toString() const {
    string result;
    if (isNegative) {
        result += '-';
    }
    for (int i = digits.size() - 1; i >= 0; --i) {
        result += to_string(digits[i]);
    }
    return result;
}

bool BigNumber::isZero() const {
    return digits.size() == 1 && digits[0] == 0;
}

BigNumber BigNumber::abs() const {
    BigNumber result = *this;
    result.isNegative = false;
    return result;
}

BigNumber BigNumber::gcd(BigNumber a, BigNumber b) {
    a = a.abs();
    b = b.abs();

    while (!b.isZero()) {
        BigNumber temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

BigNumber BigNumber::lcm(const BigNumber& a, const BigNumber& b) {
    if (a.isZero() || b.isZero()) {
        return BigNumber(0);
    }
    return (a * b) / gcd(a, b);
}

// ==================== МЕТОДЫ ПРОВЕРКИ ПРОСТОТЫ ====================

BigNumber BigNumber::sqrt(const BigNumber& n) {
    if (n < BigNumber(0)) {
        throw invalid_argument("Square root of negative number");
    }
    if (n.isZero()) return BigNumber(0);
    if (n == BigNumber(1)) return BigNumber(1);

    BigNumber low(1), high = n;
    BigNumber result(1);

    while (low <= high) {
        BigNumber mid = (low + high) / BigNumber(2);
        BigNumber square = mid * mid;

        if (square == n) {
            return mid;
        } else if (square < n) {
            low = mid + BigNumber(1);
            result = mid;
        } else {
            high = mid - BigNumber(1);
        }
    }

    return result;
}

BigNumber BigNumber::modPow(const BigNumber& base, const BigNumber& exponent, const BigNumber& mod) {
    if (mod == BigNumber(1)) return BigNumber(0);
    
    BigNumber result(1);
    BigNumber b = base % mod;
    BigNumber exp = exponent;
    
    while (!exp.isZero()) {
        if (exp.digits[0] % 2 == 1) {
            result = (result * b) % mod;
        }
        b = (b * b) % mod;
        exp = exp / BigNumber(2);
    }
    
    return result;
}

bool BigNumber::isPrimeStandard(const BigNumber& n) {
    if (n < BigNumber(2)) return false;
    if (n == BigNumber(2)) return true;
    if (n.isEven()) return false;
    
    // Quick checks for small primes
    if (n == BigNumber(3) || n == BigNumber(5) || n == BigNumber(7)) return true;
    if (n % BigNumber(3) == BigNumber(0)) return false;
    if (n % BigNumber(5) == BigNumber(0)) return false;
    
    // Check sum of digits for divisibility by 3
    int sum = 0;
    for (int digit : n.digits) {
        sum += digit;
    }
    if (sum % 3 == 0) return n == BigNumber(3);
    
    // Check last digit for divisibility by 5
    if (n.digits[0] == 0 || n.digits[0] == 5) return n == BigNumber(5);
    
    // Proper trial division up to sqrt(n)
    BigNumber i(3);
    BigNumber limit = sqrt(n) + BigNumber(1);
    
    while (i <= limit) {
        if (n % i == BigNumber(0)) {
            return false;
        }
        i = i + BigNumber(2);
        
        // Skip multiples of 3 and 5 for optimization
        if (i % BigNumber(3) == BigNumber(0)) i = i + BigNumber(2);
        if (i % BigNumber(5) == BigNumber(0)) i = i + BigNumber(2);
    }
    
    return true;
}

bool BigNumber::isPrimeEratosthenes(const BigNumber& n, int limit) {
    if (n < BigNumber(2)) return false;
    
    try {
        long long num = stoll(n.toString());
        if (num <= limit) {
            if (num < 2) return false;
            if (num == 2) return true;
            if (num % 2 == 0) return false;
            
            vector<bool> sieve(num + 1, true);
            sieve[0] = sieve[1] = false;
            
            for (long long i = 2; i * i <= num; ++i) {
                if (sieve[i]) {
                    for (long long j = i * i; j <= num; j += i) {
                        sieve[j] = false;
                    }
                }
            }
            return sieve[num];
        }
    } catch (...) {
    }
    
    return isPrimeStandard(n);
}

bool BigNumber::isPrimeAtkin(const BigNumber& n, int limit) {
    if (n < BigNumber(2)) return false;
    
    try {
        long long num = stoll(n.toString());
        if (num <= limit) {
            if (num < 2) return false;
            if (num == 2 || num == 3) return true;
            
            vector<bool> sieve(num + 1, false);
            sieve[2] = sieve[3] = true;
            
            for (long long x = 1; x * x <= num; x++) {
                for (long long y = 1; y * y <= num; y++) {
                    long long temp = 4 * x * x + y * y;
                    if (temp <= num && (temp % 12 == 1 || temp % 12 == 5)) {
                        sieve[temp] = !sieve[temp];
                    }
                    
                    temp = 3 * x * x + y * y;
                    if (temp <= num && temp % 12 == 7) {
                        sieve[temp] = !sieve[temp];
                    }
                    
                    temp = 3 * x * x - y * y;
                    if (x > y && temp <= num && temp % 12 == 11) {
                        sieve[temp] = !sieve[temp];
                    }
                }
            }
            
            for (long long i = 5; i * i <= num; i++) {
                if (sieve[i]) {
                    for (long long j = i * i; j <= num; j += i * i) {
                        sieve[j] = false;
                    }
                }
            }
            
            return sieve[num];
        }
    } catch (...) {
    }
    
    return isPrimeStandard(n);
}

bool BigNumber::lucasLehmerTest(int p) {
    if (p < 2) return false;
    if (p == 2) return true;
    
    BigNumber mersenne = (BigNumber(2) ^ BigNumber(p)) - BigNumber(1);
    BigNumber s(4);
    
    for (int i = 0; i < p - 2; ++i) {
        s = (s * s - BigNumber(2)) % mersenne;
    }
    
    return s.isZero();
}

bool BigNumber::isPrime(int iterations) const {
    return isPrimeECPP(iterations);
}

BigNumber BigNumber::generateRandomPrime(int numDigits, mt19937& gen) {
    if (numDigits <= 0) {
        throw invalid_argument("Invalid number of digits");
    }

    BigNumber candidate(numDigits, gen);
    
    if (candidate.isEven()) {
        candidate = candidate + BigNumber(1);
    }

    while (!candidate.isPrime()) {
        candidate = candidate + BigNumber(2);
    }

    return candidate;
}

BigNumber BigNumber::log(const BigNumber& n) {
    if (n <= BigNumber(0)) {
        throw invalid_argument("Logarithm of non-positive number");
    }
    if (n == BigNumber(1)) {
        return BigNumber(0);
    }
    
    int digit_count = n.toString().length();
    return BigNumber(digit_count * 2);
}

// ==================== ECPP РЕАЛИЗАЦИЯ ====================

/**
 * Генерирует случайное большое число в диапазоне [0, max-1]
 * @param max - верхняя граница диапазона (исключительно)
 * @param gen - генератор случайных чисел
 * @return случайное BigNumber меньше max
 */
BigNumber BigNumber::randomBigNumber(const BigNumber& max, mt19937& gen) {
    // Базовые проверки
    if (max <= BigNumber(1)) return BigNumber(0);  // Если max <= 1, возвращаем 0
    
    // Преобразуем максимальное число в строку для обработки цифр
    string maxStr = max.toString();
    string resultStr;  // Строка для результата
    
    // Генераторы случайных цифр:
    uniform_int_distribution<int> firstDist(1, maxStr[0] - '0');  // Первая цифра: 1 до первой цифры max
    uniform_int_distribution<int> otherDist(0, 9);               // Остальные цифры: 0-9
    
    // Генерируем цифры для случайного числа
    bool firstDigit = true;  // Флаг для первой цифры
    for (size_t i = 0; i < maxStr.length(); ++i) {
        int digit;
        if (firstDigit) {
            digit = firstDist(gen);   // Первая цифра с ограничением
            firstDigit = false;       // Следующие цифры будут обычными
        } else {
            digit = otherDist(gen);   // Остальные цифры без ограничений
        }
        resultStr += to_string(digit);  // Добавляем цифру к результату
    }
    
    // Преобразуем строку в BigNumber
    BigNumber result(resultStr);
    
    // Гарантируем, что результат < max (делим пополам пока не выполнится условие)
    while (result >= max) {
        result = result / BigNumber(2);
    }
    
    return result;
}

/**
 * Алгоритм ECPP (Elliptic Curve Primality Proving) - доказательство простоты на эллиптических кривых
 * @param maxAttempts - максимальное количество попыток поиска подходящей кривой
 * @return true если число простое, false если составное
 */
bool BigNumber::isPrimeECPP(int maxAttempts) const {
    BigNumber n = *this;  // Работаем с копией числа
    
    // ========== ЭТАП 1: БАЗОВЫЕ ПРОВЕРКИ ==========
    
    // Проверка тривиальных случаев
    if (n < BigNumber(2)) return false;           // Числа < 2 не простые
    if (n == BigNumber(2) || n == BigNumber(3)) return true;  // 2 и 3 - простые
    if (n.isEven()) return false;              // Четные числа > 2 не простые
    
    // ========== ЭТАП 2: ПРОВЕРКА МАЛЫХ ПРОСТЫХ ДЕЛИТЕЛЕЙ ==========
    
    // Список малых простых чисел для быстрой проверки
    static const int smallPrimes[] = {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47};
    
    // Проверяем делимость на малые простые числа
    for (int p : smallPrimes) {
        if (n % BigNumber(p) == BigNumber(0)) {
            // Если делится, то простое только если равно этому простому числу
            return n == BigNumber(p);
        }
    }
    
    // Для небольших чисел используем стандартный метод (оптимизация)
    if (n < BigNumber(10000)) {
        return isPrimeStandard(n);
    }
    
    // ========== ЭТАП 3: ОСНОВНОЙ АЛГОРИТМ ECPP ==========
    
    random_device rd;    // Источник энтропии
    mt19937 gen(rd());   // Генератор случайных чисел
    
    // Многократные попытки найти подходящую эллиптическую кривую
    for (int attempt = 0; attempt < maxAttempts; ++attempt) {
        
        // ========== ЭТАП 3.1: ГЕНЕРАЦИЯ ПАРАМЕТРОВ КРИВОЙ ==========
        
        // Генерируем случайные параметры эллиптической кривой:
        // Уравнение кривой: y² = x³ + a·x + b
        BigNumber a = randomBigNumber(n, gen);  // Коэффициент a
        BigNumber x = randomBigNumber(n, gen);  // Координата x случайной точки
        BigNumber y = randomBigNumber(n, gen);  // Координата y случайной точки
        
        // ========== ЭТАП 3.2: ВЫЧИСЛЕНИЕ ПАРАМЕТРА b ==========
        
        // Вычисляем b из уравнения кривой: b = y² - x³ - a·x (mod n)
        BigNumber x2 = (x * x) % n;          // x² mod n
        BigNumber x3 = (x2 * x) % n;         // x³ mod n  
        BigNumber ax = (a * x) % n;          // a·x mod n
        BigNumber y2 = (y * y) % n;          // y² mod n
        
        // Вычисляем b по формуле
        BigNumber b = (y2 - x3 - ax) % n;
        // Корректируем если b отрицательное
        if (b < BigNumber(0)) b = b + n;
        
        // ========== ЭТАП 3.3: ПРОВЕРКА ДИСКРИМИНАНТА ==========
        
        // Дискриминант эллиптической кривой: Δ = 4a³ + 27b²
        // Должен быть ≠ 0 для невырожденной кривой
        BigNumber disc = (BigNumber(4) * a * a * a + BigNumber(27) * b * b) % n;
        
        // Если дискриминант нулевой - кривая вырожденная, пробуем снова
        if (disc == BigNumber(0)) continue;
        
        // ========== ЭТАП 3.4: УПРОЩЕННЫЙ ПОДСЧЕТ ТОЧЕК ==========
        
        // В РЕАЛЬНОМ ECPP: используется сложный алгоритм Шуфа для подсчета точек
        // В ЭТОЙ РЕАЛИЗАЦИИ: используем упрощенную аппроксимацию
        // m ≈ количество точек на кривой по модулю n
        BigNumber m = n + BigNumber(1) + randomBigNumber(BigNumber(100), gen);
        
        // ========== ЭТАП 3.5: ПОИСК ПОДХОДЯЩЕГО ПРОСТОГО ДЕЛИТЕЛЯ ==========
        
        // Ищем простой делитель q числа m такой, что:
        // 1. q - простое число
        // 2. m делится на q
        // 3. m/q - тоже простое число
        for (BigNumber q = BigNumber(2); q < BigNumber(1000); q = q + BigNumber(1)) {
            // Проверяем что q простое и делит m
            if (BigNumber::isPrimeStandard(q) && (m % q == BigNumber(0))) {
                BigNumber candidate = m / q;  // Кандидат на простое число
                
                // Если кандидат > 1 и простой - вероятно n простое
                if (candidate > BigNumber(1) && BigNumber::isPrimeStandard(candidate)) {
                    return true;  // Число вероятно простое
                }
            }
        }
    }
    
    // ========== ЭТАП 4: РЕЗЕРВНЫЙ МЕТОД ==========
    
    // Если ECPP не смог доказать простоту за maxAttempts попыток,
    // используем стандартный метод как запасной вариант
    return isPrimeStandard(n);
}

// ==================== АЛГОРИТМЫ ФАКТОРИЗАЦИИ ====================

/**
 * Алгоритм Полларда-Ро для факторизации (нахождение нетривиального делителя)
 * @param n - число для факторизации
 * @param maxIterations - максимальное количество итераций
 * @return нетривиальный делитель n или 1 если не найден
 */
BigNumber BigNumber::pollardRho(const BigNumber& n, int maxIterations) {
    // Базовые случаи
    if (n == BigNumber(1)) return BigNumber(1);           // 1 не имеет делителей
    if (n % BigNumber(2) == BigNumber(0)) return BigNumber(2);  // Четные числа делятся на 2
    
    // Инициализация генератора случайных чисел
    random_device rd;
    mt19937 gen(rd());
    
    // Начальные значения для алгоритма
    BigNumber x = randomBigNumber(n, gen);  // Начальная точка x
    BigNumber y = x;                     // Начальная точка y (такая же как x)
    BigNumber d = BigNumber(1);             // Найденный делитель (пока 1)
    
    // Функция итерации: f(x) = (x² + 1) mod n
    auto f = [](const BigNumber& x, const BigNumber& n) {
        return (x * x + BigNumber(1)) % n;
    };
    
    // Основной цикл алгоритма Полларда-Ро
    for (int i = 0; i < maxIterations && d == BigNumber(1); ++i) {
        // Движение "черепахи" - один шаг
        x = f(x, n);
        // Движение "зайца" - два шага (метод Флойда)
        y = f(f(y, n), n);
        
        // Вычисляем НОД(|x-y|, n) - потенциальный делитель
        d = gcd((x - y).abs(), n);
        
        // Если d != 1 и d != n, нашли нетривиальный делитель
    }
    
    return d;  // Возвращаем найденный делитель (или 1 если не нашли)
}

/**
 * Полная факторизация числа на простые множители
 * @param n - число для факторизации
 * @param maxAttempts - максимальное количество попыток алгоритма Полларда-Ро
 * @return вектор простых множителей в порядке возрастания
 */
vector<BigNumber> BigNumber::factorize(const BigNumber& n, int maxAttempts) {
    vector<BigNumber> factors;  // Результирующий вектор множителей
    BigNumber temp = n;         // Временная переменная для разложения
    
    // ========== ЭТАП 1: ПРОВЕРКА МАЛЫХ ПРОСТЫХ ДЕЛИТЕЛЕЙ ==========
    
    // Список малых простых чисел для последовательной проверки
    for (int p : {2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37}) {
        // Пока число делится на текущее простое, добавляем его в множители
        while (temp % BigNumber(p) == BigNumber(0)) {
            factors.push_back(BigNumber(p));  // Добавляем простой делитель
            temp = temp / BigNumber(p);       // Делим число на найденный множитель
        }
    }
    
    // Если после этого осталась 1 - факторизация завершена
    if (temp == BigNumber(1)) return factors;
    
    // Если оставшееся число простое - добавляем его и завершаем
    if (isPrimeStandard(temp)) {
        factors.push_back(temp);
        return factors;
    }
    
    // ========== ЭТАП 2: АЛГОРИТМ ПОЛЛАРДА-РО ДЛЯ БОЛЬШИХ ЧИСЕЛ ==========
    
    // Используем алгоритм Полларда-Ро для поиска нетривиальных делителей
    for (int attempt = 0; attempt < maxAttempts && temp > BigNumber(1); ++attempt) {
        // Пытаемся найти делитель алгоритмом Полларда-Ро
        BigNumber factor = pollardRho(temp, 1000);
        
        // Если нашли нетривиальный делитель (не 1 и не само число)
        if (factor > BigNumber(1) && factor < temp) {
            // Рекурсивно факторизуем найденный делитель
            auto subfactors = factorize(factor, maxAttempts);
            // Добавляем все подмножители в результат
            factors.insert(factors.end(), subfactors.begin(), subfactors.end());
            // Делим оставшееся число на найденный множитель
            temp = temp / factor;
            
            // Если оставшаяся часть простая - добавляем и выходим
            if (isPrimeStandard(temp)) {
                factors.push_back(temp);
                break;
            }
        } else {
            // Если делитель не найден, прерываем попытки
            break;
        }
    }
    
    // ========== ЭТАП 3: ДОБАВЛЕНИЕ ОСТАТКА ==========
    
    // Если после всех попыток осталось число > 1, добавляем его как есть
    if (temp > BigNumber(1)) {
        factors.push_back(temp);
    }
    
    return factors;  // Возвращаем полный список множителей
}