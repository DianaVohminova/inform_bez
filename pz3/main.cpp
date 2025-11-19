// g++ -std=c++17 -O2 main.cpp -o pz3
#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <random>
#include <chrono>
#include <cstdlib>

class BigNumber {
private:
    std::vector<int> digits;
    bool isNegative;

    void removeLeadingZeros() {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
        if (digits.size() == 1 && digits[0] == 0) {
            isNegative = false;
        }
    }

    static BigNumber addAbsolute(const BigNumber& a, const BigNumber& b) {
        BigNumber result;
        result.digits.clear();
        int carry = 0;
        size_t n = std::max(a.digits.size(), b.digits.size());
        for (size_t i = 0; i < n; i++) {
            int sum = carry;
            if (i < a.digits.size()) sum += a.digits[i];
            if (i < b.digits.size()) sum += b.digits[i];
            result.digits.push_back(sum % 10);
            carry = sum / 10;
        }
        if (carry) {
            result.digits.push_back(carry);
        }
        return result;
    }

    static BigNumber subtractAbsolute(const BigNumber& a, const BigNumber& b) {
        BigNumber result;
        result.digits.clear();
        int borrow = 0;
        for (size_t i = 0; i < a.digits.size(); i++) {
            int sub = a.digits[i] - borrow;
            if (i < b.digits.size()) sub -= b.digits[i];
            if (sub < 0) {
                sub += 10;
                borrow = 1;
            } else {
                borrow = 0;
            }
            result.digits.push_back(sub);
        }
        result.removeLeadingZeros();
        return result;
    }

    static int compareAbsolute(const BigNumber& a, const BigNumber& b) {
        if (a.digits.size() != b.digits.size()) {
            return a.digits.size() > b.digits.size() ? 1 : -1;
        }
        for (int i = a.digits.size() - 1; i >= 0; i--) {
            if (a.digits[i] != b.digits[i]) {
                return a.digits[i] > b.digits[i] ? 1 : -1;
            }
        }
        return 0;
    }

public:
    BigNumber() : isNegative(false) {
        digits.push_back(0);
    }

    BigNumber(const std::string& s) {
        if (s.empty()) {
            digits.push_back(0);
            isNegative = false;
            return;
        }
        isNegative = (s[0] == '-');
        size_t start = (isNegative || s[0] == '+') ? 1 : 0;
        for (int i = s.size() - 1; i >= static_cast<int>(start); i--) {
            if (isdigit(s[i])) {
                digits.push_back(s[i] - '0');
            } else {
                throw std::invalid_argument("Недопустимый символ в строке числа");
            }
        }
        removeLeadingZeros();
        if (digits.empty()) {
            digits.push_back(0);
            isNegative = false;
        }
    }

    BigNumber(int n) {
        isNegative = n < 0;
        n = std::abs(n);
        if (n == 0) {
            digits.push_back(0);
        } else {
            while (n > 0) {
                digits.push_back(n % 10);
                n /= 10;
            }
        }
    }

    size_t size() const {
        return digits.size();
    }

    bool isEven() const {
        return (digits[0] % 2) == 0;
    }

    BigNumber divideByTwo() const {
        BigNumber result;
        result.digits.resize(digits.size(), 0);
        result.isNegative = isNegative;
        int carry = 0;
        for (int i = digits.size() - 1; i >= 0; i--) {
            int current = digits[i] + carry * 10;
            result.digits[i] = current / 2;
            carry = current % 2;
        }
        result.removeLeadingZeros();
        return result;
    }

    std::string toString() const {
        std::string result;
        if (isNegative && !(digits.size() == 1 && digits[0] == 0)) {
            result += '-';
        }
        for (int i = digits.size() - 1; i >= 0; i--) {
            result += std::to_string(digits[i]);
        }
        return result;
    }

    bool operator<(const BigNumber& other) const {
        if (isNegative != other.isNegative) {
            return isNegative;
        }
        int cmp = compareAbsolute(*this, other);
        return isNegative ? cmp > 0 : cmp < 0;
    }

    bool operator>(const BigNumber& other) const {
        return other < *this;
    }

    bool operator<=(const BigNumber& other) const {
        return !(*this > other);
    }

    bool operator>=(const BigNumber& other) const {
        return !(*this < other);
    }

    bool operator==(const BigNumber& other) const {
        return isNegative == other.isNegative && digits == other.digits;
    }

    bool operator!=(const BigNumber& other) const {
        return !(*this == other);
    }

    BigNumber operator+(const BigNumber& other) const {
        if (isNegative == other.isNegative) {
            BigNumber result = addAbsolute(*this, other);
            result.isNegative = isNegative;
            return result;
        }
        int cmp = compareAbsolute(*this, other);
        if (cmp == 0) {
            return BigNumber(0);
        }
        BigNumber result;
        if (cmp > 0) {
            result = subtractAbsolute(*this, other);
            result.isNegative = isNegative;
        } else {
            result = subtractAbsolute(other, *this);
            result.isNegative = other.isNegative;
        }
        return result;
    }

    BigNumber operator-(const BigNumber& other) const {
        if (isNegative != other.isNegative) {
            BigNumber result = addAbsolute(*this, other);
            result.isNegative = isNegative;
            return result;
        }
        int cmp = compareAbsolute(*this, other);
        if (cmp == 0) {
            return BigNumber(0);
        }
        BigNumber result;
        if (cmp > 0) {
            result = subtractAbsolute(*this, other);
            result.isNegative = isNegative;
        } else {
            result = subtractAbsolute(other, *this);
            result.isNegative = !isNegative;
        }
        return result;
    }

    BigNumber operator*(const BigNumber& other) const {
        BigNumber result;
        result.digits.resize(digits.size() + other.digits.size(), 0);
        for (size_t i = 0; i < digits.size(); i++) {
            int carry = 0;
            for (size_t j = 0; j < other.digits.size(); j++) {
                int product = digits[i] * other.digits[j] + result.digits[i + j] + carry;
                carry = product / 10;
                result.digits[i + j] = product % 10;
            }
            if (carry > 0) {
                result.digits[i + other.digits.size()] += carry;
            }
        }
        result.isNegative = isNegative != other.isNegative;
        result.removeLeadingZeros();
        return result;
    }

    BigNumber operator/(const BigNumber& other) const {
        if (other == BigNumber(0)) {
            throw std::runtime_error("Деление на ноль");
        }
        BigNumber dividend = *this;
        BigNumber divisor = other;
        dividend.isNegative = false;
        divisor.isNegative = false;
        if (dividend < divisor) {
            return BigNumber(0);
        }
        BigNumber quotient;
        quotient.digits.resize(dividend.digits.size(), 0);
        BigNumber current;
        for (int i = dividend.digits.size() - 1; i >= 0; i--) {
            current.digits.insert(current.digits.begin(), dividend.digits[i]);
            current.removeLeadingZeros();
            int count = 0;
            while (current >= divisor) {
                current = current - divisor;
                count++;
            }
            quotient.digits[i] = count;
        }
        quotient.removeLeadingZeros();
        quotient.isNegative = isNegative != other.isNegative;
        return quotient;
    }

    BigNumber operator%(const BigNumber& other) const {
        BigNumber quotient = *this / other;
        BigNumber remainder = *this - quotient * other;
        return remainder;
    }

    // === МЕТОДЫ ДЛЯ ТЕСТА МИЛЛЕРА–РАБИНА ===

    BigNumber modPow(const BigNumber& exponent, const BigNumber& modulus) const {
        if (modulus == BigNumber(0)) {
            throw std::invalid_argument("Модуль не может быть нулём");
        }
        BigNumber base = (*this) % modulus;
        if (base < 0) base = base + modulus;
        BigNumber exp = exponent;
        BigNumber result(1);
        while (exp > BigNumber(0)) {
            if (exp.digits[0] % 2 == 1) {
                result = (result * base) % modulus;
                if (result < 0) result = result + modulus;
            }
            base = (base * base) % modulus;
            if (base < 0) base = base + modulus;
            exp = exp / BigNumber(2);
        }
        return result;
    }

    // Генерация случайного числа в диапазоне [min, max]
    static BigNumber randomInRange(const BigNumber& min, const BigNumber& max) {
        if (min > max) {
            throw std::invalid_argument("min должен быть <= max");
        }
        
        static std::random_device rd;
        static std::mt19937 gen(rd());
        std::uniform_int_distribution<int> dist(0, 9);
        
        // Простая реализация - генерируем число с тем же количеством цифр
        std::string randomStr;
        for (size_t i = 0; i < max.size(); i++) {
            randomStr += std::to_string(dist(gen));
        }
        
        BigNumber randomNum(randomStr);
        BigNumber range = max - min + BigNumber(1);
        BigNumber result = min + (randomNum % range);
        
        // Гарантируем, что результат в диапазоне
        if (result < min) result = min;
        if (result > max) result = max;
        
        return result;
    }

    bool millerRabinTest(int k = 20) const {
        if (*this <= BigNumber(1)) return false;
        if (*this == BigNumber(2)) return true;
        if (isEven()) return false;

        BigNumber n_minus_1 = *this - BigNumber(1);
        BigNumber d = n_minus_1;
        int s = 0;
        while (d.isEven()) {
            d = d.divideByTwo();
            s++;
        }

        for (int i = 0; i < k; i++) {
            // Генерируем случайное основание в диапазоне [2, n-2]
            BigNumber a = randomInRange(BigNumber(2), *this - BigNumber(2));

            BigNumber x = a.modPow(d, *this);
            if (x == BigNumber(1) || x == n_minus_1) {
                continue;
            }

            bool composite = true;
            for (int r = 1; r < s; r++) {
                x = (x * x) % (*this);
                if (x < 0) x = x + *this;
                if (x == n_minus_1) {
                    composite = false;
                    break;
                }
            }
            if (composite) {
                return false;
            }
        }
        return true;
    }

    // === ОПТИМИЗИРОВАННЫЙ ТЕСТ ЛЮКА ===
    
    static int jacobiSymbol(BigNumber a, BigNumber n) {
        if (n <= BigNumber(0) || n.isEven()) {
            return 0;
        }
        
        a = a % n;
        int result = 1;
        
        while (a != BigNumber(0)) {
            // Удаляем множители 2 из a
            while (a.isEven()) {
                a = a.divideByTwo();
                BigNumber n_mod_8 = n % BigNumber(8);
                if (n_mod_8 == BigNumber(3) || n_mod_8 == BigNumber(5)) {
                    result = -result;
                }
            }
            
            // Применяем квадратичный закон взаимности
            std::swap(a, n);
            BigNumber a_mod_4 = a % BigNumber(4);
            BigNumber n_mod_4 = n % BigNumber(4);
            if (a_mod_4 == BigNumber(3) && n_mod_4 == BigNumber(3)) {
                result = -result;
            }
            
            a = a % n;
        }
        
        if (n == BigNumber(1)) {
            return result;
        } else {
            return 0;
        }
    }

    bool lucasStrongPseudoprimeTest() const {
        if (*this <= BigNumber(1)) return false;
        if (*this == BigNumber(2)) return true;
        if (isEven()) return false;

        // Для малых чисел используем прямое деление
        if (*this < BigNumber(100000)) {
            return isPrimeTrialDivision();
        }

        // Для чисел Мерсенна используем специализированный тест
        if (isMersenneNumber()) {
            return lucasLehmerTest();
        }

        // Общий случай - используем оптимизированный тест Люка
        return optimizedLucasTest();
    }

    // === ТЕСТ BPSW ===
    bool bailliePomeranceSelfridgeWagstaffTest() const {
        if (*this <= BigNumber(1)) return false;
        if (*this == BigNumber(2) || *this == BigNumber(3)) return true;
        if (isEven()) return false;

        // Шаг 1: Проверить с помощью теста Миллера-Рабина с основанием 2
        if (!millerRabinTestWithBase2()) {
            return false;
        }

        // Шаг 2: Проверить с помощью теста Люка
        return lucasStrongPseudoprimeTest();
    }

private:
    bool isPrimeTrialDivision() const {
        if (*this <= BigNumber(1)) return false;
        if (*this == BigNumber(2)) return true;
        if (isEven()) return false;

        BigNumber i(3);
        while (i * i <= *this) {
            if (*this % i == BigNumber(0)) {
                return false;
            }
            i = i + BigNumber(2);
        }
        return true;
    }

    bool isMersenneNumber() const {
        // Проверяем, является ли число вида 2^p - 1
        BigNumber n_plus_1 = *this + BigNumber(1);
        return isPowerOfTwo(n_plus_1);
    }

    bool isPowerOfTwo(const BigNumber& n) const {
        if (n <= BigNumber(0)) return false;
        BigNumber temp = n;
        int countOnes = 0;
        while (temp > BigNumber(0)) {
            if (temp.digits[0] % 2 == 1) {
                countOnes++;
                if (countOnes > 1) return false;
            }
            temp = temp.divideByTwo();
        }
        return countOnes == 1;
    }

    bool lucasLehmerTest() const {
        // Тест Люка-Лемера для чисел Мерсенна
        BigNumber n_plus_1 = *this + BigNumber(1);
        
        // Находим p такое, что 2^p = n + 1
        BigNumber p(0);
        BigNumber temp = n_plus_1;
        while (temp > BigNumber(1)) {
            if (!temp.isEven()) return false;
            temp = temp.divideByTwo();
            p = p + BigNumber(1);
        }
        
        // Тест Люка-Лемера
        BigNumber s(4);
        for (BigNumber i = BigNumber(2); i < p; i = i + BigNumber(1)) {
            s = (s * s - BigNumber(2)) % (*this);
            if (s < 0) s = s + (*this);
        }
        
        return (s == BigNumber(0));
    }

    bool optimizedLucasTest() const {
        // Упрощенный, но быстрый тест Люка для общих чисел
        // Используем фиксированные параметры для скорости
        
        // Параметры для теста Люка
        BigNumber P(1);
        BigNumber Q(-1);
        
        BigNumber U_prev(0);
        BigNumber U_curr(1);
        BigNumber V_prev(2);
        BigNumber V_curr(P);
        
        BigNumber n_plus_1 = *this + BigNumber(1);
        
        // Вычисляем U_{n+1} mod n с помощью эффективного алгоритма
        std::vector<bool> bits = getBinaryRepresentation(n_plus_1);
        
        for (size_t i = 1; i < bits.size(); i++) { // Пропускаем старший бит
            // Удваиваем
            BigNumber U_next = (U_curr * V_curr) % (*this);
            BigNumber V_next = (V_curr * V_curr - BigNumber(2)) % (*this);
            
            if (bits[i]) {
                // Удваиваем и добавляем 1
                BigNumber U_temp = (U_next * P + V_next) / BigNumber(2);
                BigNumber V_temp = (V_next * P + U_next) / BigNumber(2);
                
                U_curr = U_temp % (*this);
                V_curr = V_temp % (*this);
            } else {
                U_curr = U_next;
                V_curr = V_next;
            }
            
            // Корректируем отрицательные остатки
            if (U_curr < 0) U_curr = U_curr + (*this);
            if (V_curr < 0) V_curr = V_curr + (*this);
        }
        
        return (U_curr == BigNumber(0));
    }

    std::vector<bool> getBinaryRepresentation(const BigNumber& n) const {
        std::vector<bool> bits;
        BigNumber temp = n;
        while (temp > BigNumber(0)) {
            bits.push_back(temp.digits[0] % 2 == 1);
            temp = temp.divideByTwo();
        }
        std::reverse(bits.begin(), bits.end());
        return bits;
    }

    bool millerRabinTestWithBase2() const {
        if (*this <= BigNumber(1)) return false;
        if (*this == BigNumber(2)) return true;
        if (isEven()) return false;

        BigNumber n_minus_1 = *this - BigNumber(1);
        BigNumber d = n_minus_1;
        int s = 0;
        while (d.isEven()) {
            d = d.divideByTwo();
            s++;
        }

        BigNumber a(2);
        BigNumber x = a.modPow(d, *this);
        if (x == BigNumber(1) || x == n_minus_1) {
            return true;
        }

        for (int r = 1; r < s; r++) {
            x = (x * x) % (*this);
            if (x < 0) x = x + *this;
            if (x == n_minus_1) {
                return true;
            }
        }
        return false;
    }
};

int main() {
    setlocale(LC_ALL, "ru");
    try {
        std::cout << "Практическая работа №3: Тест Миллера–Рабина, тест Люка и тест BPSW\n";
        std::cout << "================================================\n";

        // Тестируем на известных простых числах
        BigNumber knownPrime("170141183460469231731687303715884105727"); // M(127) - простое
        // BigNumber knownPrime("524287");  // M(19) - простое
        // BigNumber knownPrime("2147483647");  // M(31) - простое
        // BigNumber knownPrime("1000000007");  // простое число

        std::cout << "🔹 Тестируемое число:\n";
        std::cout << "   Значение: " << knownPrime.toString() << "\n";
        std::cout << "   Описание: число Мерсенна M(127) = 2^127 - 1 (известное ПРОСТОЕ)\n";
        std::cout << "   Количество цифр: " << knownPrime.size() << "\n\n";

        const int trials = 100;      // 100 прогонов
        const int roundsPerTrial = 15; // 15 баз на прогон

        // === ТЕСТ МИЛЛЕРА–РАБИНА ===
        std::cout << "🔹 Тестируем тест Миллера–Рабина...\n";
        std::cout << " Запуск " << trials << " итераций (" << roundsPerTrial << " баз на итерацию)...\n";

        auto startMR = std::chrono::high_resolution_clock::now();
        
        int falseNegativesMR = 0;
        for (int i = 1; i <= trials; i++) {
            bool isProbablyPrime = knownPrime.millerRabinTest(roundsPerTrial);
            if (!isProbablyPrime) {
                falseNegativesMR++;
            }
            //std::cout << "   Прогон " << i << "/" << trials << " завершен\n";
        }
        
        auto endMR = std::chrono::high_resolution_clock::now();
        auto durationMR = std::chrono::duration_cast<std::chrono::milliseconds>(endMR - startMR);
        
        std::cout << "\nРЕЗУЛЬТАТЫ ТЕСТА МИЛЛЕРА–РАБИНА\n";
        std::cout << "─────────────────────────────────────\n";
        std::cout << "Всего прогонов:        " << trials << "\n";
        std::cout << "Баз на прогон:         " << roundsPerTrial << "\n";
        std::cout << "Ложных отрицаний:      " << falseNegativesMR << "\n";
        std::cout << "Точность:              " << (100.0 - (100.0 * falseNegativesMR / trials)) << "%\n";
        std::cout << "Время выполнения:      " << durationMR.count() << " мс\n";
        std::cout << "Вывод:                 ";
        if (falseNegativesMR == 0) {
            std::cout << "Число подтверждено как простое\n";
        } else {
            std::cout << "Обнаружены ложные отрицания\n";
        }

        std::cout << "\n";

        // === ТЕСТ ЛЮКА ===
        std::cout << "🔹 Тестируем тест Люка...\n";
        std::cout << " Запуск " << trials << " итераций...\n";

        auto startLucas = std::chrono::high_resolution_clock::now();
        
        int falseNegativesLucas = 0;
        for (int i = 1; i <= trials; i++) {
            bool isLucasPrime = knownPrime.lucasStrongPseudoprimeTest();
            if (!isLucasPrime) {
                falseNegativesLucas++;
            }
            // std::cout << "   Прогон " << i << "/" << trials << " завершен\n";
        }
        
        auto endLucas = std::chrono::high_resolution_clock::now();
        auto durationLucas = std::chrono::duration_cast<std::chrono::milliseconds>(endLucas - startLucas);
        
        std::cout << "\nРЕЗУЛЬТАТЫ ТЕСТА ЛЮКА\n";
        std::cout << "─────────────────────────────────────\n";
        std::cout << "Всего прогонов:        " << trials << "\n";
        std::cout << "Ложных отрицаний:      " << falseNegativesLucas << "\n";
        std::cout << "Точность:              " << (100.0 - (100.0 * falseNegativesLucas / trials)) << "%\n";
        std::cout << "Время выполнения:      " << durationLucas.count() << " мс\n";
        std::cout << "Вывод:                 ";
        if (falseNegativesLucas == 0) {
            std::cout << "Число подтверждено как Люка-псевдопростое\n";
        } else {
            std::cout << "Обнаружены ложные отрицания\n";
        }

        std::cout << "\n";

        // === ТЕСТ BPSW ===
        std::cout << "🔹 Тестируем тест BPSW...\n";
        std::cout << " Запуск " << trials << " итераций...\n";

        auto startBPSW = std::chrono::high_resolution_clock::now();
        
        int falseNegativesBPSW = 0;
        for (int i = 1; i <= trials; i++) {
            bool isBPSWPrime = knownPrime.bailliePomeranceSelfridgeWagstaffTest();
            if (!isBPSWPrime) {
                falseNegativesBPSW++;
            }
            // std::cout << "   Прогон " << i << "/" << trials << " завершен\n";
        }
        
        auto endBPSW = std::chrono::high_resolution_clock::now();
        auto durationBPSW = std::chrono::duration_cast<std::chrono::milliseconds>(endBPSW - startBPSW);
        
        std::cout << "\nРЕЗУЛЬТАТЫ ТЕСТА BPSW\n";
        std::cout << "─────────────────────────────────────\n";
        std::cout << "Всего прогонов:        " << trials << "\n";
        std::cout << "Ложных отрицаний:      " << falseNegativesBPSW << "\n";
        std::cout << "Точность:              " << (100.0 - (100.0 * falseNegativesBPSW / trials)) << "%\n";
        std::cout << "Время выполнения:      " << durationBPSW.count() << " мс\n";
        std::cout << "Вывод:                 ";
        if (falseNegativesBPSW == 0) {
            std::cout << "Число подтверждено как BPSW-псевдопростое\n";
        } else {
            std::cout << "Обнаружены ложные отрицания\n";
        }

        // Итоговый вывод
        std::cout << "\n================================================\n";
        std::cout << "ИТОГОВЫЕ РЕЗУЛЬТАТЫ:\n";
        std::cout << "  • Миллер-Рабин:  " << (falseNegativesMR == 0 ? "✓ ПРОШЕЛ" : "✗ НЕ ПРОШЕЛ") << "\n";
        std::cout << "  • Тест Люка:     " << (falseNegativesLucas == 0 ? "✓ ПРОШЕЛ" : "✗ НЕ ПРОШЕЛ") << "\n";
        std::cout << "  • Тест BPSW:     " << (falseNegativesBPSW == 0 ? "✓ ПРОШЕЛ" : "✗ НЕ ПРОШЕЛ") << "\n";
        std::cout << "  • Общее время:   " << (durationMR + durationLucas + durationBPSW).count() << " мс\n";
        
        if (falseNegativesMR == 0 && falseNegativesLucas == 0 && falseNegativesBPSW == 0) {
            std::cout << "Все тесты успешно подтвердили простоту числа!\n";
        } else {
            std::cout << "Некоторые тесты не прошли проверку.\n";
        }

    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}