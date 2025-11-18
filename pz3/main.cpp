#include <iostream>
#include <vector>
#include <string>
#include <algorithm>
#include <fstream>
#include <stdexcept>
#include <random>
#include <chrono>

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

    static BigNumber binaryGCD(BigNumber a, BigNumber b) {
        a.isNegative = false;
        b.isNegative = false;
        if (a == BigNumber(0)) return b;
        if (b == BigNumber(0)) return a;
        BigNumber two(2);
        BigNumber shift(0);
        while (a.isEven() && b.isEven()) {
            a = a.divideByTwo();
            b = b.divideByTwo();
            shift = shift + BigNumber(1);
        }
        while (a.isEven()) {
            a = a.divideByTwo();
        }
        do {
            while (b.isEven()) {
                b = b.divideByTwo();
            }
            if (a > b) {
                std::swap(a, b);
            }
            b = b - a;
        } while (b != BigNumber(0));
        return a * two.pow(shift);
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

    void writeToFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << toString();
            file.close();
        } else {
            throw std::runtime_error("Не удалось открыть файл: " + filename);
        }
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

    BigNumber pow(const BigNumber& exponent) const {
        if (exponent.isNegative) {
            throw std::invalid_argument("Отрицательный показатель не поддерживается");
        }
        BigNumber result(1);
        BigNumber base = *this;
        BigNumber exp = exponent;
        while (exp > BigNumber(0)) {
            if (exp.digits[0] % 2 == 1) {
                result = result * base;
            }
            base = base * base;
            exp = exp / BigNumber(2);
        }
        return result;
    }

    // === НОВЫЕ МЕТОДЫ ДЛЯ ТЕСТА МИЛЛЕРА–РАБИНА ===

    BigNumber modPow(const BigNumber& exponent, const BigNumber& modulus) const {
        if (modulus == BigNumber(0)) {
            throw std::invalid_argument("Модуль не может быть нулём");
        }
        BigNumber base = (*this) % modulus;
        BigNumber exp = exponent;
        BigNumber result(1);
        while (exp > BigNumber(0)) {
            if (exp.digits[0] % 2 == 1) {
                result = (result * base) % modulus;
            }
            base = (base * base) % modulus;
            exp = exp / BigNumber(2);
        }
        return result;
    }

    bool millerRabinTest(int k = 10) const {
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
            BigNumber a(2 + i);
            if (a >= *this - BigNumber(1)) {
                break;
            }

            BigNumber x = a.modPow(d, *this);
            if (x == BigNumber(1) || x == n_minus_1) {
                continue;
            }

            bool composite = true;
            for (int r = 1; r < s; r++) {
                x = (x * x) % (*this);
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

    // === КОНЕЦ НОВЫХ МЕТОДОВ ===

    static BigNumber gcd(BigNumber a, BigNumber b) {
        return binaryGCD(a, b);
    }

    static BigNumber lcm(const BigNumber& a, const BigNumber& b) {
        BigNumber g = gcd(a, b);
        return (a * b) / g;
    }

    static BigNumber generateRandomNumber(size_t length) {
        if (length == 0) {
            return BigNumber(0);
        }
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_int_distribution<> dis(0, 9);
        std::string numberStr;
        numberStr += std::to_string(dis(gen) % 9 + 1);
        for (size_t i = 1; i < length; i++) {
            numberStr += std::to_string(dis(gen));
        }
        return BigNumber(numberStr);
    }

    friend void swap(BigNumber& first, BigNumber& second) {
        using std::swap;
        swap(first.digits, second.digits);
        swap(first.isNegative, second.isNegative);
    }
};

int main() {
    setlocale(LC_ALL, "ru");
    try {
        std::cout << "Практическая работа №3: Тест Миллера–Рабина\n";
        std::cout << "================================================\n";

        // Известное простое число: M(127) = 2^127 - 1
        BigNumber knownPrime("170141183460469231731687303715884105727");

        std::cout << "🔹 Тестируемое число:\n";
        std::cout << "   Значение: " << knownPrime.toString() << "\n";
        std::cout << "   Описание: число Мерсенна M(127) = 2^127 - 1 (известное простое)\n";
        std::cout << "   Количество цифр: " << knownPrime.size() << "\n\n";

        const int trials = 100;
        const int roundsPerTrial = 10;

        // Открываем файл для записи лога прогонов
        std::ofstream logFile("miller_rabin_trials.log");
        if (!logFile.is_open()) {
            throw std::runtime_error("Не удалось создать файл лога: miller_rabin_trials.log");
        }

        std::cout << " Запуск " << trials << " итераций теста Миллера–Рабина...\n";
        std::cout << "   (Каждая итерация использует " << roundsPerTrial << " баз)\n";

        int falseNegatives = 0;
        for (int i = 1; i <= trials; i++) {
            bool isProbablyPrime = knownPrime.millerRabinTest(roundsPerTrial);
            if (!isProbablyPrime) {
                falseNegatives++;
                logFile << "Прогон " << i << ": FAIL\n";
            } else {
                logFile << "Прогон " << i << ": PASS\n";
            }
        }
        logFile.close();
        std::cout << "Тестирование завершено. Лог сохранён в miller_rabin_trials.log\n\n";

        // Вывод итоговой статистики
        std::cout << "РЕЗУЛЬТАТЫ ТЕСТА\n";
        std::cout << "─────────────────────────────────────\n";
        std::cout << "Всего прогонов:        " << trials << "\n";
        std::cout << "Баз на прогон:         " << roundsPerTrial << "\n";
        std::cout << "Ложных отрицаний:      " << falseNegatives << "\n";
        std::cout << "Точность:              " << (100.0 - (100.0 * falseNegatives / trials)) << "%\n";
        std::cout << "Вывод:                 ";
        if (falseNegatives == 0) {
            std::cout << "Число подтверждено как простое во всех прогонах.\n";
        } else {
            std::cout << "Обнаружены отклонения — возможна ошибка в реализации.\n";
        }

        // Сохраняем краткий итог в отдельный файл
        std::ofstream summary("miller_rabin_results.txt");
        if (summary.is_open()) {
            summary << "Число: " << knownPrime.toString() << "\n";
            summary << "Прогонов: " << trials << "\n";
            summary << "Ложных отрицаний: " << falseNegatives << "\n";
            summary << "Вывод: " << (falseNegatives == 0 ? "Простое" : "Ошибки") << "\n";
            summary.close();
        }

    } catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}