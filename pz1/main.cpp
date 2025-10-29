#include <iostream>
#include <vector>
#include <string>
#include <algorithm> // для функций max, swap
#include <fstream>
#include <stdexcept> // исключения для работы с ошибками
#include <random>
#include <chrono>

class BigNumber {
private:
    std::vector<int> digits;  // вектор для хранения цифр числа (младшие разряды в начале)
    bool isNegative;          // флаг отрицательности числа

    // Метод для удаления ведущих нулей
    void removeLeadingZeros() {
        while (digits.size() > 1 && digits.back() == 0) {
            digits.pop_back();
        }
        if (digits.size() == 1 && digits[0] == 0) {
            isNegative = false; // число 0 неотрицательно
        }
    }

    // Вспомогательный метод для сложения абсолютных значений двух чисел
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

    // Вспомогательный метод для вычитания абсолютных значений (a >= b)
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
            }
            else {
                borrow = 0;
            }
            result.digits.push_back(sub);
        }

        result.removeLeadingZeros();
        return result;
    }

    // Вспомогательный метод для сравнения абсолютных значений двух чисел
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

    // Бинарный алгоритм Евклида (алгоритм Штейна)
    // Использует следующие свойства для вычисления НОД:
    /*
    * 1. gcd(0, v) = v, gcd(u, 0) = u
    * 2. gcd(u, v) = 2 * gcd(u/2, v/2), если u,v - четные
    * 3. gcd(u, v) = gcd(u/2, v), если u-четное, v-нечетное
    * 4. gcd(u, v) = gcd(|u-v|, min(u, v)), если u, v - нечетные
    */
    static BigNumber binaryGCD(BigNumber a, BigNumber b) {
        a.isNegative = false;
        b.isNegative = false;

        // Базовые случаи: gcd(n, 0) = gcd(0, n) = n
        if (a == BigNumber(0)) return b;
        if (b == BigNumber(0)) return a;

        BigNumber two(2); // константа 2 будет часто использоваться
        BigNumber shift(0); // счетчик для учета общего деления на 2, которое нужно будет вернуть

        // Пока оба числа четные
        while (a.isEven() && b.isEven()) {
            a = a.divideByTwo(); // делим оба на 2
            b = b.divideByTwo();
            shift = shift + BigNumber(1); // увеличиваем счетчик сдвига +1
        }
        // (убираем общий множитель 2: запоминаем, сколько раз разделили на 2 в shift)

        while (a.isEven()) { // если a четное, делим на 2
            a = a.divideByTwo();
        }

        do { // если b четное, делим на 2
            while (b.isEven()) {
                b = b.divideByTwo();
            }

            // Убеждаемся что a <= b, если нет - меняем местами
            if (a > b) {
                std::swap(a, b);
            }
            // 4 свойство: gcd(u, v) = gcd(u, v-u) если u<=v и u,v нечетные
            b = b - a;
        } while (b != BigNumber(0)); // пока не получим gcd(u, 0) = u

        // 3 свойство: gcd(u, 2^j v) = gcd(u, v), если u нечетное
        // 2 не является общим делителем
        return a * two.pow(shift); 
    }

public:
    // Конструктор по умолчанию
    // Создает число 0
    BigNumber() : isNegative(false) {
        digits.push_back(0);
    }

    // Конструктор из строки
    // Преобразует строковое представление числа в BigNumber
    BigNumber(const std::string& s) {
        if (s.empty()) { // если строка пустая
            digits.push_back(0); // создаем число 0
            isNegative = false;
            return;
        }

        isNegative = (s[0] == '-'); // определяем знак числа
        // Пропускаем знак в начале строки (в отрицательных числах первый символ '-')
        size_t start = (isNegative || s[0] == '+') ? 1 : 0;

        // Обрабатываем цифры с конца строки (младшие разряды идут первыми)
        for (int i = s.size() - 1; i >= static_cast<int>(start); i--) {
            if (isdigit(s[i])) {
                // Преобразуем символ в цифру и добавляем в вектор
                digits.push_back(s[i] - '0');
            }
            else {
                throw std::invalid_argument("Недопустимый символ в строке числа");
            }
        }
        // Удаляем ведущие нули, которые могли появиться при разборе
        removeLeadingZeros();
        if (digits.empty()) {
            digits.push_back(0);
            isNegative = false;
        }
    }

    // Конструктор из целого числа
    BigNumber(int n) {
        isNegative = n < 0; // определяем знак числа
        n = std::abs(n); // берем модуль числа

        if (n == 0) {
            digits.push_back(0);
        }
        else { // если число не ноль, разбираем его по цифрам
            // (младшие разряды идут первыми в векторе)
            while (n > 0) {
                digits.push_back(n % 10);
                n /= 10;
            }
        }
    }

    // Метод для получения количества цифр
    size_t size() const {
        return digits.size();
    }

    // Проверка на четность
    bool isEven() const {
        return (digits[0] % 2) == 0;
    }

    // Деление на 2
    // Работает с вектором цифр, т.е. с младшими разрядами в начале
    BigNumber divideByTwo() const {
        BigNumber result;
        result.digits.resize(digits.size(), 0);
        result.isNegative = isNegative;

        int carry = 0; // перенос из старшего разряда
        // Идем от старших разрядов к младшим
        for (int i = digits.size() - 1; i >= 0; i--) {
            // Текущее значение: цифра + перенос*10
            int current = digits[i] + carry * 10;
            // Вычисляем цифру результата
            result.digits[i] = current / 2;
            //Сохраняем остаток от деления в перенос
            carry = current % 2;
        }
        // Удаляем ведущие нули после деления
        result.removeLeadingZeros();
        return result;
    }

    // Преобразование в строку
    std::string toString() const {
        std::string result;
        // Если число отрицательное (и не ноль), добавляем знак минус
        if (isNegative && !(digits.size() == 1 && digits[0] == 0)) {
            result += '-';
        }

        // Идем от старших разрядов к младшим (конец вектора к началу)
        for (int i = digits.size() - 1; i >= 0; i--) {
            result += std::to_string(digits[i]);
        }

        return result;
    }

    // Запись в файл
    void writeToFile(const std::string& filename) const {
        std::ofstream file(filename);
        if (file.is_open()) {
            file << toString();
            file.close();
        }
        else {
            throw std::runtime_error("Не удалось открыть файл: " + filename);
        }
    }

    // Операторы сравнения
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

    // Арифметические операции
    BigNumber operator+(const BigNumber& other) const {
        // Случай 1: оба числа имеют одинаковый знак
        if (isNegative == other.isNegative) {
            // Складываем абсолютные значения чисел
            BigNumber result = addAbsolute(*this, other);
            result.isNegative = isNegative;
            return result;
        }

        // Случай 2: числа имеют разные знаки
        // Сравниваем абсолютные значения чисел
        int cmp = compareAbsolute(*this, other);

        // Если абсолютные значения равны - возвращаем ноль
        if (cmp == 0) {
            return BigNumber(0);
        }

        BigNumber result;
        // Первое число больше по модулю
        if (cmp > 0) {
            // Вычитаем из большего числа меньшее (по модулю)
            result = subtractAbsolute(*this, other);
            result.isNegative = isNegative;
        }

        // Второе число больше по модулю
        else {
            // Вычитаем из большего числа меньшее (по модулю)
            result = subtractAbsolute(other, *this);
            result.isNegative = other.isNegative;
        }

        return result;
    }

    BigNumber operator-(const BigNumber& other) const {
        // Случай 1: числа имеют разные знаки
        if (isNegative != other.isNegative) {
            // Сложение чисел с учетом знака дает вычитание абсолютных значений
            BigNumber result = addAbsolute(*this, other);
            // Устанавливаем знак равным знаку первого числа (this)
            result.isNegative = isNegative;
            return result;
        }

        // Случай 2: числа имеют одинаковые знаки
        // Сравниваем абсолютные значения чисел
        int cmp = compareAbsolute(*this, other);

        // Если абсолютные значения равны - возвращаем ноль
        if (cmp == 0) {
            return BigNumber(0);
        }

        BigNumber result;
        // Первое число больше по модулю
        if (cmp > 0) {
            // Вычитаем из большего числа меньшее (по модулю)
            result = subtractAbsolute(*this, other);
            // Сохраняем знак исходных чисел
            result.isNegative = isNegative;
        }
        // Второе число больше по модулю
        else {
            // Вычитаем из большего числа меньшее (по модулю)
            result = subtractAbsolute(other, *this);
            // Устанавливаем противоположный знак исходным числам
            result.isNegative = !isNegative;
        }

        return result;
    }

    BigNumber operator*(const BigNumber& other) const {
        // Используем метод умножения в столбик, временная сложность O(n*m)
        // Результирующий вектор имеет размер = сумма размеров множителей
        BigNumber result;
        result.digits.resize(digits.size() + other.digits.size(), 0);

        // Выполняем умножение "в столбик"
        for (size_t i = 0; i < digits.size(); i++) {
            int carry = 0;  // перенос на следующий разряд

            // Умножаем i-ю цифру первого числа на все цифры второго числа
            for (size_t j = 0; j < other.digits.size(); j++) {
                // Вычисляем произведение цифр + значение в текущей позиции + перенос
                int product = digits[i] * other.digits[j] + result.digits[i + j] + carry;

                // Вычисляем перенос для следующей итерации
                carry = product / 10;

                // Записываем остаток от деления на 10 в текущую позицию результата
                result.digits[i + j] = product % 10;
            }

            // Если после умножения на все цифры остался перенос,
            // добавляем его в следующий разряд
            if (carry > 0) {
                result.digits[i + other.digits.size()] += carry;
            }
        }

        // Определяем знак результата: 
        // произведение отрицательно, если знаки множителей различны
        result.isNegative = isNegative != other.isNegative;

        // Удаляем ведущие нули, которые могли появиться при умножении
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

    // Возведение в степень n с использованием алгоритма быстрого возведения 
    // в степень (exponentiation by squaring), который позволяет вычислять 
    // степень за логарифмическое время.
    BigNumber pow(const BigNumber& exponent) const {
        if (exponent.isNegative) {
            throw std::invalid_argument("Отрицательный показатель не поддерживается");
        }

        BigNumber result(1); // результат
        BigNumber base = *this; // основание
        BigNumber exp = exponent; // показатель степени

        while (exp > BigNumber(0)) {
            if (exp.digits[0] % 2 == 1) {
                result = result * base;
            }
            base = base * base;
            exp = exp / BigNumber(2); // целочисленное деление
        }

        return result;
    }

    // НОД с использованием бинарного алгоритма Евклида
    static BigNumber gcd(BigNumber a, BigNumber b) {
        return binaryGCD(a, b);
    }

    // НОК двух чисел (по формуле)
    static BigNumber lcm(const BigNumber& a, const BigNumber& b) {
        BigNumber g = gcd(a, b);
        return (a * b) / g;
    }

    // Генерация случайного числа
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

    // Функция обмена
    friend void swap(BigNumber& first, BigNumber& second) {
        using std::swap;
        swap(first.digits, second.digits);
        swap(first.isNegative, second.isNegative);
    }
};

int main() {
    setlocale(LC_ALL, "ru");

    try {
        std::cout << "Генерация больших чисел..." << std::endl;
        BigNumber a = BigNumber::generateRandomNumber(10000);
        BigNumber b = BigNumber::generateRandomNumber(10000);

        std::cout << "Число a сгенерировано, количество цифр: " << a.size() << std::endl;
        std::cout << "Число b сгенерировано, количество цифр: " << b.size() << std::endl;

        a.writeToFile("a.txt");
        b.writeToFile("b.txt");

        std::cout << "Выполнение арифметических операций..." << std::endl;

        BigNumber sum = a + b;
        sum.writeToFile("sum.txt");

        BigNumber diff = a - b;
        diff.writeToFile("diff.txt");

        BigNumber product = a * b;
        product.writeToFile("product.txt");

        BigNumber quotient = a / b;
        quotient.writeToFile("quotient.txt");

        BigNumber power = a.pow(BigNumber(2));
        power.writeToFile("power.txt");

        // Вычисление НОД с замером времени
        std::cout << "Вычисление НОД..." << std::endl;
        auto start_gcd = std::chrono::high_resolution_clock::now();
        BigNumber gcd = BigNumber::gcd(a, b);
        auto end_gcd = std::chrono::high_resolution_clock::now();
        auto duration_gcd = std::chrono::duration_cast<std::chrono::milliseconds>(end_gcd - start_gcd);
        gcd.writeToFile("gcd.txt");

        // Вычисление НОК с замером времени
        std::cout << "Вычисление НОК..." << std::endl;
        auto start_lcm = std::chrono::high_resolution_clock::now();
        BigNumber lcm = BigNumber::lcm(a, b);
        auto end_lcm = std::chrono::high_resolution_clock::now();
        auto duration_lcm = std::chrono::duration_cast<std::chrono::milliseconds>(end_lcm - start_lcm);
        lcm.writeToFile("lcm.txt");

        std::cout << "Все операции завершены. Результаты сохранены в файлы." << std::endl;

        std::cout << "Размеры результатов:" << std::endl;
        std::cout << "  Сумма: " << sum.size() << " цифр" << std::endl;
        std::cout << "  Разность: " << diff.size() << " цифр" << std::endl;
        std::cout << "  Произведение: " << product.size() << " цифр" << std::endl;
        std::cout << "  Частное: " << quotient.size() << " цифр" << std::endl;
        std::cout << "  Квадрат: " << power.size() << " цифр" << std::endl;
        std::cout << "  НОД: " << gcd.size() << " цифр" << std::endl;
        std::cout << "  НОК: " << lcm.size() << " цифр" << std::endl;

        // Время выполнения операций
        std::cout << "\nВремя выполнения операций:" << std::endl;
        std::cout << "  НОД: " << duration_gcd.count()*0.001 << " с" << std::endl;
        std::cout << "  НОК: " << duration_lcm.count() * 0.001 << " с" << std::endl;

    }
    catch (const std::exception& e) {
        std::cerr << "Ошибка: " << e.what() << std::endl;
        return 1;
    }

    return 0;
}