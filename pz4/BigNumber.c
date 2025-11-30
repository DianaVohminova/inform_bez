
#include "BigNumber.h"
#include <stdio.h>
#include <ctype.h>
#include <time.h>
#include <math.h>
#include <complex.h>


// Создание BigNum длины length
BigNum CreateBigNum(int length) {
    BigNum num;
    num.length = length; // длина числа
    num.sign = 1;        // знак +

    // Выделение памяти под цифры
    num.digits = (int*)calloc(length, sizeof(int)); 
    if (!num.digits) {
        exit(1);
    }
    return num;
}

// Освобождение памяти, занятой BigNum
void FreeBigNum(BigNum *num) {
    // Если была выделена память под цифры
    if (num->digits) {
        free(num->digits); // освобождаем память
        num->digits = NULL;// обнуляем указатель
    }
    // Сбрасываем длину и знак
    num->length = 0;
    num->sign = 1;
}

// Копирование BigNum
BigNum CopyBigNum(BigNum a) {
    // Создаём новый BigNum такой же длины
    BigNum res = CreateBigNum(a.length);
    res.sign = a.sign; // копируем знак
    memcpy(res.digits, a.digits, a.length * sizeof(int));//копируем цифры
    return res;
}

// Удаление лишних нулей
void TrimBigNum(BigNum *num) {
    // Пока в старших разрядах нули и длина больше 1
    while (num->length > 1 && num->digits[num->length - 1] == 0) {
        num->length--; // уменьшаем длину
    }
    
    // Если число стало нулём
    if (num->length == 1 && num->digits[0] == 0) {
        num->sign = 1; // ноль всегда положительный
    }
}

// Генерация случайного BigNum
void GenerateRandomBigNum(BigNum *num, int length) {
    *num = CreateBigNum(length); // создаём число заданной длины
    num->sign = 1;               // знак +

    num->digits[length - 1] = rand() % 9 + 1; // старшая цифра не 0
    // Заполняем остальные цифры
    for (int i = length - 2; i >= 0; i--) {
        num->digits[i] = rand() % 10; // случайная цифра от 0 до 9
    }
}

// Преобразование строки в BigNum
BigNum MakeBigNum(const char *str) {
    BigNum res;
    res.sign = 1; // знак +

    // Если первый символ -
    if (str[0] == '-') {
        res.sign = -1; // знак -
        str++;         // пропускаем знак -
    }

    int n = strlen(str);   // длина строки
    res = CreateBigNum(n); // создаём число

    // Заполняем цифры в обратном порядке
    for (int i = 0; i < n; i++) {
        if (!isdigit(str[n - 1 - i])) continue; // пропускаем не цифры
        res.digits[i] = str[n - 1 - i] - '0';   // преобразуем символ в цифру
    }
    TrimBigNum(&res); // убираем ведущие нули
    return res;
}

// Печать BigNum в файл
void PrintBigNum(const char *filename, BigNum a) {
    FILE *f = fopen(filename, "w"); // открываем файл для записи
    if (!f) {
        printf("Ошибка: не удалось открыть файл %s\n", filename);
        exit(1);
    }

    // Если число отрицательное, печатаем знак -
    if (a.sign == -1) fprintf(f, "-");

    // Печатаем цифры от старших к младшим
    for (int i = a.length - 1; i >= 0; i--) {
        fprintf(f, "%d", a.digits[i]);
    }
    fprintf(f, "\n");
    fclose(f);
}

// Чтение BigNum из файла
BigNum ReadBigNum(const char *filename) {
    FILE *f = fopen(filename, "r"); // открываем файл для чтения
    if (!f) {
        exit(1);
    }
    char *buffer = malloc(1000005); // выделяем буфер для чтения
    fgets(buffer, 1000000, f);      // читаем строку из файла
    fclose(f);
    buffer[strcspn(buffer, "\n")] = 0; // убираем символ новой строки
    BigNum res = MakeBigNum(buffer);   // преобразуем строку в BigNum
    free(buffer);                      // освобождаем буфер
    return res;
}
    




// ------------------- АРИФМЕТИЧЕСКИЕ ОПЕРАЦИИ -------------------

// Сложение по модулю
BigNum AddAbs(BigNum a, BigNum b){
    int maxlen; // макс. длина результата сложения
    
    // Если а длиннее
    if(a.length > b.length){
        maxlen = a.length + 1; // берём длину а + 1
    } 
    else{
        maxlen = b.length + 1; // берём длину b + 1
    }

    BigNum res = CreateBigNum(maxlen);; // BigNum для результата
    res.sign = 1; // знак +

    int c = 0; // перенос между разрядами

    // Сложение в столбик
    for (int i = 0; i< maxlen; i++){
        int sum = c; // начинаем с переноса

        // Если есть цифра в а
        if (i < a.length){
            sum += a.digits[i];
        }
        // Если есть цифра в b
        if (i < b.length){
            sum += b.digits[i];
        }

        res.digits[i] = sum % 10; // записываем остаток
        c = sum / 10;             // вычисляем перенос
    }
    TrimBigNum(&res); // убираем ведущие нули

    return res;
}

// Вычитание по модулю (a >= b)
BigNum SubAbs(BigNum a, BigNum b){
    BigNum res = CreateBigNum(a.length); // результат такой же длины как а

    int f = 0; // флаг занимания

    int b_digit; // текущая цифра из b
    for (int i = 0; i < a.length; i++){
   
        // Если есть цифра в b
        if (i < b.length){
            b_digit = b.digits[i];
        }
        // Если цифры в b закончились
        else{
            b_digit = 0;
        }

        int diff = a.digits[i] - f - b_digit; // разность для текущего разряда

        // Если разность отрицательная
        if (diff < 0){ 
            diff += 10; // занимаем 10 у следующего разряда
            f = 1;      // устанавливаем флаг заёма
        }
        else{
            f = 0;      // сбрасываем флаг заёма
        }
        res.digits[i] = diff; // записываем в результат
    }
    
    TrimBigNum(&res); // убираем ведущие нули
    return res;
}

// Сравнение по модулю
int CompareAbs(BigNum a, BigNum b){
    //Сначала сравниваем длину
    if (a.length != b.length){
        if (a.length > b.length){
            return 1; // a > b
        }
        else{
            return -1;// a < b
        }
    }

    //Если длины равны, сравниваем цифры с конца (со старших разрядов)
    for (int i = a.length - 1; i >= 0; i--){
        if (a.digits[i] != b.digits[i]){
            if (a.digits[i] > b.digits[i]){
                return 1; // a > b
            }
            else{
                return -1;// a < b
            }
        }
    }

    //Если все цифры одинаковые
    return 0;
}

// Сложение с учётом знака
BigNum Add(BigNum a, BigNum b){

    //Если знаки одинаковые
    if (a.sign == b.sign){
        BigNum res = AddAbs(a, b);  // складываем по модулю
        res.sign = a.sign;
        return res;
    }

    // Если знаки разные
    int c = CompareAbs(a, b); // сравниваем модули
        
    if (c == 0){    // модули равны
        BigNum res = CreateBigNum(1);
        res.digits[0] = 0; // сумма равна 0
        return res;
    }

    else if (c > 0){ // |a| > |b|
        BigNum res = SubAbs(a, b); // a - b
        res.sign = a.sign;         // знак как у а
        return res;
    }

        else{        // |a| < |b|
            BigNum res = SubAbs(b, a); // b - a
            res.sign = b.sign;         // знак как у b
            return res;
        }
}

// Вычитание с учётом знака
BigNum Sub(BigNum a, BigNum b){
    b.sign *= -1; // меняем знак b
    return Add(a, b); // складывем a + (-b)
}

// Простое умножение 
BigNum SimpleMul(BigNum a, BigNum b) {
    // Если одно из чисел ноль
    if ((a.length == 1 && a.digits[0] == 0) || (b.length == 1 && b.digits[0] == 0)) {
        BigNum zero = CreateBigNum(1);
        zero.digits[0] = 0;
        return zero;
    }
    
    BigNum res = CreateBigNum(a.length + b.length);
    res.sign = a.sign * b.sign;
    
    
    int b_digit; // текущая цифра числа b для умножения
    // Умножение в столбик
    for (int i = 0; i < a.length; i++) {
        int carry = 0; // перенос
        for (int j = 0; j < b.length || carry; j++) { // проверка границ
            if (i + j >= res.length) continue;
            
        // Проверяем, не вышли ли за пределы числа b
        if (j < b.length) {
            b_digit = b.digits[j];  // берём j-ю цифру числа b, если индекс в пределах
        } else {
            b_digit = 0;            // иначе используем 0 
        }

        // Вычисляем произведение с учетом всех компонентов:
        long long product = res.digits[i + j] +  // текущее значение в разряде результата
                        (long long)a.digits[i] * b_digit +  // произведение i-й цифры числа a на текущую цифру числа b
                        carry; // перенос из предыдущего разряда

            res.digits[i + j] = product % 10; // записываем остаток
            carry = product / 10;             // новый перенос
        }
    }
    TrimBigNum(&res); // убираем ведущие нули
    return res;
}

// Деление BigNum на цифру (0..9)
BigNum DivByDigit(BigNum a, int d) {
    // Проверка деления на ноль
    if (d == 0) {
        exit(1);
    }
    
    // Деление на 1
    if (d == 1) {
        return CopyBigNum(a);
    }
    
    // Создаем результат такой же длины
    BigNum result = CreateBigNum(a.length);
    result.sign = a.sign;
    
    int carry = 0;
    
    // Деление в столбик (от старших разрядов к младшим)
    for (int i = a.length - 1; i >= 0; i--) {
        int current = carry * 10 + a.digits[i];
        result.digits[i] = current / d;
        carry = current % d;
    }
    
    // Убираем ведущие нули
    TrimBigNum(&result);
    
    return result;
}
// Основная функция умножения
BigNum Mul(BigNum a, BigNum b) {
    // Если одно из чисел ноль
    if ((a.length == 1 && a.digits[0] == 0) || (b.length == 1 && b.digits[0] == 0)) {
        BigNum zero = CreateBigNum(1);
        zero.digits[0] = 0;
        return zero;
    }

    // // Если числа очень маленькие, используем простое умножение
    // if (a.length <= 10 && b.length <= 10) {
    //     return SimpleMul(a, b);
    // }

    // return FFTMul(a, b); // используем быстрое умножение

    return SimpleMul(a, b);
}



// Сравнение с учётом знака
int Compare(BigNum a, BigNum b) {
    // Если знаки разные
    if (a.sign != b.sign) {
        if (a.sign > b.sign){
            return 1;   // a > b
        }
        else{
            return -1;  // a < b
        }
    }

    // Если знаки одинаковые, сравниваем модули чисел
    int abs_cmp = CompareAbs(a, b);
    // Оба полож. - чем больше модуль, тем больше число
    if (a.sign == 1){
        return abs_cmp;
    }// Оба отриц. - чем больше модуль, тем меньше число
    else{
        return -abs_cmp;
    }
}

// Деление в столбик с оценкой цифры
BigNum Div(BigNum a, BigNum b) {
    // Проверка деления на ноль
    if (b.length == 1 && b.digits[0] == 0) {
        exit(1);
    }

    // Если |a| < |b|, результат = 0
    if (CompareAbs(a, b) < 0) {
        BigNum zero = CreateBigNum(1);
        zero.digits[0] = 0;
        return zero;
    }

    int n = a.length; // длина делимого
    int m = b.length; // длина делителя

    // Результат максимум n - m + 1 цифра
    BigNum q = CreateBigNum(n - m + 1); // частное

    // Копируем делимое (чтобы вычитать)
    BigNum r = CopyBigNum(a);

    // Нормализуем делитель (старшая цифра > 0)
    int d = 10 / (b.digits[m - 1] + 1);

    if (d > 1) {
        a = MulByDigit(a, d); // умножаем а и b на d
        b = MulByDigit(b, d);
        r = MulByDigit(r, d);
    }

    // Обработка разрядов частного
    for (int k = n - m; k >= 0; k--) {
        // Оценка цифры частного
        long long r2 = (long long)r.digits[m + k] * 10 + r.digits[m + k - 1];
        int qt = (int)(r2 / b.digits[m - 1]); // предполагаемая цифра

        if (qt > 9) qt = 9; // ограничиваем 9

        // Умножаем делитель на qt и вычитаем
        BigNum tmp = MulByDigit(b, qt);
        ShiftLeftT(&tmp, k); // сдвигаем на нужную позицию

        while (CompareAbs(r, tmp) < 0) { // если переоценили
            qt--;                   // уменьшаем цифру
            FreeBigNum(&tmp);
            tmp = MulByDigit(b, qt);// пересчитываем
            ShiftLeftT(&tmp, k);
        }

        r = Sub(r, tmp); // вычитаем
        FreeBigNum(&tmp);

        q.digits[k] = qt;// записываем цифру частного
    }

    // Убираем ведущие нули
    while (q.length > 1 && q.digits[q.length - 1] == 0)
        q.length--;

    // Если знаки делимого и делителя одинаковые, частное полож.
    if (a.sign == b.sign){
        q.sign = 1;
    }
    // Если знаки разные, частное отриц.
    else{
        q.sign = -1;
    }

    return q;
}

// Сдвиг вправо на k позиций (деление на 10^k)
void ShiftRightT(BigNum *a, int k) {
    
    // Проверка на нулевой или отрицательный сдвиг
    if (k <= 0) return;
    
    // Если сдвиг больше или равен длине числа, результат 0
    if (k >= a->length) {
        free(a->digits);
        a->digits = (int*)calloc(1, sizeof(int));
        a->digits[0] = 0;
        a->length = 1;
        a->sign = 1;
        return;
    }
    
    // Новая длина после сдвига
    int newLen = a->length - k;
    
    // Сдвигаем цифры влево на k позиций
    for (int i = 0; i < newLen; i++) {
        a->digits[i] = a->digits[i + k];
    }
    
    // Обновляем длину
    a->length = newLen;
    
    // Перевыделяем память для уменьшенного массива
    a->digits = realloc(a->digits, newLen * sizeof(int));
    
    // Убираем ведущие нули
    TrimBigNum(a);
    
}

// Сдвиг влево (умножение на 10^k)
void ShiftLeftT(BigNum *a, int k) {
    // Проверка на нулевой или отрицательный сдвиг
    if (k <= 0) return;
    
    // Вычисляем новую длину числа после сдвига
    int newLen = a->length + k;
    
    // Перевыделяем память для увеличенного массива цифр
    a->digits = realloc(a->digits, newLen * sizeof(int));
    
    // Проверяем успешность выделения памяти
    if (!a->digits) {
        exit(1);
    }
    
    // Сдвигаем существующие цифры вправо на k позиций
    // Начинаем с конца, чтобы не перезаписать нужные цифры
    for (int i = newLen - 1; i >= k; i--)
        // Копируем цифру из позиции [i-k] в позицию [i]
        a->digits[i] = a->digits[i - k];
    
    // Заполняем освободившиеся младшие k разрядов нулями
    for (int i = 0; i < k; i++)
        // Устанавливаем i-ю цифру в 0
        a->digits[i] = 0;
    
    // Обновляем длину числа
    a->length = newLen;
}

// Умножение BigNum на цифру (0..9)
BigNum MulByDigit(BigNum a, int d) {
    // Умножение на 0
    if (d == 0) {
        BigNum zero = CreateBigNum(1);
        zero.digits[0] = 0;
        return zero;
    }

    BigNum res = CreateBigNum(a.length + 1); // результат
    int carry = 0; // перенос
    for (int i = 0; i < a.length; i++) {
        int prod = a.digits[i] * d + carry; // произведение + перенос
        res.digits[i] = prod % 10;          // записываем остаток
        carry = prod / 10;                  // вычисляем новый перенос
    }
    // Если остался перенос
    if (carry) res.digits[a.length] = carry;
    // Иначе уменьшаем длину
    else res.length--;
    return res;
}


// Создание BigNum из int 
BigNum CreateFromInt(int n) {
    // Если 0
    if (n == 0) {
        BigNum zero = CreateBigNum(1);  // создаём BigNum длиной 1
        zero.digits[0] = 0;             // устанавливаем единственную цифру в 0
        return zero;                    
    }
    
    int temp;
    // Если число положительное, берём его как есть
    if (n > 0) {
        temp = n;       
    } 
    // Если число отрицательное, берём его модуль
    else {
        temp = -n;      
    }
    
    // Вычисляем количество цифр в числе
    int length = 0;
    int t = temp;       // копия для подсчёта цифр
    while (t > 0) {
        length++;       // увеличиваем счётчик цифр
        t /= 10;        // убираем последнюю цифру
    }
    
    // Создаем BigNum с вычисленной длиной
    BigNum result = CreateBigNum(length);
    
    // Устанавливаем знак числа
    if (n > 0) {
        result.sign = 1;    // положительное число
    } else {
        result.sign = -1;   // отрицательное число
    }
    
    // Заполняем массив цифр (в обратном порядке)
    for (int i = 0; i < length; i++) {
        result.digits[i] = temp % 10;  // берём последнюю цифру
        temp /= 10;                    // убираем последнюю цифру из временной переменной
    }
    
    return result;
}

// Возведение в степень 
BigNum Pow(BigNum a, int n) {
    if (n == 0) {
        return CreateFromInt(1);
    }
    if (n == 1) {
        return CopyBigNum(a);
    }

    BigNum result = CreateFromInt(1);
    BigNum base = CopyBigNum(a);
    int exponent = n;
    
    while (exponent > 0) {
        if (exponent & 1) {
            BigNum tmp = Mul(result, base);
            FreeBigNum(&result);
            result = tmp;
        }
        
        exponent >>= 1;
        if (exponent > 0) {
            BigNum tmp = Mul(base, base);
            FreeBigNum(&base);
            base = tmp;
        }
    }
    
    FreeBigNum(&base);
    return result;
}



//--------НОД и НОК--------------------

// BigNum NOD(BigNum a, BigNum b) {
//     a.sign = 1; b.sign = 1;
    
//     if (IsZero(a)) return CopyBigNum(b);
//     if (IsZero(b)) return CopyBigNum(a);
    
//     BigNum temp;
//     int shift = 0;
    
//     // Убираем общие степени двойки
//     while (IsEven(a) && IsEven(b)) {
//         ShiftRight(&a);
//         ShiftRight(&b);
//         shift++;
//     }
    
//     do {
//         // Убираем степени двойки из каждого числа
//         while (IsEven(a)) ShiftRight(&a);
//         while (IsEven(b)) ShiftRight(&b);
        
//         // Сравниваем и вычитаем
//         int cmp = CompareAbs(a, b);
//         if (cmp >= 0) {
//             temp = SubAbs(a, b);
//             FreeBigNum(&a);
//             a = temp;
//         } else {
//             temp = SubAbs(b, a);
//             FreeBigNum(&b);
//             b = temp;
//         }
//     } while (!IsZero(a) && !IsZero(b));
    
//     // Возвращаем ненулевое число
//     BigNum result = IsZero(a) ? CopyBigNum(b) : CopyBigNum(a);
    
//     // Восстанавливаем степени двойки
//     for (int i = 0; i < shift; i++) {
//         ShiftLeft(&result);
//     }
    
//     if (!IsZero(a)) FreeBigNum(&a);
//     if (!IsZero(b)) FreeBigNum(&b);
    
//     return result;
// }


//  НОД
BigNum NOD(BigNum a, BigNum b) {
    
    // Тривиальные случаи
    if (IsZero(a)) return CopyBigNum(b);
    if (IsZero(b)) return CopyBigNum(a);
    
    // Для маленьких a используем быстрый алгоритм
    if (a.length <= 2) {
        int a_val = 0;
        for (int i = a.length - 1; i >= 0; i--) a_val = a_val * 10 + a.digits[i];
        
        // Вычисляем НОД через остатки
        int remainder = 0;
        for (int i = b.length - 1; i >= 0; i--) {
            remainder = (remainder * 10 + b.digits[i]) % a_val;
        }
        
        // Алгоритм Евклида для маленьких чисел
        int x = a_val, y = remainder;
        while (y != 0) {
            int temp = y;
            y = x % y;
            x = temp;
        }
        
        return CreateFromInt(x);
    }
    
    // Для больших чисел используем бинарный алгоритм
    BigNum x = CopyBigNum(a);
    BigNum y = CopyBigNum(b);
    x.sign = 1;
    y.sign = 1;
    
    int shift = 0;
    int iterations = 0;
    const int MAX_ITERATIONS = 20;
    
    // Убираем общие степени двойки
    while (IsEven(x) && IsEven(y)) {
        ShiftRight(&x);
        ShiftRight(&y);
        shift++;
    }
    
    while (!IsZero(y) && iterations < MAX_ITERATIONS) {
        
        while (IsEven(x)) ShiftRight(&x);
        while (IsEven(y)) ShiftRight(&y);
        
        if (Compare(x, y) >= 0) {
            BigNum temp = Sub(x, y);
            FreeBigNum(&x);
            x = temp;
            ShiftRight(&x);
        } else {
            BigNum temp = Sub(y, x);
            FreeBigNum(&y);
            y = temp;
            ShiftRight(&y);
        }
        
        iterations++;
    }
    
    BigNum result = IsZero(x) ? CopyBigNum(y) : CopyBigNum(x);
    
    // Восстанавливаем степени двойки
    for (int i = 0; i < shift; i++) {
        ShiftLeft(&result);
    }
    
    // for (int i = result.length - 1; i >= 0; i--) {
    //     printf("%d", result.digits[i]);
    // }
    // printf("\n");
    
    FreeBigNum(&x);
    FreeBigNum(&y);
    
    return result;
}
// НОК 

BigNum NOK(BigNum a, BigNum b) {
    
    // Метод: НОК(a,b) = a × (b / НОД(a,b))
    // Это требует только одного деления меньшего масштаба
    
    BigNum g = NOD(a, b);
    
    // Делим b на НОД 
    BigNum b_div_g = Div(b, g);
    
    // Умножаем a на результат
    BigNum nok = Mul(a, b_div_g);
    
    FreeBigNum(&g);
    FreeBigNum(&b_div_g);
    return nok;
}

// Проверка на ноль
int IsZero(BigNum a) {
    TrimBigNum(&a); // убираем ведущие нули
    return (a.length == 1 && a.digits[0] == 0); // проверяем на 0
}

// Проверка на чётность

int IsEven(BigNum n) {
    int result = (n.digits[0] % 2 == 0);
    return result;
}

// Оптимизированный сдвиг вправо (деление на 2)
void ShiftRight(BigNum *a) {

    int carry = 0;
    for (int i = a->length - 1; i >= 0; i--) {
        int cur = carry * 10 + a->digits[i];  // текущее значение с переносом
        a->digits[i] = cur / 2;         // делим на 2
        carry = cur % 2;                // остаток
    }
    TrimBigNum(a);                      // убираем ведущие нули
}

// Оптимизированный сдвиг влево (умножение на 2)
void ShiftLeft(BigNum *a) {
    int carry = 0;
    int new_length = a->length + 1; // возможное увеличение длины
    a->digits = realloc(a->digits, new_length * sizeof(int)); // перевыделяем память
    
    for (int i = 0; i < a->length; i++) {
        int cur = a->digits[i] * 2 + carry; // умножаем на 2 + перенос
        a->digits[i] = cur % 10; // записываем остаток
        carry = cur / 10; // вычисляем перенос
    }
    
    if (carry > 0) { // если остался перенос
        a->digits[a->length] = carry;
        a->length = new_length;
    }
    
    TrimBigNum(a);
}



// -------- БЫСТРОЕ ПРЕОБРАЗОВАНИЕ ФУРЬЕ (FFT) ----------

// Рекурсивная реализация алгоритма FFT 
static void fft(long double complex *a, int n, int invert) {
    // Этап 1: Бит-реверсивная перестановка
    for (int i = 1, j = 0; i < n; i++) {
        int bit = n >> 1;  // cтарший бит для реверса
        // Находим следующий индекс в бит-реверсивном порядке
        for (; j & bit; bit >>= 1) j ^= bit;
        j ^= bit;
        // Меняем местами элементы, если нужно
        if (i < j) {
            long double complex tmp = a[i]; 
            a[i] = a[j]; 
            a[j] = tmp;
        }
    }

    // Этап 2: Бабочковые операции
    for (int len = 2; len <= n; len <<= 1) {
        // Вычисляем угол для поворачивающего множителя (корня из единицы)
        long double ang;

        // Определяем направление вращения в комплексной плоскости:
        if (invert) {
            ang = 2 * M_PI / len * (-1);  // Обратное преобразование - отрицательный угол
        } else {
            ang = 2 * M_PI / len * 1;     // Прямое преобразование - положительный угол
        }
        long double complex wlen = cosl(ang) + I * sinl(ang);  // комплексная экспонента
                // Обрабатываем все блоки длины len
        for (int i = 0; i < n; i += len) {
            long double complex w = 1.0L + 0.0L*I;  // текущий поворачивающий множитель
            // Бабочковая операция для половинок блока
            for (int j = 0; j < len/2; j++) {
                // Берем элементы из двух половин
                long double complex u = a[i+j];              // элемент из первой половины
                long double complex v = a[i+j+len/2] * w;    // элемент из второй половины с поворотом
                
                // Выполняем бабочковую операцию
                a[i+j] = u + v;          // сумма
                a[i+j+len/2] = u - v;    // разность
                
                w *= wlen;  // умножаем на корень для следующей итерации
            }
        }
    }

    // Для обратного преобразования нормируем результат
    if (invert) {
        for (int i = 0; i < n; i++) a[i] /= n;
    }
}


// --------- СВЕРТКА С ПОМОЩЬЮ FFT ----------

// Вычисление свертки двух массивов с использованием FFT
static long long * convolution(const int *a, int na, const int *b, int nb, int *out_len) {
    // Определяем необходимый размер для FFT (степень двойки)
    int need = na + nb;  // максимальная длина свертки
    int n = 1;
    while (n < need) n <<= 1;  // находим ближайшую степень двойки >= need

    // Выделяем память под комплексные представления массивов
    long double complex *fa = calloc(n, sizeof(long double complex));
    long double complex *fb = calloc(n, sizeof(long double complex));

    // Копируем исходные данные в комплексные массивы
    for (int i = 0; i < na; i++) fa[i] = a[i];
    for (int i = 0; i < nb; i++) fb[i] = b[i];

    // Выполняем прямое FFT для обоих массивов
    fft(fa, n, 0);  // прямое преобразование для a
    fft(fb, n, 0);  // прямое преобразование для b
    
    // Поэлементно перемножаем результаты в частотной области
    for (int i = 0; i < n; i++) fa[i] *= fb[i];
    // Выполняем обратное FFT для получения свертки
    fft(fa, n, 1);  // обратное преобразование

    // Преобразуем комплексный результат в целочисленный
    long long *res = calloc(n, sizeof(long long));
    for (int i = 0; i < n; i++) res[i] = llroundl(creall(fa[i]));  // округляем вещественную часть

    // Освобождаем временные массивы
    free(fa); free(fb);
    
    // Возвращаем длину результата через указатель
    *out_len = n;
    return res;
}

// -------- УПАКОВКА/РАСПАКОВКА ДЛЯ FFT --------

// Упаковка BigNum в блоки для эффективного умножения через FFT
static int * pack_blocks(const BigNum *a, int *blocks_len) {
    // Вычисляем количество блоков (каждый блок содержит BASE_DIGS цифр)
    int blocks = (a->length + BASE_DIGS - 1) / BASE_DIGS;
    
    // Выделяем память под массив блоков
    int *arr = calloc(blocks, sizeof(int));
    
    // Заполняем блоки цифрами из BigNum
    for (int i = 0; i < blocks; i++) {
        int val = 0, pow10 = 1;  // инициализация значения блока и степени 10
        
        // Собираем BASE_DIGS цифр в один блок
        for (int j = 0; j < BASE_DIGS; j++) {
            int idx = i*BASE_DIGS + j;  // абсолютный индекс цифры в BigNum
            
            // Если цифра существует, добавляем ее к значению блока
            if (idx < a->length) val += a->digits[idx] * pow10;
            pow10 *= 10;  // увеличиваем степень для следующей цифры
        }
        arr[i] = val;  // сохраняем значение блока
    }
    
    // Возвращаем длину через указатель
    *blocks_len = blocks;
    return arr;
}
// Распаковка блоков обратно в BigNum после умножения
static BigNum unpack_blocks_to_BigNum(long long *blocks, int blocks_len) {
    // Убираем ведущие нулевые блоки
    while (blocks_len > 1 && blocks[blocks_len-1] == 0) blocks_len--;

    // Обрабатываем переносы между блоками
    long long carry = 0;
    for (int i = 0; i < blocks_len; i++) {
        long long cur = blocks[i] + carry;  // текущее значение + перенос
        
        if (cur >= 0) {
            // Положительное значение - нормализуем по основанию FFT_BASE
            blocks[i] = cur % FFT_BASE;  // оставляем остаток
            carry = cur / FFT_BASE;      // переносим целую часть
        } else {
            // Отрицательное значение - корректируем заем
            long long k = (llabs(cur) + FFT_BASE - 1) / FFT_BASE;  // вычисляем корректирующий множитель
            blocks[i] = cur + k * FFT_BASE;  // добавляем основание
            carry = -k;             // учитываем заем
        }
    }
    // Обрабатываем оставшийся перенос (расширяем массив при необходимости)
    while (carry > 0) {
        blocks = realloc(blocks, (blocks_len+1)*sizeof(long long));  // расширяем массив
        blocks[blocks_len++] = carry % FFT_BASE;  // добавляем новый блок
        carry /= FFT_BASE;             // продолжаем обработку переноса
    }

    // Вычисляем общее количество цифр в результате
    int total_digits = blocks_len * BASE_DIGS;  // максимально возможное количество
    long long last = blocks[blocks_len-1];      // последний (старший) блок
    int extra = 0;
    
    // Подсчитываем реальное количество цифр в последнем блоке
    while (last > 0) { extra++; last /= 10; }
    if (extra == 0) extra = 1;  // гарантируем хотя бы одну цифру
    
    // Корректируем общее количество цифр
    total_digits = (blocks_len-1)*BASE_DIGS + extra;

    // Создаем BigNum для результата
    BigNum res = CreateBigNum(total_digits);
    int pos = 0;  // позиция в массиве цифр результата
    // Распаковываем блоки в цифры
    for (int i = 0; i < blocks_len; i++) {
        long long v = blocks[i];  // текущий блок
        
        // Разбираем блок на отдельные цифры
        for (int j = 0; j < BASE_DIGS; j++) {
            if (pos >= total_digits) break;  // проверка выхода за границы
            res.digits[pos++] = v % 10;      // берём младшую цифру
            v /= 10;                         // сдвигаемся к следующей цифре
        }
    }
    
    // Убираем возможные ведущие нули
    TrimBigNum(&res);
    return res;
}

// --------- ОСНОВНАЯ ФУНКЦИЯ УМНОЖЕНИЯ ЧЕРЕЗ FFT ---------

// Умножение больших чисел с использованием быстрого преобразования Фурье
BigNum FFTMul(BigNum a, BigNum b) {
    // Проверка тривиальных случаев - умножение на ноль
    if ((a.length == 1 && a.digits[0] == 0) || (b.length == 1 && b.digits[0] == 0)) {
        BigNum zero = CreateBigNum(1);
        zero.digits[0] = 0;
        return zero;
    }

    // Упаковываем числа в блоки для FFT
    int na, nb;  // Длины упакованных массивов
    int *A = pack_blocks(&a, &na);  // упаковка первого числа
    int *B = pack_blocks(&b, &nb);  // упаковка второго числа

    // Вычисляем свертку через FFT
    int conv_len;
    long long *conv = convolution(A, na, B, nb, &conv_len);

    // Освобождаем временные упакованные массивы
    free(A); free(B);
    // Распаковываем результат свертки обратно в BigNum
    BigNum result = unpack_blocks_to_BigNum(conv, conv_len);

    // Освобождаем массив свертки
    free(conv);
    return result;
}