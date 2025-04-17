#include "PAdic.h"
#include <vector>

PAdic::PAdic(int p, int prec): Base(p), Precision(prec) {
    if (p <= 1) throw std::invalid_argument("Prime must be > 1");
}
std::vector<int> PAdic::double_to_padic(double x){
    std::vector<int> digits;
    if (x < 0 || x >= 1) throw std::invalid_argument("x must be in [0,1)");
    float temp = x;
    for (int i = 0; i < Precision; ++i) {
        temp *= Base;
        int digit = static_cast<int>(temp);
        digits.push_back(digit);
        temp -= digit;
    }
    return digits;
}

float PAdic::padic_to_double(std::vector<int>* digits) {
    double result = 0.0;
    double weight = 1.0 / Base;
    for (int digit : *digits) {
        result += digit * weight;
        weight /= Base;
    }
    return result;
}

float PAdic::add(double* firstNb, double* secondNb){
    std::vector<int> pDecompFirst = double_to_padic(*firstNb);
    std::vector<int> pDecompSecond = double_to_padic(*secondNb);

    std::vector<int> result;
    int carry = 0;
    for (int i = 0; i < std::max(pDecompFirst.size(), pDecompSecond.size()); ++i) {
        int sum = carry;
        if (i < pDecompFirst.size()) sum += pDecompFirst[i];
        if (i < pDecompSecond.size()) sum += pDecompSecond[i];

        result.push_back(sum % Base);
        carry = sum / Base;
    }
    return padic_to_double(&result);
}

/*std::vector<int> PAdic::rotate(float* number, float* angle) {
    std::vector<int> a_digits = float_to_padic(*number);
    std::vector<int> b_digits = float_to_padic(*angle);
    std::vector<int> result(Precision, 0);

    // Convolution multiplication
    for (int i = 0; i < Precision; ++i) {
        for (int j = 0; j < Precision; ++j) {
            if (i + j >= Precision) continue;
            result[i + j] += a_digits[i] * b_digits[j];
        }
    }

    // Carry propagation
    int carry = 0;
    for (int k = 0; k < Precision; ++k) {
        result[k] += carry;
        carry = result[k] / Base;
        result[k] %= Base;
    }
    return result;
}*/
