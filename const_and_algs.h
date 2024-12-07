#ifndef CONST_AND_ALGS_H
#define CONST_AND_ALGS_H

const double g = 9.81;                 // Ускорение свободного падения, м/с²
const double m = 98000;                // Масса ЛА, кг
const double S = 201;                  // Площадь крыла, м²
const double rho0 = 1.225;             // Плотность воздуха на уровне моря, кг/м³
const double H_start = 300;            // Начальная высота, м
const double H_end = 5000;             // Конечная высота, м
const double V_start = 83.33;          // Начальная скорость, м/с
const double V_end = 222.22;           // Конечная скорость, м/с
const double H_step = 100;             // Шаг высоты, м
const double V_step = 10;              // Шаг скорости, м/с
const double P_total = 103000 * 3;     // Суммарная тяга двигателей, Н
const double m_initial = 98000;        // Масса ЛА, кг
const double Cp = 0.86;                // Удельный расход топлива, кг/Н·ч

// Функция для расчета плотности воздуха
double air_density(double h) {
    if (h < 11000) {
        return rho0 * pow(1 - 0.0065 * h / 288.15, 4.256);
    } else {
        return rho0 * 0.22336 * exp(-0.0001577 * (h - 11000));
    }
}

// Аппроксимация зависимости Cx от Cy (условно, из поляры)
double calculate_Cx(double Cy) {
    return 0.02 + 0.04 * Cy * Cy;
}

#endif //CONST_AND_ALGS_H
