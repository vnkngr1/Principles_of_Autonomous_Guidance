#ifndef TASK1_H
#define TASK1_H

#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

// Константы
const double g = 9.81;                 // Ускорение свободного падения, м/с²
const double m = 98000;                // Масса ЛА, кг
const double S = 201;                  // Площадь крыла, м²
const double rho0 = 1.225;             // Плотность воздуха на уровне моря, кг/м³
const double H_start = 300;            // Начальная высота, м
const double V_start = 83.33;          // Начальная скорость, м/с
const double V_end = 222.22;           // Конечная скорость, м/с
const double V_step = 10;              // Шаг скорости, м/с
const double P_total = 103000 * 3;     // Суммарная тяга двигателей, Н

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

int alg_task1() {
    // Векторы для хранения данных
    std::vector<double> heights;
    std::vector<double> velocities;

    // Инициализация начальных значений
    double H_current = H_start;
    heights.push_back(H_current);

    for (double V = V_start; V <= V_end; V += V_step) {
        // Плотность воздуха
        double rho = air_density(H_current);

        // Динамическое давление
        double q = 0.5 * rho * V * V;

        // Коэффициент подъемной силы
        double Cy = m * g / (q * S);

        // Сопротивление
        double Cx = calculate_Cx(Cy);
        double D = Cx * q * S;

        // Угол наклона траектории
        double theta = asin((P_total - D) / (m * g));

        // Изменение высоты
        double delta_H = V * sin(theta) * V_step / (V + V_step);

        // Обновление высоты
        H_current += delta_H;
        heights.push_back(H_current);
        velocities.push_back(V);
    }

    // Вывод результатов в файл
    std::ofstream outFile("trajectory_data.csv");
    if (outFile.is_open()) {
        outFile << "Velocity (m/s),Height (m)\n";
        for (size_t i = 0; i < velocities.size(); ++i) {
            outFile << velocities[i] << "," << heights[i] << "\n";
        }
        outFile.close();
        std::cout << "Данные траектории сохранены в файл trajectory_data.csv" << std::endl;
    } else {
        std::cerr << "Не удалось открыть файл для записи!" << std::endl;
    }

    return 0;
}

#endif //TASK1_H
