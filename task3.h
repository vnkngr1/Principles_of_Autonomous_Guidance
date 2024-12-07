#ifndef TASK3_H
#define TASK3_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>

#include "const_and_algs.h"

// Функция для расчета времени перехода между высотами
double calculate_time(double H1, double H2, double V, double Cy) {
    double rho = air_density(H1);
    double q = 0.5 * rho * V * V;
    double theta = asin(Cy / (q * S / (m * g)));
    if (sin(theta) <= 0) return std::numeric_limits<double>::infinity();
    return (H2 - H1) / (V * sin(theta));
}

int alg_task3() {
    // Создаем сетки скоростей и высот
    int num_heights = static_cast<int>((H_end - H_start) / H_step) + 1;
    int num_speeds = static_cast<int>((V_end - V_start) / V_step) + 1;

    std::vector<std::vector<double>> time_grid(num_heights, std::vector<double>(num_speeds, std::numeric_limits<double>::infinity()));
    time_grid[0] = std::vector<double>(num_speeds, 0.0); // Время на начальной высоте = 0

    // Динамическое программирование
    for (int j = 1; j < num_heights; ++j) {
        double H_current = H_start + (j - 1) * H_step;
        double H_next = H_start + j * H_step;

        for (int i = 0; i < num_speeds; ++i) {
            double V = V_start + i * V_step;
            double rho = air_density(H_current);
            double q = 0.5 * rho * V * V;
            double Cy = m * g / (q * S);

            for (int k = 0; k < num_speeds; ++k) {
                double V_next = V_start + k * V_step;
                double delta_t = calculate_time(H_current, H_next, V, Cy);

                if (time_grid[j - 1][i] + delta_t < time_grid[j][k]) {
                    time_grid[j][k] = time_grid[j - 1][i] + delta_t;
                }
            }
        }
    }

    // Поиск минимального времени для конечной высоты
    double min_time = *std::min_element(time_grid[num_heights - 1].begin(), time_grid[num_heights - 1].end());
    std::cout << "Минимальное время полета: " << min_time / 60.0 << " минут." << std::endl;

    return 0;
}

#endif //TASK3_H
