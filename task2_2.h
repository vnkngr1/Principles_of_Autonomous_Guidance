#ifndef TASK2_2_H
#define TASK2_2_H

#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <fstream>

#include "const_and_algs.h"

int alg_task2_2() {
    // Сетка высот и скоростей
    int num_heights = static_cast<int>((H_end - H_start) / H_step) + 1;
    int num_speeds = static_cast<int>((V_end - V_start) / V_step) + 1;

    // Массив для хранения минимального времени полета
    std::vector<std::vector<double>> time_grid(num_heights, std::vector<double>(num_speeds, std::numeric_limits<double>::infinity()));
    time_grid[0] = std::vector<double>(num_speeds, 0.0); // На начальной высоте время = 0

    // Динамическое программирование
    for (int j = 1; j < num_heights; ++j) {
        double H_current = H_start + (j - 1) * H_step;
        double H_next = H_start + j * H_step;

        for (int i = 0; i < num_speeds; ++i) {
            double V = V_start + i * V_step;
            double rho = air_density(H_current);
            double q = 0.5 * rho * V * V;
            double Cy = m * g / (q * S);
            double Cx = calculate_Cx(Cy);
            double D = Cx * q * S;

            for (int k = 0; k < num_speeds; ++k) {
                double V_next = V_start + k * V_step;
                double theta = asin((P_total - D) / (m * g));
                if (sin(theta) <= 0) continue; // Невозможный угол

                // Расчет времени перехода
                double delta_t = (H_next - H_current) / (V * sin(theta));

                // Обновление минимального времени
                if (time_grid[j - 1][i] + delta_t < time_grid[j][k]) {
                    time_grid[j][k] = time_grid[j - 1][i] + delta_t;
                }
            }
        }
    }

    // Поиск минимального времени полета на конечной высоте
    double min_time = *std::min_element(time_grid[num_heights - 1].begin(), time_grid[num_heights - 1].end());
    std::cout << "Минимальное время полета: " << min_time / 60.0 << " минут." << std::endl;

    // Вывод данных траектории в файл
    std::ofstream outFile("time_optimized_trajectory.csv");
    if (outFile.is_open()) {
        outFile << "Height (m),Velocity (m/s),Time (s)\n";
        for (int j = 0; j < num_heights; ++j) {
            for (int i = 0; i < num_speeds; ++i) {
                outFile << H_start + j * H_step << "," << V_start + i * V_step << "," << time_grid[j][i] << "\n";
            }
        }
        outFile.close();
        std::cout << "Данные траектории сохранены в файл time_optimized_trajectory.csv" << std::endl;
    } else {
        std::cerr << "Не удалось открыть файл для записи!" << std::endl;
    }

    return 0;
}

#endif //TASK2_2_H
