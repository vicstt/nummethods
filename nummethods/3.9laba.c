#include <stdio.h>
#include <math.h>

// Структура для хранения узлов интерполяции
typedef struct {
    double x;
    double y;
} Point;

// Функция вычисления интерполяционного многочлена Лагранжа
double lagrange_interpolation(Point points[], int n, double x) {
    double result = 0.0;
    
    for (int i = 0; i < n; i++) {
        double term = points[i].y;
        for (int j = 0; j < n; j++) {
            if (j != i) {
                term *= (x - points[j].x) / (points[i].x - points[j].x);
            }
        }
        result += term;
    }
    
    return result;
}

// Функция для выбора ближайших узлов интерполяции
void select_nearest_points(Point all_points[], int total_count, Point selected_points[], 
                          int selected_count, double x) {
    // Создаем временный массив для сортировки по расстоянию
    Point temp[total_count];
    for (int i = 0; i < total_count; i++) {
        temp[i] = all_points[i];
    }
    
    // Сортируем точки по расстоянию до x (пузырьковая сортировка)
    for (int i = 0; i < total_count - 1; i++) {
        for (int j = 0; j < total_count - i - 1; j++) {
            if (fabs(temp[j].x - x) > fabs(temp[j + 1].x - x)) {
                Point tmp = temp[j];
                temp[j] = temp[j + 1];
                temp[j + 1] = tmp;
            }
        }
    }
    
    // Выбираем нужное количество ближайших точек
    for (int i = 0; i < selected_count; i++) {
        selected_points[i] = temp[i];
    }
}

int main() {
    // Исходные данные из таблицы
    Point all_points[] = {
        {-9.73, -0.7654},
        {-7.87, -0.1332},
        {-6.01, 1.5269},
        {-4.15, 0.8327},
        {-2.29, -0.9364},
        {-0.43, -1.0428},
        {1.43, 0.4213},
        {3.29, 0.9816},
        {5.15, 0.2153}
    };
    
    int total_points = sizeof(all_points) / sizeof(all_points[0]);
    double x_target = -3.548;
    
    printf("Интерполяция функции в точке x = %.3f\n\n", x_target);
    printf("Исходные данные:\n");
    printf("i\tx_i\ty_i\n");
    for (int i = 0; i < total_points; i++) {
        printf("%d\t%.2f\t%.4f\n", i, all_points[i].x, all_points[i].y);
    }
    printf("\n");
    
    // Интерполяция многочленом Лагранжа 2-й степени
    Point points_degree2[3];
    select_nearest_points(all_points, total_points, points_degree2, 3, x_target);
    
    printf("Узлы для многочлена 2-й степени:\n");
    printf("i\tx_i\ty_i\n");
    for (int i = 0; i < 3; i++) {
        printf("%d\t%.2f\t%.4f\n", i, points_degree2[i].x, points_degree2[i].y);
    }
    
    double result_degree2 = lagrange_interpolation(points_degree2, 3, x_target);
    printf("L2(%.3f) = %.6f\n\n", x_target, result_degree2);
    
    // Проверка в узловой точке
    double check_point = points_degree2[1].x;
    double check_value = lagrange_interpolation(points_degree2, 3, check_point);
    printf("Проверка в узловой точке x = %.2f:\n", check_point);
    printf("L2(%.2f) = %.6f\n", check_point, check_value);
    printf("Фактическое значение y = %.4f\n", points_degree2[1].y);
    printf("Погрешность: %.6f\n\n", fabs(check_value - points_degree2[1].y));
    
    // Интерполяция многочленом Лагранжа 3-й степени
    Point points_degree3[4];
    select_nearest_points(all_points, total_points, points_degree3, 4, x_target);
    
    printf("Узлы для многочлена 3-й степени:\n");
    printf("i\tx_i\ty_i\n");
    for (int i = 0; i < 4; i++) {
        printf("%d\t%.2f\t%.4f\n", i, points_degree3[i].x, points_degree3[i].y);
    }
    
    double result_degree3 = lagrange_interpolation(points_degree3, 4, x_target);
    printf("L3(%.3f) = %.6f\n\n", x_target, result_degree3);
    
    // Проверка в узловой точке
    check_point = points_degree3[1].x;
    check_value = lagrange_interpolation(points_degree3, 4, check_point);
    printf("Проверка в узловой точке x = %.2f:\n", check_point);
    printf("L3(%.2f) = %.6f\n", check_point, check_value);
    printf("Фактическое значение y = %.4f\n", points_degree3[1].y);
    printf("Погрешность: %.6f\n\n", fabs(check_value - points_degree3[1].y));
    
    // Сравнение результатов
    printf("Сравнение результатов интерполяции:\n");
    printf("L2(%.3f) = %.6f\n", x_target, result_degree2);
    printf("L3(%.3f) = %.6f\n", x_target, result_degree3);
    printf("Разность: %.6f\n", fabs(result_degree3 - result_degree2));
    
    return 0;
}