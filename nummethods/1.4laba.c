#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>

#define N 4
#define EPSILON 0.0001
#define MAX_ITERATIONS 1000

void print_matrix(double A[N][N], const char* title) {
    printf("%s:\n", title);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%8.2f ", A[i][j]);
        }
        printf("\n");
    }
}

void print_vector(double v[N], const char* title) {
    printf("%s: ", title);
    for (int i = 0; i < N; i++) {
        printf("%8.2f ", v[i]);
    }
    printf("\n");
}

void print_system(double A[N][N], double b[N]) {
    printf("Система уравнений:\n");
    for (int i = 0; i < N; i++) {
        printf("%.0fx₁", A[i][0]);
        for (int j = 1; j < N; j++) {
            if (A[i][j] >= 0) {
                printf(" + %.0fx%d", A[i][j], j + 1);
            } else {
                printf(" - %.0fx%d", -A[i][j], j + 1);
            }
        }
        printf(" = %.0f\n", b[i]);
    }
}

bool check_convergence(double A[N][N]) {
    printf("\nПроверка условия сходимости:\n");

    bool row_dominance = true;
    for (int i = 0; i < N; i++) {
        double diagonal = fabs(A[i][i]);
        double sum_others = 0.0;
        for (int j = 0; j < N; j++) {
            if (j != i) {
                sum_others += fabs(A[i][j]);
            }
        }

        if (diagonal <= sum_others) {
            row_dominance = false;
            printf("Условие диагонального преобладания по строке %d не выполняется: |%.1f| <= %.1f\n", 
                   i + 1, diagonal, sum_others);
        } else {
            printf("Условие диагонального преобладания по строке %d выполняется: |%.1f| > %.1f\n", 
                   i + 1, diagonal, sum_others);
        }
    }

    bool col_dominance = true;
    for (int j = 0; j < N; j++) {
        double diagonal = fabs(A[j][j]);
        double sum_others = 0.0;
        for (int i = 0; i < N; i++) {
            if (i != j) {
                sum_others += fabs(A[i][j]);
            }
        }

        if (diagonal <= sum_others) {
            col_dominance = false;
            printf("Условие диагонального преобладания по столбцу %d не выполняется: |%.1f| <= %.1f\n", 
                   j + 1, diagonal, sum_others);
        } else {
            printf("Условие диагонального преобладания по столбцу %d выполняется: |%.1f| > %.1f\n", 
                   j + 1, diagonal, sum_others);
        }
    }

    bool converges = row_dominance || col_dominance;

    if (converges) {
        printf("\nДостаточное условие сходимости выполняется (по строкам или столбцам)!\n");
        printf("Метод простых итераций сходится!\n");
    } else {
        printf("\nНи диагональное преобладание по строкам, ни по столбцам не выполняется!\n");
        printf("Метод простых итераций может не сходиться!\n");
    }

    return converges;
}

void simple_iteration_method(double A[N][N], double b[N], double solution[N], 
                           int *iterations, double *final_error) {
    double C[N][N], d[N];
    
    // Преобразование системы к виду x = Cx + d
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i == j) {
                C[i][j] = 0.0;
            } else {
                C[i][j] = -A[i][j] / A[i][i];
            }
        }
        d[i] = b[i] / A[i][i];
    }
    
    double x_prev[N], x_new[N];
    for (int i = 0; i < N; i++) {
        x_prev[i] = d[i];
    }
    
    *iterations = 0;
    
    printf("\nИтерационный процесс:\n");
    printf("%-8s %-12s %-12s %-12s %-12s %-12s\n", 
           "Итерация", "x1", "x2", "x3", "x4", "Погрешность");
    printf("----------------------------------------------------------------------------\n");
    
    while (true) {
        (*iterations)++;
        
        // Вычисление нового приближения
        for (int i = 0; i < N; i++) {
            double sum_cx = 0.0;
            for (int j = 0; j < N; j++) {
                sum_cx += C[i][j] * x_prev[j];
            }
            x_new[i] = sum_cx + d[i];
        }
        
        // Вычисление погрешности
        double error = 0.0;
        for (int i = 0; i < N; i++) {
            double diff = fabs(x_new[i] - x_prev[i]);
            if (diff > error) {
                error = diff;
            }
        }
        
        if (*iterations <= 10 || *iterations % 10 == 0 || error < EPSILON) {
            printf("%-8d %-12.6f %-12.6f %-12.6f %-12.6f %-12.6f\n",
                   *iterations, x_new[0], x_new[1], x_new[2], x_new[3], error);
        }
        
        // Проверка условия остановки
        if (error < EPSILON) {
            break;
        }
        
        if (*iterations >= MAX_ITERATIONS) {
            printf("Достигнуто максимальное количество итераций!\n");
            break;
        }
        
        // Копирование нового приближения для следующей итерации
        for (int i = 0; i < N; i++) {
            x_prev[i] = x_new[i];
        }
    }
    
    for (int i = 0; i < N; i++) {
        solution[i] = x_new[i];
    }
    *final_error = *final_error; 
}

double check_solution(double A[N][N], double b[N], double x[N]) {
    printf("\nПроверка решения:\n");
    double max_error = 0.0;
    
    for (int i = 0; i < N; i++) {
        double result = 0.0;
        for (int j = 0; j < N; j++) {
            result += A[i][j] * x[j];
        }
        double error = fabs(result - b[i]);
        if (error > max_error) {
            max_error = error;
        }
        printf("Уравнение %d: %4f (должно быть %2f), погрешность: %4f\n",
               i + 1, result, b[i], error);
    }
    
    return max_error;
}

void rearrange_system(double A[N][N], double b[N], double A_new[N][N], double b_new[N]) {
    printf("\nПроверка сходимости перестановкой уравнений\n");
    
    // Копирование исходной системы
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A_new[i][j] = A[i][j];
        }
        b_new[i] = b[i];
    }
    
    // Частичный выбор главного элемента по строкам
    for (int i = 0; i < N; i++) {
        double max_val = fabs(A_new[i][i]);
        int max_row = i;
        
        for (int j = i + 1; j < N; j++) {
            if (fabs(A_new[j][i]) > max_val) {
                max_val = fabs(A_new[j][i]);
                max_row = j;
            }
        }
        
        if (max_row != i) {
            // Перестановка строк
            for (int k = 0; k < N; k++) {
                double temp = A_new[i][k];
                A_new[i][k] = A_new[max_row][k];
                A_new[max_row][k] = temp;
            }
            double temp_b = b_new[i];
            b_new[i] = b_new[max_row];
            b_new[max_row] = temp_b;
        }
    }
}

int main() {
    double A[N][N] = {
        {6.0, -5.0, -19.0, -7.0},
        {-3.0, -7.0, 4.0, 14.0},
        {-3.0, 12.0, -2.0, 5.0},
        {18.0, -4.0, 9.0, -3.0}
    };
    
    double b[N] = {-39.0, 37.0, 33.0, 43.0};
    
    print_system(A, b);
    
    bool converges = check_convergence(A);
    
    double solution[N];
    int iterations;
    double final_error;
    
    if (converges) {
        simple_iteration_method(A, b, solution, &iterations, &final_error);
        
        printf("\nРезультат:\n");
        printf("Количество итераций: %d\n", iterations);
        printf("Финальная погрешность: %.8f\n", final_error);
        for (int i = 0; i < N; i++) {
            printf("x%d = %.8f\n", i + 1, solution[i]);
        }
        
        check_solution(A, b, solution);
    } else {
        double A_rearranged[N][N], b_rearranged[N];
        rearrange_system(A, b, A_rearranged, b_rearranged);
        
        printf("\nПереставленная система:\n");
        print_system(A_rearranged, b_rearranged);
        
        printf("\nПроверка сходимости переставленной системы:\n");
        bool converges_rearranged = check_convergence(A_rearranged);
        
        if (converges_rearranged) {
            simple_iteration_method(A_rearranged, b_rearranged, solution, &iterations, &final_error);
            
            printf("\nРезультат после %d итераций:\n", iterations);
            printf("Текущая погрешность: %.8f\n", final_error);
            for (int i = 0; i < N; i++) {
                printf("x%d = %.8f\n", i + 1, solution[i]);
            }
            
            check_solution(A, b, solution);
        } else {
            printf("\nДаже после перестановки система не удовлетворяет достаточному условию сходимости.\n");
            simple_iteration_method(A_rearranged, b_rearranged, solution, &iterations, &final_error);
            
            printf("\nРезультат после %d итераций:\n", iterations);
            printf("Текущая погрешность: %.8f\n", final_error);
            for (int i = 0; i < N; i++) {
                printf("x%d = %.8f\n", i + 1, solution[i]);
            }
            
            check_solution(A, b, solution);
        }
    }
    
    return 0;
}