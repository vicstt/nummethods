#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 8

void print_vector(double vector[], int n, const char* title) {
    printf("%s: ", title);
    for (int i = 0; i < n; i++) {
        printf("%8.3f ", vector[i]);
    }
    printf("\n");
}

int check_stability_conditions(double a[], double b[], double c[], int n) {
    int conditions_met = 1;
    
    printf("Проверка условий устойчивости метода прогонки:\n");
    
    for (int i = 0; i < n; i++) {
        if (i == 0) {
            int condition = fabs(b[0]) >= fabs(c[0]);
            printf("Условие 1 (i=0): |%.1f| >= |%.1f| - %s\n", 
                   b[0], c[0], condition ? "выполнено" : "не выполнено");
            if (!condition) conditions_met = 0;
        }
        else if (i == n - 1) {
            int condition = fabs(b[n - 1]) >= fabs(a[n - 2]);
            printf("Условие 1 (i=%d): |%.1f| >= |%.1f| - %s\n", 
                   n-1, b[n-1], a[n-2], condition ? "выполнено" : "не выполнено");
            if (!condition) conditions_met = 0;
        }
        else {
            int condition = fabs(b[i]) >= fabs(a[i - 1]) + fabs(c[i]);
            printf("Условие 1 (i=%d): |%.1f| >= |%.1f| + |%.1f| - %s\n", 
                   i, b[i], a[i-1], c[i], condition ? "выполнено" : "не выполнено");
            if (!condition) conditions_met = 0;
        }
    }
    
    double condition2 = fabs(c[0] / b[0]);
    int cond2_ok = condition2 < 1;
    printf("Условие 2: |%.1f/%.1f| = %.3f < 1 - %s\n", 
           c[0], b[0], condition2, cond2_ok ? "выполнено" : "не выполнено");
    if (!cond2_ok) conditions_met = 0;
    
    double condition3 = fabs(a[n - 2] / b[n - 1]);
    int cond3_ok = condition3 < 1;
    printf("Условие 3: |%.1f/%.1f| = %.3f < 1 - %s\n", 
           a[n-2], b[n-1], condition3, cond3_ok ? "выполнено" : "не выполнено");
    if (!cond3_ok) conditions_met = 0;
    
    if (conditions_met) {
        printf("Все условия устойчивости выполнены\n");
    } else {
        printf("Не все условия устойчивости выполнены\n");
    }
    
    return conditions_met;
}

void thomas_solve(double a[], double b[], double c[], double d[], double x[], int n) {
    double b_copy[N], d_copy[N];
    for (int i = 0; i < n; i++) {
        b_copy[i] = b[i];
        d_copy[i] = d[i];
    }
    
    for (int i = 1; i < n; i++) {
        double m = a[i - 1] / b_copy[i - 1];
        b_copy[i] = b_copy[i] - m * c[i - 1];
        d_copy[i] = d_copy[i] - m * d_copy[i - 1];
    }
    
    x[n - 1] = d_copy[n - 1] / b_copy[n - 1];
    
    for (int i = n - 2; i >= 0; i--) {
        x[i] = (d_copy[i] - c[i] * x[i + 1]) / b_copy[i];
    }
}

double determinant(double a[], double b[], double c[], int n) {
    double det_prev = 1.0;
    double det_current = b[0];
    
    for (int i = 1; i < n; i++) {
        double det_next = b[i] * det_current - a[i - 1] * c[i - 1] * det_prev;
        det_prev = det_current;
        det_current = det_next;
    }
    
    return det_current;
}

void check_solution(double matrix[N][N], double right_side[], double solution[], int n) {
    printf("\nПроверка решения:\n");
    for (int i = 0; i < n; i++) {
        double calculated = 0.0;
        for (int j = 0; j < n; j++) {
            calculated += matrix[i][j] * solution[j];
        }
        printf("Уравнение %d: вычислено %8.3f, должно быть %8.3f\n", 
               i + 1, calculated, right_side[i]);
    }
}

void solve_system() {
    double matrix[N][N] = {
        {7, 3, 0, 0, 0, 0, 0, 0},
        {5, 10, -4, 0, 0, 0, 0, 0},
        {0, 3, 11, 6, 0, 0, 0, 0},
        {0, 0, 5, -11, -5, 0, 0, 0},
        {0, 0, 0, 2, 13, -9, 0, 0},
        {0, 0, 0, 0, 3, 11, 8, 0},
        {0, 0, 0, 0, 0, 3, 8, -2},
        {0, 0, 0, 0, 0, 0, 2, 5}
    };
    
    double right_side[N] = {27, 31, 59, -67, 19, 15, 26, 30};
    
    double a[N-1];  
    double b[N];    
    double c[N-1]; 
    
    for (int i = 0; i < N; i++) {
        b[i] = matrix[i][i];
        if (i > 0) {
            a[i-1] = matrix[i][i-1];
        }
        if (i < N-1) {
            c[i] = matrix[i][i+1];
        }
    }
    
    int stability_ok = check_stability_conditions(a, b, c, N);
    
    double solution[N];
    thomas_solve(a, b, c, right_side, solution, N);
    
    double det_value = determinant(a, b, c, N);
    
    printf("\nРешение системы методом прогонки:\n");
    for (int i = 0; i < N; i++) {
        printf("x%d = %8.6f\n", i + 1, solution[i]);
    }
    
    printf("\nОпределитель матрицы: %.6f\n", det_value);
    
    check_solution(matrix, right_side, solution, N);
}

int main() {
    solve_system();
    return 0;
}