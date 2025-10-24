#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPSILON 0.0001
#define MAX_ITER 1000
#define N 2


void f(double x, double y, double *f1, double *f2) {
    *f1 = 2.0 * sin(x) - 2.0 * cos(y) + 1.0;
    *f2 = x*x + x*y + y*y + 2.0*x - 1.0;
}

void jacobian(double x, double y, double J[2][2]) {
    J[0][0] = 2.0 * cos(x);      // df1/dx
    J[0][1] = 2.0 * sin(y);      // df1/dy
    J[1][0] = 2.0*x + y + 2.0;   // df2/dx
    J[1][1] = x + 2.0*y;         // df2/dy
}

void phi_functions_1(double x, double y, double *phi1, double *phi2) {
    *phi1 = x + 0.002 * (2.0 * sin(x) - 2.0 * cos(y) + 1.0);
    *phi2 = y + 0.002 * (x*x + x*y + y*y + 2.0*x - 1.0);
}

void phi_functions_2(double x, double y, double *phi1, double *phi2) {
    *phi1 = x + 0.05 * (2.0 * sin(x) - 2.0 * cos(y) + 1.0);
    *phi2 = sqrt(1.0 - x*x - x*y - 2.0*x);  
}

void phi_functions_3(double x, double y, double *phi1, double *phi2) {
    *phi1 = asin((2.0 * cos(y) - 1.0) / 2.0); 
    *phi2 = y - 0.1 * (x*x + x*y + y*y + 2.0*x - 1.0);
}

void phi_functions_4(double x, double y, double *phi1, double *phi2) {
    *phi1 = x - 0.1 * (2.0 * sin(x) - 2.0 * cos(y) + 1.0);
    *phi2 = (-x - sqrt(x*x - 4.0*(x*x + 2.0*x - 1.0))) / 2.0; 
}

// Якобиан итерирующих функций
void jacobian_phi(double x, double y, void (*phi_func)(double, double, double*, double*), double J_phi[2][2]) {
    double h = 1e-6;
    double phi1, phi2;
    
    // ∂φ₁/∂x
    phi_func(x + h, y, &phi1, &phi2);
    double phi1_x_plus = phi1;
    phi_func(x - h, y, &phi1, &phi2);
    double phi1_x_minus = phi1;
    J_phi[0][0] = (phi1_x_plus - phi1_x_minus) / (2.0 * h);
    
    // ∂φ₁/∂y
    phi_func(x, y + h, &phi1, &phi2);
    double phi1_y_plus = phi1;
    phi_func(x, y - h, &phi1, &phi2);
    double phi1_y_minus = phi1;
    J_phi[0][1] = (phi1_y_plus - phi1_y_minus) / (2.0 * h);
    
    // ∂φ₂/∂x
    phi_func(x + h, y, &phi1, &phi2);
    double phi2_x_plus = phi2;
    phi_func(x - h, y, &phi1, &phi2);
    double phi2_x_minus = phi2;
    J_phi[1][0] = (phi2_x_plus - phi2_x_minus) / (2.0 * h);
    
    // ∂φ₂/∂y
    phi_func(x, y + h, &phi1, &phi2);
    double phi2_y_plus = phi2;
    phi_func(x, y - h, &phi1, &phi2);
    double phi2_y_minus = phi2;
    J_phi[1][1] = (phi2_y_plus - phi2_y_minus) / (2.0 * h);
}

double matrix_norm(double A[2][2]) {
    double norm1 = fabs(A[0][0]) + fabs(A[0][1]);
    double norm2 = fabs(A[1][0]) + fabs(A[1][1]);
    return (norm1 > norm2) ? norm1 : norm2;
}

// LU-декомпозиция без выбора главного элемента
int lu_decomposition(double A[N][N], double LU[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            LU[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        if (fabs(LU[i][i]) < 1e-12) {
            return 0;
        }

        for (int k = i + 1; k < N; k++) {
            double factor = LU[k][i] / LU[i][i];
            LU[k][i] = factor; 

            for (int j = i + 1; j < N; j++) {
                LU[k][j] -= factor * LU[i][j];
            }
        }
    }
    return 1;
}

// Прямая подстановка (решение Ly = b)
void forward_substitution(double LU[N][N], double b[N], double y[N]) {
    y[0] = b[0];
    
    for (int i = 1; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];
        }
    }
}

// Обратная подстановка (решение Ux = y)
void backward_substitution(double LU[N][N], double y[N], double x[N]) {
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= LU[i][j] * x[j]; 
        }
        x[i] /= LU[i][i]; 
    }
}

// Решение системы Ax = b через LU-разложение
int solve_system_lu(double A[N][N], double b[N], double x[N]) {
    double LU[N][N];
    
    if (!lu_decomposition(A, LU)) {
        printf("Ошибка: матрица вырождена\n");
        return 0;
    }
    
    // Решаем систему A * x = b
    double y[N];
    forward_substitution(LU, b, y);    
    backward_substitution(LU, y, x); 
    
    return 1;
}

int inverse_matrix_lu(double A[N][N], double A_inv[N][N]) {
    double LU[N][N];
    
    if (!lu_decomposition(A, LU)) {
        printf("Ошибка: матрица вырождена, обратной не существует\n");
        return 0;
    }
    
    for (int j = 0; j < N; j++) {
        double b[N] = {0.0, 0.0};
        b[j] = 1.0;
        
        double y[N], x[N];
        forward_substitution(LU, b, y);  
        backward_substitution(LU, y, x);     
        
        for (int i = 0; i < N; i++) {
            A_inv[i][j] = x[i];
        }
    }
    
    return 1;
}

// Умножение матрицы на вектор
void matrix_vector_mult(double A[N][N], double v[N], double result[N]) {
    for (int i = 0; i < N; i++) {
        result[i] = 0.0;
        for (int j = 0; j < N; j++) {
            result[i] += A[i][j] * v[j];
        }
    }
}

// Проверка условия сходимости для метода Ньютона
double check_newton_convergence(double x, double y) {
    double J[N][N];
    jacobian(x, y, J);
    
    // Проверяем, что якобиан не вырожден
    double det = J[0][0]*J[1][1] - J[0][1]*J[1][0];
    return fabs(det);
}

int newton_method(double x0, double y0, double *x_root, double *y_root) {
    printf("\n1. Метод Ньютона:\n");

    double x = x0, y = y0;
    int iter = 0;

    double det = check_newton_convergence(x0, y0);
    printf("Определитель якобиана в начальной точке: %.6f\n", det);
    if (det < 1e-12) {
        printf("Якобиан близок к вырожденному - сходимость не гарантирована!\n");
    }

    do {
        double f1, f2;
        f(x, y, &f1, &f2);
        
        // Критерий остановки
        if (fabs(f1) < EPSILON && fabs(f2) < EPSILON) {
            *x_root = x;
            *y_root = y;
            printf("Корень найден: (%.8f, %.8f)\n", x, y);
            printf("Значения функций: f1=%.8f, f2=%.8f\n", f1, f2);
            printf("Количество итераций: %d\n", iter);
            return iter;
        }
        
        // Вычисляем якобиан
        double J[N][N];
        jacobian(x, y, J);
        
        double J_inv[N][N];
        if (!inverse_matrix_lu(J, J_inv)) {
            printf("Ошибка. Матрица вырождена)\n");
            return -1;
        }
        
        // delta = -J^{-1} * f
        double f_vec[N] = {f1, f2};
        double delta[N];
        matrix_vector_mult(J_inv, f_vec, delta);
        
        delta[0] = -delta[0];
        delta[1] = -delta[1];
        
        double x_new = x + delta[0];
        double y_new = y + delta[1];
        
        if (sqrt(delta[0]*delta[0] + delta[1]*delta[1]) < EPSILON) {
            *x_root = x_new;
            *y_root = y_new;
            printf("Корень найден: (%.8f, %.8f)\n", x_new, y_new);
            printf("Значения функций: f1=%.8f, f2=%.8f\n", f1, f2);
            printf("Количество итераций: %d\n", iter + 1);
            return iter + 1;
        }
        
        x = x_new;
        y = y_new;
        iter++;
        
        if (iter >= MAX_ITER) {
            printf("Превышено максимальное число итераций\n");
            *x_root = x;
            *y_root = y;
            return iter;
        }
        
    } while (1);
}

int simple_iteration_method(double x0, double y0, void (*phi_func)(double, double, double*, double*), 
                           double *x_root, double *y_root) {
    printf("\n2. Метод простой итерации:\n");
    double x = x0, y = y0;
    int iter = 0;

    // Проверка условия сходимости
    double J_phi[2][2];
    jacobian_phi(x0, y0, phi_func, J_phi);
    double q = matrix_norm(J_phi);
    printf("Норма якобиана итерирующих функций: %.6f\n", q);
    if (q >= 1.0) {
        printf("Условие сходимости не выполняется (q >= 1)\n");
    }

    do {
        double phi1, phi2;
        phi_func(x, y, &phi1, &phi2);
        
        double x_new = phi1;
        double y_new = phi2;
        
        double dx = x_new - x;
        double dy = y_new - y;
        
        if (sqrt(dx*dx + dy*dy) < EPSILON) {
            *x_root = x_new; *y_root = y_new;
            
            double f1, f2;
            f(x_new, y_new, &f1, &f2);
            printf("Корень найден: (%.6f, %.6f)\n", x_new, y_new);
            printf("Значения функций: f1=%.8f, f2=%.8f\n", f1, f2);
            printf("Количество итераций: %d\n", iter + 1);
            return iter + 1;
        }
        
        x = x_new; y = y_new;
        iter++;
        
        if (iter >= MAX_ITER) {
            printf("Достигнуто максимальное число итераций\n");
            *x_root = x; *y_root = y;
            return iter;
        }
    } while (1);
}

int seidel_method(double x0, double y0, void (*phi_func)(double, double, double*, double*), 
                 double *x_root, double *y_root) {
    printf("\n3. Метод Зейделя:\n");
    double x = x0, y = y0;
    int iter = 0;

    // Проверка условия сходимости
    double J_phi[2][2];
    jacobian_phi(x0, y0, phi_func, J_phi);
    double q = matrix_norm(J_phi);
    printf("Норма якобиана итерирующих функций: %.6f\n", q);
    if (q >= 1.0) {
        printf("Условие сходимости не выполняется (q >= 1)\n");
    }

    do {
        double phi1, phi2;
        
        // Сначала обновляем x
        phi_func(x, y, &phi1, &phi2);
        double x_new = phi1;
        
        // Затем обновляем y с новым x
        phi_func(x_new, y, &phi1, &phi2);
        double y_new = phi2;
        
        double dx = x_new - x;
        double dy = y_new - y;
        
        if (sqrt(dx*dx + dy*dy) < EPSILON) {
            *x_root = x_new; *y_root = y_new;
            
            double f1, f2;
            f(x_new, y_new, &f1, &f2);
            printf("Корень найден: (%.6f, %.6f)\n", x_new, y_new);
            printf("Значения функций: f1=%.8f, f2=%.8f\n", f1, f2);
            printf("Количество итераций: %d\n", iter + 1);
            return iter + 1;
        }
        
        x = x_new; y = y_new;
        iter++;
        
        if (iter >= MAX_ITER) {
            printf("Достигнуто максимальное число итераций\n");
            *x_root = x; *y_root = y;
            return iter;
        }
    } while (1);
}

int main() {
    printf("Решение системы нелинейных уравнений:\n");
    printf("  f1(x, y) = 2*sin(x) - 2*cos(y) + 1 = 0\n");
    printf("  f2(x, y) = x^2 + x*y + y^2 + 2*x - 1 = 0\n\n");

    struct {
        double x0, y0;
        void (*phi_func)(double, double, double*, double*);
        const char *description;
    } 
    
    test_cases[] = {
        {-3.0, 1.1, phi_functions_1},
        {-0.1, 1.2, phi_functions_2},
        {0.4, -0.4, phi_functions_3,},
        {0.0, -0.9, phi_functions_4}
    };
    
    int n_cases = 4;

    for (int i = 0; i < n_cases; i++) {
        printf("Начальное приближение: (%.4f, %.4f)\n", 
               test_cases[i].x0, test_cases[i].y0);

        double x_newton, y_newton;
        double x_simple, y_simple;
        double x_seidel, y_seidel;

        newton_method(test_cases[i].x0, test_cases[i].y0, &x_newton, &y_newton);
        simple_iteration_method(test_cases[i].x0, test_cases[i].y0, test_cases[i].phi_func, &x_simple, &y_simple);
        seidel_method(test_cases[i].x0, test_cases[i].y0, test_cases[i].phi_func, &x_seidel, &y_seidel);

        printf("\n");
    }

    return 0;
}