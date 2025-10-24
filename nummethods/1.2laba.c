#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5
#define EPS 1e-10

void print_matrix(double matrix[N][N], const char* title) {
    printf("\n%s:\n", title);
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%10.6f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_vector(double vector[N], const char* title) {
    printf("\n%s: ", title);
    for (int i = 0; i < N; i++) {
        printf("%10.6f ", vector[i]);
    }
    printf("\n");
}

int lu_decomposition(double A[N][N], double LU[N][N]) {

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            LU[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < N; i++) {

        if (fabs(LU[i][i]) < EPS) {
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
    y[0] = b[0];  // L[0][0] = 1
    
    for (int i = 1; i < N; i++) {
        y[i] = b[i];
        for (int j = 0; j < i; j++) {
            y[i] -= LU[i][j] * y[j];  // LU[i][j] содержит L[i][j]
        }
    }
}

// Обратная подстановка (решение Ux = y)
void backward_substitution(double LU[N][N], double y[N], double x[N]) {
    for (int i = N - 1; i >= 0; i--) {
        x[i] = y[i];
        for (int j = i + 1; j < N; j++) {
            x[i] -= LU[i][j] * x[j];  // LU[i][j] содержит U[i][j]
        }
        x[i] /= LU[i][i];  // Делим на диагональный элемент U
    }
}

// Решение системы Ax = b через LU-разложение
int solve_lu(double A[N][N], double b[N], double x[N]) {
    double LU[N][N];
    
    // Выполняем LU-разложение
    if (!lu_decomposition(A, LU)) {
        printf("Ошибка: матрица вырождена\n");
        return 0;
    }
    
    // Решаем систему A * x = b
    double y[N];
    forward_substitution(LU, b, y);    // Ly = b
    backward_substitution(LU, y, x);   // Ux = y
    
    return 1;
}

int inverse_matrix_lu(double A[N][N], double A_inv[N][N]) {
    double LU[N][N];
    
    // Выполняем LU-разложение
    if (!lu_decomposition(A, LU)) {
        printf("Ошибка: матрица вырождена, обратной не существует\n");
        return 0;
    }
    
    // Решаем A * A_inv = I для каждого столбца единичной матрицы
    for (int j = 0; j < N; j++) {
        // j-ый столбец единичной матрицы
        double b[N] = {0.0, 0.0};
        b[j] = 1.0;
        
        double y[N], x[N];
        forward_substitution(LU, b, y);      // Ly = b
        backward_substitution(LU, y, x);     // Ux = y
        
        // Записываем результат в j-ый столбец обратной матрицы
        for (int i = 0; i < N; i++) {
            A_inv[i][j] = x[i];
        }
    }
    
    return 1;
}

double determinant_from_lu(double U[N][N]) {
    double det = 1.0;
    for (int i = 0; i < N; i++) {
        det *= U[i][i];
    }
    return det;
}

void extract_L_and_U(double LU[N][N], double L[N][N], double U[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            if (i > j) {
                L[i][j] = LU[i][j];  
                U[i][j] = 0.0;
            } else if (i < j) {
                U[i][j] = LU[i][j];  
                L[i][j] = 0.0;
            } else { 
                L[i][j] = 1.0;  
                U[i][j] = LU[i][j];
            }
        }
    }
}

int main() {

    double A_ext[N][N + 1] = {
        {1.0, -5.0, -4.0, 5.0, 3.0, 29.0},
        {-7.0, 3.0, 1.0, -4.0, 2.0, 17.0},
        {-2.0, -3.0, 4.0, 2.0, -5.0, 20.0},
        {2.0, -5.0, -3.0, 7.0, 3.0, 17.0},
        {-6.0, 3.0, -5.0, 1.0, -2.0, 24.0}
    };

    double A[N][N], b[N], LU[N][N], x[N], A_inv[N][N];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A[i][j] = A_ext[i][j];
        }
        b[i] = A_ext[i][N];
    }

    printf("Расширенная матрица [A | b]:\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j <= N; j++) {
            printf("%10.6f ", A_ext[i][j]);
        }
        printf("\n");
    }

    if (!lu_decomposition(A, LU)) {
        printf("Ошибка: LU-разложение невозможно (матрица вырожденная)\n");
        return 1;
    }

    print_matrix(LU, "Матрица LU");

    double L[N][N], U[N][N];
    extract_L_and_U(LU, L, U);

    print_matrix(L, "Матрица L (нижняя треугольная)");
    print_matrix(U, "Матрица U (верхняя треугольная)");

    printf("\nПроверка LU-разложения (L * U - A):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double product = 0.0;
            for (int k = 0; k < N; k++) {
                product += L[i][k] * U[k][j];
            }
            double error = product - A[i][j];
            printf("%10.6f ", error);
        }
        printf("\n");
    }

    if (!solve_lu(A, b, x)) {
        printf("Ошибка: Система не имеет единственного решения\n");
        return 1;
    }

    printf("\nРешение системы:\n");
    for (int i = 0; i < N; i++) {
        printf("x%d = %.6f\n", i + 1, x[i]);
    }

    printf("\nПроверка решения (A*x - b):\n");
    double residuals[N] = {0};
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            residuals[i] += A[i][j] * x[j];
        }
        residuals[i] -= b[i];
        printf("Уравнение %d: %.10f\n", i + 1, residuals[i]);
    }

    double det = determinant_from_lu(U);
    printf("\nОпределитель матрицы A: %.6f\n", det);

    if (!inverse_matrix_lu(A, A_inv)) {
        printf("Ошибка: Обратная матрица не существует\n");
        return 1;
    }

    print_matrix(A_inv, "Обратная матрица A⁻¹");

    printf("\nПроверка обратной матрицы (A * A⁻¹):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double product = 0.0;
            for (int k = 0; k < N; k++) {
                product += A[i][k] * A_inv[k][j];
            }
            printf("%10.6f ", product);
        }
        printf("\n");
    }

    printf("\nПроверка обратной матрицы (A⁻¹ * A):\n");
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            double product = 0.0;
            for (int k = 0; k < N; k++) {
                product += A_inv[i][k] * A[k][j];
            }
            printf("%10.6f ", product);
        }
        printf("\n");
    }

    return 0;
}