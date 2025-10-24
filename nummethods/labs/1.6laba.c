#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <complex.h> 

#define N 5
#define EPSILON 0.001
#define MIN_TOLERANCE 1e-12

void print_matrix(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%10.6f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void matrix_copy(double src[N][N], double dest[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            dest[i][j] = src[i][j];
        }
    }
}

void matrix_multiply(double A[N][N], double B[N][N], double C[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            C[i][j] = 0.0;
            for (int k = 0; k < N; k++) {
                C[i][j] += A[i][k] * B[k][j];
            }
        }
    }
}

void matrix_transpose(double A[N][N], double result[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            result[i][j] = A[j][i];
        }
    }
}

double vector_norm(double v[N]) {
    double norm_sq = 0.0;
    for (int i = 0; i < N; i++) {
        norm_sq += v[i] * v[i];
    }
    return sqrt(norm_sq);
}

void householder_reflection(double v[N], int n, double H[N][N]) {
    double v_norm = vector_norm(v);
    double u[N];
    
    for (int i = 0; i < n; i++) {
        u[i] = v[i];
    }
    
    // Первый столбец
    if (v[0] >= 0) {
        u[0] += v_norm;
    } else {
        u[0] -= v_norm;
    }
    
    double u_norm = vector_norm(u);
    
    // Проверка на нулевой вектор
    if (fabs(u_norm) < MIN_TOLERANCE) {
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                H[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        return;
    }
    
    // Нормализуем u
    for (int i = 0; i < n; i++) {
        u[i] /= u_norm;
    }
    
    // Заполняем остальные элементы нулями (для полной матрицы N×N)
    for (int i = n; i < N; i++) {
        u[i] = 0.0;
    }
    
    // Создаем матрицу Хаусхолдера: H = I - 2*u*u^T
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            H[i][j] = -2.0 * u[i] * u[j];
        }
        H[i][i] += 1.0;
    }
}

void qr_decomposition(double A[N][N], double Q[N][N], double R[N][N]) {
    double Q_temp[N][N], R_temp[N][N], H_k[N][N], H_k_T[N][N];
    
    // Инициализация Q как единичной матрицы
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Q[i][j] = (i == j) ? 1.0 : 0.0;
        }
    }
    
    matrix_copy(A, R);
    
    for (int k = 0; k < N; k++) {
        // Создаем подвектор x для преобразования Хаусхолдера
        double x[N];
        int sub_size = N - k;
        
        for (int i = 0; i < N; i++) {
            if (i < k) {
                x[i] = 0.0;
            } else {
                x[i] = R[i][k];
            }
        }
        
        // Проверка на нулевой подвектор
        double x_sub[N];
        for (int i = 0; i < sub_size; i++) {
            x_sub[i] = x[k + i];
        }
        
        if (vector_norm(x_sub) < MIN_TOLERANCE) {
            continue;
        }
        
        // Создаем преобразование Хаусхолдера для подвектора
        double H_k_small[N][N];
        householder_reflection(x_sub, sub_size, H_k_small);
        
        // Расширяем преобразование до полной матрицы
        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                H_k[i][j] = (i == j) ? 1.0 : 0.0;
            }
        }
        
        for (int i = k; i < N; i++) {
            for (int j = k; j < N; j++) {
                H_k[i][j] = H_k_small[i - k][j - k];
            }
        }
        
        // Применяем преобразование к R: R = H_k * R
        matrix_multiply(H_k, R, R_temp);
        matrix_copy(R_temp, R);
        
        // Обновляем Q: Q = Q * H_k^T
        matrix_transpose(H_k, H_k_T);
        matrix_multiply(Q, H_k_T, Q_temp);
        matrix_copy(Q_temp, Q);
    }
}

void extract_complex_eigenvalues(double B[2][2], double complex *lambda1, double complex *lambda2) {
    double a = B[0][0];
    double b = B[0][1];
    double c = B[1][0];
    double d = B[1][1];

    double trace = a + d;
    double det = a * d - b * c;

    double discriminant = trace * trace - 4 * det;

    if (discriminant < 0) {
        double real_part = trace / 2.0;
        double imag_part = sqrt(-discriminant) / 2.0;
        *lambda1 = real_part + I * imag_part;
        *lambda2 = real_part - I * imag_part;
    } else {
        double sqrt_disc = sqrt(discriminant);
        *lambda1 = (trace + sqrt_disc) / 2.0;
        *lambda2 = (trace - sqrt_disc) / 2.0;
    }
}

void qr_algorithm(double A[N][N], double eigenvalues[N], int *iterations) {
    double Ak[N][N], Ak_next[N][N], Q[N][N], R[N][N];
    
    matrix_copy(A, Ak);
    *iterations = 0;
    
    printf("QR-алгоритм для нахождения собственных значений\n");
    
    for (*iterations = 1; *iterations < 50; (*iterations)++) {
        qr_decomposition(Ak, Q, R);
        matrix_multiply(R, Q, Ak_next);
        matrix_copy(Ak_next, Ak);
        
        // Проверка сходимости
        double max_off_diag = 0.0;
        for (int i = 1; i < N; i++) {
            for (int j = 0; j < i; j++) {
                if (fabs(Ak[i][j]) > max_off_diag) {
                    max_off_diag = fabs(Ak[i][j]);
                }
            }
        }
        
        if (max_off_diag < EPSILON) {
            break;
        }
    }
        
    // Извлекаем собственные значения из диагонали
    for (int i = 0; i < N; i++) {
        eigenvalues[i] = Ak[i][i];
    }

    printf("\nИтоговая матрица:\n");
    print_matrix(Ak);

    for (int i = 0; i < N; i++) {
        if (i < N - 1 && fabs(Ak[i + 1][i]) > EPSILON) {
            double B[2][2] = {
                {Ak[i][i], Ak[i][i+1]},
                {Ak[i+1][i], Ak[i+1][i+1]}
            };
            double complex lambda1, lambda2;
            extract_complex_eigenvalues(B, &lambda1, &lambda2);
            printf("Комплексные собственные значения: %.6f + %.6fi, %.6f %.6fi\n",
                   creal(lambda1), cimag(lambda1), creal(lambda2), cimag(lambda2));
            i++;
        } else {
            printf("Вещественное собственное значение: %.8f\n", Ak[i][i]);
        }
    }
}

int main() {
    double A[N][N] = {
        {1.0, -5.0, -4.0, 7.0, 6.0},
        {-4.0, 3.0, 5.0, 2.0, 2.0},
        {0.0, 3.0, 9.0, 7.0, 5.0},
        {7.0, 5.0, -4.0, 1.0, 3.0},
        {-9.0, 3.0, -5.0, 1.0, 0.0}
    };
    
    printf("Исходная матрица:\n");
    print_matrix(A);
    
    printf("\n");
    
    double eigenvalues[N];
    int iterations;
    
    qr_algorithm(A, eigenvalues, &iterations);
    
    printf("\nКоличество итераций: %d\n", iterations);
    
    // Проверка: след матрицы должен равняться сумме собственных значений
    double trace = 0.0;
    for (int i = 0; i < N; i++) {
        trace += A[i][i];
    }
    
    double sum_eigenvalues = 0.0;
    for (int i = 0; i < N; i++) {
        sum_eigenvalues += eigenvalues[i];
    }
    
    printf("\nПроверка:\n");
    printf("След матрицы (сумма диагональных элементов) = %.6f\n", trace);
    printf("Сумма собственных значений = %.6f\n", sum_eigenvalues);
    
    return 0;
}
