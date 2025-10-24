#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define N 5
#define EPS 1e-10

void print_matrix(double matrix[N][N]) {
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            printf("%10.6f ", matrix[i][j]);
        }
        printf("\n");
    }
}

void print_vector(double vector[N]) {
    for (int i = 0; i < N; i++) {
        printf("%10.6f ", vector[i]);
    }
    printf("\n");
}

int gauss_elimination(double A[N][N], double b[N], double x[N]) {
    double Ab[N][N+1];

    // Создаем расширенную матрицу [A|b]
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            Ab[i][j] = A[i][j];
        }
        Ab[i][N] = b[i];
    }

    // Прямой ход метода Гаусса 
    for (int i = 0; i < N; i++) {
        if (fabs(Ab[i][i]) < EPS) {
            return 0; 
        }

        double pivot = Ab[i][i];
        for (int j = i; j <= N; j++) {
            Ab[i][j] /= pivot;
        }

        for (int k = i + 1; k < N; k++) {
            double factor = Ab[k][i]; 
            for (int j = i; j <= N; j++) {
                Ab[k][j] -= factor * Ab[i][j];
            }
        }
    }

    // Обратный ход 
    for (int i = N - 1; i >= 0; i--) {
        x[i] = Ab[i][N]; 
        for (int j = i + 1; j < N; j++) {
            x[i] -= Ab[i][j] * x[j]; 
        }
    }

    return 1; 
}

double determinant(double A[N][N]) {
    double det = 1.0;
    double A_copy[N][N];

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            A_copy[i][j] = A[i][j];
        }
    }

    for (int i = 0; i < N; i++) {
        if (fabs(A_copy[i][i]) < EPS) {
            return 0.0;
        }

        det *= A_copy[i][i];

        double pivot = A_copy[i][i];
        for (int j = i; j < N; j++) {
            A_copy[i][j] /= pivot;
        }

        for (int k = i + 1; k < N; k++) {
            double factor = A_copy[k][i];
            for (int j = i + 1; j < N; j++) {
                A_copy[k][j] -= factor * A_copy[i][j];
            }
        }
    }

    return det;
}

int inverse_matrix(double A[N][N], double A_inv[N][N]) {
    double b[N];      
    double x[N]; 

    for (int j = 0; j < N; j++) {
        for (int i = 0; i < N; i++) {
            b[i] = (i == j) ? 1.0 : 0.0;
        }

        if (!gauss_elimination(A, b, x)) {
            return 0; 
        }

        for (int i = 0; i < N; i++) {
            A_inv[i][j] = x[i];
        }
    }

    return 1; 
}

int main() {
    double A[N][N] = {
            {1.0, -5.0, -4.0, 5.0, 3.0},
            {-7.0, 3.0, 1.0, -4.0, 2.0},
            {-2.0, -3.0, 4.0, 2.0, -5.0},
            {2.0, -5.0, -3.0, 7.0, 3.0},
            {-6.0, 3.0, -5.0, 1.0, -2.0}
    };

    double b[N] = {29.0, 17.0, 20.0, 17.0, 24.0};
    double x[N];
    double A_inv[N][N];

    printf("Матрица A:\n");
    print_matrix(A);

    printf("\nВектор b:\n");
    print_vector(b);

    if (gauss_elimination(A, b, x)) {
        printf("\nРешение системы:\n");
        for (int i = 0; i < N; i++) {
            printf("x%d = %.6f\n", i + 1, x[i]);
        }

        printf("\nПроверка решения (A*x, b):\n");
        for (int i = 0; i < N; i++) {
            double calculated = 0.0;
            for (int j = 0; j < N; j++) {
                calculated += A[i][j] * x[j];
            }
            printf("Уравнение %d: вычислено %.6f, должно быть %.6f\n",
                   i + 1, calculated, b[i]);
        }
    } else {
        printf("\nСистема не имеет единственного решения!\n");
    }

    double det = determinant(A);
    printf("\nОпределитель матрицы A: %.6f\n", det);

    if (inverse_matrix(A, A_inv)) {
        printf("\nОбратная матрица A⁻¹:\n");
        print_matrix(A_inv);

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
    } else {
        printf("\nОбратная матрица не существует (матрица вырожденная)!\n");
    }

    return 0;
}