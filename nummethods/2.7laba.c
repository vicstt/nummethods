#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#define EPSILON 0.0001

double f(double x) {
    return 0.25 * sinh(x * x) - 5 * x * x - 25 * x + 21;
}

double f_prime(double x) {
    return 0.25 * 2 * x * cosh(x * x) - 10 * x - 25;
}

double f_double_prime(double x) {
    return 0.25 * (2 * cosh(x * x) + 4 * x * x * sinh(x * x)) - 10;
}

double phi(double x, double lambda) {
    return x + lambda * f(x);
}

// Проверка достаточного условия сходимости: |f(x)*f''(x)| < (f'(x))^2
int check_convergence_condition(double a, double b) {

    int n_check = 1000;
    int condition_violated = 0; 

    for (int i = 0; i <= n_check; i++) {
        double x = a + i * (b - a) / n_check;
        double fx = f(x);
        double fp = f_prime(x);
        double fpp = f_double_prime(x);

        if (fabs(fp) < 1e-8) {
            printf("Ошибка: f'(x) ≈ 0 в точке x=%.6f. Условие не применимо.\n", x);
            return 0;
        }

        double left = fabs(fx * fpp);
        double right = fp * fp;

        if (left >= right) {
            if (!condition_violated) {
                condition_violated = 1;
            }
        }
    }

    if (!condition_violated) {
        printf("Достаточное условие сходимости выполнено.\n");
    } else {
        printf("Достаточное условие сходимости нарушено.\n");
    }

    return 1;
}

// Функция для выбора начального приближения по правилу: f(x0)*f''(x0) > 0
int choose_initial_point(double a, double b, double *x0) {
    double fa = f(a);
    double fb = f(b);
    double fpp_a = f_double_prime(a);
    double fpp_b = f_double_prime(b);
    
    int case_a = (fa * fpp_a > 0);
    int case_b = (fb * fpp_b > 0);
    
    if (case_a && case_b) {
        if (fabs(fa) < fabs(fb)) {
            *x0 = a;
        } else {
            *x0 = b;
        }
        return 3;
    }
    else if (case_a) {
        *x0 = a;
        return 1;
    }
    else if (case_b) {
        *x0 = b;
        return 2;
    }
    else {
        printf("Метод может не сходиться!\n");
        *x0 = (a + b) / 2.0;
        return 4;
    }
}

// Метод простой итерации
int simple_iteration(double a, double b, double lambda, double *root) {
    printf("\n1. Метод простой итерации\n");
    printf("Интервал: [%.4f, %.4f], lambda = %.3f\n", a, b, lambda);

    if (f(a) * f(b) >= 0) {
        printf("Ошибка: f(a)*f(b) >= 0. Условие сходимости не выполнено.\n");
        return -1;
    }

    // Проверка |phi'(x)| <= q < 1 на [a,b]
    double max_phi_prime = 0.0;
    int n_check = 1000; 
    for (int i = 0; i <= n_check; i++) {
        double x = a + i * (b - a) / n_check;
        double phi_prime = 1.0 + lambda * f_prime(x);
        if (fabs(phi_prime) > max_phi_prime) max_phi_prime = fabs(phi_prime);
    }

    if (max_phi_prime >= 1.0) {
        printf("Условие сходимости не выполняется (q >= 1)\n");
    } else {
        printf("Условие сходимости выполнено.\n");
    }

    double x_prev = (a + b) / 2.0;
    double x_curr;
    int iter = 0;

    do {
        x_curr = phi(x_prev, lambda);
        iter++;

        if (iter > 1000) {
            printf("Превышено максимальное число итераций (%d).\n", iter);
            return -2;
        }

        if (fabs(x_curr - x_prev) < EPSILON) {
            *root = x_curr;
            printf("Корень найден: x = %.8f\n", *root);
            printf("Количество итераций: %d\n", iter);
            return iter;
        }

        x_prev = x_curr;
    } while (1);
}

// Метод Ньютона
int newton_method(double a, double b, double *root) {
    printf("\n2. Метод Ньютона\n");
    printf("Интервал: [%.4f, %.4f]\n", a, b);

    if (f(a) * f(b) >= 0) {
        printf("Ошибка: f(a)*f(b) >= 0. Условие сходимости не выполнено.\n");
    }

    // Проверка f'(x) != 0 на [a,b]
    int n_check = 1000;
    for (int i = 0; i <= n_check; i++) {
        double x = a + i * (b - a) / n_check;
        if (fabs(f_prime(x)) < 1e-8) {
            printf("Ошибка: f'(x) ≈ 0 в точке x=%.6f. Условие сходимости не выполнено.\n", x);
            return -2;
        }
    }

    // Проверка постоянства знака f''(x) и достаточного условия сходимости
    double fpp_a = f_double_prime(a);
    double fpp_b = f_double_prime(b);
    if (fpp_a * fpp_b < 0) {
        printf("f''(x) меняет знак на [a,b]. Условие сходимости не выполнено.\n");
    }

    if (!check_convergence_condition(a, b)) {
        return -2;
    }

    double x0;
    int case_type = choose_initial_point(a, b, &x0);
    double x_prev = x0;
    double x_curr;
    int iter = 0;

    do {
        double fp = f_prime(x_prev);
        if (fabs(fp) < 1e-8) {
            printf("Ошибка: производная близка к нулю на итерации %d.\n", iter);
            return -3;
        }

        x_curr = x_prev - f(x_prev) / fp;
        iter++;

        if (iter > 1000) {
            printf("Превышено максимальное число итераций (%d).\n", iter);
            return -4;
        }

        if (fabs(x_curr - x_prev) < EPSILON) {
            *root = x_curr;
            printf("Корень найден: x = %.8f\n", *root);
            printf("Количество итераций: %d\n", iter);
            return iter;
        }

        x_prev = x_curr;
    } while (1);
}

// Метод секущих
int secant_method(double a, double b, double *root) {
    printf("\n3. Метод секущих\n");
    printf("Интервал: [%.4f, %.4f]\n", a, b);

    // Проверка f(a)*f(b) < 0
    if (f(a) * f(b) >= 0) {
        printf("Ошибка: f(a)*f(b) >= 0. Условие сходимости не выполнено.\n");
    }

    // Проверка f'(x) != 0 на [a,b]
    int n_check = 1000;
    for (int i = 0; i <= n_check; i++) {
        double x = a + i * (b - a) / n_check;
        if (fabs(f_prime(x)) < 1e-8) {
            printf("Предупреждение: f'(x) ≈ 0 в точке x=%.6f. Условие сходимости не выполнено.\n", x);
        }
    }


    // Проверка постоянства знака f''(x)
    double fpp_a = f_double_prime(a);
    double fpp_b = f_double_prime(b);
    if (fpp_a * fpp_b < 0) {
        printf("Предупреждение: f''(x) меняет знак на [a,b]. Условие сходимости не выполнено.\n");
    }

    if (!check_convergence_condition(a, b)) {
        return -2;
    }

    double x0;
    int case_type = choose_initial_point(a, b, &x0);

    double fp_x0 = f_prime(x0);
    if (fabs(fp_x0) < 1e-8) {
        printf("Ошибка: производная в x0 близка к нулю.\n");
        return -2;
    }

    double x1 = x0 - f(x0) / fp_x0;

    double x_prev2 = x0; // x(k-1)
    double x_prev1 = x1; // x(k)
    double x_curr;
    int iter = 1;

    do {
        double denom = f(x_prev1) - f(x_prev2);
        if (fabs(denom) < 1e-8) {
            printf("Ошибка: знаменатель близок к нулю на итерации %d.\n", iter);
            return -3;
        }

        x_curr = x_prev1 - f(x_prev1) * (x_prev1 - x_prev2) / denom;
        iter++;

        if (iter > 1000) {
            printf("Превышено максимальное число итераций (%d).\n", iter);
            return -4;
        }

        if (fabs(x_curr - x_prev1) < EPSILON) {
            *root = x_curr;
            printf("Корень найден: x = %.8f\n", *root);
            printf("Количество итераций: %d\n", iter);
            return iter;
        }

        x_prev2 = x_prev1;
        x_prev1 = x_curr;
    } while (1);
}

// Метод хорд
int chord_method(double a, double b, double *root) {
    printf("\n4. Метод хорд\n");
    printf("Интервал: [%.4f, %.4f]\n", a, b);

    // Проверка f(a)*f(b) < 0
    if (f(a) * f(b) >= 0) {
        printf("Ошибка: f(a)*f(b) >= 0. Условие сходимости не выполнено.\n");
        return -1;
    }

    if (!check_convergence_condition(a, b)) {
        return -2;
    }

    double z; 
    int case_type = choose_initial_point(a, b, &z);

    // Определяем неподвижный конец
    double fixed_end, moving_point;
    if (f(a) * f_double_prime(a) > 0) {
        fixed_end = a;      
        moving_point = b;   
    } else {
        fixed_end = b;      
        moving_point = a; 
    }

    double x_prev = moving_point;
    double x_curr;
    int iter = 0;

    do {
        double f_prev = f(x_prev);
        double f_fixed = f(fixed_end);
        double denom = f_prev - f_fixed;
        if (fabs(denom) < 1e-8) {
            printf("Ошибка: знаменатель близок к нулю на итерации %d.\n", iter);
            return -2;
        }

        x_curr = x_prev - f_prev * (x_prev - fixed_end) / denom;
        iter++;

        if (iter > 1000) {
            printf("Превышено максимальное число итераций (%d).\n", iter);
            return -3;
        }

        if (fabs(x_curr - x_prev) < EPSILON && fabs(f(x_curr)) < EPSILON) {
            *root = x_curr;
            printf("Корень найден: x = %.8f\n", *root);
            printf("Количество итераций: %d\n", iter);
            return iter;
        }

        x_prev = x_curr;
    } while (1);
}

// Метод дихотомии
int bisection_method(double a, double b, double *root) {
    printf("\n5. Метод дихотомии\n");
    printf("Интервал: [%.4f, %.4f]\n", a, b);

    // Проверка f(a)*f(b) < 0
    if (f(a) * f(b) >= 0) {
        printf("Ошибка: f(a)*f(b) >= 0. Условие сходимости не выполнено.\n");
        return -1;
    }

    double x_prev = a;
    double x_curr = b;
    int iter = 0;

    do {
        double c = (x_prev + x_curr) / 2.0;
        double fc = f(c);

        if (fc == 0.0) {
            *root = c;
            printf("Корень найден точно: x = %.8f\n", *root);
            printf("Количество итераций: %d\n", iter + 1);
            return iter + 1;
        }

        if (f(x_prev) * fc < 0) {
            x_curr = c;
        } else {
            x_prev = c;
        }

        iter++;

        if (iter > 1000) {
            printf("Превышено максимальное число итераций (%d).\n", iter);
            return -2;
        }

        if (fabs(x_curr - x_prev) < EPSILON) {
            *root = (x_prev + x_curr) / 2.0;
            printf("Корень найден: x = %.8f\n", *root);
            printf("Количество итераций: %d\n", iter);
            return iter;
        }
    } while (1);
}

int main() {
    printf("Решение уравнения: 0.25 * sinh(x^2) - 5*x^2 - 25*x + 21 = 0\n");

    printf("\nКорень 1: в интервале [0, 1]\n");

    double root1;
    int iter1;

    // Метод простой итерации 
    iter1 = simple_iteration(0.0, 1.0, 0.01, &root1);
    if (iter1 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root1, f(root1));
    }

    // Метод Ньютона
    iter1 = newton_method(0.0, 1.0, &root1);
    if (iter1 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root1, f(root1));
    }

    // Метод секущих
    iter1 = secant_method(0.0, 1.0, &root1);
    if (iter1 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root1, f(root1));
    }

    // Метод хорд
    iter1 = chord_method(0.0, 1.0, &root1);
    if (iter1 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root1, f(root1));
    }

    // Метод дихотомии
    iter1 = bisection_method(0.0, 1.0, &root1);
    if (iter1 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root1, f(root1));
    }

    // Второй корень: в интервале [2, 3]
    printf("\n\nКорень 2: в интервале [2, 3]\n");

    double root2;
    int iter2;

    // Метод простой итерации
    iter2 = simple_iteration(2.0, 3.0, -0.005, &root2);
    if (iter2 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root2, f(root2));
    }

    // Метод Ньютона
    iter2 = newton_method(2.0, 3.0, &root2);
    if (iter2 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root2, f(root2));
    }

    // Метод секущих
    iter2 = secant_method(2.0, 3.0, &root2);
    if (iter2 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root2, f(root2));
    }

    // Метод хорд
    iter2 = chord_method(2.0, 3.0, &root2);
    if (iter2 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root2, f(root2));
    }

    // Метод дихотомии
    iter2 = bisection_method(2.0, 3.0, &root2);
    if (iter2 > 0) {
        printf("Значение функции в корне: f(%.8f) = %.8f\n", root2, f(root2));
    }

    return 0;
}