import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# Задаем систему уравнений
def equations(vars):
    x, y = vars
    eq1 = 2 * np.sin(x) - 2 * np.cos(y) + 1  # 2sinx - 2cosy = -1
    eq2 = x**2 + x*y + y**2 + 2*x - 1       # x² + xy + y² + 2x - 1 = 0
    return [eq1, eq2]

# Создаем более детальную сетку значений x и y
x = np.linspace(-4, 4, 800)
y = np.linspace(-3, 3, 800)
X, Y = np.meshgrid(x, y)

# Вычисляем уравнения системы
Z1 = 2*np.sin(X) - 2*np.cos(Y) + 1  # 2sin(x) - 2cos(y) = -1
Z2 = X**2 + X*Y + Y**2 + 2*X - 1    # x² + xy + y² + 2x - 1 = 0

# Создаем график
plt.figure(figsize=(14, 10))

# Первое уравнение: 2sin(x) - 2cos(y) = -1
# Используем contour для плавных линий
contour1 = plt.contour(X, Y, Z1, levels=[0], colors='red', linewidths=2.5)

# Второе уравнение: x² + xy + y² + 2x - 1 = 0
contour2 = plt.contour(X, Y, Z2, levels=[0], colors='blue', linewidths=2.5)

# Создаем кастомные элементы для легенды
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='red', lw=3, label='2sin(x) - 2cos(y) = -1'),
    Line2D([0], [0], color='blue', lw=3, label='x² + xy + y² + 2x - 1 = 0'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green',
           markersize=10, markeredgecolor='black', label='Найденные корни')
]

# Находим начальные приближения (точки пересечения)
initial_guesses = [
    [-2, -2], [-2, 0], [-2, 2],
    [0, -2], [0, 0], [0, 2],
    [2, -2], [2, 0], [2, 2]
]

solutions = []
for guess in initial_guesses:
    try:
        solution, info, ier, msg = fsolve(equations, guess, full_output=True)
        if ier == 1:  # Решение найдено успешно
            # Проверяем, что решение действительно удовлетворяет уравнениям
            if abs(equations(solution)[0]) < 1e-8 and abs(equations(solution)[1]) < 1e-8:
                # Проверяем на уникальность
                is_unique = True
                for sol in solutions:
                    if np.linalg.norm(np.array(sol) - np.array(solution)) < 0.1:
                        is_unique = False
                        break
                if is_unique:
                    solutions.append(solution)
                    # Отмечаем точку на графике
                    plt.plot(solution[0], solution[1], 'go', markersize=10,
                            markeredgecolor='black', zorder=5)
                    plt.text(solution[0] + 0.15, solution[1] + 0.15,
                            f'({solution[0]:.2f}, {solution[1]:.2f})',
                            fontsize=11,
                            bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.8),
                            zorder=6)
    except Exception as e:
        continue

# Дополнительная визуализация - заливка для лучшего понимания
plt.contourf(X, Y, Z1, levels=[-10, 0], colors=['red'], alpha=0.1)
plt.contourf(X, Y, Z1, levels=[0, 10], colors=['lightcoral'], alpha=0.1)
plt.contourf(X, Y, Z2, levels=[-10, 0], colors=['blue'], alpha=0.1)
plt.contourf(X, Y, Z2, levels=[0, 10], colors=['lightblue'], alpha=0.1)

plt.xlabel('x', fontsize=12)
plt.ylabel('y', fontsize=12)
plt.title('График системы нелинейных уравнений\nТочки пересечения - начальные приближения для методов',
          fontsize=14, pad=20)
plt.grid(True, alpha=0.3)
plt.axhline(0, color='black', linewidth=0.8)
plt.axvline(0, color='black', linewidth=0.8)
plt.legend(handles=legend_elements, loc='upper right', fontsize=11)
plt.axis('equal')
plt.xlim(-3, 4)
plt.ylim(-2, 2)

plt.text(0.02, 0.98, f'Найдено корней: {len(solutions)}',
         transform=plt.gca().transAxes, fontsize=12,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
         verticalalignment='top')

plt.show()
