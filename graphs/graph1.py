import numpy as np
import matplotlib.pyplot as plt

def f(x):
    return 0.25 * np.sinh(x**2) - 5*x**2 - 25*x + 21

# Создаем сетку
x = np.linspace(-1, 3, 1000)
y = f(x)

root_intervals = []
for i in range(len(x) - 1):
    if y[i] * y[i + 1] < 0:  # смена знака = корень между этими точками
        a, b = x[i], x[i + 1]
        root_intervals.append((a, b))

# Находим приближенные решения (середины интервалов)
solutions = []
for a, b in root_intervals:
    solution = (a + b) / 2
    solutions.append(round(solution, 2))  # округляем до 2 знаков

# Строим график
plt.figure(figsize=(12, 8))
plt.plot(x, y, 'b-', linewidth=3, label=r'$0.25 \cdot \sinh(x^2) - 5x^2 - 25x + 21 = 0$')
plt.axhline(y=0, color='black', linestyle='-', linewidth=1)

# Сетка и оформление
plt.grid(True, alpha=0.3)
plt.xlabel('x', fontsize=12)
plt.ylabel('f(x)', fontsize=12)
plt.title('График уравнения sh(x²) - 5x² - 25x + 21 = 0', fontsize=14)

from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], color='blue', lw=3, label='sinh(x^2) - 5x^2 - 25x + 21 = 0'),
    Line2D([0], [0], marker='o', color='w', markerfacecolor='green',
           markersize=10, markeredgecolor='black', label='Найденные корни')
]

# Определяем интервалы корней
root_intervals = [
    (0, 1),
    (2, 3)
]

# Отмечаем интервалы
colors = ['hotpink', 'mediumorchid']
for i, (a, b) in enumerate(root_intervals):
    plt.axvspan(a, b, alpha=0.2, color=colors[i])

# Находим и отмечаем решения
for i, sol in enumerate(solutions):
    plt.plot(sol, 0, 'ro', markersize=10, markeredgecolor='black', zorder=5)
    plt.text(sol + 0.1, 0.5, f'({sol}, 0)', fontsize=11,
             bbox=dict(boxstyle="round,pad=0.3", facecolor="yellow", alpha=0.8))

plt.text(0.02, 0.98, 'Найдено корней: 2',
         transform=plt.gca().transAxes, fontsize=12,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8),
         verticalalignment='top')

plt.text(0.02, 0.90, 'Интервалы корней:\n[0, 1]\n[2, 3]',
         transform=plt.gca().transAxes, fontsize=11,
         bbox=dict(boxstyle="round,pad=0.3", facecolor="lightblue", alpha=0.8),
         verticalalignment='top')

plt.ylim(-45, 45)
plt.tight_layout()
plt.show()
