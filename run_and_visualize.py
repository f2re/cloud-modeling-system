"Скрипт для проведения моделирования и создания 3D анимации.

Этот скрипт выполняет следующие шаги:
1. Инициализирует и запускает CMS модель на заданный промежуток времени.
2. Сохраняет 3D поля скорости и концентрации реагента на каждом n-ом шаге.
3. "Выпускает" виртуальные частицы и рассчитывает их траектории на основе полей скорости.
4. Создает и сохраняет 3D анимацию, показывающую эволюцию поля концентрации
   и движение частиц, с вращением камеры для наглядности.
"
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits.mplot3d import Axes3D
from scipy.interpolate import RegularGridInterpolator
import os

from cms.model import CMSModel
from cms.config import GridConfig, PhysicsConfig

# --- 1. Настройка параметров моделирования и визуализации ---

print("Шаг 1: Настройка параметров...")

# Параметры сетки (уменьшенная для скорости демонстрации)
 g_config = GridConfig(nx=30, ny=30, nz=30, dx=100, dy=100, dz=100)
 p_config = PhysicsConfig()

# Параметры симуляции
 sim_duration_minutes = 10 # Увеличено время моделирования
 dt_s = 2.0                # Используем большой шаг благодаря IMEX
 num_steps = int(sim_duration_minutes * 60 / dt_s)
 save_interval = 5         # Сохранять каждый 5-й кадр

# Параметры визуализации
 N_PARTICLES = 50          # Количество трассируемых частиц
 OUTPUT_FILENAME = "seeding_animation.mp4"

# --- 2. Проведение моделирования ---

print(f"Шаг 2: Запуск моделирования ({num_steps} шагов)...")

model = CMSModel(g_config, p_config)

# Начальные условия: небольшой сдвиг ветра для интересного движения
 model.u += 2.0  # м/с
 model.v += 1.0  # м/с

# Точка засева в центре
 mid_x, mid_y, mid_z = g_config.nx // 2, g_config.ny // 2, g_config.nz // 2
 model.c_agi[mid_x, mid_y, mid_z] = 1.0 # Начальная концентрация реагента

# Списки для хранения истории полей
 history_c_agi = []
 history_u, history_v, history_w = [], [], []

for step in range(num_steps):
    if step % 10 == 0:
        print(f"  Шаг {step}/{num_steps}...")
     model.step(dt=dt_s)
    
     if step % save_interval == 0:
        # Сохраняем внутреннюю часть сетки без "призрачных" ячеек
         history_c_agi.append(model.grid.get_inner(model.c_agi).copy())
         history_u.append(model.grid.get_inner(model.u).copy())
         history_v.append(model.grid.get_inner(model.v).copy())
         history_w.append(model.grid.get_inner(model.w).copy())

num_frames = len(history_c_agi)
print(f"Моделирование завершено. Сохранено {num_frames} кадров.")

# --- 3. Расчет траекторий частиц ---

print("Шаг 3: Расчет траекторий частиц...")

# Начальные позиции частиц вблизи точки засева
 np.random.seed(42)
 particles = np.random.rand(N_PARTICLES, 3) * 2 * g_config.dx +
            np.array([mid_x*g_config.dx, mid_y*g_config.dy, mid_z*g_config.dz]) - g_config.dx

 particle_trajectories = np.zeros((num_frames, N_PARTICLES, 3))
 particle_trajectories[0, :, :] = particles

# Оси координат для интерполяции
 x = model.grid.x
 y = model.grid.y
 z = model.grid.z

for i in range(1, num_frames):
    # Создаем интерполяторы для полей скорости на предыдущем шаге
     interp_u = RegularGridInterpolator((x, y, z), history_u[i-1])
     interp_v = RegularGridInterpolator((x, y, z), history_v[i-1])
     interp_w = RegularGridInterpolator((x, y, z), history_w[i-1])
    
    # Получаем скорости в точках нахождения частиц
     u_vel = interp_u(particles)
     v_vel = interp_v(particles)
     w_vel = interp_w(particles)
    
    # Обновляем позиции частиц (простая схема Эйлера)
     dt_anim = dt_s * save_interval
     particles[:, 0] += u_vel * dt_anim
     particles[:, 1] += v_vel * dt_anim
     particles[:, 2] += w_vel * dt_anim
    
    # Записываем новые позиции
     particle_trajectories[i, :, :] = particles

print("Траектории рассчитаны.")

# --- 4. Создание 3D анимации ---

print("Шаг 4: Создание 3D анимации (это может занять несколько минут)...")

fig = plt.figure(figsize=(12, 9))
 ax = fig.add_subplot(111, projection='3d')

# Координаты сетки для отображения
 X, Y, Z = np.meshgrid(x, y, z, indexing='ij')

# Порог для отображения поля концентрации (чтобы не загромождать график)
 c_threshold = 0.01

# Инициализация объектов графика
# Поле концентрации (полупрозрачные точки)
 c_agi_plot = ax.scatter([], [], [], c=[], s=[], alpha=0.3, cmap='viridis')
# Частицы (яркие точки)
 particles_plot = ax.scatter([], [], [], c='red', s=40, marker='o', depthshade=False, label='Частицы')

def init():
    """Инициализация фона анимации."""
     ax.set_xlim(0, g_config.nx * g_config.dx)
     ax.set_ylim(0, g_config.ny * g_config.dy)
     ax.set_zlim(0, g_config.nz * g_config.dz)
     ax.set_xlabel('X (м)')
     ax.set_ylabel('Y (м)')
     ax.set_zlabel('Z (м)')
     ax.set_title('3D Визуализация поля концентрации и траекторий частиц')
     ax.legend()
     return c_agi_plot, particles_plot

def update(frame):
    """Обновление кадра анимации."""
     print(f"  Рендеринг кадра {frame+1}/{num_frames}...")
    
    # --- Обновление поля концентрации ---
     c_data = history_c_agi[frame]
     mask = c_data > c_threshold
    
    # Обновляем данные scatter plot'а
     c_agi_plot._offsets3d = (X[mask], Y[mask], Z[mask])
    
    # Размер и цвет точек пропорциональны концентрации
     sizes = c_data[mask] * 2000
     colors = c_data[mask]
     c_agi_plot.set_sizes(sizes)
     c_agi_plot.set_array(colors)
    
    # --- Обновление позиций частиц ---
     particle_pos = particle_trajectories[frame]
     particles_plot._offsets3d = (particle_pos[:, 0], particle_pos[:, 1], particle_pos[:, 2])

    # --- Вращение камеры ---
    # Вращаем по оси Z (azimuth)
     ax.view_init(elev=20., azim=frame * 360 / num_frames)
    
    # Обновление заголовка
     time_minutes = frame * dt_anim / 60
     ax.set_title(f'Время: {time_minutes:.1f} мин. Кадр: {frame+1}/{num_frames}')
    
     return c_agi_plot, particles_plot

# Создание и сохранение анимации
 ani = FuncAnimation(fig, update, frames=num_frames, init_func=init, blit=False)

# Сохранение в файл. Требует наличия ffmpeg.
# dpi - качество, fps - кадры в секунду
try:
    ani.save(OUTPUT_FILENAME, writer='ffmpeg', fps=10, dpi=150, progress_callback=lambda i, n: print(f"Сохранение кадра {i+1}/{n}"))
    print(f"\nАнимация успешно сохранена в файл: {os.path.abspath(OUTPUT_FILENAME)}")
except Exception as e:
    print(f"\nОшибка при сохранении анимации: {e}")
    print("Убедитесь, что у вас установлен `ffmpeg`. Попробуйте `conda install ffmpeg` или `sudo apt-get install ffmpeg`.")

plt.close(fig) # Закрываем фигуру, чтобы она не отображалась в блокноте
