import numpy as np
from cms.core.grid import Grid
from cms.config import PhysicsConfig

class Turbulence:
    """
    Модуль для расчета турбулентной диффузии и дисперсии.
    Использует модель крупных вихрей (Large Eddy Simulation, LES)
    с параметризацией Смагоринского для вычисления подсеточной вязкости.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 4.2, 6.5
    """
    def __init__(self, grid: Grid, config: PhysicsConfig):
        """
        Инициализирует модуль турбулентности.
        
        Args:
            grid: Объект расчетной сетки.
            config: Конфигурация с физическими константами.
        """
        self.grid = grid
        self.config = config
        self.cs_sq = 0.18**2  # Квадрат константы Смагоринского

    def compute_eddy_viscosity(self, u: np.ndarray, v: np.ndarray, w: np.ndarray) -> np.ndarray:
        """
        Вычисляет подсеточную турбулентную вязкость (nu_t) по модели Смагоринского.
        
        Для обеспечения численной стабильности добавлены ограничители (клиппинг)
        на градиенты скорости и на итоговое значение вязкости.
        
        Ref: IMPLEMENTATION_GUIDE.md, Eq 4.2, Section 6.5
        Формула: nu_t = (Cs * Delta)^2 * |S|, где |S| - норма тензора скоростей деформации.

        Args:
            u, v, w: 3D поля компонент скорости.

        Returns:
            3D поле турбулентной вязкости (м^2/с).
        """
        # Ширина фильтра, равная среднему размеру ячейки сетки
        delta = (self.grid.dx * self.grid.dy * self.grid.dz)**(1/3)
        
        # --- Стабилизация: Ограничение градиентов скорости ---
        # Это предотвращает "взрыв" квадратичных членов при резких градиентах,
        # например, на фронтах или в сильных сдвиговых потоках.
        max_grad = 1.0  # Максимальный градиент скорости (1/с), физически обоснованный предел
        du_dx = np.clip(np.gradient(u, self.grid.dx, axis=0), -max_grad, max_grad)
        dv_dy = np.clip(np.gradient(v, self.grid.dy, axis=1), -max_grad, max_grad)
        dw_dz = np.clip(np.gradient(w, self.grid.dz, axis=2), -max_grad, max_grad)
        
        # Упрощенная норма тензора скоростей деформации (только диагональные члены)
        s_sq = 2 * (du_dx**2 + dv_dy**2 + dw_dz**2)
        
        # Вычисление турбулентной вязкости
        nu_t = (self.cs_sq * delta**2) * np.sqrt(s_sq)
        
        # --- Стабилизация: Ограничение итоговой вязкости ---
        # Предотвращает получение нереалистично больших значений nu_t.
        nu_t_max = 1000.0  # Физически реалистичный максимум для атмосферы (м^2/с)
        nu_t = np.clip(nu_t, 0.0, nu_t_max)
        
        return nu_t

    def compute_diffusion(self, field: np.ndarray, nu_t: np.ndarray) -> np.ndarray:
        """
        Вычисляет тенденцию для скалярного поля из-за турбулентной диффузии.
        
        Физический смысл: сглаживание поля со скоростью, пропорциональной
        турбулентной вязкости. Математически это дивергенция диффузионного потока.
        
        Формула: div(nu_t * grad(field))

        Args:
            field: 3D скалярное поле (например, qv, theta).
            nu_t: 3D поле турбулентной вязкости.

        Returns:
            3D поле тенденции (скорости изменения) поля `field`.
        """
        # 1. Вычислить градиент поля
        grad_x = np.gradient(field, self.grid.dx, axis=0)
        grad_y = np.gradient(field, self.grid.dy, axis=1)
        grad_z = np.gradient(field, self.grid.dz, axis=2)
        
        # 2. Вычислить диффузионные потоки (K-теория)
        flux_x = nu_t * grad_x
        flux_y = nu_t * grad_y
        flux_z = nu_t * grad_z
        
        # 3. Вычислить дивергенцию потоков
        div_flux = (np.gradient(flux_x, self.grid.dx, axis=0) +
                    np.gradient(flux_y, self.grid.dy, axis=1) +
                    np.gradient(flux_z, self.grid.dz, axis=2))
                    
        return div_flux