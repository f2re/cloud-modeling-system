import numpy as np
from cms.core.grid import Grid
from cms.dynamics.advection import WENO5
from cms.config import PhysicsConfig

class NavierStokesSolver:
    """
    Решатель уравнений Навье-Стокса для сжимаемой жидкости.
    
    Этот модуль отвечает за вычисление тенденций для полей скорости (u, v, w),
    плотности (rho) и потенциальной температуры (theta), обусловленных
    динамическими процессами: адвекцией, градиентом давления и плавучестью.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 4.1
    """
    def __init__(self, grid: Grid, config: PhysicsConfig, use_gpu: bool = False):
        """
        Инициализирует решатель.
        
        Args:
            grid: Объект расчетной сетки.
            config: Конфигурация с физическими константами.
            use_gpu: Флаг использования GPU (передается в модуль адвекции).
        """
        self.grid = grid
        self.config = config
        self.weno = WENO5(grid, use_gpu=use_gpu)

    def compute_tendencies(self, 
                           u: np.ndarray, v: np.ndarray, w: np.ndarray, 
                           rho: np.ndarray, theta: np.ndarray, 
                           p: np.ndarray) -> tuple:
        """
        Вычисляет тенденции (правые части уравнений) для динамических переменных.
        
        Формула (упрощенно):
        d(phi)/dt = -div(phi * U) + S
        где S - источниковые члены (градиент давления, плавучесть).

        Args:
            u, v, w: 3D поля компонент скорости.
            rho: 3D поле плотности.
            theta: 3D поле потенциальной температуры.
            p: 3D поле давления.

        Returns:
            Кортеж с 3D полями тенденций для (u, v, w, rho, theta).
        """
        # --- 1. Адвективные члены ---
        # Вычисляются с помощью схемы WENO5 для минимизации численной диффузии.
        # Физически это означает перенос самих себя потоком.
        du_dt = self.weno.advect(u, u, v, w)
        dv_dt = self.weno.advect(v, u, v, w)
        dw_dt = self.weno.advect(w, u, v, w)
        drho_dt = self.weno.advect(rho, u, v, w)
        dtheta_dt = self.weno.advect(theta, u, v, w)

        # --- 2. Члены градиента давления ---
        # Сила, действующая на жидкость из-за разницы давлений.
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 4.1 (второй член)
        # grad(p) / rho
        dp_dx = np.gradient(p, self.grid.dx, axis=0)
        dp_dy = np.gradient(p, self.grid.dy, axis=1)
        dp_dz = np.gradient(p, self.grid.dz, axis=2)

        du_dt -= dp_dx / rho
        dv_dt -= dp_dy / rho
        dw_dt -= dp_dz / rho

        # --- 3. Сила плавучести (только для вертикальной скорости w) ---
        # Возникает из-за разницы потенциальной температуры относительно фонового состояния.
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 4.1 (третий член)
        # g * (theta_v - theta_v0) / theta_v0
        
        # УПРОЩЕНИЕ: В текущей версии виртуальная потенциальная температура (theta_v)
        # аппроксимируется обычной потенциальной температурой (theta).
        theta0 = 300.0  # Фоновая (референсная) температура, K
        buoyancy = self.config.g * (theta - theta0) / theta0
        dw_dt += buoyancy

        return du_dt, dv_dt, dw_dt, drho_dt, dtheta_dt