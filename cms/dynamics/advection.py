import numpy as np
from typing import Tuple
from cms.core.grid import Grid

# Попытка импортировать быстрые ядра, скомпилированные с Numba
try:
    from cms.dynamics.advection_numba import (
        weno5_reconstruct_x, 
        weno5_reconstruct_y, 
        weno5_reconstruct_z, 
        compute_advection_flux
    )
    HAS_NUMBA = True
except ImportError:
    HAS_NUMBA = False

# Попытка импортировать быстрые ядра, написанные для GPU с CUDA
try:
    from numba import cuda
    from cms.dynamics.advection_cuda import (
        weno5_reconstruct_x_kernel,
        weno5_reconstruct_y_kernel,
        weno5_reconstruct_z_kernel,
        compute_advection_flux_kernel
    )
    HAS_CUDA = cuda.is_available()
except ImportError:
    HAS_CUDA = False

class WENO5:
    """
    Класс-обертка для реализации адвекции по схеме WENO 5-го порядка.
    
    WENO (Weighted Essentially Non-Oscillatory) - это современная численная схема
    высокого порядка точности, предназначенная для решения гиперболических уравнений
    (к которым относятся уравнения адвекции). Ее главное преимущество - способность
    сохранять резкие градиенты в полях (например, на границе облака) без
    возникновения численных осцилляций.
    
    Этот класс может использовать разные "бэкенды" для вычислений:
    - CUDA (если доступен GPU и установлен пакет numba.cuda)
    - Numba (компиляция Python-кода в быстрый машинный код для CPU)
    - Pure NumPy (самый медленный, используется как запасной вариант)
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 6.1 & 6.3
    """
    def __init__(self, grid: Grid, use_gpu: bool = False):
        """
        Инициализирует схему адвекции.
        
        Args:
            grid: Объект расчетной сетки.
            use_gpu: Флаг, указывающий на попытку использовать GPU.
        """
        self.grid = grid
        self.use_gpu = use_gpu and HAS_CUDA
        self.use_numba = HAS_NUMBA and not self.use_gpu
        
        if use_gpu and not HAS_CUDA:
            print("ВНИМАНИЕ: Запрошено использование GPU, но CUDA недоступен. Расчет будет на CPU.")

    def _reconstruct_weno5(self, f: np.ndarray, axis: int) -> np.ndarray:
        """
        (Только для NumPy) Выполняет реконструкцию WENO5 на поле f вдоль оси axis.
        Возвращает реконструированные значения на гранях ячеек.
        Эта версия медленная и используется только если Numba и CUDA недоступны.
        """
        # Сдвинутые версии поля для построения 5-точечного шаблона
        f_m2 = np.roll(f, 2, axis=axis)
        f_m1 = np.roll(f, 1, axis=axis)
        f_0  = f
        f_p1 = np.roll(f, -1, axis=axis)
        f_p2 = np.roll(f, -2, axis=axis)

        # Индикаторы гладкости (beta)
        b0 = (13/12) * (f_m2 - 2*f_m1 + f_0)**2 + (1/4) * (f_m2 - 4*f_m1 + 3*f_0)**2
        b1 = (13/12) * (f_m1 - 2*f_0 + f_p1)**2 + (1/4) * (f_m1 - f_p1)**2
        b2 = (13/12) * (f_0 - 2*f_p1 + f_p2)**2 + (1/4) * (3*f_0 - 4*f_p1 + f_p2)**2

        # Идеальные веса для гладкого поля
        d0, d1, d2 = 0.1, 0.6, 0.3

        # Нелинейные веса alpha
        eps = 1e-6 # Малая константа для избежания деления на ноль
        a0 = d0 / (eps + b0)**2
        a1 = d1 / (eps + b1)**2
        a2 = d2 / (eps + b2)**2

        # Финальные веса w
        w_sum = a0 + a1 + a2
        w0, w1, w2 = a0/w_sum, a1/w_sum, a2/w_sum

        # Полиномы-кандидаты 3-го порядка
        q0 = (1/3)*f_m2 - (7/6)*f_m1 + (11/6)*f_0
        q1 = -(1/6)*f_m1 + (5/6)*f_0 + (1/3)*f_p1
        q2 = (1/3)*f_0 + (5/6)*f_p1 - (1/6)*f_p2

        f_reconstructed = w0*q0 + w1*q1 + w2*q2
        return f_reconstructed

    def advect(self, q: np.ndarray, u: np.ndarray, v: np.ndarray, w: np.ndarray) -> np.ndarray:
        """
        Вычисляет тенденцию адвекции для скалярного поля q.
        Формула: -div(q * U), где U - вектор скорости (u, v, w).
        """
        dq_dt = np.zeros_like(q)
        nx, ny, nz = q.shape
        
        # --- Бэкенд на GPU (CUDA) ---
        if self.use_gpu:
            d_q = cuda.to_device(q)
            d_u = cuda.to_device(u)
            d_v = cuda.to_device(v)
            d_w = cuda.to_device(w)
            d_dq_dt = cuda.to_device(dq_dt)
            d_q_recons = cuda.device_array_like(d_q)
            
            threadsperblock = (8, 8, 8)
            blockspergrid_x = (nx + 7) // 8
            blockspergrid_y = (ny + 7) // 8
            blockspergrid_z = (nz + 7) // 8
            blockspergrid = (blockspergrid_x, blockspergrid_y, blockspergrid_z)
            
            weno5_reconstruct_x_kernel[blockspergrid, threadsperblock](d_q, d_q_recons, nx, ny, nz)
            compute_advection_flux_kernel[blockspergrid, threadsperblock](d_q_recons, d_u, d_dq_dt, 0, self.grid.dx, nx, ny, nz)
            
            weno5_reconstruct_y_kernel[blockspergrid, threadsperblock](d_q, d_q_recons, nx, ny, nz)
            compute_advection_flux_kernel[blockspergrid, threadsperblock](d_q_recons, d_v, d_dq_dt, 1, self.grid.dy, nx, ny, nz)
            
            weno5_reconstruct_z_kernel[blockspergrid, threadsperblock](d_q, d_q_recons, nx, ny, nz)
            compute_advection_flux_kernel[blockspergrid, threadsperblock](d_q_recons, d_w, d_dq_dt, 2, self.grid.dz, nx, ny, nz)
            
            d_dq_dt.copy_to_host(dq_dt)
            return dq_dt
        
        # --- Бэкенд на CPU (Numba) ---
        if self.use_numba:
            q_recons = np.empty_like(q)
            
            weno5_reconstruct_x(q, q_recons, nx, ny, nz)
            compute_advection_flux(q_recons, u, dq_dt, 0, self.grid.dx, nx, ny, nz)
            
            weno5_reconstruct_y(q, q_recons, nx, ny, nz)
            compute_advection_flux(q_recons, v, dq_dt, 1, self.grid.dy, nx, ny, nz)
            
            weno5_reconstruct_z(q, q_recons, nx, ny, nz)
            compute_advection_flux(q_recons, w, dq_dt, 2, self.grid.dz, nx, ny, nz)
            
            return dq_dt

        # --- Запасной вариант на Pure NumPy (медленно) ---
        for axis, vel in enumerate([u, v, w]):
            q_star = self._reconstruct_weno5(q, axis)
            dx = [self.grid.dx, self.grid.dy, self.grid.dz][axis]
            flux = vel * q_star
            f_p = flux
            f_m = np.roll(flux, 1, axis=axis)
            dq_dt -= (f_p - f_m) / dx
            
        return dq_dt