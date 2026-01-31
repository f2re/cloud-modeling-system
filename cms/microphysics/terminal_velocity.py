import numpy as np
from scipy.special import gamma as gamma_func

def compute_terminal_velocity_double_moment(qx, nx, rho, rho_x, 
                                             a_v, b_v, alpha=2.0):
    """
    Вычисляет масс-взвешенную терминальную скорость для гидрометеоров,
    описываемых гамма-распределением.
    
    Ref: IMPLEMENTATION_GUIDE раздел 3.4
         Morrison et al. (2005)
         
    Args:
        qx (np.ndarray): Массовая концентрация (mixing ratio) [кг/кг].
        nx (np.ndarray): Числовая концентрация [м⁻³].
        rho (np.ndarray): Плотность воздуха [кг/м⁻³].
        rho_x (float): Плотность вещества гидрометеора (воды, льда) [кг/м⁻³].
        a_v (float): Параметр степенного закона для скорости v(D) = a*D^b.
        b_v (float): Параметр степенного закона для скорости v(D) = a*D^b.
        alpha (float, optional): Параметр формы для гамма-распределения. Defaults to 2.0.

    Returns:
        np.ndarray: Масс-взвешенная терминальная скорость [м/с].
    """
    # Избегаем деления на ноль
    nx_safe = np.maximum(nx, 1e-12)
    qx_safe = np.maximum(qx, 0.0)
    
    # Средняя масса частицы
    # m_x = rho * qx / nx -> Эта формула неверна, т.к. qx - безразмерный mixing ratio
    # Правильно: qx * rho - массовая концентрация в кг/м^3
    mass_concentration_x = qx_safe * rho
    m_x = mass_concentration_x / nx_safe
    
    # Средний диаметр
    D_x = np.zeros_like(qx)
    mask = m_x > 1e-18
    D_x[mask] = (6.0 * m_x[mask] / (np.pi * rho_x))**(1.0/3.0)
    
    # Корректировка на плотность воздуха (для падения в разреженной атмосфере)
    rho_corr = np.power(rho / 1.225, -0.5) # (rho_0/rho)^0.5
    
    # Масс-взвешенная скорость
    v_tx = a_v * np.power(D_x, b_v) * rho_corr
    
    # Корректировка на гамма-распределение
    # множитель = Gamma( (4+b)/alpha ) / Gamma( 4/alpha )
    # где 4 = 3 (для объема) + 1
    if alpha > 0:
        factor = gamma_func((4.0 + b_v) / alpha) / gamma_func(4.0 / alpha)
        v_tx *= factor
    
    return v_tx
