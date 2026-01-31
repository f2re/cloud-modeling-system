import numpy as np
from cms.config import G

def compute_terminal_velocity(r, rho_particle, rho_air, eta_air=1.81e-5):
    """
    Вычисляет терминальную скорость оседания сферических частиц по закону Стокса.

    Эта реализация предполагает ламинарный режим (малые числа Рейнольдса),
    что справедливо для мелких аэрозольных частиц. Уравнения 5.17-5.19
    из монографии могут содержать более сложные зависимости.

    Args:
        r (float or np.ndarray): Радиус частиц (м).
        rho_particle (float): Плотность материала частиц (кг/м^3).
        rho_air (float): Плотность воздуха (кг/м^3).
        eta_air (float, optional): Динамическая вязкость воздуха (Па·с).
                                   Значение по умолчанию для ~15°C.

    Returns:
        float or np.ndarray: Терминальная скорость оседания (м/с).
                             Положительное значение означает движение вниз.
    """
    # Закон Стокса для терминальной скорости
    # v_t = (2 * g * r^2 * (rho_p - rho_a)) / (9 * eta)
    
    # Убедимся, что r - это numpy массив для векторных операций
    r = np.asarray(r)

    numerator = 2 * G * r**2 * (rho_particle - rho_air)
    denominator = 9 * eta_air

    # Скорость может быть нулевой, если плотности равны или радиус 0
    if denominator == 0:
        return np.zeros_like(r)

    v_terminal = numerator / denominator

    return v_terminal
