import numpy as np

def compute_wind_profile(z, z0, u_star, L, kappa=0.4, gamma=16, beta=5):
    """
    Вычисляет логарифмический профиль ветра с учетом стратификации атмосферы
    по теории подобия Монина-Обухова (уравнения 5.11-5.14 из монографии).

    Ref: Businger-Dyer relations.

    Args:
        z (float or np.ndarray): Высота(ы) над поверхностью (м).
        z0 (float): Параметр шероховатости (м).
        u_star (float): Динамическая скорость (скорость трения) (м/с).
        L (float): Длина Монина-Обухова (м).
        kappa (float, optional): Постоянная фон Кармана. Defaults to 0.4.
        gamma (float, optional): Константа для неустойчивых условий. Defaults to 16.
        beta (float, optional): Константа для устойчивых условий. Defaults to 5.

    Returns:
        float or np.ndarray: Скорость ветра на высоте z (м/с).
    """
    # Убедимся, что z - это numpy массив для векторных операций
    z = np.asarray(z)
    
    # Инициализируем psi_m как массив нулей той же формы, что и z
    psi_m = np.zeros_like(z, dtype=float)
    
    zeta = z / L

    # Векторизованная обработка условий
    # Неустойчивые условия (L < 0 -> zeta < 0)
    unstable_mask = zeta < 0
    if np.any(unstable_mask):
        zeta_unstable = zeta[unstable_mask]
        x = (1 - gamma * zeta_unstable)**0.25
        psi_m[unstable_mask] = (2 * np.log((1 + x) / 2) + 
                              np.log((1 + x**2) / 2) - 
                              2 * np.arctan(x) + np.pi / 2)

    # Устойчивые условия (L > 0 -> zeta > 0)
    stable_mask = zeta > 0
    if np.any(stable_mask):
        zeta_stable = zeta[stable_mask]
        psi_m[stable_mask] = -beta * zeta_stable

    # Основное уравнение профиля ветра
    # (предполагаем, что высота смещения d=0)
    # np.log(z / z0) вызовет ошибку, если z=0. z0 должно быть > 0.
    # Добавим небольшое смещение для z, чтобы избежать log(0), если z - массив, начинающийся с 0.
    wind_speed = (u_star / kappa) * (np.log(np.maximum(z, z0) / z0) - psi_m)
    
    # Ветер не может быть отрицательным
    return np.maximum(0, wind_speed)
