import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig
from cms.microphysics.terminal_velocity import compute_terminal_velocity_double_moment

class WarmMicrophysics:
    """
    Модуль, описывающий микрофизические процессы в теплых облаках (без ледяной фазы).
    Реализована двухмоментная схема, прогнозирующая как массовую долю (q),
    так и числовую концентрацию (N) для облачных и дождевых капель.
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 3.5.1
    """
    def __init__(self, config: PhysicsConfig):
        """
        Инициализирует модуль с физическими константами.
        
        Args:
            config: Датакласс PhysicsConfig.
        """
        self.config = config
        # Параметры для степенного закона терминальной скорости дождевых капель v(D) = a * D^b
        self.a_v_rain = 846.0
        self.b_v_rain = 0.8

    def compute_rates(self, 
                      qc: np.ndarray, qr: np.ndarray, 
                      nc: np.ndarray, nr: np.ndarray, 
                      rho: np.ndarray) -> Tuple[np.ndarray, ...]:
        """
        Вычисляет скорости изменения (тенденции) для всех переменных теплой микрофизики.
        
        Args:
            qc, qr: Массовые доли облачных и дождевых капель (кг/кг).
            nc, nr: Числовые концентрации облачных и дождевых капель (1/м^3).
            rho: Плотность воздуха (кг/м^3).

        Returns:
            Кортеж с тенденциями (dqc, dqr, dnc, dnr).
        """
        # --- 1. Автоконверсия (превращение облачных капель в дождевые) ---
        # Формула основана на схеме Morrison.
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.1
        # (dqr/dt)_auto = 1350 * qc^2.47 * nc^-1.79 / rho^1.47
        auto_rate = np.zeros_like(qc)
        mask = qc > 1e-9 # Процесс активен только при наличии облачных капель

        # Для численной стабильности при работе с большими степенями используется
        # логарифмическое представление. Ref: IMPLEMENTATION_GUIDE.md, Section 6.4.1
        log_qc = np.log(np.maximum(qc[mask], 1e-15))
        log_nc = np.log(np.maximum(nc[mask], 1e-3))
        log_rho = np.log(np.maximum(rho[mask], 0.1))

        log_auto = (np.log(1350.0) +
                    2.47 * log_qc -
                    1.79 * log_nc -
                    1.47 * log_rho)
        
        # Ограничение результата для предотвращения экстремальных значений
        log_auto_safe = np.clip(log_auto, -50, 10)
        auto_rate[mask] = np.exp(log_auto_safe)
        
        dqr_auto = auto_rate   # Увеличение массы дождя
        dqc_auto = -auto_rate  # Уменьшение массы облаков

        # --- 2. Аккреция (сбор облачных капель дождевыми) ---
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.1
        # (dqr/dt)_accr = (pi/4) * E_cr * nr * qc * Dr^2 * vr
        e_cr = 1.0  # Коэффициент столкновения (упрощенно)
        
        # Защита от отрицательных значений и деления на ноль
        qr_safe = np.maximum(qr, 0.0)
        nr_safe = np.maximum(nr, 1e-12)
        
        # Средний диаметр дождевых капель (Dr)
        mass_conc_r = qr_safe * rho
        m_r = mass_conc_r / nr_safe
        dr = np.power(6.0 * m_r / (np.pi * self.config.rho_w), 1.0/3.0)
        dr[m_r < 1e-18] = 0.0

        # Терминальная скорость дождевых капель (vr)
        vr = compute_terminal_velocity_double_moment(
            qr, nr, rho, self.config.rho_w,
            self.a_v_rain, self.b_v_rain
        )
        
        accr_rate = (np.pi / 4.0) * e_cr * nr * qc * dr**2 * vr
        
        dqr_accr = accr_rate   # Увеличение массы дождя
        dqc_accr = -accr_rate  # Уменьшение массы облаков

        # --- 3. Расчет тенденций для числовой концентрации ---
        # Этот блок был исправлен для обеспечения сохранения массы и стабильности
        
        # Масса одной новой дождевой капли (определяется минимальным радиусом из конфига)
        m_rain_new = (4.0/3.0) * np.pi * self.config.r_rain_min**3 * self.config.rho_w
        
        # Средняя масса одной облачной капли
        qc_safe = np.maximum(qc, 1e-12)
        nc_safe_tend = np.maximum(nc, 1e-3)
        m_c = (qc_safe * rho) / nc_safe_tend
        m_c_inv = 1.0 / np.maximum(m_c, 1e-15)  # Используем обратную массу для избежания деления

        # Прирост числа дождевых капель от автоконверсии.
        # dN/dt = (dq/dt * rho) / m_particle, где dq/dt - тенденция массовой доли (кг/кг/с)
        dnr_auto = (dqr_auto * rho) / m_rain_new
        
        # Убыль числа облачных капель от автоконверсии и аккреции
        dnc_auto = (dqc_auto * rho) * m_c_inv
        dnc_accr = (dqc_accr * rho) * m_c_inv
        
        # --- 4. Самосбор дождевых капель (уменьшение их числа при слиянии) ---
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.5.1
        # (dnr/dt)_self = -5.78 * nr^2 * dr^3
        dnr_self = -5.78 * nr_safe**2 * dr**3
        
        # --- 5. Объединение всех тенденций ---
        dqc = dqc_auto + dqc_accr
        dqr = dqr_auto + dqr_accr
        dnc = dnc_auto + dnc_accr
        dnr = dnr_auto + dnr_self

        return dqc, dqr, dnc, dnr