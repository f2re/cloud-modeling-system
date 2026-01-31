import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig
from cms.microphysics.terminal_velocity import compute_terminal_velocity_double_moment

class WarmMicrophysics:
    """
    Double-moment warm rain microphysics (Kessler/Morrison style).
    Ref: IMPLEMENTATION_GUIDE.md Section 3.5.1
    """
    def __init__(self, config: PhysicsConfig):
        self.config = config
        # Parameters for rain terminal velocity v(D) = a * D^b
        self.a_v_rain = 846.0
        self.b_v_rain = 0.8

    def compute_rates(self, 
                      qc: np.ndarray, qr: np.ndarray, 
                      nc: np.ndarray, nr: np.ndarray, 
                      rho: np.ndarray) -> Tuple[np.ndarray, ...]:
        """
        Computes transformation rates for warm microphysics.
        Returns (dqc, dqr, dnc, dnr) tendencies.
        """
        # 1. Autoconversion (cloud -> rain) using stable logarithmic formulation
        # Eq 3.5.1: (dqr/dt)_auto = 1350 * qc^2.47 * nc^-1.79 / rho^1.47
        auto_rate = np.zeros_like(qc)
        mask = qc > 1e-9

        # Используем логарифмическое представление для стабильности
        log_qc = np.log(np.maximum(qc[mask], 1e-15))
        log_nc = np.log(np.maximum(nc[mask], 1e-3))
        log_rho = np.log(np.maximum(rho[mask], 0.1))

        log_auto = (np.log(1350.0) +
                    2.47 * log_qc -
                    1.79 * log_nc -
                    1.47 * log_rho)
        
        # Ограничиваем результат для предотвращения экстремальных значений
        log_auto_safe = np.clip(log_auto, -50, 10)
        auto_rate[mask] = np.exp(log_auto_safe)
        
        dqr_auto = auto_rate
        dqc_auto = -auto_rate
        
        # 2. Accretion (cloud by rain)
        # Eq 3.5.1: (dqr/dt)_accr = (pi/4) * E_cr * nr * qc * Dr^2 * vr
        e_cr = 1.0
        
        # Safety: Clamp inputs
        qr_safe = np.maximum(qr, 0.0)
        nr_safe = np.maximum(nr, 1e-12)
        
        # Mean diameter of rain drops
        mass_conc_r = qr_safe * rho
        m_r = mass_conc_r / nr_safe
        dr = np.power(6.0 * m_r / (np.pi * self.config.rho_w), 1.0/3.0)
        dr[m_r < 1e-18] = 0.0

        # Terminal velocity of rain drops
        vr = compute_terminal_velocity_double_moment(
            qr, nr, rho, self.config.rho_w,
            self.a_v_rain, self.b_v_rain
        )
        
        accr_rate = (np.pi / 4.0) * e_cr * nr * qc * dr**2 * vr
        
        dqr_accr = accr_rate
        dqc_accr = -accr_rate

        # 3. Рассчитать тенденции для числовой концентрации
        
        # Масса одной новой дождевой капли (из минимального радиуса)
        m_rain_new = (4.0/3.0) * np.pi * self.config.r_rain_min**3 * self.config.rho_w
        
        # Средняя масса облачной капли
        qc_safe = np.maximum(qc, 1e-12)
        nc_safe_tend = np.maximum(nc, 1e-3)
        m_c = (qc_safe * rho) / nc_safe_tend
        m_c_inv = 1.0 / np.maximum(m_c, 1e-15) # Обратная масса для избежания деления

        # Прирост числа дождевых капель от автоконверсии
        # dN/dt = (dq/dt * rho) / m_particle
        dnr_auto = (dqr_auto * rho) / m_rain_new
        
        # Убыль числа облачных капель от автоконверсии и аккреции
        dnc_auto = (dqc_auto * rho) * m_c_inv
        dnc_accr = (dqc_accr * rho) * m_c_inv
        
        # 4. Self-collection (уменьшение числа дождевых капель)
        # Eq 3.5.1: (dnr/dt)_self = -5.78 * nr^2 * dr^3
        dnr_self = -5.78 * nr_safe**2 * dr**3
        
        # 5. Объединить тенденции
        dqc = dqc_auto + dqc_accr
        dqr = dqr_auto + dqr_accr
        dnc = dnc_auto + dnc_accr
        dnr = dnr_auto + dnr_self

        return dqc, dqr, dnc, dnr
