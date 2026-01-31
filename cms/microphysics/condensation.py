import numpy as np
from cms.config import PhysicsConfig

class CondensationEvaporation:
    """
    Явные процессы конденсации и испарения (уравнения 5.32-5.35 монографии)
    Ref: IMPLEMENTATION_GUIDE.md Section 3 (Implicitly required)
    """
    def __init__(self, config: PhysicsConfig):
        self.config = config
        self.Rv = 461.5  # Газовая постоянная водяного пара, Дж/(кг*К)
        self.D_v = 2.26e-5  # Коэффициент диффузии водяного пара в воздухе, м^2/с

    def compute_saturation_vapor_pressure(self, T: np.ndarray) -> np.ndarray:
        """
        Рассчитывает давление насыщенного пара над водой.
        Формула Тетенса.
        T in Kelvin.
        Returns pressure in Pa.
        """
        # Добавлена защита от численного переполнения
        T_c = np.clip(T - 273.15, -80.0, 80.0)
        return 611.2 * np.exp(17.67 * T_c / (T_c + 243.5))

    def compute_condensation_rate(self, qv: np.ndarray, qc: np.ndarray, T: np.ndarray, p: np.ndarray, rho: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """
        Рассчитывает скорость конденсации/испарения с использованием
        устойчивой схемы релаксации к насыщению.
        Возвращает кортеж (dqc_dt, dqv_dt).
        """
        epsilon = self.config.rd / self.Rv # 0.622
        
        # 1. Рассчитать равновесную массовую долю водяного пара (q_sat)
        e_sat = self.compute_saturation_vapor_pressure(T)
        q_sat = epsilon * e_sat / np.maximum(p - e_sat, 1e-5)

        # 2. Определить отклонение от насыщения
        delta_q = qv - q_sat
        
        # 3. Применить схему релаксации
        # Положительное delta_q -> пересыщение (конденсация)
        # Отрицательное delta_q -> недосыщение (испарение)
        tau = 10.0 # Время релаксации к равновесию [с]
        
        # Скорость изменения [кг/кг / с]
        dqv_dt = -delta_q / tau
        
        # 4. Ограничить скорость, чтобы избежать отрицательных значений
        # Не испарять больше облачной воды, чем есть (qc)
        dqv_dt = np.where(dqv_dt > 0, np.minimum(dqv_dt, qc / tau), dqv_dt)
        # Не конденсировать больше водяного пара, чем есть (qv)
        dqv_dt = np.where(dqv_dt < 0, np.maximum(dqv_dt, -qv / tau), dqv_dt)
        
        dqc_dt = -dqv_dt
        
        return dqc_dt, dqv_dt