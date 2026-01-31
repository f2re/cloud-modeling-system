import numpy as np
from typing import Tuple
from cms.config import PhysicsConfig

class Activation:
    """
    Модуль, описывающий взаимодействие аэрозолей с облаками,
    в частности, процесс активации ядер конденсации облаков (CCN).
    
    Ref: IMPLEMENTATION_GUIDE.md, Section 3.6
    """
    def __init__(self, config: PhysicsConfig):
        """
        Инициализирует модуль с физическими константами.
        
        Args:
            config: Датакласс PhysicsConfig.
        """
        self.config = config
        # Константы для континентального типа аэрозоля
        # Ref: IMPLEMENTATION_GUIDE.md, Section 3.6
        self.k1 = 0.7
        self.k2 = 1.5
        
    def compute_activation(self, 
                           w: np.ndarray, 
                           n_ccn: np.ndarray, 
                           nc: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
        """
        Вычисляет скорость активации капель (образования новых облачных капель на CCN).
        
        Процесс зависит от максимального пересыщения (S_max), которое, в свою очередь,
        зависит от вертикальной скорости и концентрации CCN.

        Args:
            w: 3D поле вертикальной скорости (м/с).
            n_ccn: 3D поле концентрации ядер конденсации (1/м^3).
            nc: 3D поле концентрации облачных капель (1/м^3).

        Returns:
            Кортеж с тенденциями для концентрации облачных капель и CCN (dnc_dt, dn_ccn_dt).
        """
        # 1. Вычисление максимального пересыщения (S_max)
        # Активация происходит только в восходящих потоках (w > 0)
        w_up = np.maximum(w, 0.0)
        
        # Защита от деления на ноль для n_ccn
        n_ccn_safe = np.maximum(n_ccn, 1e6) # Минимум 1 ядро/см^3
        
        # Коэффициенты из параметризации Твоми, аппроксимирующей S_max
        c1 = 100.0 # Масштабирующий коэффициент
        
        # S_max пропорционально w^(3/4) * N_ccn^(-1/2)
        s_max = c1 * (w_up**0.75) * (n_ccn_safe**-0.5)
        
        # 2. Вычисление числа активированных капель (N_act)
        # Ref: IMPLEMENTATION_GUIDE.md, Eq 3.6
        # N_c_act = C2 * N_ccn^k1 * S_max^k2
        c2 = 1.0 # Масштабирующий коэффициент
        n_act = c2 * (n_ccn_safe**self.k1) * (s_max**self.k2)
        
        # Количество активированных капель не может превышать доступное количество CCN
        n_act = np.minimum(n_act, n_ccn_safe)
        
        # 3. Вычисление скорости активации
        # Новые капли активируются только если целевое число n_act больше текущего nc.
        # Используется схема релаксации к целевому значению.
        tau = 1.0 # Характерное время активации (с)
        
        activation_tendency = np.maximum(n_act - nc, 0.0) / tau
        
        # dnc_dt - скорость появления новых облачных капель
        dnc_dt = activation_tendency
        # dn_ccn_dt - скорость истощения аэрозольных частиц
        dn_ccn_dt = -activation_tendency
        
        return dnc_dt, dn_ccn_dt